/*! \brief Functors for various methods
 * \author PICKSC
 * \file   clock.cpp
 *
 * Includes spatial advection, electric field advection, and electric field update routines
 *
 */
//--------------------------------------------------------------
//  Standard libraries
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <mpi.h>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

//  Declerations
#include "state.h"
#include "input.h"
#include "clock.h"

//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Stepper::Stepper(double starttime, double __dt, double abs_tol, double rel_tol, size_t _maxfails): 
    current_time(starttime), dt_next(0.5*__dt), _dt(0.5*__dt),
    atol(abs_tol), rtol(rel_tol), 
    acceptability(0.), err_val(0.), 
    failed_steps(0), max_failures(_maxfails),overall_check(true), _success(false),
    Nbc(Input::List().BoundaryCells), world_rank(0), world_size(1)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        acceptabilitylist = new double[world_size];
    }
//--------------------------------------------------------------
Stepper:: ~Stepper(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] acceptabilitylist;
}
//--------------------------------------------------------------
Stepper& Stepper::operator++() 
{
    dt_next = min(dt_next,1.01*_dt);
    dt_next = min(dt_next,Input::List().dt);
    failed_steps = 0;
    current_time += _dt;
    _dt = dt_next; 
    _success = false;
    return *this;
}

//--------------------------------------------------------------
//  Collect all of the terms
void Stepper::update_dt(const State1D& Y_old, const State1D& Ystar, State1D& Y_new){
//--------------------------------------------------------------

    // err_val =  check_temperature(Ystar,Y);  

    /// Ex is checked for convergence
    err_val =  check_flds(Ystar,Y_new,acceptability);

    /// Error is determined
    acceptability = err_val/(atol + rtol*acceptability);
    if (acceptability > 1) overall_check = false;
    else overall_check = true;

    // std::cout << "\n acc = " << acceptability;
    // std::cout << "\n dt = " << _dt;

    /// Error is shared so that global timestep can be determined on rank 0
    MPI_Gather(&acceptability, 1, MPI_DOUBLE, acceptabilitylist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (world_rank == 0)
    {   

        /// Determine whether to proceed
        for (size_t iprocess(0); iprocess < world_size; ++iprocess)
        {
            overall_check = ((!(acceptabilitylist[iprocess] > 1)) && overall_check);
            acceptability = max(acceptabilitylist[iprocess],acceptability);
        }

        /// Update timestep and if failed, add an iteration
        if (overall_check > 0)
        {
            dt_next = 0.9*_dt/pow(acceptability,0.2);   
        }
        else
        {
            dt_next = 0.9*_dt/pow(acceptability,0.25);
            
            ++failed_steps;
            if (failed_steps > max_failures) 
            {
                fprintf(stderr, "Time Stepper failed to converge within %d steps \n", max_failures);
                MPI_Finalize();
                exit(1);
            }
        }
    }

    /// Share success and new timestep
    MPI_Bcast(&overall_check, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_next, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /// If failed, restore old state and updated time step.
    /// Success time-step is updated at the end of outer loop.
    if (!overall_check)
    {
        _dt = dt_next;
        Y_new = Y_old;
    }
    else
    {
        _success = true;
    }

    // return (overall_check > 0);
}
//--------------------------------------------------------------
//  Collect all of the terms
void Stepper::update_dt(const State2D& Y_old, const State2D& Ystar, State2D& Y_new){
//--------------------------------------------------------------

    // err_val =  check_temperature(Ystar,Y);  

    /// Ex is checked for convergence
    err_val =  check_flds(Ystar,Y_new,acceptability);

    /// Error is determined
    acceptability = err_val/(atol + rtol*acceptability);
    if (acceptability > 1) overall_check = false;
    else overall_check = true;


    /// Error is shared so that global timestep can be determined on rank 0
    MPI_Gather(&acceptability, 1, MPI_DOUBLE, acceptabilitylist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (world_rank == 0)
    {   

        /// Determine whether to proceed
        for (size_t iprocess(0); iprocess < world_size; ++iprocess)
        {
            overall_check = ((!(acceptabilitylist[iprocess] > 1)) && overall_check);
            acceptability = max(acceptabilitylist[iprocess],acceptability);
        }

        /// Update timestep and if failed, add an iteration
        if (overall_check > 0)
        {
            dt_next = 0.95*_dt/pow(acceptability,0.2);
        }
        else
        {
            dt_next = 0.9*_dt/pow(acceptability,0.25);
            
            ++failed_steps;
            if (failed_steps > max_failures) 
            {
                fprintf(stderr, "Time Stepper failed to converge within %d steps \n", max_failures);
                MPI_Finalize();
                exit(1);
            }
        }
    }

    /// Share success and new timestep
    MPI_Bcast(&overall_check, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt_next, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /// If failed, restore old state and updated time step.
    /// Success time-step is updated at the end of outer loop.
    if (!overall_check)
    {
        _dt = dt_next;
        Y_new = Y_old;
    }
    else
    {
        _success = true;
    }
}

//--------------------------------------------------------------
//  Collect all of the terms
double Stepper::check_temperature(const State1D& Ystar, const State1D& Y){
//--------------------------------------------------------------
    // // Add up temperature over whole grid
    // for (size_t s(0); s < Ystar.Species(); ++s)
    // {
    //     tmp_err = Ystar.DF(s).getpressure() - Y.DF(s).getpressure();
    //     tmp_err *= tmp_err;
    //     return tmp_err.sum();
    // }
    
}



//--------------------------------------------------------------
//  Collect all of the terms
double Stepper::check_flds(const State1D& Ystar, const State1D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    valarray<double> tmp_err(Y.FLD(0).numx());

    for (size_t i(0); i < 1; ++i)
    {
        for (size_t ix(Nbc); ix < tmp_err.size() - Nbc; ++ix)
        {
            tmp_err[ix] = Ystar.FLD(i)(ix).real()-Y.FLD(i)(ix).real();
            maxval = max(maxval,abs(Y.FLD(i)(ix).real()));
            maxval = max(maxval,abs(Ystar.FLD(i)(ix).real()));
        }

        tmp_err = abs(tmp_err);

        return tmp_err.sum();
    }
}
//--------------------------------------------------------------
//  Collect all of the terms
double Stepper::check_flds(const State2D& Ystar, const State2D& Y, double& maxval){
//--------------------------------------------------------------
    maxval = 0.;
    valarray<double> tmp_err(Y.FLD(0).numx()*Y.FLD(0).numy());
    size_t ixy;
    for (size_t i(0); i < 1; ++i)
    {
        ixy = 0;
        for (size_t ix(Nbc); ix < Y.FLD(i).numx() - Nbc; ++ix)
        {
            for (size_t iy(Nbc); iy < Y.FLD(i).numy() - Nbc; ++iy)
            {
                tmp_err[ixy] = Ystar.FLD(i)(ix,iy).real()-Y.FLD(i)(ix,iy).real();
                maxval = max(maxval,abs(Y.FLD(i)(ix,iy).real()));
                maxval = max(maxval,abs(Ystar.FLD(i)(ix,iy).real()));
                ++ixy;
            }
        }

        tmp_err = abs(tmp_err);

        return tmp_err.sum();
    }
}


//--------------------------------------------------------------
//  Collect all of the terms
double Stepper::check_js(const State1D& Ystar, const State1D& Y){
//--------------------------------------------------------------
    
    // for (size_t s(0); s < Ystar.Species(); ++s)
    // {
    //     for (size_t i(0); i < 3; ++i)
    //     {
    //         tmp_err = Ystar.DF(s).getcurrent(i) - Y.DF(s).getcurrent(i);
    //         tmp_err *= tmp_err;

    //         // if (tmp_err.sum() < total_tol) return true;        
    //         return tmp_err.sum();
    //     }
    // }

}
