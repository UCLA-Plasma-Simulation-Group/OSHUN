/*! \brief Main Loop
 * \author PICKSC
 * \date   September 1, 2016
 * \file   main.cpp
 *
 * This file includes
 * - the main loop.
 * - the definition for the clock class
 * - the definition of periodic boundaries
 *
 */

// Standard libraries
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdarg.h>


#include <map>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <ctime>


// My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

#include "external/exprtk.hpp"
#include "external/spline.h"
#include "external/highfive/H5DataSet.hpp"

// Misc Declerations
#include "input.h"
#include "state.h"
#include "clock.h"
#include "formulary.h"
#include "setup.h"
#include "fluid.h"
#include "vlasov.h"
#include "collisions.h"
#include "functors.h"
#include "parallel.h"
#include "implicitE.h"
#include "particletracker.h"
#include "export.h"



//--------------------------------------------------------------
//**************************************************************

//----------------------------------------------------------------------------------------------------------------------------
    double startmessages(){

        Formulary formulas;
        // double Tref = Input::List().pth_ref*Input::List().pth_ref;
        double ND = 1.72e9*sqrt(pow(formulas.T0,3.0)/formulas.n);
        double nuei_wp = sqrt(2.0/Formulary::pi)/9.0/ND*
        formulas.LOGei(1.0,pow(Input::List().pth_ref,2.0),formulas.Zeta)*formulas.Zeta;

        double nu_ei = nuei_wp * formulas.wp;
        double dt_out(Input::List().t_stop / Input::List().n_outsteps);
        double CLF(Input::List().dt);
        size_t n_outsteps(Input::List().n_outsteps);
        double numh  = size_t(static_cast<int>(dt_out/CLF))+1;
        double h     = dt_out/static_cast<double>(numh);


        std::cout << "\n\n";
        std::cout<< "------- ---------------- */*/*/*/*/* ---------------- -------\n";
        std::cout<< "------- ,-----.  ,---.  ,--.  ,--.,--. ,--.,--.  ,--. -------\n";
        std::cout<< "-------'  .-.  ''   .-' |  '--'  ||  | |  ||  ,'.|  | -------\n";
        std::cout<< "-------|  | |  |`.  `-. |  .--.  ||  | |  ||  |' '  | -------\n";
        std::cout<< "-------'  '-'  '.-'    ||  |  |  |'  '-'  '|  | `   | -------\n";
        std::cout<< "------- `-----' `-----' `--'  `--' `-----' `--'  `--' -------\n";
        std::cout<< "------- ---------------- */*/*/*/*/* ---------------- -------\n";



        std::cout << "\n\n";
        std::cout << "--------------- OSHUN Beta - 1/2D+3P --------------- \n";
        std::cout << "    Particle-in-Cell and Kinetic Simulation Center   \n";
        std::cout << "------------------- UCLA - 2017 -------------------- \n";
        
        std::cout << "----------------- Reference Units ------------------ \n";
        std::cout << "normalizing density                   = " << formulas.n << " / cc \n";
        std::cout << "plasma frequency                      = " << formulas.wp << " Hz \n";
        std::cout << "background ions - Z                   = " << formulas.Zeta << " \n";
        std::cout << "normalizing magnetic field            = " << formulas.B0 << " T \n";
        std::cout << "Reference Temperature                 = " << formulas.T0 << " eV \n";
        std::cout << "Corresponding thermal velocity        = " << Input::List().pth_ref*Formulary::cL  << " m/s \n";

        std::cout << "\n";
        std::cout << "Corresponding e-i log Lambda          = " << formulas.LOGei(1.0,pow(Input::List().pth_ref,2.0),formulas.Zeta) << " \n";
        std::cout << "Corresponding e-i collision frequency = " << nu_ei << " Hz \n";
        std::cout << "Corresponding tau_e-i / tau_p         = " << 1.0/nuei_wp << " \n";
        std::cout << "Corresponding e-i mean free path      = " << Input::List().pth_ref*Formulary::cL/nu_ei *1e6<< " microns \n";

        std::cout << "\n";
        std::cout << "Corresponding e-e log Lambda          = " << formulas.LOGee(1.0,pow(Input::List().pth_ref,2.0)) << " \n";
        std::cout << "Corresponding e-e collision frequency = " << 1.0/formulas.Tau_e(1.0,pow(Input::List().pth_ref,2.0)) << " Hz \n";
            // std::cout << "Corresponding thermal mean free path =" << Input::List().pth_ref/cL / nu_ei  << "\n";
            // std::cout << "Corresponding thermal velocity  c =" << pow(Input::List().pth_ref*cL,2.0)*0.5*me << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "--------------  Simulation parameters -------------- \n";
        std::cout << "skin depth (normalizing distance)     = " << formulas.skindepth*1e6 << " microns \n";
        std::cout << "plasma period (normalizing time)      = " << 1.0/formulas.wp*1e15 << " fs \n";
        std::cout << "Time step     (normalizing time)      = " << h << " plasma periods \n";
        std::cout << "\n";
        std::cout << "\n";



        return 1.0/formulas.wp*1e15;
    }

//**************************************************************
int main(int argc, char** argv) {

    MPI_Init(&argc,&argv);
    time_t tstart, tend;
    tstart = omp_get_wtime();

    omp_set_num_threads(Input::List().ompthreads);
    

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    /////////     MAIN LOOP                                                      ///////
    /////////       CONTAINS AN IF STATEMENT FOR EXPLICIT OR IMPLICIT E SOLVER  ///////
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    if (Input::List().dim == 1)   /// 1-D
    {
        ///  Initiate the Parallel Environment and decompose the Computational Domain
        std::cout << "\nInitializing parallel environment ...";
        Parallel_Environment_1D PE;

        std::cout << "     done \n\n";
        
        if (PE.RANK() == 0) {
            Export_Files::Folders();
        }
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    ///  Set up the grid
    ///    Moves all of the relevant data from the input deck into a single container
    ///
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
        if (!PE.RANK()) std::cout << "\nInitializing grid ...";
        Grid_Info grid(Input::List().ls, Input::List().ms,
                        Input::List().xminLocal, Input::List().xmaxLocal, Input::List().NxLocal,
                        Input::List().xminGlobal, Input::List().xmaxGlobal, Input::List().NxGlobal,
                        // Input::List().pmax, Input::List().numps,
                        Input::List().dp,
                        Input::List().dpx,Input::List().dpy,Input::List().dpz);
                        // Input::List().Npx, Input::List().Npy, Input::List().Npz);
        if (!PE.RANK()) std::cout << "     done \n";
        
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    ///  CLOCK
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
        int tout_start;
        if (Input::List().isthisarestart) 
        {
            tout_start = Input::List().restart_time;
        }
        else
        {
            tout_start = 0;
        }
        size_t t_out(tout_start+1);
        double start_time(0.);
    
        double dt_out(Input::List().t_stop / (Input::List().n_outsteps+1));
        double dt_dist_out(Input::List().t_stop / (Input::List().n_distoutsteps+1));
        double dt_big_dist_out(Input::List().t_stop / (Input::List().n_bigdistoutsteps+1));
        double dt_restart(Input::List().t_stop / (Input::List().n_restarts+1));
        
        double next_out(dt_out);
        double next_dist_out(dt_dist_out);
        double next_big_dist_out(dt_big_dist_out);

        double next_restart(dt_restart);

    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    ///  INITIALIZATION
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    
        if (!PE.RANK()) std::cout << "Initializing restart environment ...";
        Export_Files::Restart_Facility Re(PE.RANK());
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Initializing state variable ...";
        State1D Y( grid.axis.Nx(0), Input::List().ls, Input::List().ms, 
            Input::List().dp, 
            Input::List().qs, Input::List().mass, 
            Input::List().hydromass, Input::List().hydrocharge, 
            Input::List().numparticles, Input::List().particlemass, Input::List().particlecharge);
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Initializing plasma profile ...";
        Setup_Y::initialize(Y, grid);
        if (!PE.RANK()) std::cout << "     done \n";
        

        if (!PE.RANK()) std::cout << "Initializing collision module ...";
        collisions_1D collide(Y);
        if (!PE.RANK()) std::cout << "     done \n";    
        
        if (!PE.RANK()) std::cout << "Initializing hydro module ...";
        Hydro_Functor         HydroFunc(grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
        if (!PE.RANK()) std::cout << "     done \n";
        
        if (!PE.RANK()) std::cout << "Initializing particle tracker ...";
        Particle_Pusher       Particle_Push(Input::List().par_xpos, Input::List().par_px, Input::List().par_py,  Input::List().par_pz,
            Input::List().xminLocalnobnd[0], Input::List().xmaxLocalnobnd[0], Input::List().NxLocalnobnd[0], Y.particles());
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (Input::List().isthisarestart){
            if (!PE.RANK()) std::cout << "Reading restart files ...";
            Re.Read(PE.RANK(),tout_start,Y,start_time);
            if (!PE.RANK()) std::cout << "     done \n";
        }
        
        if (!PE.RANK()) std::cout << "Initializing output module ...";
        Output_Data::Output_Preprocessor  output( grid, Input::List().oTags);
        if (!PE.RANK()) std::cout << "     done \n";
        
        if (!PE.RANK()) std::cout << "Output #0 ...";    
        output( Y, grid, tout_start, start_time, Input::List().dt, PE );
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Distribution function output #0 ...";    
        output.distdump(Y, grid, tout_start, start_time, Input::List().dt, PE );
        output.bigdistdump(Y, grid, tout_start, start_time, Input::List().dt, PE );
        if (!PE.RANK()) std::cout << "     done \n";
    
        double plasmaperiod;
        if (!(PE.RANK())){
            plasmaperiod = startmessages();
        }
    
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
    //  ITERATION LOOP
    // --------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------------------------------------------------------------------------------------------------------
        if (Input::List().implicit_E) 
        {
            if (!PE.RANK())
            { 
                std::cout << "Starting Semi-Implicit, 1D OSHUN\n";
            }
            
            Algorithms::RK2<State1D> RK(Y);
            
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // IMPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_implicitE_p1 impE_p1_Functor(Input::List().ls, Input::List().ms, 
                                                            Input::List().dp,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            VlasovFunctor1D_implicitE_p2 impE_p2_Functor(Input::List().ls, Input::List().ms,                                                         
                                                            Input::List().dp,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            
            // --------------------------------------------------------------------------------------------------------------------------------
            using Electric_Field_Methods::Efield_Method;
            Electric_Field_Methods::Implicit_E_Field eim(grid.axis);

            if (!Input::List().collisions) {
                if (!PE.RANK())
                    std::cout << "\n Need collisions for implicit E field solver. \n Exiting. \n";
                exit(0);
            }
            
            
            // Algorithms::RKCK54<State1D> RK54(Y);
            // Algorithms::RK4<State1D> RK(Y);
            // State1D Y_star(Y), Y_old(Y);
            // bool success(false);

            Stepper step(start_time,Input::List().dt,Input::List().abs_tol,Input::List().rel_tol,Input::List().max_fails);

            for(step; step.time() < Input::List().t_stop; ++step)
            {
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, step.time());

                // Y_old = Y;

                // while(!success)
                // {
                //     RK54(Y_star,Y,step.dt(),&impE_p1_Functor);
                //     success = step.update_dt(Y_old,Y_star, Y);
                // }
                
                Y = RK(Y, step.dt(), &impE_p1_Functor);                                                             /// Vlasov - Updates the distribution function: Spatial Advection and B Field "action".
                PE.Neighbor_ImplicitE_Communications(Y);                                                            /// Boundaries
                eim.advance(&RK, Y, collide,&impE_p2_Functor, step.dt());                                           /// Finds new electric field
                Y = RK(Y, step.dt(), &impE_p2_Functor);         

                if (Input::List().collisions)
                    collide.advance(Y,step.time(),step.dt());                                               ///  Fokker-Planck   //

                // if (Input::List().hydromotion)
                //     Y = RK(Y, step.dt(), &HydroFunc);                                      /// Hydro Motion

                if (Input::List().trav_wave) Setup_Y::applytravelingwave(grid, Y, step.time(), step.dt());

                PE.Neighbor_Communications(Y);                                         ///  Boundaries      //
                // --------------------------------------------------------------------------------------------------------------------------------
                // --------------------------------------------------------------------------------------------------------------------------------
                /// Output
                // --------------------------------------------------------------------------------------------------------------------------------
                if (step.time() > next_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
                    output.distdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_dist_out += dt_dist_out;
                }
                
                if (step.time() > next_big_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
                    output.bigdistdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_big_dist_out += dt_big_dist_out;
                }

                if (step.time() > next_restart)
                {
                    if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";
                    Re.Write(PE.RANK(), t_out, Y, step.time());
                    next_restart += dt_restart;
                }

                if (step.time() > next_out)
                {
                    if (!(PE.RANK()))
                    {
                        cout << "\n dt = " << step.dt();
                        cout << " , Output #" << t_out;
                        
                    }

                    output(Y, grid, t_out, step.time(), step.dt(), PE);
                    Y.checknan();

                    next_out += dt_out;
                    ++t_out;
                }

                // success = false;                    
            }
        } 
        else 
        {
            if (!PE.RANK())
            { 
                std::cout << "Starting Fully-Explicit, 1D OSHUN\n";
            }
            
            // Algorithms::RK4<State1D> RK(Y);

            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // EXPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_explicitE rkF(Input::List().ls, Input::List().ms, 
                                                            Input::List().dp,
                                          grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            // --------------------------------------------------------------------------------------------------------------------------------

            // Algorithms::RKCK54<State1D> RK54(Y);
            Algorithms::RKBS54<State1D> RK54(Y);
            // Algorithms::RKT54<State1D> RK54(Y);
            // Algorithms::RK4<State1D> RK(Y);
            State1D Y_star(Y), Y_old(Y);

            Stepper step(start_time,Input::List().dt,Input::List().abs_tol,Input::List().rel_tol,Input::List().max_fails);

            for(step; step.time() < Input::List().t_stop; ++step)
            {
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, step.time());

                Y_old = Y;

                // RK(Y,step.dt(),&rkF);

                while(!step.success())
                {
                    RK54(Y_star,Y,step.dt(),&rkF);
                    step.update_dt(Y_old,Y_star, Y);
                }

                if (Input::List().collisions)
                    collide.advance(Y,step.time(),step.dt());                                           ///  Fokker-Planck   //

                // if (Input::List().hydromotion)
                    // Y = RK(Y, step.dt(), &HydroFunc);                                                   /// Hydro Motion

                if (Input::List().trav_wave) 
                {                    
                    Setup_Y::applytravelingwave(grid, Y, step.time(), step.dt());
                }

                PE.Neighbor_Communications(Y);                                         ///  Boundaries      //
                // --------------------------------------------------------------------------------------------------------------------------------
                // --------------------------------------------------------------------------------------------------------------------------------
                /// Output
                // --------------------------------------------------------------------------------------------------------------------------------
                if (step.time() > next_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
                    output.distdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_dist_out += dt_dist_out;
                }
                
                if (step.time() > next_big_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
                    output.bigdistdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_big_dist_out += dt_big_dist_out;
                }

                if (step.time() > next_restart)
                {
                    if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";
                    Re.Write(PE.RANK(), t_out, Y, step.time());
                    next_restart += dt_restart;
                }

                if (step.time() > next_out)
                {
                    if (!(PE.RANK()))
                    {
                        cout << "\n dt = " << step.dt();
                        cout << " , Output #" << t_out;
                        
                    }

                    output(Y, grid, t_out, step.time(), step.dt(), PE);
                    Y.checknan();

                    next_out += dt_out;
                    ++t_out;
                }
            }
        }
        tend = omp_get_wtime();
        if (!(PE.RANK())){
            cout << "Simulation took "<< difftime(tend, tstart) <<" second(s)."<< endl;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    /////////     MAIN LOOP FOR 2D CODE                                         ///////
    /////////       CONTAINS AN IF STATEMENT FOR EXPLICIT OR IMPLICIT E SOLVER  ///////
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////    
    else
    {
        ///  Initiate the Parallel Environment and decompose the Computational Domain
        std::cout << "\nInitializing parallel environment ...";
        Parallel_Environment_2D PE;
        std::cout << "     done \n\n";
        
        if (PE.RANK() == 0) {
            Export_Files::Folders();
        }
    
        ///  Set up the grid
        ///    Moves all of the relevant data from the input deck into a single container
        if (!PE.RANK()) std::cout << "\nInitializing grid ...";
        Grid_Info grid(Input::List().ls, Input::List().ms,
                        Input::List().xminLocal, Input::List().xmaxLocal, Input::List().NxLocal,
                        Input::List().xminGlobal, Input::List().xmaxGlobal, Input::List().NxGlobal,
                        Input::List().dp,
                        Input::List().dpx,Input::List().dpy,Input::List().dpz);
                        
        if (!PE.RANK()) std::cout << "     done \n";
        
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        ///  CLOCK
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        int tout_start;
        if (Input::List().isthisarestart) tout_start = Input::List().restart_time;
        else  tout_start = 0;
        size_t t_out(tout_start+1);
        double start_time(0.);
    
        double dt_out(Input::List().t_stop / (Input::List().n_outsteps+1));
        double dt_dist_out(Input::List().t_stop / (Input::List().n_distoutsteps+1));
        double dt_big_dist_out(Input::List().t_stop / (Input::List().n_bigdistoutsteps+1));
        double dt_restart(Input::List().t_stop / (Input::List().n_restarts+1));
        
        double next_out(dt_out);
        double next_dist_out(dt_dist_out);
        double next_big_dist_out(dt_big_dist_out);

        double next_restart(dt_restart);

        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        ///  INITIALIZATION
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        if (!PE.RANK()) std::cout << "Initializing restart environment ...";
        Export_Files::Restart_Facility Re(PE.RANK());
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Initializing state variable ...";
        State2D Y( grid.axis.Nx(0), grid.axis.Nx(1), Input::List().ls, Input::List().ms, 
            Input::List().dp, 
            Input::List().qs, Input::List().mass, 
            Input::List().hydromass, Input::List().hydrocharge);
        if (!PE.RANK()) std::cout << "     done \n";

    
        if (!PE.RANK()) std::cout << "Initializing plasma profile ...";
        Setup_Y::initialize(Y, grid);
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Initializing collision module ...";
        collisions_2D collide(Y);
        if (!PE.RANK()) std::cout << "     done \n";    
    
        // if (!PE.RANK()) std::cout << "Initializing hydro module ...";
        // Hydro_Functor         HydroFunc(grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
        // if (!PE.RANK()) std::cout << "     done \n";
        
        //         if (!PE.RANK()) std::cout << "Initializing particle tracker ...";
        //         Particle_Pusher       Particle_Push(Input::List().par_xpos, Input::List().par_px, Input::List().par_py,  Input::List().par_pz,
        //             Input::List().xminLocalnobnd[0], Input::List().xmaxLocalnobnd[0], Input::List().NxLocalnobnd[0], Y.particles());
        //         if (!PE.RANK()) std::cout << "     done \n";
        
        if (Input::List().isthisarestart){
            if (!PE.RANK()) std::cout << "Reading restart files ...";
            Re.Read(PE.RANK(),tout_start,Y,start_time);
            if (!PE.RANK()) std::cout << "     done \n";
        }
        
        if (!PE.RANK()) std::cout << "Initializing output module ...";
        Output_Data::Output_Preprocessor  output( grid, Input::List().oTags);
        if (!PE.RANK()) std::cout << "     done \n";
                     
        if (!PE.RANK()) std::cout << "Output #0 ...";    
        output( Y, grid, tout_start, start_time, Input::List().dt, PE );
        if (!PE.RANK()) std::cout << "     done \n";
    
        if (!PE.RANK()) std::cout << "Distribution function output #0 ...";    
        output.distdump(Y, grid, tout_start, start_time, Input::List().dt, PE );
        output.bigdistdump(Y, grid, tout_start, start_time, Input::List().dt, PE );
        if (!PE.RANK()) std::cout << "     done \n";
    
        double plasmaperiod;
        if (!(PE.RANK())){
            plasmaperiod = startmessages();
        }
    
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        //  ITERATION LOOP
        // --------------------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------------------
        if (Input::List().implicit_E) 
        {
            if (!PE.RANK())
            { 
                std::cout << "Starting Semi-Implicit, 2D OSHUN\n";
            }
            Algorithms::RK2<State2D> RK(Y);
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // IMPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor2D_implicitE_p1 impE_p1_Functor(Input::List().ls, Input::List().ms, 
                                                        // Input::List().pmax, Input::List().numps,
                                                            Input::List().dp,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0),
                                                         grid.axis.xmin(1), grid.axis.xmax(1), grid.axis.Nx(1));
            VlasovFunctor2D_implicitE_p2 impE_p2_Functor(Input::List().ls, Input::List().ms, 
                                                            // Input::List().pmax, Input::List().numps,
                                                            Input::List().dp,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0),
                                                         grid.axis.xmin(1), grid.axis.xmax(1), grid.axis.Nx(1));
            // --------------------------------------------------------------------------------------------------------------------------------
            using Electric_Field_Methods::Efield_Method;
            Electric_Field_Methods::Implicit_E_Field eim(grid.axis);

            if (!Input::List().collisions) {
                if (!PE.RANK())
                    std::cout << "\n Need collisions for implicit E field solver. \n Exiting. \n";
                exit(0);
            }

            // Algorithms::RKCK54<State1D> RK54(Y);
            // Algorithms::RK4<State1D> RK(Y);
            // State1D Y_star(Y), Y_old(Y);
            // bool success(false);

            Stepper step(start_time,Input::List().dt,Input::List().abs_tol,Input::List().rel_tol,Input::List().max_fails);

            for(step; step.time() < Input::List().t_stop; ++step)
            {
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, step.time());

                // Y_old = Y;

                // while(!success)
                // {
                //     RK54(Y_star,Y,step.dt(),&impE_p1_Functor);
                //     success = step.update_dt(Y_old,Y_star, Y);
                // }
                
                Y = RK(Y, step.dt(), &impE_p1_Functor);                                                             /// Vlasov - Updates the distribution function: Spatial Advection and B Field "action".
                PE.Neighbor_ImplicitE_Communications(Y);                                                            /// Boundaries
                eim.advance(&RK, Y, collide,&impE_p2_Functor, step.dt());                                           /// Finds new electric field
                Y = RK(Y, step.dt(), &impE_p2_Functor);         

                if (Input::List().collisions)
                    collide.advance(Y,step.time(),step.dt());                                               ///  Fokker-Planck   //

                // if (Input::List().hydromotion)
                    // Y = RK(Y, step.dt(), &HydroFunc);                                      /// Hydro Motion

                if (Input::List().trav_wave) Setup_Y::applytravelingwave(grid, Y, step.time(), step.dt());

                PE.Neighbor_Communications(Y);                                         ///  Boundaries      //
                // --------------------------------------------------------------------------------------------------------------------------------
                // --------------------------------------------------------------------------------------------------------------------------------
                /// Output
                // --------------------------------------------------------------------------------------------------------------------------------
                if (step.time() > next_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
                    output.distdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_dist_out += dt_dist_out;
                }
                
                if (step.time() > next_big_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
                    output.bigdistdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_big_dist_out += dt_big_dist_out;
                }

                if (step.time() > next_restart)
                {
                    if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";
                    Re.Write(PE.RANK(), t_out, Y, step.time());
                    next_restart += dt_restart;
                }

                if (step.time() > next_out)
                {
                    if (!(PE.RANK()))
                    {
                        cout << "\n dt = " << step.dt();
                        cout << " , Output #" << t_out;
                        
                    }

                    output(Y, grid, t_out, step.time(), step.dt(), PE);
                    Y.checknan();

                    next_out += dt_out;
                    ++t_out;
                }

                // success = false;                    
            }
        } 
        else 
        {
            if (!PE.RANK())
            { 
                std::cout << "Starting Fully-Explicit, 2D OSHUN\n";
            }
            Algorithms::RK4<State2D> RK(Y);
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // EXPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            // std::cout << "\n 10 \n";
            VlasovFunctor2D_explicitE rkF(Input::List().ls, Input::List().ms, 
                                                            Input::List().dp,
                                          grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0),
                                          grid.axis.xmin(1), grid.axis.xmax(1), grid.axis.Nx(1));

            Algorithms::RKCK54<State2D> RK54(Y);
            // Algorithms::RK4<State1D> RK(Y);
            State2D Y_star(Y), Y_old(Y);

            Stepper step(start_time,Input::List().dt,Input::List().abs_tol,Input::List().rel_tol,Input::List().max_fails);

            for(step; step.time() < Input::List().t_stop; ++step)
            {
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, step.time());

                Y_old = Y;

                while(!step.success())
                {
                    RK54(Y_star,Y,step.dt(),&rkF);
                    step.update_dt(Y_old,Y_star, Y);
                }

                if (Input::List().collisions)
                    collide.advance(Y,step.time(),step.dt());                                           ///  Fokker-Planck   //

                // if (Input::List().hydromotion)
                    // Y = RK(Y, step.dt(), &HydroFunc);                                                   /// Hydro Motion

                if (Input::List().trav_wave) Setup_Y::applytravelingwave(grid, Y, step.time(), step.dt());

                PE.Neighbor_Communications(Y);                                         ///  Boundaries      //
                // --------------------------------------------------------------------------------------------------------------------------------
                // --------------------------------------------------------------------------------------------------------------------------------
                /// Output
                // --------------------------------------------------------------------------------------------------------------------------------
                if (step.time() > next_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Dist Output #" << t_out << "\n";
                    output.distdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_dist_out += dt_dist_out;
                }
                
                if (step.time() > next_big_dist_out)
                {
                    if (!(PE.RANK())) cout << " \n Big Dist Output #" << t_out << "\n";
                    output.bigdistdump(Y, grid, t_out, step.time(), step.dt(), PE);
                    next_big_dist_out += dt_big_dist_out;
                }

                if (step.time() > next_restart)
                {
                    if (!(PE.RANK())) cout << " \n Restart Output #" << t_out << "\n";
                    Re.Write(PE.RANK(), t_out, Y, step.time());
                    next_restart += dt_restart;
                }

                if (step.time() > next_out)
                {
                    if (!(PE.RANK()))
                    {
                        cout << "\n dt = " << step.dt();
                        cout << " , Output #" << t_out;
                        
                    }

                    output(Y, grid, t_out, step.time(), step.dt(), PE);
                    Y.checknan();

                    next_out += dt_out;
                    ++t_out;
                }                
            }
        }
        tend = omp_get_wtime();
        if (!(PE.RANK())){
            cout << "Simulation took "<< difftime(tend, tstart) <<" second(s)."<< endl;
        }
    }

    MPI_Finalize();
	return 0;
}
//**************************************************************
//---------------------------------------------------------------------------
