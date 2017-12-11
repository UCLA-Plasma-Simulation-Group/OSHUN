/*!\brief  Initialization routines - Definitions
* \author  PICKSC
 * \date   2017
 * \file   setup.cpp
 *
 * In here are the initialization routines for the State object.
 *
 */

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <cfloat>

#include <math.h>
// #include <boost/math/special_functions/bessel.hpp>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "external/exprtk.hpp"

//  Declerations
#include "input.h"
#include "state.h"
#include "formulary.h"
#include "parser.h"
#include "setup.h"

//**************************************************************
//**************************************************************
//--------------------------------------------------------------
void Setup_Y::applyexternalfields(Grid_Info &grid, State1D& Y, double time)
{

    valarray<double> Ex_profile( grid.axis.Nx(0));
    valarray<double> Ey_profile( grid.axis.Nx(0));
    valarray<double> Ez_profile( grid.axis.Nx(0));
    valarray<double> Bx_profile( grid.axis.Nx(0));
    valarray<double> By_profile( grid.axis.Nx(0));
    valarray<double> Bz_profile( grid.axis.Nx(0));

    double ex_time_coeff(0.0);
    double ey_time_coeff(0.0);
    double ez_time_coeff(0.0);
    double bx_time_coeff(0.0);
    double by_time_coeff(0.0);
    double bz_time_coeff(0.0);

    Parser::parseprofile(grid.axis.x(0), Input::List().ex_profile_str, Ex_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().ey_profile_str, Ey_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().ez_profile_str, Ez_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().bx_profile_str, Bx_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().by_profile_str, By_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().bz_profile_str, Bz_profile);

    Parser::parseprofile(time, Input::List().ex_time_profile_str, ex_time_coeff);
    Parser::parseprofile(time, Input::List().ey_time_profile_str, ey_time_coeff);
    Parser::parseprofile(time, Input::List().ez_time_profile_str, ez_time_coeff);
    Parser::parseprofile(time, Input::List().bx_time_profile_str, bx_time_coeff);
    Parser::parseprofile(time, Input::List().by_time_profile_str, by_time_coeff);
    Parser::parseprofile(time, Input::List().bz_time_profile_str, bz_time_coeff);

    for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    {
        if (ex_time_coeff > 1e-20) Y.EMF().Ex()(ix) = Ex_profile[ix]*ex_time_coeff;
        if (ey_time_coeff > 1e-20) Y.EMF().Ey()(ix) = Ey_profile[ix]*ey_time_coeff;
        if (ez_time_coeff > 1e-20) Y.EMF().Ez()(ix) = Ez_profile[ix]*ez_time_coeff;
        if (bx_time_coeff > 1e-20) Y.EMF().Bx()(ix) = Bx_profile[ix]*bx_time_coeff;
        if (by_time_coeff > 1e-20) Y.EMF().By()(ix) = By_profile[ix]*by_time_coeff;
        if (bz_time_coeff > 1e-20) Y.EMF().Bz()(ix) = Bz_profile[ix]*bz_time_coeff;

    }
}
//---------------------------------------------------------------------------
//--------------------------------------------------------------
void Setup_Y::applyexternalfields(Grid_Info &grid, State2D& Y, double time)
{

    Array2D<double> Ex_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> Ey_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> Ez_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> Bx_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> By_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> Bz_profile( grid.axis.Nx(0), grid.axis.Nx(1));

    double ex_time_coeff(0.0);
    double ey_time_coeff(0.0);
    double ez_time_coeff(0.0);
    double bx_time_coeff(0.0);
    double by_time_coeff(0.0);
    double bz_time_coeff(0.0);

    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().ex_profile_str, Ex_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().ey_profile_str, Ey_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().ez_profile_str, Ez_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().bx_profile_str, Bx_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().by_profile_str, By_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().bz_profile_str, Bz_profile);

    Parser::parseprofile(time, Input::List().ex_time_profile_str, ex_time_coeff);
    Parser::parseprofile(time, Input::List().ey_time_profile_str, ey_time_coeff);
    Parser::parseprofile(time, Input::List().ez_time_profile_str, ez_time_coeff);
    Parser::parseprofile(time, Input::List().bx_time_profile_str, bx_time_coeff);
    Parser::parseprofile(time, Input::List().by_time_profile_str, by_time_coeff);
    Parser::parseprofile(time, Input::List().bz_time_profile_str, bz_time_coeff);

    for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    {
        for (size_t iy(0);iy<Y.SH(0,0,0).numy();++iy)
        {

            if (ex_time_coeff > 1e-20) Y.EMF().Ex()(ix,iy) = Ex_profile(ix,iy)*ex_time_coeff;
            if (ey_time_coeff > 1e-20) Y.EMF().Ey()(ix,iy) = Ey_profile(ix,iy)*ey_time_coeff;
            if (ez_time_coeff > 1e-20) Y.EMF().Ez()(ix,iy) = Ez_profile(ix,iy)*ez_time_coeff;
            if (bx_time_coeff > 1e-20) Y.EMF().Bx()(ix,iy) = Bx_profile(ix,iy)*bx_time_coeff;
            if (by_time_coeff > 1e-20) Y.EMF().By()(ix,iy) = By_profile(ix,iy)*by_time_coeff;
            if (bz_time_coeff > 1e-20) Y.EMF().Bz()(ix,iy) = Bz_profile(ix,iy)*bz_time_coeff;
        }

    }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void Setup_Y::applytravelingwave(Grid_Info &grid, State1D& Y, double time, double stepsize)
{

    for (size_t n(0); n < Input::List().num_waves; ++n)
    {

        valarray<double> Ex_profile( grid.axis.Nx(0));
        valarray<double> Ey_profile( grid.axis.Nx(0));
        valarray<double> Ez_profile( grid.axis.Nx(0));
        valarray<double> Bx_profile( grid.axis.Nx(0));
        valarray<double> By_profile( grid.axis.Nx(0));
        valarray<double> Bz_profile( grid.axis.Nx(0));

        /// Parser::parse the strings
        Parser::parseprofile(grid.axis.x(0), time, Input::List().ex_wave_profile_str[n], Ex_profile);
        Parser::parseprofile(grid.axis.x(0), time, Input::List().ey_wave_profile_str[n], Ey_profile);
        Parser::parseprofile(grid.axis.x(0), time, Input::List().ez_wave_profile_str[n], Ez_profile);
        Parser::parseprofile(grid.axis.x(0), time, Input::List().bx_wave_profile_str[n], Bx_profile);
        Parser::parseprofile(grid.axis.x(0), time, Input::List().by_wave_profile_str[n], By_profile);
        Parser::parseprofile(grid.axis.x(0), time, Input::List().bz_wave_profile_str[n], Bz_profile);


        double time_coeff(-1.0);

        double pulse_start(Input::List().trav_wave_center[n] -
                        Input::List().trav_wave_flat[n]*0.5 -
                        Input::List().trav_wave_rise[n] );

        double pulse_end(Input::List().trav_wave_center[n] +
                        Input::List().trav_wave_flat[n]*0.5 +
                        Input::List().trav_wave_fall[n] );

        double normalized_time(0.0);
        
        if (time >= pulse_start && time <= pulse_end)
        {
        
            if (time < Input::List().trav_wave_center[n] - 0.5*Input::List().trav_wave_flat[n])
            {
                normalized_time = (time - pulse_start)/ Input::List().trav_wave_rise[n];

                time_coeff  = 6.0 * pow(normalized_time,5.0);
                time_coeff -= 15.0 * pow(normalized_time,4.0);
                time_coeff += 10.0 * pow(normalized_time,3.0);
            }
            else if (time >= Input::List().trav_wave_center[n] - 0.5*Input::List().trav_wave_flat[n] &&
                time <= Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n])
            {
                time_coeff = 1.0;
            }
            else if (time > Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n])
            {
                normalized_time = (time - (Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n]))
                / Input::List().trav_wave_fall[n];

                time_coeff  = 15.0 * pow(normalized_time,4.0);
                time_coeff -= 6.0 * pow(normalized_time,5.0);
                time_coeff -= 10.0 * pow(normalized_time,3.0);
                time_coeff += 1.0;
            }

            for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
            {
                Y.EMF().Ex()(ix) += Ex_profile[ix]*time_coeff*stepsize;
                Y.EMF().Ey()(ix) += Ey_profile[ix]*time_coeff*stepsize;
                Y.EMF().Ez()(ix) += Ez_profile[ix]*time_coeff*stepsize;
                Y.EMF().Bx()(ix) += Bx_profile[ix]*time_coeff*stepsize;
                Y.EMF().By()(ix) += By_profile[ix]*time_coeff*stepsize;
                Y.EMF().Bz()(ix) += Bz_profile[ix]*time_coeff*stepsize;
            }
        }
    }
}
//---------------------------------------------------------------------------
void Setup_Y::applytravelingwave(Grid_Info &grid, State2D& Y, double time, double stepsize)
{

    for (size_t n(0); n < Input::List().num_waves; ++n)
    {

        Array2D<double> Ex_profile( grid.axis.Nx(0), grid.axis.Nx(1));
        Array2D<double> Ey_profile( grid.axis.Nx(0), grid.axis.Nx(1));
        Array2D<double> Ez_profile( grid.axis.Nx(0), grid.axis.Nx(1));
        Array2D<double> Bx_profile( grid.axis.Nx(0), grid.axis.Nx(1));
        Array2D<double> By_profile( grid.axis.Nx(0), grid.axis.Nx(1));
        Array2D<double> Bz_profile( grid.axis.Nx(0), grid.axis.Nx(1));

        /// Parser::parse the strings
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().ex_wave_profile_str[n], Ex_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().ey_wave_profile_str[n], Ey_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().ez_wave_profile_str[n], Ez_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().bx_wave_profile_str[n], Bx_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().by_wave_profile_str[n], By_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), time, Input::List().bz_wave_profile_str[n], Bz_profile);


        double time_coeff(-1.0);

        double pulse_start(Input::List().trav_wave_center[n] -
                        Input::List().trav_wave_flat[n]*0.5 -
                        Input::List().trav_wave_rise[n] );

        double pulse_end(Input::List().trav_wave_center[n] +
                        Input::List().trav_wave_flat[n]*0.5 +
                        Input::List().trav_wave_fall[n] );

        double normalized_time(0.0);
        
        if (time >= pulse_start && time <= pulse_end)
        {
        
            if (time < Input::List().trav_wave_center[n] - 0.5*Input::List().trav_wave_flat[n])
            {
                normalized_time = (time - pulse_start)/ Input::List().trav_wave_rise[n];

                time_coeff  = 6.0 * pow(normalized_time,5.0);
                time_coeff -= 15.0 * pow(normalized_time,4.0);
                time_coeff += 10.0 * pow(normalized_time,3.0);
            }
            else if (time >= Input::List().trav_wave_center[n] - 0.5*Input::List().trav_wave_flat[n] &&
                time <= Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n])
            {
                time_coeff = 1.0;
            }
            else if (time > Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n])
            {
                normalized_time = (time - (Input::List().trav_wave_center[n] + 0.5*Input::List().trav_wave_flat[n]))
                / Input::List().trav_wave_fall[n];

                time_coeff  = 15.0 * pow(normalized_time,4.0);
                time_coeff -= 6.0 * pow(normalized_time,5.0);
                time_coeff -= 10.0 * pow(normalized_time,3.0);
                time_coeff += 1.0;
            }

            for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
            {
                for (size_t iy(0);iy<Y.SH(0,0,0).numy();++iy)
                {
                    Y.EMF().Ex()(ix,iy) += Ex_profile(ix,iy)*time_coeff*stepsize;
                    Y.EMF().Ey()(ix,iy) += Ey_profile(ix,iy)*time_coeff*stepsize;
                    Y.EMF().Ez()(ix,iy) += Ez_profile(ix,iy)*time_coeff*stepsize;
                    Y.EMF().Bx()(ix,iy) += Bx_profile(ix,iy)*time_coeff*stepsize;
                    Y.EMF().By()(ix,iy) += By_profile(ix,iy)*time_coeff*stepsize;
                    Y.EMF().Bz()(ix,iy) += Bz_profile(ix,iy)*time_coeff*stepsize;
                }
            }
        }
    }
}
//**************************************************************
//**************************************************************
//--------------------------------------------------------------
void Setup_Y::initialize(State1D &Y, Grid_Info &grid){
//--------------------------------------------------------------
//   Initialization of the harmonics and the fields
//--------------------------------------------------------------

    Y = static_cast<complex<double> >(0.0);


    valarray<double> dens_profile( grid.axis.Nx(0));
    valarray<double> temp_profile( grid.axis.Nx(0));
    valarray<double> f10x_profile( grid.axis.Nx(0));
    valarray<double> f20x_profile( grid.axis.Nx(0));
    valarray<double> pedestal_profile( grid.axis.Nx(0));

    valarray<double> hydro_dens_profile( grid.axis.Nx(0));
    valarray<double> hydro_temp_profile( grid.axis.Nx(0));
    valarray<double> hydro_vel_profile( grid.axis.Nx(0));
    valarray<double> hydro_Z_profile( grid.axis.Nx(0));
    
    for (size_t s(0); s < Input::List().qs.size(); ++s) {

        temp_profile = 0.0;
        
        Parser::parseprofile(grid.axis.x(0), Input::List().dens_profile_str[s], dens_profile);
        // std::cout << "\n11\n";
        Parser::parseprofile(grid.axis.x(0), Input::List().temp_profile_str[s], temp_profile);
        // std::cout << "\n12\n";
        Parser::parseprofile(grid.axis.x(0), Input::List().f10x_profile_str[s], f10x_profile);
        // std::cout << "\n13\n";
        // Parser::parseprofile(grid.axis.x(0), Input::List().f20x_profile_str[s], f20x_profile);
        // Parser::parseprofile(grid.axis.x(0), Input::List().f_pedestal[s], pedestal_profile);
        // std::cout << "\n14\n";
        init_f0(s, Y.SH(s,0,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, Y.DF(s).mass(), pedestal_profile);
        // std::cout << "\n15\n";
        if (Input::List().init_f1) init_f1(s, Y.SH(s,1,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, f10x_profile, Y.SH(s,0,0), Y.DF(s).mass());
        
        // if (Input::List().init_f2) init_f2(s, Y.SH(s,2,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, f20x_profile, Y.DF(s).mass());

    }
    
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      Setup the electromagnetic fields
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Y.EMF() = static_cast<complex<double> >(0.0);

    
    Parser::parseprofile(grid.axis.x(0), Input::List().hydro_dens_profile_str, hydro_dens_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().hydro_temp_profile_str, hydro_temp_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().hydro_vel_profile_str, hydro_vel_profile);
    Parser::parseprofile(grid.axis.x(0), Input::List().hydro_Z_profile_str, hydro_Z_profile);

    
    Y.HYDRO().vxarray()             = 0.0;
    Y.HYDRO().vyarray()             = 0.0;
    Y.HYDRO().vzarray()             = 0.0;
    Y.HYDRO().densityarray()        = 1.0;
    Y.HYDRO().temperaturearray()    = 1.0;
    Y.HYDRO().Zarray()              = 1.0;


    for (size_t i(0); i < temp_profile.size(); ++i)
    {

        Y.HYDRO().density(i)        = hydro_dens_profile[i];
        Y.HYDRO().temperature(i)    = hydro_temp_profile[i];
        Y.HYDRO().vx(i)             = hydro_vel_profile[i];
        Y.HYDRO().Z(i)              = hydro_Z_profile[i];
//            Y.HYDRO().vy(i)       = hydro_vel_profile[i];
//            Y.HYDRO().vz(i)       = hydro_vel_profile[i];
        // std::cout << "velocity(" << grid.axis.x(0)[i] << ")=" << Y.HYDRO().velocity(i) << "\n";
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Setup_Y::initialize(State2D &Y, Grid_Info &grid){
//--------------------------------------------------------------
//   Initialization of the harmonics and the fields
//--------------------------------------------------------------

    Y = static_cast<complex<double> >(0.0);


    Array2D<double> dens_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> temp_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> f10x_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> f20x_profile( grid.axis.Nx(0), grid.axis.Nx(1));
    Array2D<double> pedestal_profile( grid.axis.Nx(0), grid.axis.Nx(1));

    Array2D<double> hydro_dens_profile( grid.axis.Nx(0),grid.axis.Nx(1));
    Array2D<double> hydro_temp_profile( grid.axis.Nx(0),grid.axis.Nx(1));
    Array2D<double> hydro_vel_profile( grid.axis.Nx(0),grid.axis.Nx(1));
    Array2D<double> hydro_Z_profile( grid.axis.Nx(0),grid.axis.Nx(1));



    for (size_t s(0); s < Input::List().qs.size(); ++s) {

        temp_profile = 0.0;
        
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().dens_profile_str[s], dens_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().temp_profile_str[s], temp_profile);
        Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().f10x_profile_str[s], f10x_profile);
        // Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().f_pedestal[s], pedestal_profile);
        
        init_f0(s, Y.SH(s,0,0), grid.axis.p(s), grid.axis.x(0), grid.axis.x(1), dens_profile, temp_profile, Y.DF(s).mass(), pedestal_profile);

        // if (Input::List().init_f1) init_f1(s, Y.SH(s,1,0), grid.axis.p(s), grid.axis.x(0),  grid.axis.x(1), dens_profile, temp_profile, f10x_profile, Y.SH(s,0,0), Y.DF(s).mass());
        
        // if (Input::List().init_f2) init_f2(s, Y.SH(s,2,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, f20x_profile, Y.DF(s).mass());

    }

    // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    // {
    //     for (size_t iy(0);iy<Y.SH(0,0,0).numy();++iy)
    //     {

    //         std::cout << "T[" << ix << "," << iy << "] = " << temp_profile(ix,iy) << "\n";
    //     }
    // }
    
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      Setup the electromagnetic fields
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Y.EMF() = static_cast<complex<double> >(0.0);

    
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().hydro_dens_profile_str, hydro_dens_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().hydro_temp_profile_str, hydro_temp_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().hydro_vel_profile_str, hydro_vel_profile);
    Parser::parseprofile(grid.axis.x(0), grid.axis.x(1), Input::List().hydro_Z_profile_str, hydro_Z_profile);

    
    Y.HYDRO().vxarray()             = 0.0;
    Y.HYDRO().vyarray()             = 0.0;
    Y.HYDRO().vzarray()             = 0.0;
    Y.HYDRO().densityarray()        = 1.0;
    Y.HYDRO().temperaturearray()    = 1.0;
    Y.HYDRO().Zarray()              = 1.0;


    for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    {
        for (size_t iy(0);iy<Y.SH(0,0,0).numy();++iy)
        {
            Y.HYDRO().density(ix,iy)        = hydro_dens_profile(ix,iy);
            Y.HYDRO().temperature(ix,iy)    = hydro_temp_profile(ix,iy);
            Y.HYDRO().Z(ix,iy)              = hydro_Z_profile(ix,iy);
        }
    }

}
//--------------------------------------------------------------
//**************************************************************
/**
 * @brief      Initializes super-Maxwellian or Maxwell-Juttner according to density and temperature and m (super-maxwellian-ness, m=2 for Maxwellian)
 */

//-----------------------------------------------------------------------------
void Setup_Y:: init_f0(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
                       valarray<double>& density, valarray<double>& temperature, const double mass, const valarray<double>& pedestal){
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    double alpha, coeff, coefftemp, coefftemp_relativistic;
    double m = Input::List().super_gaussian_m;

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/alpha/alpha/alpha/tgamma(3.0/m);
    coeff *= sqrt(M_PI)/4.0;

    for (int j(0); j < h.numx(); ++j){

        // std::cout << "\n";

        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]*mass,1.5);
        
        // coefftemp_relativistic = density[j]/(4.0*M_PI*temperature[j]*pow(mass,3.)*boost::math::cyl_bessel_k(2.,1./temperature[j]));

        // std::cout << "\n\n coefftemp_relativistic = " << coefftemp_relativistic << "\n\n";

        
        for (int k(0); k < h.nump(); ++k){

            // New formulation for temperature distribution and super-Gaussians
            h(k,j) = coefftemp*exp(-1.0*pow((p[k])/alpha/sqrt(2.0*temperature[j]*mass),m));
            h(k,j)+= pedestal[j];
            // std::cout << "f0[" << p[k] << "] = " << exp(-1.0*pow((p[k])/alpha/sqrt(2.0*temperature[j]*mass),m)) << "\n";
            // 
            // Maxwell-Jutner distribution
            if (Input::List().relativity) 
                h(k,j) = coefftemp_relativistic*exp(-sqrt(1.0+p[k]*p[k])/temperature[j]);

            // std::cout << "\n\n h(" << k << "," << j << ") = " << h(k,j) << "\n\n";
        }
        // exit(1);    

    }

}
//----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void Setup_Y:: init_f0(size_t s, SHarmonic2D& h, const valarray<double>& p, const valarray<double>& x, const valarray<double>& y,
                       Array2D<double>& density, Array2D<double>& temperature, const double mass, const Array2D<double>& pedestal){
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    double alpha, coeff, coefftemp, coefftemp_relativistic;
    double m = Input::List().super_gaussian_m;

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/alpha/alpha/alpha/tgamma(3.0/m);
    coeff *= sqrt(M_PI)/4.0;

    for (int ix(0); ix < h.numx(); ++ix)
    {
        for (int iy(0); iy < h.numy(); ++iy)
        {

            // std::cout << "\n";

            coefftemp = coeff*density(ix,iy)/pow(2.0*M_PI*temperature(ix,iy)*mass,1.5);
            
            // coefftemp_relativistic = density(ix,iy)/(4.0*M_PI*temperature(ix,iy)*pow(mass,3.)*boost::math::cyl_bessel_k(2.,1./temperature(ix,iy)));

            // std::cout << "\n\n coefftemp_relativistic = " << coefftemp_relativistic << "\n\n";

            
            for (int k(0); k < h.nump(); ++k){

                // New formulation for temperature distribution and super-Gaussians
                h(k,ix,iy) = coefftemp*exp(-1.0*pow((p[k])/alpha/sqrt(2.0*temperature(ix,iy)*mass),m));
                h(k,ix,iy)+= pedestal(ix,iy);
                // std::cout << "f0[" << ix << "," << iy << "," << p[k] << "] = " << h(k,ix,iy) << "\n";
                // 
                // Maxwell-Jutner distribution
                // if (Input::List().relativity) 
                //     h(k,ix,iy) = coefftemp_relativistic*exp(-sqrt(1.0+p[k]*p[k])/temperature(ix,iy));

                // std::cout << "\n\n h(" << k << "," << j << ") = " << h(k,j) << "\n\n";
            }
              

        }
        // exit(1);  
    }

}
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
/**
 * @brief      To initialize f1 by using the steady state linearized vlasov equation. This menas using
 *              collision term, v grad f, and E df0/dv term in the VFP equation. This gives f1 = -v^3 * ( v df0/dx - E(x) df/dv)
 *              Assuming no df0/dx, it gives the result below.
 */
//----------------------------------------------------------------------------------------------------------------------------
void Setup_Y:: init_f1(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
                       valarray<double>& density, valarray<double>& temperature, valarray<double>& f10x, const SHarmonic1D& f0, const double mass){
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
    double alpha, coeff, coefftemp;
    double m = Input::List().super_gaussian_m;

    SHarmonic1D df0(f0); df0.Dp();

    double idp = -0.5/(p[1]-p[0]);

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/alpha/alpha/alpha/tgamma(3.0/m);
    coeff *= sqrt(M_PI)/4.0;


    for (int j(0); j < h.numx(); ++j){

        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]*mass,1.5)*f10x[j];
        for (int k(0); k < h.nump(); ++k){
            // New formulation for temperature distribution and super-Gaussians
            h(k,j) = idp*df0(k,j)*pow(p[k],3.0)*f10x[j];
        }

    }

}
//----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void Setup_Y:: init_f2(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
                       valarray<double>& density, valarray<double>& temperature, valarray<double>& f20x, const double mass){
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    double alpha, coeff, coefftemp;
    double m = Input::List().super_gaussian_m;

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/alpha/alpha/alpha/tgamma(3.0/m);
    coeff *= sqrt(M_PI)/4.0;


    for (int j(0); j < h.numx(); ++j){

        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]/mass,1.5);
        for (int k(0); k < h.nump(); ++k){
            // New formulation for temperature distribution and super-Gaussians
            h(k,j)= 2.0*coefftemp*exp(-pow((p[k]-3.0)/alpha/sqrt(2.0*temperature[j]/mass),m));

        }

    }

}
//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
