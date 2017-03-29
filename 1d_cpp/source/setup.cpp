/*!\brief  Initialization routines - Definitions
* \author PICKSC
 * \date   August 9, 2016
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
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "exprtk.hpp"

//  Declerations
#include "input.h"
#include "state.h"
#include "formulary.h"
#include "setup.h"
//  	#include "parallel.h"


//**************************************************************
//**************************************************************
// INITIALIZATION
//**************************************************************
//**************************************************************
// void Setup_Y::prepare(Grid_Info &grid, State1D& Y)
// {
//     parseprofile(grid.axis.x(0), Input::List().Ex_profile_str, Input::List().Ex_profile);
//     parseprofile(grid.axis.x(0), Input::List().Ey_profile_str, Input::List().Ey_profile);
//     parseprofile(grid.axis.x(0), Input::List().Ez_profile_str, Input::List().Ez_profile);
//     parseprofile(grid.axis.x(0), Input::List().Bx_profile_str, Input::List().Bx_profile);
//     parseprofile(grid.axis.x(0), Input::List().By_profile_str, Input::List().By_profile);
//     parseprofile(grid.axis.x(0), Input::List().Bz_profile_str, Input::List().Bz_profile);
// }
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

    parseprofile(grid.axis.x(0), Input::List().ex_profile_str, Ex_profile);
    parseprofile(grid.axis.x(0), Input::List().ey_profile_str, Ey_profile);
    parseprofile(grid.axis.x(0), Input::List().ez_profile_str, Ez_profile);
    parseprofile(grid.axis.x(0), Input::List().bx_profile_str, Bx_profile);
    parseprofile(grid.axis.x(0), Input::List().by_profile_str, By_profile);
    parseprofile(grid.axis.x(0), Input::List().bz_profile_str, Bz_profile);

    parseprofile(time, Input::List().ex_time_profile_str, ex_time_coeff);
    parseprofile(time, Input::List().ey_time_profile_str, ey_time_coeff);
    parseprofile(time, Input::List().ez_time_profile_str, ez_time_coeff);
    parseprofile(time, Input::List().bx_time_profile_str, bx_time_coeff);
    parseprofile(time, Input::List().by_time_profile_str, by_time_coeff);
    parseprofile(time, Input::List().bz_time_profile_str, bz_time_coeff);

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
//---------------------------------------------------------------------------
void Setup_Y::applytravelingwave(Grid_Info &grid, State1D& Y, double time)
{

    for (size_t n(0); n < Input::List().num_waves; ++n)
    {

        valarray<double> Ex_profile( grid.axis.Nx(0));
        valarray<double> Ey_profile( grid.axis.Nx(0));
        valarray<double> Ez_profile( grid.axis.Nx(0));
        valarray<double> Bx_profile( grid.axis.Nx(0));
        valarray<double> By_profile( grid.axis.Nx(0));
        valarray<double> Bz_profile( grid.axis.Nx(0));


        /// Parse the strings   
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().ex_wave_profile_str[n], Ex_profile);
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().ey_wave_profile_str[n], Ey_profile);
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().ez_wave_profile_str[n], Ez_profile);
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().bx_wave_profile_str[n], Bx_profile);
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().by_wave_profile_str[n], By_profile);
        parsetwovariableprofile(grid.axis.x(0), time, Input::List().bz_wave_profile_str[n], Bz_profile);

        double time_coeff(-1.0); 

        // parseprofile(time, Input::List().wave_time_envelope_str[n], time_coeff);

        double pulse_start(Input::List().trav_wave_center[n] - 
                        Input::List().trav_wave_flat[n]*0.5 - 
                        Input::List().trav_wave_rise[n] );

        double pulse_end(Input::List().trav_wave_center[n] +
                        Input::List().trav_wave_flat[n]*0.5 + 
                        Input::List().trav_wave_fall[n] );

        double normalized_time(0.0);

        // std::cout << "\n pulse_start, end = " << pulse_start << "," << pulse_end;


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
        }

        if (time_coeff > 0.0)
        {
            for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
            {
                Y.EMF().Ex()(ix) += Ex_profile[ix]*time_coeff;
                Y.EMF().Ey()(ix) += Ey_profile[ix]*time_coeff;
                Y.EMF().Ez()(ix) += Ez_profile[ix]*time_coeff;
                Y.EMF().Bx()(ix) += Bx_profile[ix]*time_coeff;
                Y.EMF().By()(ix) += By_profile[ix]*time_coeff;
                Y.EMF().Bz()(ix) += Bz_profile[ix]*time_coeff;
            }    
        }
    }
}
//**************************************************************
//--------------------------------------------------------------
void Setup_Y::initialize(State1D &Y, Grid_Info &grid){
//--------------------------------------------------------------
//   Initialization of the harmonics and the fields
//--------------------------------------------------------------

    Y = 0.0;

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

        Y.SH(s,0,0) = 0.0;

        temp_profile = 0.0;

        parseprofile(grid.axis.x(0), Input::List().dens_profile_str[s], dens_profile);
        parseprofile(grid.axis.x(0), Input::List().temp_profile_str[s], temp_profile);
        
        parseprofile(grid.axis.x(0), Input::List().f10x_profile_str[s], f10x_profile);
        // parseprofile(grid.axis.x(0), Input::List().f20x_profile_str[s], f20x_profile);
        parseprofile(grid.axis.x(0), Input::List().f_pedestal[s], pedestal_profile);

        init_f0(s, Y.SH(s,0,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, Y.DF(s).mass(), pedestal_profile);

        if (Input::List().init_f1) init_f1(s, Y.SH(s,1,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, f10x_profile, Y.SH(s,0,0), Y.DF(s).mass());

        // if (Input::List().init_f2) init_f2(s, Y.SH(s,2,0), grid.axis.p(s), grid.axis.x(0), dens_profile, temp_profile, f20x_profile, Y.DF(s).mass());

    }

//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      Setup the electromagnetic fields
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Y.EMF() = static_cast<complex<double> >(0.0);

    parseprofile(grid.axis.x(0), Input::List().hydro_dens_profile_str, hydro_dens_profile);
    parseprofile(grid.axis.x(0), Input::List().hydro_temp_profile_str, hydro_temp_profile);
    parseprofile(grid.axis.x(0), Input::List().hydro_vel_profile_str, hydro_vel_profile);
    parseprofile(grid.axis.x(0), Input::List().hydro_Z_profile_str, hydro_Z_profile);


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
//**************************************************************


//-----------------------------------------------------------------------------
void Setup_Y:: init_f0(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
                       valarray<double>& density, valarray<double>& temperature, const double mass, const valarray<double>& pedestal){
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
    double alpha, coeff, coefftemp;
    double m = Input::List().super_gaussian_m;

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/alpha/alpha/alpha/tgamma(3.0/m);
    coeff *= sqrt(M_PI)/4.0;

//    std::cout << "\n alpha = " << alpha << ", coeff = "  << coeff << "\n";
    for (int j(0); j < h.numx(); ++j){

        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]*mass,1.5);
        for (int k(0); k < h.nump(); ++k){
            // New formulation for temperature distribution and super-Gaussians
            h(k,j) = coefftemp*exp(-1.0*pow((p[k])/alpha/sqrt(2.0*temperature[j]*mass),m));
            h(k,j)+= pedestal[j];

            // h(k,j)= coefftemp*exp(-pow((p[k]-3.0)/alpha/sqrt(2.0*temperature[j]/mass),m));
            // if (p[k] < 3.0) h(k,j) += 1e-2;
//            if (j==3) std :: cout << "p = " << p[k] << ", h = " << h(k,j).real() << "\n";

        }

        // h(0,j) = 1.425577e2;
        // h(1,j) = 1.497164e1;
        // h(2,j) = 1.527841e-1;
        // h(3,j) = 1.515014e-4;
        // h(4,j) = 1.459772e-8;
        // h(5,j) = 1.366732e-13;
        // h(6,j) = 1.243401e-19;
        // h(7,j) = 1.099181e-26;



    }

}
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
void Setup_Y:: init_f1(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
                       valarray<double>& density, valarray<double>& temperature, valarray<double>& f10x, const SHarmonic1D& f0, const double mass){
//----------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
    double alpha, coeff, coefftemp;
    double m = Input::List().super_gaussian_m;

    SHarmonic1D df0(f0); df0.Dp();

    double idp = -0.5/(p[1]-p[0]);

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    coeff = m/4.0*sqrt(M_PI)*alpha/alpha/alpha/tgamma(3.0/m);

    for (int j(0); j < h.numx(); ++j){

        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]*mass,1.5)*f10x[j];
        for (int k(0); k < h.nump(); ++k){
            // New formulation for temperature distribution and super-Gaussians
//            h(k,j)  = coefftemp*exp(-pow((p[k])/alpha/sqrt(2.0*temperature[j]*mass),m))*pow(p[k],3.0)*(p[k]/mass/temperature[j]);
            h(k,j) = idp*df0(k,j)*pow(p[k],3.0)*f10x[j];

            // h(k,j)= coefftemp*exp(-pow((p[k]-3.0)/alpha/sqrt(2.0*temperature[j]/mass),m));
            // if (p[k] < 3.0) h(k,j) += 1e-2;
//            if (j==3) std :: cout << "p = " << p[k] << ", h = " << h(k,j) << "\n";

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
    //Matrix2D<double> coeff(pt);
    double m = 2.0;

    alpha = sqrt(3.0*tgamma(3.0/m)/2.0/tgamma(5.0/m));
    // coeff = m/4.0/M_PI/alpha/alpha/alpha/tgamma(3.0/m);
    coeff = 1.0;

    for (int j(0); j < h.numx(); ++j){
        // std::cout << "temp = " << temperature[j] << "\n";
        // coefftemp = coeff*density(s,j)/pow(temperature(s,j),1.5);
        coefftemp = coeff*density[j]/pow(2.0*M_PI*temperature[j]/mass,1.5);
        for (int k(0); k < h.nump(); ++k){
            // New formulation for temperature distribution and super-Gaussians


            // h(k,j)= coefftemp*exp(-pow(p[k]/alpha/sqrt(temperature(s,j)),m));
            h(k,j)= 2.0*coefftemp*exp(-pow((p[k]-3.0)/alpha/sqrt(2.0*temperature[j]/mass),m));
            // if (j==3) std :: cout << "p = " << p[k] << ", h = " << h(k,j) << "\n";
            // h(i,j,k) *= coeff(j,k);
            // std::cout << h(k,j) << "\n";
        }


        // Old formulation for momentum distribution
        // alpha(j,k) = 1.0/(2.0*pt(j,k)*pt(j,k));
        // coeff(j,k) = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt(j,k)*pt(j,k)*pt(j,k));
        // std::cout << "pt=" << pt(j,k) << "\n";
    }

}
//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
void Setup_Y::checkparse(parser_t& parser, std::string& expression_str, expression_t& expression){
//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------


    if (!parser.compile(expression_str,expression))
    {
        printf("Error: %s\tExpression: %s\n",
               parser.error().c_str(),
               expression_str.c_str());

        for (std::size_t i = 0; i < parser.error_count(); ++i)
        {
            error_t error = parser.get_error(i);
            printf("Error: %02d Position: %02d Type: [%s] Msg: %s Expr: %s\n",
                   static_cast<int>(i),
                   static_cast<int>(error.token.position),
                   exprtk::parser_error::to_str(error.mode).c_str(),
                   error.diagnostic.c_str(),
                   expression_str.c_str());
        }

        exit(1);
    }


}
//----------------------------------------------------------------------------------------------------------------------------    
/**
 * @brief      { function_description }
 *
 * @param      str_profile  The string profile
 * @param      profile      The profile
 */
//----------------------------------------------------------------------------------------------------------------------------
void Setup_Y::parseprofile( const valarray<double>& grid, std::string& str_profile, valarray<double>& profile){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        profile = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr;
        symbol_table.add_variable("x",xstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);

        for (size_t i(0); i < profile.size(); ++i) {
            xstr = grid[i];
            // std::cout<< "i= " << i << ",xstr = " << xstr << "\n";

            // result =
            // Temperature_map(s,i) = expression_Temperature.value();
            profile[i] = expression.value();
            // std::cout  << i << " , " << grid[i] << " , " << profile[i] << "\n";
            // printf("Result: %10.5f\n",expression.value());
        }
    }

    if (str_profile[0] == 'p'){
        std::string pcw_pairs = str_profile.substr(posL+1,posR-(posL+1));
        vector<size_t> positions_comma, positions_semicolon; // holds all the positions that sub occurs within str
        vector<double> loc,val;
        size_t pos = pcw_pairs.find(',', 0);
        while(pos != string::npos)
        {
            positions_comma.push_back(pos);
            pos = pcw_pairs.find(',',pos+1);
        }

        pos = pcw_pairs.find(';', 0);
        while(pos != string::npos)
        {
            positions_semicolon.push_back(pos);
            pos = pcw_pairs.find(';',pos+1);
        }

        if (positions_comma.size() != positions_semicolon.size()){
            std::cout << "Error in piecewise input, the values must be supplied in pairs with ',' in between the two values in each pair, and ending each pair with a ';'. E.g. {0.1,1.0;0.2,1.5;} \n";

            exit(1);
        }


        loc.push_back(std::stod(pcw_pairs.substr(0,positions_comma[0]-1)));
        val.push_back(std::stod(pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)));

        // std::cout << "loc0 = " << loc[0] << endl;
        // std::cout << "val0 = " << val[0] << endl;

        for (size_t i(1); i < positions_comma.size(); ++i) {

            loc.push_back(std::stod(pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)));
            val.push_back(std::stod(pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)));

            // std::cout << "loci = " << loc[i] << endl;
            // std::cout << "vali = " << val[i] << endl;

        }


        double tmp_x, xx;
        size_t t(0);
        for (size_t i(0); i < grid.size(); ++i) {
            t = 1;
            xx = grid[i];
            if (grid[i] < loc[0])               xx += loc[loc.size()-1] - loc[0];
            if (grid[i] > loc[loc.size()-1])  xx += loc[0]-loc[loc.size()-1] ;

            while ((xx > loc[t]) && (t < loc.size()-1)) ++t;
            tmp_x = val[t-1] + (val[t]-val[t-1])/(loc[t]-loc[t-1])
                               * (xx - loc[t-1]);
            // std::cout << i << "     " << x[i] << "      " << tmp_x << std::endl;
            //             for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x;
            profile[i] = tmp_x;
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------------    
/**
 * @brief      { function_description }
 *
 * @param      str_profile  The string profile
 * @param      profile      The profile
 */
//----------------------------------------------------------------------------------------------------------------------------
void Setup_Y::parseprofile(const double& input, std::string& str_profile, double& output){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        output = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr;
        symbol_table.add_variable("t",xstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);

        // for (size_t i(0); i < profile.size(); ++i) {
        xstr = input;
        // std::cout<< "i= " << i << ",xstr = " << xstr << "\n";

        // result = 
        // Temperature_map(s,i) = expression_Temperature.value();
        output = expression.value();
        // std::cout  << " \n 10: " << input << " , " << output << "\n";
        // exit(1);
        // printf("Result: %10.5f\n",expression.value());

    }

    //     if (str_profile[0] == 'p'){
    //         std::string pcw_pairs = str_profile.substr(posL+1,posR-(posL+1));
    //         vector<size_t> positions_comma, positions_semicolon; // holds all the positions that sub occurs within str
    //         vector<double> loc,val;
    //         size_t pos = pcw_pairs.find(',', 0);
    //         while(pos != string::npos)
    //         {
    //             positions_comma.push_back(pos);
    //             pos = pcw_pairs.find(',',pos+1);
    //         }

    //         pos = pcw_pairs.find(';', 0);
    //         while(pos != string::npos)
    //         {
    //             positions_semicolon.push_back(pos);
    //             pos = pcw_pairs.find(';',pos+1);
    //         }

    //         if (positions_comma.size() != positions_semicolon.size()){
    //             std::cout << "Error in piecewise input, the values must be supplied in pairs with ',' in between the two values in each pair, and ending each pair with a ';'. E.g. {0.1,1.0;0.2,1.5;} \n";

    //             exit(1);
    //         }


    //         loc.push_back(std::stod(pcw_pairs.substr(0,positions_comma[0]-1)));
    //         val.push_back(std::stod(pcw_pairs.substr(positions_comma[0]+1,positions_semicolon[0]-positions_comma[0]-1)));

    //         // std::cout << "loc0 = " << loc[0] << endl;
    //         // std::cout << "val0 = " << val[0] << endl;

    //         for (size_t i(1); i < positions_comma.size(); ++i) {

    //                 loc.push_back(std::stod(pcw_pairs.substr(positions_semicolon[i-1]+1,positions_comma[i]-positions_semicolon[i-1]-1)));
    //                 val.push_back(std::stod(pcw_pairs.substr(positions_comma[i]+1,positions_semicolon[i]-positions_comma[i]-1)));

    //                 // std::cout << "loci = " << loc[i] << endl;
    //                 // std::cout << "vali = " << val[i] << endl;

    //         }


    //         double tmp_x, xx;
    //         size_t t(0);
    //         for (size_t i(0); i < grid.size(); ++i) {
    //             t = 1;
    //             xx = grid[i];
    //             if (grid[i] < loc[0])               xx += loc[loc.size()-1] - loc[0];
    //             if (grid[i] > loc[loc.size()-1])  xx += loc[0]-loc[loc.size()-1] ;

    //             while ((xx > loc[t]) && (t < loc.size()-1)) ++t; 
    //             tmp_x = val[t-1] + (val[t]-val[t-1])/(loc[t]-loc[t-1]) 
    //                                 * (xx - loc[t-1]);
    //             // std::cout << i << "     " << x[i] << "      " << tmp_x << std::endl;
    // //             for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x; 
    //             profile[i] = tmp_x; 
    //         }
    //     } 
}
//----------------------------------------------------------------------------------------------------------------------------    
/**
 * @brief      { function_description }
 *
 * @param      str_profile  The string profile
 * @param      profile      The profile
 */
//----------------------------------------------------------------------------------------------------------------------------
void Setup_Y::parsetwovariableprofile(const valarray<double>& grid, const double& input, std::string& str_profile, valarray<double>& output){

    /// Find curly brackets
    std::size_t posL = str_profile.find("{");
    std::size_t posR = str_profile.find("}");

    /// Uniform profile, convert string to double and set.
    if (str_profile[0] == 'c')
    {
        // std::cout<< "expression_str is  " << str_profile.substr(posL+1,posR-(posL+1)) << "\n";
        output = std::stod(str_profile.substr(posL+1,posR-(posL+1)));
    }

    /// Profile defined by function
    if (str_profile[0] == 'f')
    {
        /// Grab everything inside curly brackets
        std::string expression_str = str_profile.substr(posL+1,posR-(posL+1));

        symbol_table_t symbol_table;
        symbol_table.add_constants();
        double xstr, tstr;
        symbol_table.add_variable("x",xstr);
        symbol_table.add_variable("t",tstr);

        expression_t expression;

        expression.register_symbol_table(symbol_table);

        parser_t parser;

        checkparse(parser, expression_str, expression);
        tstr = input;

        for (size_t i(0); i < output.size(); ++i) {
            xstr = grid[i];
            
            output[i] = expression.value();            
        }
        
        

    }

}

//----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
