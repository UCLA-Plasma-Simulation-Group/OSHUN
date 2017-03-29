/*! \brief Collisions - Definitions
 * \author PICKSC
 * \date   September 1, 2016
 * \file   collisions.cpp
 * 
 * Contains:
 * 1) the explicit energy-conserving algorithm for effect of 
 * electron-electron collisions on f_00
 * 2) implicit algorithm used for e-e + e-i collisions on f_lm
 * And all their containers.
 */


//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "exprtk.hpp"

//  Declarations
#include "input.h"
#include "state.h"
#include "setup.h"
#include "formulary.h"
#include "nmethods.h"
#include "collisions.h"




//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//*********************************************************************************************
self_f00_implicit_step::self_f00_implicit_step(const size_t &nump, const double &pmax, const double &_mass, const double& _deltat, bool& _ib):
        vr(Algorithms::MakeCAxis(0.0,pmax,nump)), C_RB(0.0,nump+1), D_RB(0.0,nump+1), delta_CC(0.0,nump+1),
        mass(_mass), dt(_deltat),
        dvr(0.0,nump), vrh(0.0,nump), dtoverv2(nump),
        p2dp(0.0,nump),p2dpm1(0.0,nump),phdp(0.0,nump),phdpm1(0.0,nump), p4dp(0.0,nump), laser_Inv_Uav6(0.0,nump),
        I4_Lnee(0.0), ib(_ib)
{
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// Collisions
    //classical electron radius
    double re(2.8179402894e-13);   //< classical electron radius
    double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re)); //< sqrt(4*pi*n_{e0}*re).
    c_kpre = re*kp; //< (4*pi*re)^1.5*sqrt(n_{e0})/3.

    // Determine vr
    for (size_t i(0); i < vr.size(); ++i) {
        vr[i] /= mass;
    }

    for (size_t i(0); i < vr.size()-1; ++i) {
        vrh[i] = 0.5*(vr[i+1] - vr[i]);
        dvr[i] = (vr[i+1] - vr[i]);
    }
    vrh[vrh.size()-1] = 0.0;
    dvr[vrh.size()-1] = dvr[vrh.size()-2];

    for (size_t i(1); i < vr.size(); ++i) {
        p2dp[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);
        p2dpm1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);
        p4dp[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);
//        p4dpm1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);
        phdp[i]   = 0.5 * vrh[i]           * (vrh[i]-vrh[i-1]);
        phdpm1[i] = 0.5 * vrh[i-1]         * (vrh[i]-vrh[i-1]);
    }


    for (size_t i(0); i < vr.size(); ++i) {
        dtoverv2[i] = dt / vr[i] / vr[i] / dvr[i];
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// Laser
    for (size_t i(0); i < vr.size(); ++i) {
        laser_Inv_Uav6[i] = pow( 2.0/(vr[i+1]+vr[i]), 6);
    }

    double omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Input::List().lambda_0));
    double omega_p(5.64 * 1.0e+4*sqrt(Input::List().density_np));
    vw_coeff_cube = omega_p/omega_0 * c_kpre;
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// Cooling

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void self_f00_implicit_step::update_C_Rosenbluth(valarray<double> &fin) {
    /// Remember that C is defined on the boundaries of the velocity grid
    /// Therefore, C[0] is C_{1/2} aka C(v=0)
    /// and, C[1] is C_{3/2} aka C(v[1st point, 0th index C style])

    C_RB[0] = 0.0;
    I4_Lnee = 0.0;
    for (int n(1); n < C_RB.size(); ++n) {
//        C_RB[n]  = p2dp[n-1] * fin[n-1] + p2dpm1[n-1] * fin[n - 2];
        C_RB[n]  = vr[n - 1] * vr[n - 1] * dvr[n - 1] * fin[n - 1];


        C_RB[n] += C_RB[n - 1];


        I4_Lnee  += p4dp[n - 1] * fin[n - 1];

    }
    C_RB *= 4.0 * M_PI;
    I4_Lnee *= 4.0 * M_PI;
}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
double self_f00_implicit_step::update_D_Rosenbluth(const size_t& k, valarray<double>& fin, const double& delta) {

    double answer(0.0);

    valarray<double> innersum(fin.size() - 1);  /// Only needs to be np-1

    int n(innersum.size()-1);                /// Initialize at last point, which is the sum from np-1 to np-1
    /// Now n = np - 1

    innersum[n]  = (1.0 - delta) * fin[n + 1] + delta * fin[n];         /// So innersum = f_{np} * (1-delta_{np-1/2})
    innersum[n] *= (vr[n + 1] * vr[n + 1] - vr[n] * vr[n]);             ///             + f_{np-1} * (delta_{np-1/2})


    --n;

    for (n; n > -1; --n) {
        innersum[n] = (1.0 - delta) * fin[n + 1] + delta * fin[n];
        innersum[n] *= (vr[n + 1] * vr[n + 1] - vr[n] * vr[n]);

        innersum[n] += innersum[n + 1];

//        std::cout << "\n innersum[" << n << "] = " << innersum[n] << "\n";

    }

    /// Using indexing from Kingham2004. v is defined from 1 to nv, obviously
    /// So k = 1 is the first velocity cell and distribution function point.
    /// Therefore, all the l's should be referenced to with l - 1
    for (size_t l(1); l < k + 1; ++l) {
        answer += vr[l - 1] * vr[l - 1] * dvr[l - 1] * innersum[l - 1];
    }

    answer *= 4.0 * M_PI;
    answer *= 1.0 / (vr[k - 1] + vr[k]);

    return answer;
}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void self_f00_implicit_step::update_D_and_delta(valarray<double>& fin){
    size_t iterations(0);
    bool iteration_check(0);
    double D(0.0);
    double Dold(10.0);
    double delta(0.5);


    /// Remember that D and delta are defined on the boundaries of the velocity grid
    /// Therefore, D[0] = D_{1/2} = D(v=0)
    /// and, D[1] D = D_{3/2} D(v[1st point, 0th index C-style])


    D_RB[0]             = 0.0;
    D_RB[D_RB.size()-1] = 0.0;

    delta_CC[0]             = 0.5;
    delta_CC[D_RB.size()-1] = 0.0;

    for (size_t k(1); k < fin.size(); ++k){
        delta = 0.5;
        D = 0.0;
        Dold = 10.0;
        iterations = 0;
        iteration_check = 0;
        do
        {
            D = update_D_Rosenbluth(k,fin,delta);
//            std::cout << ", D[" << k << "] = " << D;
            if (fabs(D-Dold) < Input::List().RB_D_tolerance*(1.0+fabs(D+Dold))) iteration_check = 1;

            delta = calc_delta_ChangCooper(k, C_RB[k], D);
//            std::cout << ", delta[" << k << "] = " << delta;
            ++iterations;

            if (iterations > Input::List().RB_D_itmax) iteration_check = 1;
            Dold = D;

        } while(!iteration_check);
//        std::cout << "\n Chang-cooper iterations = " << iterations << "\n";


        D_RB[k]     = D;
        delta_CC[k] = delta;



//        std::cout << "\n C[" << k << "] = " << C_RB[k];
//        std::cout << ", D[" << k << "] = " << D;
//        std::cout << ", delta[" << k << "] = " << delta;


    }
//        exit(1);

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
double self_f00_implicit_step::calc_delta_ChangCooper(const size_t& k, const double& C, const double& D){
    double answer(0.5);
    double W = 0.5*(dvr[k-1]+dvr[k])*C/D;

    if (W > 1.0e-8){
        answer = 1.0/W - 1.0/(exp(W)-1.0);
    }
    return answer;
}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void self_f00_implicit_step::update_D_inversebremsstrahlung(const double& Z0, const double& heating_coefficient, const double& vos){

    double vw_cube;// = vw_coeff_cube; * ZLn_ei * 4.0*M_PI*I2;;

    double temperature = I4_Lnee/3.0/C_RB[C_RB.size()-1];
    vw_cube  = Z0*formulas.Zeta*formulas.LOGei(C_RB[C_RB.size()-1],temperature,Z0*formulas.Zeta);  ///< ZLogLambda
    vw_cube *= vw_coeff_cube * C_RB[C_RB.size()-1];

    double ZLnee = Z0*formulas.Zeta*formulas.LOGee(C_RB[C_RB.size()-1],temperature);

    double g, b0, nueff, xsi;

//    xsi     = 3.84+(142.59-65.48*vos/sqrt(temperature)/(27.3*vos/sqrt(temperature)+vos*vos/temperature);


    for (size_t ip(1); ip < D_RB.size()-1; ++ip){
//        b0      = (vos*vos+5.0*vr[ip]*vr[ip])/(5.0*vr[ip]*vr[ip]);
//        nueff   = ZLnee/pow(vr[ip]*vr[ip]+vos*vos/xsi,1.5);
//        g       = 1.0 - b0*nueff*vw_cube/(1.0+b0*b0*nuOeff*vw_cube);



        g = laser_Inv_Uav6[ip] * vw_cube * vw_cube;
        g = 1.0 / (1.0 + g);

        D_RB[ip] += heating_coefficient * g / vr[ip];
    }

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void self_f00_implicit_step::takestep(valarray<double>  &fin, valarray<double> &fh, const double& Z0, const double& vos, const double& cooling) {

    double collisional_coefficient;
    double heating_coefficient;

    Array2D<double> LHS(fin.size(),fin.size());

    ///  Calculate Rosenbluth and Chang-Cooper quantities
    update_C_Rosenbluth(fin);   /// Also fills in I4_Lnee (the temperature for the Lnee calculation)
    update_D_and_delta(fin);    /// And takes care of boundaries

    /// Normalizing quantities (Inspired by previous collision routines and OSHUN notes by M. Tzoufras)
    collisional_coefficient  = formulas.LOGee(C_RB[C_RB.size()-1],I4_Lnee/3.0/C_RB[C_RB.size()-1]);
    collisional_coefficient *= 4.0*M_PI/3.0*c_kpre;



    heating_coefficient  = formulas.Zeta*Z0*formulas.LOGei(C_RB[C_RB.size()-1],I4_Lnee/3.0/C_RB[C_RB.size()-1],Z0*formulas.Zeta);  ///< ZLogLambda
    heating_coefficient *= c_kpre / 6.0 * pow(vos,2.0) * C_RB[C_RB.size()-1];
    heating_coefficient /= collisional_coefficient;

    if (ib) update_D_inversebremsstrahlung(Z0, heating_coefficient, vos);

    /// Fill in matrix

    size_t ip(0);

    /// Boundaries by hand -- This operates on f(0)
    LHS(ip, ip + 1) = - dtoverv2[ip] * collisional_coefficient
                      * (C_RB[ip + 1] * (1.0 - delta_CC[ip + 1])
                         + D_RB[ip + 1] / dvr[ip + 1]);

    LHS(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                            * (C_RB[ip + 1] * delta_CC[ip + 1] - D_RB[ip + 1] / dvr[ip + 1]
                               - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);


    for (ip = 1; ip < fin.size() - 1; ++ip){
        LHS(ip, ip + 1) = - dtoverv2[ip] * collisional_coefficient
                          * (C_RB[ip + 1] * (1.0 - delta_CC[ip + 1])
                             + D_RB[ip + 1] / dvr[ip + 1]);

        LHS(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                                * (C_RB[ip + 1] * delta_CC[ip + 1] - D_RB[ip + 1] / dvr[ip + 1]
                                   - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);

        LHS(ip, ip - 1) = dtoverv2[ip] * collisional_coefficient
                          * (C_RB[ip] * delta_CC[ip]
                             - D_RB[ip] / dvr[ip]);

    }

    ip = fin.size() - 1;

    LHS(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                            * (C_RB[ip + 1] * delta_CC[ip + 1]
                               - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);

    LHS(ip, ip - 1) = dtoverv2[ip] * collisional_coefficient
                      * (C_RB[ip] * delta_CC[ip]
                         - D_RB[ip] / dvr[ip]);

//    std::cout << "\n\n LHS = \n";
//    for (size_t i(0); i < LHS.dim1(); ++i) {
//        std::cout << "i = " << i << " :::: ";
//        for (size_t j(0); j < LHS.dim2(); ++j) {
//            std::cout << LHS(i, j) << "   ";
//        }
//        std::cout << "\n";
//    }


    Thomas_Tridiagonal(LHS,fin,fh);

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void self_f00_implicit_step::getleftside(valarray<double>  &fin, const double& Z0, const double& vos, const double& cooling,  Array2D<double> &LHStemp) {

    double collisional_coefficient;
    double heating_coefficient;

//    Array2D<double> LHS(fin.size(),fin.size());

    ///  Calculate Rosenbluth and Chang-Cooper quantities
    update_C_Rosenbluth(fin);   /// Also fills in I4_Lnee (the temperature for the Lnee calculation)
    update_D_and_delta(fin);    /// And takes care of boundaries

    /// Normalizing quantities (Inspired by previous collision routines and OSHUN notes by M. Tzoufras)
    collisional_coefficient  = formulas.LOGee(C_RB[C_RB.size()-1],I4_Lnee/3.0/C_RB[C_RB.size()-1]);
    collisional_coefficient *= 4.0*M_PI/3.0*c_kpre;



    heating_coefficient  = formulas.Zeta*Z0*formulas.LOGei(C_RB[C_RB.size()-1],I4_Lnee/3.0/C_RB[C_RB.size()-1],Z0*formulas.Zeta);  ///< ZLogLambda
    heating_coefficient *= c_kpre / 6.0 * pow(vos,2.0) * C_RB[C_RB.size()-1];
    heating_coefficient /= collisional_coefficient;

    if (ib) update_D_inversebremsstrahlung(Z0, heating_coefficient,vos);

    /// Fill in matrix

    size_t ip(0);

    /// Boundaries by hand -- This operates on f(0)
    LHStemp(ip, ip + 1) = - dtoverv2[ip] * collisional_coefficient
                      * (C_RB[ip + 1] * (1.0 - delta_CC[ip + 1])
                         + D_RB[ip + 1] / dvr[ip + 1]);

    LHStemp(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                            * (C_RB[ip + 1] * delta_CC[ip + 1] - D_RB[ip + 1] / dvr[ip + 1]
                               - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);


    for (ip = 1; ip < fin.size() - 1; ++ip){
        LHStemp(ip, ip + 1) = - dtoverv2[ip] * collisional_coefficient
                          * (C_RB[ip + 1] * (1.0 - delta_CC[ip + 1])
                             + D_RB[ip + 1] / dvr[ip + 1]);

        LHStemp(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                                * (C_RB[ip + 1] * delta_CC[ip + 1] - D_RB[ip + 1] / dvr[ip + 1]
                                   - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);

        LHStemp(ip, ip - 1) = dtoverv2[ip] * collisional_coefficient
                          * (C_RB[ip] * delta_CC[ip]
                             - D_RB[ip] / dvr[ip]);

    }

    ip = fin.size() - 1;

    LHStemp(ip    , ip) = 1.0 - dtoverv2[ip] * collisional_coefficient
                            * (C_RB[ip + 1] * delta_CC[ip + 1]
                               - C_RB[ip] * (1.0 - delta_CC[ip]) - D_RB[ip] / dvr[ip]);

    LHStemp(ip, ip - 1) = dtoverv2[ip] * collisional_coefficient
                      * (C_RB[ip] * delta_CC[ip]
                         - D_RB[ip] / dvr[ip]);

//    std::cout << "\n\n LHS = \n";
//    for (size_t i(0); i < LHS.dim1(); ++i) {
//        std::cout << "i = " << i << " :::: ";
//        for (size_t j(0); j < LHS.dim2(); ++j) {
//            std::cout << LHS(i, j) << "   ";
//        }
//        std::cout << "\n";
//    }

//
//    Thomas_Tridiagonal(LHS,fin,fh);

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//*********************************************************************************************
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
self_f00_implicit_collisions::self_f00_implicit_collisions(const DistFunc1D& DFin, const double& deltat)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        :   fin(0.0, DFin(0,0).nump()), fout(0.0, DFin(0,0).nump()),
            xgrid(Algorithms::MakeCAxis(Input::List().xminLocal[0],Input::List().xmaxLocal[0],Input::List().NxLocal[0])),
            IB_heating(Input::List().IB_heating), MX_cooling(Input::List().MX_cooling),
            heatingprofile(0.0,DFin(0,0).numx()), coolingprofile(0.0,DFin(0,0).numx()),
            ib(((DFin.q() == -1) && (DFin.mass() == 1))),
            collide(DFin(0,0).nump(),DFin.pmax(),DFin.mass(), deltat, ib)
{

    Nbc = Input::List().BoundaryCells;
    szx = Input::List().NxLocal[0];  // size of useful x axis
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
void self_f00_implicit_collisions::loop(SHarmonic1D& f00, valarray<double>& Zarray, const double time, SHarmonic1D& f00h){

    //-------------------------------------------------------------------
    //  This loop scans all the locations in configuration space
    //  and calls the implicit Chang-Cooper/Langdon/Epperlein algorithm
    //-------------------------------------------------------------------
    //      Initialization
    //-------------------------------------------------------------------

    int ixtimesnump;
    Array2D<double> LHS(f00.numx()*f00.nump(),f00.numx()*f00.nump());
    valarray<double> RHSv(0.0,f00.numx()*f00.nump());
    valarray<double> solution(0.0,f00.numx()*f00.nump());

    Array2D<double> LHStemp(f00.nump(),f00.nump());

    double timecoeff;
    if (IB_heating && ib){

        /// Get time and heating profile
        /// Ray-trace would go here
        Setup_Y::parseprofile(xgrid, Input::List().intensity_profile_str, heatingprofile);
        Setup_Y::parseprofile(time, Input::List().intensity_time_profile_str, timecoeff);

        /// Make vos(x,t)
        heatingprofile *= (Input::List().lambda_0 * sqrt(7.3e-19*Input::List().I_0))*timecoeff;

    }


    for (size_t ix(0); ix < szx; ++ix){
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Copy data for a specific location in space to valarray
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t ip(0); ip < fin.size(); ++ip){
            fin[ip] = (f00(ip,ix)).real();
//            std::cout << "fin[" << ip << "] = " << fin[ip] << "\n";
        }
//
//        collide.takestep(fin,fout,Zarray[ix],heatingprofile[ix],coolingprofile[ix]);

        collide.getleftside(fin,Zarray[ix],heatingprofile[ix],coolingprofile[ix],LHStemp);

        ixtimesnump = ix * f00.nump();

        for (size_t ip(0); ip < fin.size(); ++ip){
            if (ip == 0)
            {
                LHS(ip + ixtimesnump, ip + 1 + ixtimesnump) = LHStemp(ip, ip + 1);
            }
            else if (ip == f00.nump() - 1)
            {
                LHS(ip + ixtimesnump, ip - 1 + ixtimesnump) = LHStemp(ip, ip - 1);
            }
            else {
                LHS(ip + ixtimesnump, ip - 1 + ixtimesnump) = LHStemp(ip, ip - 1);
                LHS(ip + ixtimesnump, ip + 1 + ixtimesnump) = LHStemp(ip, ip + 1);
            }

            LHS(ip + ixtimesnump, ip + ixtimesnump) = LHStemp(ip, ip);

            RHSv[ip + ixtimesnump] = fin[ip];
//            std::cout << "fin[" << ip << "] = " << fin[ip] << "\n";
        }


        // Return updated data to the harmonic
        /*for (size_t ip(0); ip < fin.size(); ++ip){
            f00h(ip,ix) = fout[ip];
//            std::cout << "fout[" << ip << "] = " << fout[ip] << "\n";
        }*/
    }

    Thomas_Tridiagonal(LHS,RHSv,solution);

    /// Unpack
    for (size_t ix(0); ix < szx; ++ix){
        ixtimesnump = ix * f00.nump();
        for (size_t ip(0); ip < fin.size(); ++ip){
            f00h(ip,ix) = solution[ip + ixtimesnump];
        }
    }


    //-------------------------------------------------------------------

}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
/**
 * @brief      Constructor that needs a distribution function input.
 *
 * @param      fslope  The distribution function that is input and transformed.
*/
self_f00_explicit_step::self_f00_explicit_step(const size_t& nump, const double& pmax, const double& _mass)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
// : fh(0.0, nump), ///< Dummy distribution
// :fh(fslope), ///< Dummy distribution
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///       Velocity axis
        :
        mass(_mass),
//        vr(Algorithms::MakeAxis(pmax/(2.0*nump-1.0),pmax,nump)),
        vr(Algorithms::MakeCAxis(0.0,pmax,nump)),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        U4(0.0, nump), ///< Constants for Integrals
        U4m1(0.0, nump), ///< Constants for Integrals
        U2(0.0, nump), ///< Constants for Integrals
        U2m1(0.0, nump), ///< Constants for Integrals
        U1(0.0, nump), ///< Constants for Integrals
        U1m1(0.0, nump), ///< Constants for Integrals
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        J1(0.0, nump), ///<      Integrals
        I2(0.0, nump), ///<      Integrals
        I4(0.0, nump), ///<      Integrals
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        U3(0.0, nump), ///<      Integrals
        Qn(0.0, nump), ///<      Integrals
        Pn(0.0, nump)  ///<      Integrals
///<

{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //classical electron radius
    double re(2.8179402894e-13);   //< classical electron radius
    double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re)); //< sqrt(4*pi*n_{e0}*re).
    c_kpre = 4.0*M_PI/3.0*re*kp; //< (4*pi*re)^1.5*sqrt(n_{e0})/3.
    NB = Input::List().NB_algorithms;

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Determine vr
    for (size_t i(0); i < vr.size(); ++i) {
        vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));
    }

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Determine U4, U4m1, U2, U2m1, U1, U1m1
    for (size_t i(1); i < U4.size(); ++i) {
        U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);
        U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);
        U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);
        U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);
        U1[i]   = 0.5 * vr[i]            * (vr[i]-vr[i-1]);
        U1m1[i] = 0.5 * vr[i-1]          * (vr[i]-vr[i-1]);
    }

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Determine U3
    for (size_t i(0); i < U3.size(); ++i) {
        U3[i] = pow(vr[i],3);
    }
    // Determine Qn
    Qn[0] = 1.0 / ((vr[0]*vr[0]*vr[1])/2.0);
    for (size_t i(1); i < Qn.size()-1; ++i) {
        Qn[i] = 1.0 / (vr[i]*vr[i]*(vr[i+1]-vr[i-1])/2.0);
    }
    // Determine Pn
    for (size_t i(0); i < Pn.size()-1; ++i) {
        Pn[i] = 1.0 / ((vr[i+1]-vr[i])/2.0*(vr[i+1]+vr[i]));


    }

}
//--------------------------------------------------------------

/// Calculate G for low momentum cells
///
/// @param[in]  n     Number of low momentum cells being treated differently.
/// @param[in]  fin   Input distribution function.
///
/// @return     Returns G value for n cells.
///
double self_f00_explicit_step::G(const int& n, const valarray<double>& fin) {
//-------------------------------------------------------------------
    double i2s, i4s;
    double f00( (fin[0] - fin[1]*(vr[0]*vr[0])/(vr[1]*vr[1]))/
                (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1])) );

    i2s = f00*pow(vr[n],3)/3.0 + (fin[1]-f00)*pow(vr[n],5)/(vr[1]*vr[1])*0.2;
    i4s = f00*pow(vr[n],5)*0.2 + (fin[1]-f00)*pow(vr[n],7)/(vr[1]*vr[1]*7.0);

    return fin[n]*i4s + (pow(vr[n],3)*fin[n]-3.0*i2s) * J1[n];
}
//-------------------------------------------------------------------

/// Compute collision integrals and advance to next step
///
/// @param[in]  fin   Input distribution function
///
/// @return     Output distribution function
void self_f00_explicit_step::takestep(const valarray<double>& fin, valarray<double>& fh) {
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------
    double Ln_ee;

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate the integrals in a standard manner
    I4[0] = 0;
    for (int n(1); n < I4.size(); ++n) {
        I4[n]  = U4[n]*fin[n]+U4m1[n]*fin[n-1];
        I4[n] += I4[n-1];
    }

    I2[0] = 0;
    for (int n(1); n < I2.size(); ++n) {
        I2[n]  = U2[n]*fin[n]+U2m1[n]*fin[n-1];
        I2[n] += I2[n-1];
    }

    J1[J1.size()-1] = 0;
    for (int n(J1.size()-2); n > -1; --n) {
        J1[n]  = U1[n+1]*fin[n+1]+U1m1[n+1]*fin[n];
        J1[n] += J1[n+1];

    }


//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Use I2 and I4 to calculate the Coulomb logarithm now, the
//     arrays I2 I4 will be modified later and it will be too late.

    // if (mass == 100.0){
    // std::cout<< "\n\n\n temp = " << mass*I4[I4.size()-1]/3.0/I2[I4.size()-1] << "\n\n\n";
    // }
    Ln_ee = formulas.LOGee(4.0*M_PI*I2[I2.size()-1],I4[I4.size()-1]/3.0/I2[I4.size()-1]);

    // // I2_temperature = I2[I2.size()-1];
    // for (int k(0); k < I4.size(); ++k){
    //    I4[k] /=  vr[k]*vr[k];
    // }


//     <><><><><><><><><><><><><><><><><><>
//     Tony's Energy Conserving Algorithm 
//     <><><><><><><><><><><><><><><><><><>

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate G using the standard integrals
//     J1 needs to be used again later so it is not modified!
    for (int n(0); n < I4.size(); ++n) {
        I2[n] *= J1[n];          // J1(n) * I2(n)
        I2[n] *= -3.0;           // -3 * J1(n) * I2(n)
        I4[n] += U3[n] * J1[n];  // I4(n) + u_n^3 * J1(n)
        I4[n] *= fin[n];         // fn * I4(n) + u_n^3 * fn * J1(n)
        I4[n] += I2[n];          // Gn = fn * I4(n) + u_n^3 * fn * J1(n) - 3 * J1(n) * I2(n)
        //

    }

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate G assuming a parabolic f(v << vt)
    for (int n(0); n < NB; ++n) {
        I4[n] = G(n,fin);

    }

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Find -DG
    fh[0]  = (-1.0)*I4[0];
    for (int n(0); n < I4.size()-1; ++n) I4[n] -= I4[n+1];

//     Find -DG/(vDv)
    fh[0] *= 2.0/ (vr[0]*vr[0]);
    for (int n(0); n < I4.size()-1; ++n) I4[n] *= Pn[n];

//     Find DDG/(vDv)
    fh[0] -= I4[0];
    for (int n(0); n < I4.size()-1; ++n) I4[n] -= I4[n+1];

//     Find DDG/(v^3*DDv)
    fh[0] *= Qn[0];
    for (int n(0); n < I4.size()-1; ++n) fh[n+1] = Qn[n+1]*I4[n];

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Normalize
    fh *=  c_kpre *  Ln_ee;
}
//-------------------------------------------------------------------
//*******************************************************************

//*******************************************************************
//*******************************************************************
//  Runge Kutta loop for the electron-electron collisions
//*******************************************************************

//-------------------------------------------------------------------
//*******************************************************************
//------------------------------------------------------------------------------
/// @brief      Constructor for RK4 method on f00
///
/// @param      fin         Input distribution
/// @param[in]  tout_start  Hmm...
///
self_f00_RKfunctor::self_f00_RKfunctor(const size_t& nump, const double& pmax, const double& mass)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        :collide(nump,pmax,mass){}
//--------------------------------------------------------------

/**
 * { item_description }
 */
void self_f00_RKfunctor::operator()(const valarray<double>& fin, valarray<double>& fslope) {
    collide.takestep(fin,fslope);
}

void self_f00_RKfunctor::operator()(const valarray<double>& fin, valarray<double>& fslope, size_t dir) {}
void self_f00_RKfunctor::operator()(const valarray<double>& fin, const valarray<double>& f2in, valarray<double>& fslope) {}
//*******************************************************************
//*******************************************************************
//   Definition for the Explicit Integration method for the 
//   energy relaxation of the distribution function
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
self_f00_explicit_collisions::self_f00_explicit_collisions(const DistFunc1D& DFin, const double& deltat)//, int tout_start) //DistFunc1D& DFin,
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        : fin(0.0, DFin(0,0).nump()),
          RK(fin),
          rkf00(DFin(0,0).nump(), DFin.pmax(), DFin.mass()),
          num_h(size_t(deltat/Input::List().small_dt)+1)
{
    h = deltat/static_cast<double>(num_h);

    Nbc = Input::List().BoundaryCells;
    szx = Input::List().NxLocal[0];  // size of useful x axis
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
void self_f00_explicit_collisions::loop(SHarmonic1D& f00,SHarmonic1D& f00h){

    //-------------------------------------------------------------------
    //  This loop scans all the locations in configuration space
    //  and calls the RK4 for explicit integration at each location
    //-------------------------------------------------------------------
    //      Initialization

    for (size_t ix(0); ix < szx; ++ix){
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // Copy data for a specific location in space to valarray
        for (size_t ip(0); ip < fin.size(); ++ip){
            fin[ip] = (f00(ip,ix)).real();
        }

        // std::cout << "\n\n h = " << num_h << ", \n \n";
        // Time loop: Update the valarray
        for (size_t h_step(0); h_step < num_h; ++h_step){
            fin = RK(fin,h,&rkf00);
        }

        // Return updated data to the harmonic
        for (int ip(0); ip < fin.size(); ++ip){
            f00h(ip,ix) = fin[ip];
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    }


}

//*******************************************************************
//*******************************************************************
//  Definition of collisions with high order harmonics FLM
//*******************************************************************
//*******************************************************************


//*******************************************************************
//--------------------------------------------------------------
self_flm_implicit_step::self_flm_implicit_step(double pmax, size_t nump, double _mass)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        :
//       Pre-Calculated Constants

//        vr(Algorithms::MakeAxis(pmax/(2.0*nump-1.0),pmax,nump)),
        vr(Algorithms::MakeCAxis(0.0,pmax,nump)),
//       Constants for Integrals
        U4(0.0,  nump),
        U4m1(0.0,nump),
        U2(0.0,  nump),
        U2m1(0.0,nump),
        U1(0.0,  nump),
        U1m1(0.0,nump),
//       Integrals
        J1m(0.0,  nump),
        I0(0.0,  nump),
        I2(0.0,  nump),
        Scattering_Term(0.0, nump),
        df0(0.0,  nump),
        ddf0(0.0,  nump),
//       Matrices
        Alpha_Tri(nump, nump),
        if_tridiagonal(Input::List().if_tridiagonal),
        mass(_mass)
//       Make matrix vk/vn
{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    double re(2.8179402894e-13);           //classical electron radius
    double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re));

    kpre = re*kp;

    // Determine vr
    for (size_t i(0); i < vr.size(); ++i) {
//        vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));
        vr[i] /= mass;
    }

//    U4[0] = pow(vr[0],4)          * (vr[1]-vr[0]);
//    U2[0] = vr[0]*vr[0]           * (vr[1]-vr[0]);
//    U1[0] = vr[0]                 * (vr[1]-vr[0]);

    // Determine U4, U4m1, U2, U2m1, U1, U1m1 for the integrals
    for (size_t i(1); i < U4.size(); ++i) {
        U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);
        U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);
        U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);
        U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);
        U1[i]   = 0.5 * vr[i]            * (vr[i]-vr[i-1]);
        U1m1[i] = 0.5 * vr[i-1]          * (vr[i]-vr[i-1]);

//        U4[i]   = pow(vr[i],4)          * (vr[i]-vr[i-1]);
//        U2[i]   = vr[i]*vr[i]           * (vr[i]-vr[i-1]);
//        U1[i]   = vr[i]                 * (vr[i]-vr[i-1]);

    }
}
//--------------------------------------------------------------

//------------------------------------------------------------------------------
/// @brief      Resets coefficients and integrals to use in the matrix solve.
///
/// @param[in]  fin      Input distribution function
/// @param[in]  Delta_t  timestep
///
void  self_flm_implicit_step::reset_coeff(const valarray<double>& fin, double Zvalue, const double Delta_t) {
//-------------------------------------------------------------------
//  Reset the coefficients based on f_0^0 
//-------------------------------------------------------------------

//     Calculate Dt
    Dt = Delta_t;

    double idp = 1.0/(vr[2] - vr[1]);


//     INTEGRALS
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     Temperature integral I2 = 4*pi / (v^2) * int_0^v f(u)*u^4du
//     Density integral I0 = 4*pi*int_0^v f(u)*u^2du 

//    I2[0] = U4[0]*fin[0];
//    I0[0] = U2[0]*fin[0];
    I2[0] = 0;
    I0[0] = 0;

    for (int k(1); k < I2.size(); ++k) {
        I2[k]  = U4[k]*fin[k]+U4m1[k]*fin[k-1];
        I2[k] += I2[k-1];

        I0[k]  = U2[k]*fin[k]+U2m1[k]*fin[k-1];
        I0[k] += I0[k-1];

//        I2[k]  = U4[k]*fin[k];
//        I2[k] += I2[k-1];
//
//        I0[k]  = U2[k]*fin[k];
//        I0[k] += I0[k-1];
    }
    I0_density = 4.0*M_PI*I0[I0.size()-1];


    I2_temperature = I2[I2.size()-1]/3.0/I0[I0.size()-1];


    //     Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
    J1m[J1m.size()-1] = 0;
    for (int k(J1m.size()-2); k > -1; --k) {
        J1m[k]  = U1[k+1]*fin[k+1]+U1m1[k+1]*fin[k];
        J1m[k] += J1m[k+1];
    }


    /*J1m[J1m.size()-1] = fin[J1m.size()-1]*U1[U1.size()-1];

    *//*std::cout << "\nI2[" << U1.size()-1 << "] = " << I2[U1.size()-1] << endl;
    std::cout << "I0[" << U1.size()-1 << "] = " << I0[U1.size()-1] << endl;
    std::cout << "J1[" << U1.size()-1 << "] = " << J1m[U1.size()-1] << endl;
    std::cout << "\n\n";*//*


    for (int k(J1m.size()-2); k > -1; --k) {
        J1m[k]  = U1[k]*fin[k];//+U1m1[k+1]*fin[k];
        J1m[k] += J1m[k+1];

        *//*std::cout << "\nI2[" << k << "] = " << I2[k] << endl;
        std::cout << "I0[" << k << "] = " << I0[k] << endl;
        std::cout << "J1[" << k << "] = " << J1m[k] << endl;
        std::cout << "\n\n";*//*
    }*/


    for (int k(0); k < J1m.size(); ++k){
        I2[k]  /=  vr[k] * vr[k] / 4.0 / M_PI;
        I0[k]  *= 4.0 * M_PI;
        J1m[k] *= 4.0 * M_PI * vr[k];
    }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//     COULOMB LOGARITHMS
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    _LOGee  = formulas.LOGee(I0_density,I2_temperature);
    _ZLOGei = (formulas.Zeta*Zvalue)*formulas.LOGei(I0_density,I2_temperature,(formulas.Zeta*Zvalue));
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


//     BASIC INTEGRALS FOR THE TRIDIAGONAL PART
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     Temporary arrays
    valarray<double>  TriI1(fin), TriI2(fin);

    for (int i(0); i < TriI1.size(); ++i){
        TriI1[i] = I0[i] + (2.0*J1m[i] - I2[i]) / 3.0 ;   // (-I2 + 2*J_{-1} + 3*I0) / 3
        TriI2[i] = ( I2[i] + J1m[i] ) / 3.0 ;             // ( I2 + J_{-1} ) / 3
    }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


//     SCATTERING TERM
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Scattering_Term    = TriI1;
    Scattering_Term   *= _LOGee;                        // Electron-electron contribution
    Scattering_Term[0] = 0.0;
    Scattering_Term   += _ZLOGei * I0_density;          // Ion-electron contribution

    for (int i(0); i < Scattering_Term.size(); ++i){    // Multiply by 1/v^3
        Scattering_Term[i] /= pow(vr[i],3);
    }
    Scattering_Term *=  kpre * Dt;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


//     MAKE TRIDIAGONAL ARRAY
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Alpha_Tri = 0.0;


    double IvDnDm1, IvDnDp1, Ivsq2Dn;
    Alpha_Tri(0,0) = 8.0 * M_PI * fin[0];                                         //  8*pi*f0[0]
    for (int i(1); i < TriI1.size()-1; ++i){
        IvDnDm1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]   - vr[i-1]));       //  ( v*D_n*D_{n-1/2} )^(-1)
        IvDnDp1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1] - vr[i]  ));       //  ( v*D_n*D_{n+1/2} )^(-1)
        Ivsq2Dn = 1.0 / (vr[i] * vr[i]                 * (vr[i+1] - vr[i-1]));       //  ( v^2 * 2*D_n )^(-1)
        Alpha_Tri(i, i  ) = 8.0 * M_PI * fin[i] - TriI2[i] * (IvDnDm1 + IvDnDp1);
        Alpha_Tri(i, i-1) = TriI2[i] * IvDnDm1 - TriI1[i] * Ivsq2Dn;
        Alpha_Tri(i, i+1) = TriI2[i] * IvDnDp1 + TriI1[i] * Ivsq2Dn;
    }

    /*Alpha_Tri(0,0) = 8.0 * M_PI * fin[0];                                         //  8*pi*f0[0]
    for (int i(1); i < TriI1.size()-1; ++i){
        Alpha_Tri(i, i  ) = - 2.0 * TriI2[i] / vr[i] * idp * idp;
        Alpha_Tri(i, i-1) = TriI2[i] / vr[i] * idp * idp - TriI1[i] / 2.0 / vr[i] / vr[i] * idp;
        Alpha_Tri(i, i+1) = TriI2[i] / vr[i] * idp * idp + TriI1[i] / 2.0 / vr[i] / vr[i] * idp;
    }*/

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    Alpha_Tri *=  (-1.0) * _LOGee * kpre * Dt;         // (-1) because the matrix moves to the LHS in the equation
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/*

//     CALCULATE DERIVATIVES 
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    df0[0]  = (fin[1] - fin[0]);
    df0[0] /= (vr[1] - vr[0]);
//     Evaluate the derivative
    for (int n(1); n < fin.size()-1; ++n) {
        df0[n]  = fin[n+1] - fin[n-1];
        df0[n] /= vr[n+1] - vr[n-1];
    }
    df0[df0.size()-1]  = fin[df0.size()-1] - fin[df0.size()-2];
    df0[df0.size()-1] /= vr[df0.size()-1] - vr[df0.size()-2];

    ddf0[0]  = fin[0] - 2.0*fin[1] + fin[2];
    ddf0[0] /= vr[1] - vr[0];
    ddf0[0] /= vr[1] - vr[0];
//     Evaluate the second derivative
//        -df/dv_{n-1/2}
    for (int n(1); n < fin.size()-1; ++n) {
        ddf0[n]  = fin[n+1] - 2.0*fin[n] + fin[n-1];
        ddf0[n] /= vr[n] - vr[n-1];
        ddf0[n] /= vr[n] - vr[n-1];
//        ddf0[n]  = (fin[n-1]-fin[n]) ;
//        ddf0[n] /= vr[n] - vr[n-1] ;
    }

    ddf0[ddf0.size()-1] = fin[ddf0.size()-1] - 2.0*fin[ddf0.size()-2] + fin[ddf0.size()-3];
    ddf0[ddf0.size()-1] /= vr[vr.size()-1] - vr[vr.size()-2];
    ddf0[ddf0.size()-1] /= vr[vr.size()-1] - vr[vr.size()-2];

//        D(df/dv)/Dv
//    for (int n(1); n < fin.size()-1; ++n) {
//        ddf0[n] -= ddf0[n+1];
//        ddf0[n]  /= 0.5*(vr[n+1]-vr[n-1]);
//    }

//     Calculate zeroth cell
//    double f00 = ( fin[0] - ( (vr[0]*vr[0])/(vr[1]*vr[1]) ) *fin[1] )
//                 / (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]));
//    ddf0[0] = 2.0 * (fin[1] - f00) / (vr[1]*vr[1]);
//    df0[0] = ddf0[0] * vr[0];

//     Calculate 1/(2v)*(d^2f)/(dv^2),  1/v^2*df/dv
    for (int n(0); n < fin.size()-1; ++n) {
        df0[n]  /= vr[n]*vr[n];
        ddf0[n] /= 2.0  *vr[n];
    }
*/


    //     Evaluate the derivative
    for (int n(1); n < fin.size()-1; ++n) {
        df0[n]  = fin[n+1]-fin[n-1];
        df0[n] /= vr[n+1]-vr[n-1];
    }

//     Evaluate the second derivative
//        -df/dv_{n-1/2}
    for (int n(1); n < fin.size(); ++n) {
        ddf0[n]  = (fin[n-1]-fin[n]) ;
        ddf0[n] /= vr[n] - vr[n-1] ;
    }
//        D(df/dv)/Dv
    for (int n(1); n < fin.size()-1; ++n) {
        ddf0[n] -= ddf0[n+1];
        ddf0[n]  /= 0.5*(vr[n+1]-vr[n-1]);
    }

//     Calculate zeroth cell
    double f00 = ( fin[0] - ( (vr[0]*vr[0])/(vr[1]*vr[1]) ) *fin[1] )
                 / (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]));
    ddf0[0] = 2.0 * (fin[1] - f00) / (vr[1]*vr[1]);
    df0[0] = ddf0[0] * vr[0];

//     Calculate 1/(2v)*(d^2f)/(dv^2),  1/v^2*df/dv
    for (int n(0); n < fin.size()-1; ++n) {
        df0[n]  /= vr[n]*vr[n];
        ddf0[n] /= 2.0  *vr[n];
    }

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

}
//-------------------------------------------------------------------

//------------------------------------------------------------------------------
/// @brief      Perform a matrix solve to calculate effect of collisions on f >= 1
///
/// @param      fin   Input distribution function
/// @param[in]  el    Number of elements in matrix (?)
///
void  self_flm_implicit_step::advance(valarray<complex<double> >& fin, const int el) {
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------
    Array2D<double> Alpha(Alpha_Tri);
    valarray<complex<double> > fout(fin);


//      ZEROTH CELL FOR TRIDIAGONAL ARRAY
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (el > 1) {
        Alpha(0,0) = 0.0;
    }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    if ( !(if_tridiagonal) ) {

//         CONSTRUCT COEFFICIENTS and then full array
//         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        double LL(el);
        double A1(         (LL+1.0)*(LL+2.0) / ((2.0*LL+1.0)*(2.0*LL+3.0)) );
        double A2( (-1.0) *(LL-1.0)* LL      / ((2.0*LL+1.0)*(2.0*LL-1.0)) );
        double B1( (-1.0) *( 0.5 *LL*(LL+1.0) +(LL+1.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0)) );
        double B2( (       (-0.5)*LL*(LL+1.0) +(LL+2.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0)) );
        double B3(         ( 0.5 *LL*(LL+1.0) +(LL-1.0) ) / ((2.0*LL+1.0)*(2.0*LL-1.0)) );
        double B4(         ( 0.5 *LL*(LL+1.0) - LL      ) / ((2.0*LL+1.0)*(2.0*LL-1.0)) );

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        for (int i(0); i < Alpha.dim1()-1; ++i){
            double t1( A1*ddf0[i] + B1*df0[i] );
            t1 *= (-1.0) * _LOGee * kpre * Dt;
            double t2( A1*ddf0[i] + B2*df0[i] );
            t2 *= (-1.0) * _LOGee * kpre * Dt;
            double t3( A2*ddf0[i] + B3*df0[i] );
            t3 *= (-1.0) * _LOGee * kpre * Dt;
            double t4( A2*ddf0[i] + B4*df0[i] );
            t4 *= (-1.0) * _LOGee * kpre * Dt;

            Alpha(i,0) += t1 * ( 2.0*M_PI*pow(vr[0]/vr[i],el+2)*vr[0]*vr[0]*(vr[1]-vr[0]) );
            Alpha(i,0) += t3 * ( 2.0*M_PI*pow(vr[0]/vr[i],el)  *vr[0]*vr[0]*(vr[1]-vr[0]) );

            for (int j(1); j < i; ++j){
                Alpha(i,j) += t1 * ( 2.0*M_PI*pow(vr[j]/vr[i],el+2)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) );
                Alpha(i,j) += t3 * ( 2.0*M_PI*pow(vr[j]/vr[i],el)  *vr[j]*vr[j]*(vr[j+1]-vr[j-1]) );
            }

            Alpha(i,i) += t1 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) );
            Alpha(i,i) += t3 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) );

            Alpha(i,i) += t2 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) );
            Alpha(i,i) += t4 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) );

            for (int j(i+1); j < Alpha.dim2()-1; ++j){
                Alpha(i,j) += t2 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el-1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) );
                Alpha(i,j) += t4 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el+1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) );
            }
        }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    }

//      INCLUDE SCATTERING TERM
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    double ll1(static_cast<double>(el));
    ll1 *= (-0.5)*(ll1 + 1.0);

    for (int i(0); i < Alpha.dim1(); ++i){
        Alpha(i,i) += 1.0 - ll1 * Scattering_Term[i];
    }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




     /// SOLVE A * Fout  = Fin
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    for (int i(0); i < fin.size(); ++i){   // Estimate Fout = Fin /  A(i,i)
//        fout[i] /= Alpha(i,i);
//    }

    if ( if_tridiagonal ) {
        if ( !(Thomas_Tridiagonal(Alpha, fin, fout)) ) {  // Invert A * fout = fin
            cout << "WARNING: Matrix is not diagonally dominant" << endl;
        }
    }
    else {
        if ( !(Gauss_Seidel(Alpha, fin, fout)) ) {  // Invert A * fout = fin
            cout << "WARNING: Matrix is not diagonally dominant" << endl;
        }
    }





    fin = fout;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

}
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//*******************************************************************
//   Definition for the Anisotropic Collisions
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
self_flm_implicit_collisions::self_flm_implicit_collisions(const DistFunc1D &DFin, const double& deltat)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        : f00(0.0, DFin(0).nump()),
          fc(0.0,DFin(0).nump()),
          l0(DFin.l0()),
          m0(DFin.m0()),//Yin.SH(0,0).m0())
          implicit_step(DFin.pmax(),DFin(0).nump(),DFin.mass()),
          Dt(deltat)    //
{


    Nbc = Input::List().BoundaryCells;
    szx = Input::List().NxLocal[0];

    if (m0 == 0) {
        f1_m_upperlimit = 1;
    }
    else { 
        f1_m_upperlimit = 2;
    }
        

}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
void self_flm_implicit_collisions::advancef1(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh){
//-------------------------------------------------------------------
//  This is the collision calculation for the harmonics f10, f11 
//-------------------------------------------------------------------

// For every location in space within the domain of this node
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 		  
    for (size_t ix(0); ix < szx; ++ix){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
        // "00" harmonic --> Valarray
        for (size_t ip(0); ip < fc.size(); ++ip){
            f00[ip] = (DF(0,0)(ip,ix)).real();
        }
        // Reset the integrals and coefficients
        implicit_step.reset_coeff(f00, Zarray[ix], Dt);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
        // Loop over the harmonics for this (x,y)
        for(int m(0); m < f1_m_upperlimit; ++m){

            // This harmonic --> Valarray
            for (int ip(0); ip < DF(1,m).nump(); ++ip) {
                fc[ip] = DF(1,m)(ip,ix);
            }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      
            // Take an implicit step
            implicit_step.advance(fc, 1);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      

            //  Valarray --> This harmonic
            for (int ip(0); ip < DF(1,m).nump(); ++ip) {
                DFh(1,m)(ip,ix) = fc[ip];
            }
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    }
}
//-------------------------------------------------------------------



//-------------------------------------------------------------------
void self_flm_implicit_collisions::advanceflm(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh)
{
//-------------------------------------------------------------------
//  This is the calculation for the high order harmonics 
//    To be specific, this routine does l=2 to l0
//    To be specific, this routine differs from 'flm_loop1D' aboe only in one place: the loop over m all m from 0 to ((m0 < l)? m0:l)
//-------------------------------------------------------------------

    // For every location in space within the domain of this node

    for (size_t ix(0); ix < szx; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        // "00" harmonic --> Valarray
        for (size_t ip(0); ip < fc.size(); ++ip){
            f00[ip] = (DF(0,0)(ip,ix)).real();
        }
        // Reset the integrals and coefficients
        // --> Assume you did that when you called f1_loop
        implicit_step.reset_coeff(f00, Zarray[ix], Dt);
        // Loop over the harmonics for this (x,y)
        for(int l(2); l < l0+1 ; ++l){
            for(int m(0); m < ((m0 < l)? m0:l)+1; ++m){
//                std::cout << "(l,m) = " << l << "," << m << "\n";
                // This harmonic --> Valarray
                for (size_t ip(0); ip < fc.size(); ++ip){
                    fc[ip] = (DF(l,m))(ip,ix);
                }

                // Take an implicit step
                implicit_step.advance(fc, l);

                //  Valarray --> This harmonic
                for (size_t ip(0); ip < fc.size(); ++ip){
                    DFh(l,m)(ip,ix) = fc[ip];
                }
            }
        }
    }

    // return;  

}
//-------------------------------------------------------------------


////*******************************************************************
//-------------------------------------------------------------------
self_collisions::self_collisions(const DistFunc1D& DFin, const double& deltat)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        : self_f00_exp_collisions(DFin, deltat),
          self_f00_imp_collisions(DFin, deltat),
          self_flm_collisions(DFin, deltat){}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
void self_collisions::advancef00(SHarmonic1D& f00, valarray<double>& Zarray, const double time,SHarmonic1D& f00h){
//-------------------------------------------------------------------
    
    if (Input::List().f00_implicitorexplicit){
        if (Input::List().f00_implicitorexplicit == 2) {
            self_f00_imp_collisions.loop(f00,Zarray,time,f00h);
        }
        else if (Input::List().f00_implicitorexplicit == 1) {
            self_f00_exp_collisions.loop(f00,f00h);
        }   
    }
}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
void self_collisions::advanceflm(DistFunc1D& DFin, valarray<double>& Zarray, DistFunc1D& DFh){
//-------------------------------------------------------------------
    self_flm_collisions.advanceflm(DFin,Zarray,DFh);
}

//-------------------------------------------------------------------
void self_collisions::advancef1(DistFunc1D& DFin,  valarray<double>& Zarray, DistFunc1D& DFh){
//-------------------------------------------------------------------
    self_flm_collisions.advancef1(DFin,Zarray,DFh);
}
////*******************************************************************
//-------------------------------------------------------------------
collisions::collisions(const State1D& Yin, const double& deltat):Yh(Yin)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
{
    for(size_t s(0); s < Yin.Species(); ++s){
        self_coll.push_back( self_collisions(Yin.DF(s),deltat)  );


        if (Yin.Species() > 1){
            for (size_t sind(0); sind < Yin.Species(); ++sind){
//                if (s!=sind) unself_f00_coll.push_back(interspecies_f00_explicit_collisions(Yin.DF(s),Yin.DF(sind),deltat));
            }



        }

    }
    // exit(1);
}
//-------------------------------------------------------------------
void collisions::advance(State1D& Yin, const Clock& W)
//-------------------------------------------------------------------
{
    Yh = complex<double>(0.0,0.0);

    if (Input::List().f00_implicitorexplicit > 0) advancef0(Yin,W,Yh);
    
    if (Input::List().flm_collisions){
        advancef1(Yin,Yh);
        advanceflm(Yin,Yh);
    }    

    for (size_t s(0); s < Yin.Species(); ++s){
        for (size_t i(0); i < Yin.DF(s).dim(); ++i){
            Yin.DF(s)(i) = Yh.DF(s)(i);
        }
    }
}
//-------------------------------------------------------------------
void collisions::advancef0(State1D& Yin, const Clock& W, State1D& Yh)
//-------------------------------------------------------------------
{
    size_t sdummy(0);


    for(size_t s(0); s < Yin.Species(); ++s)
    {
        self_coll[s].advancef00(Yin.DF(s)(0,0),Yin.HYDRO().Zarray(),W.time(),Yh.DF(s)(0,0));

        if (Yin.Species() > 1)
        {
            for (size_t sind(0); sind < Yin.Species(); ++sind){
                // std::cout << "\n\n sdummy = " << sdummy << ",sind = " << sind << ", s = " << s << "\n\n";
                if (s!=sind)  {
                    // std::cout << "\n\n12\n\n";
                    // unself_f00_coll[sdummy].rkloop(Yin.SH(s,0,0),Yin.SH(sind,0,0)); ++sdummy;

                }
            }
        }


        // Yin.DF(s).checknan();

    }

}


//-------------------------------------------------------------------
void collisions::advancef1(State1D& Yin, State1D& Yh)
//-------------------------------------------------------------------
{
    size_t sdummy(0);


    for(size_t s(0); s < Yin.Species(); ++s)
    {
        self_coll[s].advancef1(Yin.DF(s), Yin.HYDRO().Zarray(), Yh.DF(s) );

        if (Yin.Species() > 1)
        {
            for (size_t sind(0); sind < Yin.Species(); ++sind){

                if (s!=sind)  {


                }
            }
        }

    }

}
//-------------------------------------------------------------------
void collisions::advanceflm(State1D& Yin, State1D& Yh)
//-------------------------------------------------------------------
{
    size_t sdummy(0);


    for(size_t s(0); s < Yin.Species(); ++s)
    {
        self_coll[s].advanceflm(Yin.DF(s), Yin.HYDRO().Zarray(), Yh.DF(s) );

        if (Yin.Species() > 1)
        {
            for (size_t sind(0); sind < Yin.Species(); ++sind){

                if (s!=sind)  {


                }
            }
        }


        // Yin.DF(s).checknan();

    }

}


//-------------------------------------------------------------------
//-------------------------------------------------------------------
vector<self_collisions> collisions::self(){

    return self_coll;

}
// t = tnew;

//-------------------------------------------------------------------



//*******************************************************************