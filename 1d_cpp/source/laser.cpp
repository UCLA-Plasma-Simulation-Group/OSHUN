/*!\brief  Laser source - Definitions
 * \author PICKSC
 * \date   August 9, 2016
 * \file   laser.cpp
 *
 * The laser package contains:
 * - Inverse bremsstrahlung heating
 * - External field driver (TBD)
 * - Phenomenological heat source (TBD)
 * - Ray tracing (TBD)
 * - Antenna (TBD)
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
#include "laser.h"


void Setup_Y::checkparse(parser_t& parser, std::string& expression_str, expression_t& expression);
void Setup_Y::parseprofile( const valarray<double>& grid, std::string& str_profile, valarray<double>& profile);

//*******************************************************************
//*******************************************************************
// Inverse Bremsstrahlung operator
//*******************************************************************
//*******************************************************************

//*******************************************************************
//--------------------------------------------------------------
IB_f00::IB_f00(valarray<double>& fslope, double pmax)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : fh(fslope),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Velocity
//         vr(Algorithms::MakeAxis(pmax/(2.0*fslope.size()-1.0),pmax,fslope.size())),
          vr(Algorithms::MakeCAxis(0.0,pmax,fslope.size())),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Constants for Integrals
          U4(0.0, fslope.size()),
          U4m1(0.0, fslope.size()),
          U2(0.0, fslope.size()),
          U2m1(0.0, fslope.size()),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Inv_Uav6(0.0, fslope.size()),
          gn(0.0, fslope.size()),
          Qn(0.0, fslope.size()),
          Pn(0.0, fslope.size())
{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //classical electron radius
    double re(2.8179402894e-13);
    double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re));
    double omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Input::List().lambda_0));
    double omega_p(5.64 * 1.0e+4*sqrt(Input::List().density_np));


    vw_coeff_cube = omega_p/omega_0 * kp*re;
    Qn_coeff = kp*re/6.0;

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
    }

    // Determine <vr>
    Inv_Uav6_nm1 = pow( 2.0/vr[0], 6);
    for (size_t i(0); i < vr.size(); ++i) {
        Inv_Uav6[i] = pow( 2.0/(vr[i+1]+vr[i]), 6);
    }

    // Determine Qn
    Qn[0] = 1.0 / ((vr[0]*vr[0]*vr[1])/2.0);
    for (size_t i(1); i < Qn.size()-1; ++i) {
        Qn[i] = 1.0 / (vr[i]*vr[i]*(vr[i+1]-vr[i-1])/2.0);
    }
    // Determine Pn
    Pnm1 = 2.0/(vr[0]*vr[0]);
    for (size_t i(0); i < Pn.size()-1; ++i) {
        Pn[i] = 1.0 / ((vr[i+1]-vr[i])/2.0*(vr[i+1]+vr[i]));
    }

}
//--------------------------------------------------------------

//-------------------------------------------------------------------
valarray<double>& IB_f00::Getslope(const valarray<double>& fin, const double Zval ,const double vos) {
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    double ZLn_ei;
    double I4, I2;
    double vw_cube, qn_c;
    // double p0overp1_sq(Input::List().pr(0)/ Input::List().pr(1));
    double p0overp1_sq(vr[1]/ vr[0]);
    p0overp1_sq *= p0overp1_sq;
    /*debug*/ // cout << p0overp1_sq<<"\n";

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double f00((fin[0]-fin[1]*p0overp1_sq)/(1.0-p0overp1_sq));

    double b0(0.0);
    double nueff(0.0);
    double xsi(0.0);
    /*debug*/ // cout << "(p0/p1)^2 =" << p0overp1_sq << ",   f[1] = " << fin[1] << ",   f[0] = " << fin[0]<< ",   f[-1] = " << f00 <<"\n";
    /*debug*/ // cout << "P[-1] = " << Pnm1 << ",    Inv_U_av-1 = " << Inv_Uav6_nm1 << "\n";


//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate the integrals in a standard manner
    I4 = 0;
    for (int n(1); n < U4.size(); ++n) {
        I4 += U4[n]*fin[n]+U4m1[n]*fin[n-1];
    }

    I2 = 0;
    for (int n(1); n < U2.size(); ++n) {
        I2 += U2[n]*fin[n]+U2m1[n]*fin[n-1];
    }

    ZLn_ei = (formulas.Zeta*Zval)*formulas.LOGei(4.0*M_PI*I2,I4/3.0/I2,formulas.Zeta*Zval);
    vw_cube = vw_coeff_cube * ZLn_ei * 4.0*M_PI*I2;
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xsi     = 3.84+(142.59-65.48*vos/sqrt(I2))/(27.3*vos/sqrt(I2)+vos*vos/I2);
    for (size_t i(0); i < gn.size(); ++i) { //last location initialized with "1"




        gn[i]   = Inv_Uav6[i]*vw_cube*vw_cube;
        gn[i] = 1.0/(1.0+gn[i]);



        /*debug*/ //  cout << "g["<< vr[i]<<"] = " << gn[i]<< "\n";
        // // // // // // // // // // // // // // // // // // // //
        // // // // // // // // // // // // // // // // // // // //
        /// Weng et. al. PRE 80 056406 2009
        ///
        /// Implementation untested
        ///
        // b0      = (vos*vos+5.0*vr[i]*vr[i])/5.0*vr[i]*vr[i];
        // nueff   = 1.0/pow((vr[i+1]*vr[i])/2.0+vos*vos/xsi,3.0);
        // gn[i]   = 1.0 - b0*nueff*vw_cube/(1.0+b0*b0*nueff*vw_cube);

        //
        //
        //
        //
    }

    /*debug*/// cout << "vos = "<< vos << ",    n/np = " << 4.0*M_PI*I2<<  ",    ZLOG = " << ZLn_ei <<",     vw = "<< pow(vw_cube,1.0/3.0) << "\n";

    qn_c = Qn_coeff*vos*vos*ZLn_ei*4.0*M_PI*I2; // where 4.0*M_PI*I2 is ne/np

    fh[0] = qn_c * Qn[0] * (  Pn[0]   * gn[0]   * (fin[1]-fin[0])
                              - Pnm1 /(1.0 + Inv_Uav6_nm1*vw_cube*vw_cube) * (fin[0]-f00));
    /*debug*/ // cout << "fh["<< vr[0]<<"] = " << fh[0]<< "\n";
    for (size_t i(1); i < fh.size()-1; ++i) {
        fh[i] = qn_c * Qn[i] * (  Pn[i]   * gn[i]   * (fin[i+1]-fin[i])
                                  - Pn[i-1] * gn[i-1] * (fin[i]-fin[i-1]));
        /*debug*/ // cout << "fh["<< vr[i]<<"] = " << fh[i]<< "\n"; 
    }
    return fh;
}
//-------------------------------------------------------------------
//*******************************************************************


//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//*******************************************************************
//*******************************************************************
//  Runge Kutta loop for the electron-electron collisions
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
RK4_IB::RK4_IB(valarray<double>& fin, double pmax, int tout_start)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        :t(static_cast<double>(tout_start) *
           Input::List().t_stop / Input::List().n_outsteps),  // Initialize time = 0
         Tout(t + Input::List().t_stop / Input::List().n_outsteps),                                   // The next output time
         fh(0.0,fin.size()),
         f0(0.0,fin.size()),
         f1(0.0,fin.size()),    // 3 local "State" Variables
         f(fin),                                 // Initialize a reference to "State" Y
         Inversebremsstrahlung(fh, pmax){              // Initialize "Actions"

}
//--------------------------------------------------------------

//  Output time
double& RK4_IB::tout() {return Tout;}

//  Output time
size_t& RK4_IB::numh() {return num_h;}

//  Output time
double& RK4_IB::th() {return h;}

//  Real time of the simulation (can only be modified in the RK)
double& RK4_IB::time()  {return t;}

//  Call Advection Actions 
valarray<double>& RK4_IB::F(const valarray<double>& fin, const double vosc, const double Zval) {
    return Inversebremsstrahlung.Getslope(fin,vosc,Zval);
}

//--------------------------------------------------------------
RK4_IB& RK4_IB::advance(const double vosc, const double Zval){
//--------------------------------------------------------------
//  Take a step using RK4
//--------------------------------------------------------------

//      Initialization
    f0 = f; f1 = f;

//      Step 1
    F(f1,vosc,Zval);                          // fh = F(f1)
    fh *= (0.5*h);   f1 += fh;      // f1 = f1 + (h/2)*fh
    fh *= (1.0/3.0); f  += fh;      // f  = f  + (h/6)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
    F(f1,vosc,Zval);     f1  = f0;            // fh = F(f1)
    fh *= (0.5*h);   f1 += fh;      // f1 = f0 + (h/2)*fh
    fh *= (2.0/3.0); f  += fh;      // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
    F(f1,vosc,Zval);                          // fh = F(f1)
    fh *= h;          f0 += fh;     // f1 = f0 + h*Yh
    fh *= (1.0/3.0);  f  += fh;     // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 4
    F(f0,vosc,Zval);                          // fh = F(f0)
    fh *= (h/6.0);    f += fh;      // f  = f  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    t += h;

    return *this;

}
//-------------------------------------------------------------------
//*******************************************************************





//**************************************************************
//**************************************************************
// Inverse Brehmsstrahlung heating
//**************************************************************
//**************************************************************
//--------------------------------------------------------------
// InverseBremsstrahlung::InverseBremsstrahlung(DistFunc1D& DFin, int tout_start, const valarray<double>& grid)
InverseBremsstrahlung::InverseBremsstrahlung(double pmax, size_t nump, size_t numx, int tout_start,const valarray<double>& grid) //DistFunc1D& DFin,
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        :t(static_cast<double>(tout_start) *
           Input::List().t_stop / Input::List().n_outsteps),         // Initialize time = tout_start * dt_out
         tout(t + Input::List().t_stop / Input::List().n_outsteps),  // Initialize tout = time + dt_out 
         fc(0.0, nump),
         rk4_ib(fc, pmax, tout_start),
         vr(Algorithms::MakeCAxis(0.0,pmax,nump)),
//        vr(Algorithms::MakeAxis(pmax/(2.0*nump-1.0),pmax,nump)),
//       ------------------------
         InvBremsstrahlung(Input::List().inverse_bremsstrahlung),
//       ------------------------
         omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Input::List().lambda_0)),
         omega_p(5.64 * 1.0e+4*sqrt(Input::List().density_np)),
         w0overwp(omega_0/omega_p),
//       -----
//       -----
         U2(0.0, nump),
         U2m1(0.0, nump),
//       -----
         vos(Input::List().lambda_0 * sqrt(7.3e-19*Input::List().I_0)),
         IL_xprofile(numx)

{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       The intensity profile is: I0(t)*e^(-2*(r/w0)^2) and absorption takes place over
//       depth L so that I0_r(t) ~ I0(t)*e^(-((x-x0)/L)^2). The idea is that if you 
//       integrate over x you get I0(t)

    // using Inputdata::IN;
    Nbc = Input::List().BoundaryCells;
    szx = Input::List().NxLocal[0];      // size of useful x axis

    for (size_t i(1); i < U2.size(); ++i) {
        U2[i]   = 0.5 * vr[i]*vr[i]          * (vr[i]-vr[i-1]);
        U2m1[i] = 0.5 * vr[i-1]*vr[i-1]      * (vr[i]-vr[i-1]);
    }

    // Construct the intensity profile map
    Setup_Y::parseprofile(grid, Input::List().intensity_profile_str, IL_xprofile);
    if ( abs(InvBremsstrahlung) > 1 ) {
        InvBremsstrahlung = 0;
    }

}
//--------------------------------------------------------------

//--------------------------------------------------------------
int InverseBremsstrahlung:: IBSOURCE()  const {return InvBremsstrahlung;}
//--------------------------------------------------------------

//--------------------------------------------------------------
void InverseBremsstrahlung::loop(SHarmonic1D& f0, valarray<double>& Zarray, const double& tnew){
//--------------------------------------------------------------
//  Calculate the effect of the laser source
//--------------------------------------------------------------

    //cout << "omega_0 = " << omega_0 <<"\n";
    //cout << "omega_p = " << omega_p <<"\n";
    //cout << "np/nc = " << pow(omega_p / omega_0,2) <<"\n";
    // cout << "vos = " << vos <<"\n";
    // cout << "polarization = " << POL() <<"\n";

    // double dt = tnew -t;
    // double Envelope_Value = t_Profile(tnew)/dt; // This is Envelope(t)*dt/dt = Envelope(t)
    /*debug*/ // cout << "vos("<<tnew<<") = " << Envelope_Value << "\n";

    // debug: cout << "Envelope = " << Envelope_Value << /*",  Carrier Signal = " << Carrier_Signal(tnew) <<*/"\n";

    /*switch (POL()) {
        case 1: {
           for (size_t i(0); i < Y.SH(0,0).numx(); ++i){
                for (size_t j(0); j < Y.SH(0,0).numy(); ++j){
                    LaserFields.Ex()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j]
                            * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                }
            }
            break;
        }
        case 2: {
            for (size_t i(0); i < Y.SH(0,0).numx(); ++i){
                for (size_t j(0); j < Y.SH(0,0).numy(); ++j){
                    LaserFields.Ey()(i,j) = vos  * IL_xprofile[i] * IL_yprofile[j]
                             * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                }
            }
            break;
        }
        case 3: {
            for (size_t i(0); i < Y.SH(0,0).numx(); ++i){
                for (size_t j(0); j < Y.SH(0,0).numy(); ++j){
                    LaserFields.Ez()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j]
                            * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                }
            }
            break;
        }
        default:
            break;
    }*/

    // Array2D< double > Vos(Y.SH(0,0).numx(),Y.SH(0,0).numy());
    valarray< double > Vos(0.0,f0.numx());//,Y.SH(0,0).numy());

    // for (size_t iy(0); iy < Y.SH(0,0).numy(); ++iy){
    for (size_t ix(0); ix < f0.numx(); ++ix){
        Vos[ix] = vos * sqrt(IL_xprofile[ix]);// * IL_yprofile[iy]);
    }
    // }

//      Initialization
    rk4_ib.tout()  = tnew;
    rk4_ib.numh()  = static_cast<size_t>(static_cast<int>((tnew-t)/(Input::List().smaller_dt)))+1;
    rk4_ib.th()    = (tnew-t)/static_cast<double>(rk4_ib.numh());


    // for (size_t iy(0); iy < szy-2*Nbc; ++iy){
    for (size_t ix(0); ix < szx-2*Nbc; ++ix){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        // Copy data for a specific location in space to valarray
        for (size_t ip(0); ip < fc.size(); ++ip){
            // fc[ip] = (Y.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
            fc[ip] = f0(ip,ix+Nbc).real();
        }
        double I2_before(0.0);
        for (int n(1); n < U2.size(); ++n) {
            I2_before+= U2[n]*fc[n]+U2m1[n]*fc[n-1];
        }

        // Time loop: Update the valarray
        rk4_ib.time() = t;
        for (size_t h_step(0); h_step < rk4_ib.numh(); ++h_step){
            rk4_ib.advance(Vos[ix+Nbc], Zarray[ix+Nbc]);//,iy+Nbc));
            /*debug*/ // cout << "Vos("<<ix<<","<<iy<<") = "<< Vos(ix+Nbc,iy+Nbc) << "\n";
        }
        double I2_after(0.0);
        for (int n(1); n < U2.size(); ++n) {
            I2_after += U2[n]*fc[n]+U2m1[n]*fc[n-1];
        }
        fc *= I2_before/I2_after;
        // Return updated data to the harmonic
        for (int ip(0); ip < fc.size(); ++ip){
            // Y.SH(0,0)(ip,ix+Nbc,iy+Nbc) = fc[ip];
            f0(ip,ix+Nbc) = fc[ip];
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    }
    // }
    t = tnew;

}
//--------------------------------------------------------------

//--------------------------------------------------------------
//     double InverseBremsstrahlung::t_Profile(const double& t2){ 
// //--------------------------------------------------------------
// //  Calculate the effect of the laser profile
// //--------------------------------------------------------------
//         double t;

//         if ( t2 < rise_time ) {
//             t = t2 / rise_time;
//             if ( polynomial_Ivst ) {
//                 return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t);
//             }
//             if ( linear_Ivst ) {
//                 return t * (t2-t);
//             }
//             return 0;
//         } else {
//             if ( t2 < rise_time + flat_time) { 
//                 return t2-t;
//             } else {
//                 if ( t2 < rise_time + flat_time + fall_time) { 
//                     t = (rise_time + flat_time + fall_time - t2) / fall_time;
//                     if ( polynomial_Ivst ) {
//                         return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t);
//                     }
//                     if ( linear_Ivst ) {
//                         return t *(t2-t);
//                     }
//                     return 0;
//                  } 
//                  else return -1.0;
//              } 
//          }
//     }
//--------------------------------------------------------------
//**************************************************************


