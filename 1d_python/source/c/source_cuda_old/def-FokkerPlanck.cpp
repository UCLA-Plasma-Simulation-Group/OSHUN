///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	    June 7th 201
//	Last Modified:	Dec  2nd 2011
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the 
//   Backward Euler method for the collisions
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


// Standard libraries 


#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <math.h>
#include <stdio.h>
#include <float.h>

// My libraries 
#include "matrices.h"
#include "Eigen/Dense"
using namespace Eigen;

// Interface 
#include "decl-nmethods.h"
#include "decl-input.h"
#include "decl-state.h"
#include "decl-fokkerplanck.h"


//extern int cheapo_debugging;
//extern int debug_flag;

 int cheapo_debugging = 0;
 int debug_flag = 0;
#
//#define DB_STUB 1
#

//*******************************************************************
//  Coulomb Logarithms from NRL Plasma Formulary 2009, p 34
//*******************************************************************

//*******************************************************************
//-------------------------------------------------------------------
     double LOGee(double ne, double Te) {
//-------------------------------------------------------------------
//   Calculate the Coulomb logarithm for electron-electron collisions
//-------------------------------------------------------------------
//      Note: that the results here assume the distribution functions
//      are nonrelativistic, as is the case for the rest of the F-P 
//      part of the code.

        static double lnee;
        
//      if the density is positive
        if (ne > 0.000000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= Inputdata::IN().list().density_np;

            Te = log(Te); 
            ne = log(ne);

            lnee = 23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0));

            if (lnee > 2.0) return lnee;
        }
        // Default minimum "2"

        return 2.0; 
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
     double ZLOGei(double ne, double Te) {
//-------------------------------------------------------------------
//   Calculate the Coulomb logarithm * Zeta for electron-ion collisions
//-------------------------------------------------------------------
//      Note: that the results here assume the distribution functions
//      are nonrelativistic, as is the case for the rest of the F-P 
//      part of the code.

        static double lnei;

//      if the density is positive
        if (ne > 0.0000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= Inputdata::IN().list().density_np;

            Te = log(Te); 
            ne = log(ne);

            lnei = 24.0 - 0.5*ne + Te;
            if (lnei > 2.0) return lnei * Inputdata::IN().list().Zeta;
        }
 
        return  2.0 * Inputdata::IN().list().Zeta;
    }
//-------------------------------------------------------------------
//*******************************************************************

//*******************************************************************
//*******************************************************************
// Tony Bell's Energy & Number Conserving Algorithm
//*******************************************************************
//*******************************************************************

//*******************************************************************
//--------------------------------------------------------------
    EE_isotropic_conserva::EE_isotropic_conserva(valarray<double>& fslope)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : fh(fslope),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Velocity
         vr(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Constants for Integrals
         U4(0.0, Inputdata::IN().list().nump), 
         U4m1(0.0, Inputdata::IN().list().nump), 
         U2(0.0, Inputdata::IN().list().nump), 
         U2m1(0.0, Inputdata::IN().list().nump), 
         U1(0.0, Inputdata::IN().list().nump), 
         U1m1(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Integrals
         J1(0.0, Inputdata::IN().list().nump), 
         I2(0.0, Inputdata::IN().list().nump), 
         I4(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         U3(0.0, Inputdata::IN().list().nump), 
         Qn(0.0, Inputdata::IN().list().nump), 
         Pn(0.0, Inputdata::IN().list().nump)
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         //classical electron radius
         double re(2.8179402894e-13);           
         double kp(sqrt(4.0*M_PI*(Inputdata::IN().list().density_np)*re));
         c_kpre = 4.0*M_PI/3.0*re*kp;
         NB = Inputdata::IN().list().NB_algorithms;

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         // Determine vr
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             vr[i] = Inputdata::IN().inp().pr(i);
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

//-------------------------------------------------------------------
    double EE_isotropic_conserva::G(const int& n, const valarray<double>& fin) {
//-------------------------------------------------------------------
        double i2s, i4s;
        double f00( (fin[0] - fin[1]*(vr[0]*vr[0])/(vr[1]*vr[1]))/ 
                       (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1])) );
        
        i2s = f00*pow(vr[n],3)/3.0 + (fin[1]-f00)*pow(vr[n],5)/(vr[1]*vr[1])*0.2;
        i4s = f00*pow(vr[n],5)*0.2 + (fin[1]-f00)*pow(vr[n],7)/(vr[1]*vr[1]*7.0);

        return fin[n]*i4s + (pow(vr[n],3)*fin[n]-3.0*i2s) * J1[n];
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    valarray<double>& EE_isotropic_conserva::Getslope(const valarray<double>& fin) {
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

       Ln_ee = LOGee(4.0*M_PI*I2[I2.size()-1],4.0*M_PI*I4[I4.size()-1]);

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

       return fh;
    }
//-------------------------------------------------------------------
//*******************************************************************


//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//*******************************************************************
//*******************************************************************
// Michail Tzoufras' Non-Conserving Algorithm
// ------------------------------------------
// This algorithm is similar to the one by Mcdonald & Rosenbluth
// and it is ineffective, it is not currently in use. This algorithm
// may be effective for very few collision times, but over long
// time scales it cannot be expected to work well.
//*******************************************************************
//*******************************************************************

//*******************************************************************
//--------------------------------------------------------------
    EE_non_const::EE_non_const(valarray<double>& fslope)        // NOT IN USE
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : fh(fslope),
         df(fh),
         ddf(fh),
//       Pre-Calculated Constants
         vr(0.0, Inputdata::IN().list().nump), 
         Inv_Dvh(0.0, Inputdata::IN().list().nump), 
         Inv_Dvc(0.0, Inputdata::IN().list().nump), 
         Inv_v(0.0, Inputdata::IN().list().nump), 
         Inv_v2(0.0, Inputdata::IN().list().nump), 
         Inv_v4(0.0, Inputdata::IN().list().nump), 
//       Constants for Integrals
         U4(0.0, Inputdata::IN().list().nump), 
         U4m1(0.0, Inputdata::IN().list().nump), 
         U2(0.0, Inputdata::IN().list().nump), 
         U2m1(0.0, Inputdata::IN().list().nump), 
         U1(0.0, Inputdata::IN().list().nump), 
         U1m1(0.0, Inputdata::IN().list().nump), 
//       Integrals
         I1(0.0, Inputdata::IN().list().nump), 
         I2(0.0, Inputdata::IN().list().nump), 
         I4(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         U3(0.0, Inputdata::IN().list().nump), 
         Qn(0.0, Inputdata::IN().list().nump), 
         Pn(0.0, Inputdata::IN().list().nump)
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
         double re(2.8179402894e-13);           //classical electron radius
         double kp(sqrt(4.0*M_PI*(Inputdata::IN().list().density_np)*re));
         c_kpre = 4.0*M_PI/3.0*re*kp;

         NB = Inputdata::IN().list().NB_algorithms;
         if (NB < 2) NB = 1;
         if (NB > vr.size()-1) NB = vr.size()-1;

         // Determine vr
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             vr[i] = Inputdata::IN().inp().pr(i);
             vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));         
         }

         f0_c1 = 1.0 / (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]));
         f0_c2 = (vr[0]*vr[0])/(vr[1]*vr[1])*f0_c1;
         fpp0_c  = 2.0/(vr[1]*vr[1]);

         // Determine v^(-1), v^(-2), v^(-4)
         for (size_t i(0); i < U4.size(); ++i) {
             Inv_v[i]  = 1.0 / vr[i];
             Inv_v2[i] = Inv_v[i]  * Inv_v[i];
             Inv_v4[i] = Inv_v2[i] * Inv_v2[i];
         }

         // Determine Dv_n+1/2 and put it at n
         for (size_t i(0); i < U4.size(); ++i) {
             Inv_Dvh[i] = 1.0 / (vr[i+1]-vr[i]);
         }

         // Determine 2*Dv_n
         for (size_t i(1); i < U4.size()-1; ++i) {
             Inv_Dvc[i] = 1.0 / (vr[i+1]-vr[i-1]);
         }

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
         for (size_t i(1); i < Qn.size()-1; ++i) {
             Qn[i] = 1.0 / (vr[i]*vr[i]*(vr[i+1]-vr[i-1])/2.0);
         }        
         // Determine Pn 
         for (size_t i(0); i < Pn.size()-1; ++i) {
             Pn[i] = 1.0 / ((vr[i+1]-vr[i])/2.0*(vr[i+1]+vr[i]));
         }    

    }
//--------------------------------------------------------------


//-------------------------------------------------------------------
    valarray<double>&  EE_non_const::Getslope(const valarray<double>& fin) { // NOT IN USE
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------
       double Ln_ee;

//     Evaluate the integrals
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

//     Use I2 and I4 to calculate the Coulomb logarithm
       Ln_ee = LOGee(4.0*M_PI*I2[I2.size()-1],4.0*M_PI*I4[I4.size()-1]);

       I1[I1.size()-1] = 0;
       for (int n(I1.size()-2); n > -1; --n) {
           I1[n]  = U1[n+1]*fin[n+1]+U1m1[n+1]*fin[n]; 
           I1[n] += I1[n+1];
       }
 
//     Michail's Algorithm
//     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//     Evaluate the derivative
       for (int n(1); n < NB/*I4.size()-1*/; ++n) {
           df[n]  = fin[n+1]-fin[n-1];
           df[n] *= Inv_Dvc[n];
       }

//     Evaluate the second derivative
//        -df/dv
       for (int n(1); n < NB+1/*I4.size()-1*/; ++n) {
           ddf[n]  = fin[n-1]-fin[n];
           ddf[n] *= Inv_Dvh[n-1];
       }
//        D(df/dv)/Dv
       for (int n(1); n < NB/*I4.size()-1*/; ++n) {
           ddf[n] -= ddf[n+1];
           ddf[n]  *= 2.0*Inv_Dvc[n];
       }

//     fh from 1 to NB...
       for (int n(1); n < NB/*I4.size()-1*/; ++n) {
           fh[n]  = I4[n] * (Inv_v4[n]*(vr[n]*ddf[n]-df[n]));
           fh[n] += I2[n] * (3.0*Inv_v2[n]*df[n]);
           fh[n] += I1[n] * (2.0*Inv_v[n]*df[n]+ddf[n]);
           fh[n] += 3.0*fin[n]*fin[n];
       }

//     Calculate zeroth cell
       f00 = f0_c1*fin[0]-f0_c2*fin[1];
       fpp00 = fpp0_c * (fin[1] - f00);
       fh[0] = 3.0 * (f00*f00+fpp00*I1[0]);
//     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//     Tony's Energy Conserving Algorithm from NB to maximum
//     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
       for (int n(NB-1); n < I4.size(); ++n) {
           I2[n] *= I1[n];       // I1(n) * I2(n)
           I2[n] *= -3.0;        // -3 * I1(n) * I2(n)
           I1[n] *= U3[n];       // u_n^3 * I1(n)
           I4[n] += I1[n];       // I4(n) + u_n^3 * I1(n) 
           I4[n] *= fin[n];      // fn * I4(n) + u_n^3 * fn * I1(n)
           I4[n] += I2[n];       // Gn = fn * I4(n) + u_n^3 * fn * I1(n) - 3 * I1(n) * I2(n)
       }

       for (int n(NB-1); n < I4.size()-1; ++n) I4[n] -= I4[n+1];
       for (int n(NB-1); n < I4.size()-1; ++n) I4[n] *= Pn[n];
       for (int n(NB-1); n < I4.size()-2; ++n) I4[n] -= I4[n+1];
       for (int n(NB-1); n < I4.size()-2; ++n) fh[n+1] = Qn[n+1]*I4[n];

//     Calculate maximum cell
       fh[fh.size()-1] = 0.0;
//     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//     Normalize
       fh *=  c_kpre *  Ln_ee;

       return fh;
    }
//-------------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><



//*******************************************************************
//*******************************************************************
//  Runge Kutta loop for the electron-electron collisions
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
    RungeKutta4_ee::RungeKutta4_ee(valarray<double>& fin, int tout_start)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),  // Initialize time = 0
        Tout(t + Inputdata::IN().cont().dt_out),                                   // The next output time
        fh(0.0,Inputdata::IN().list().nump), 
        f0(0.0,Inputdata::IN().list().nump), 
        f1(0.0,Inputdata::IN().list().nump),    // 3 local "State" Variables
        f(fin),                                 // Initialize a reference to "State" Y
        Collide(fh){                            // Initialize "Actions"
    
    }
//--------------------------------------------------------------

//  Output time
    double& RungeKutta4_ee::tout() {return Tout;}

//  Output time
    size_t& RungeKutta4_ee::numh() {return num_h;}

//  Output time
    double& RungeKutta4_ee::th() {return h;}

//  Real time of the simulation (can only be modified in the RK)
    double& RungeKutta4_ee::time()  {return t;}

//  Call Advection Actions 
    valarray<double>& RungeKutta4_ee::F(const valarray<double>& fin) { 
        return Collide.Getslope(fin);
    }

//--------------------------------------------------------------
    RungeKutta4_ee& RungeKutta4_ee::advance(){    
//--------------------------------------------------------------
//  Take a step using RK4
//--------------------------------------------------------------
		//ScopeTimer("e*e collsion step");
		ScopeTimer ___t("e*e collsion step");
//      Initialization
        f0 = f; f1 = f; 

//      Step 1
        F(f1);                          // fh = F(f1)
        fh *= (0.5*h);   f1 += fh;      // f1 = f1 + (h/2)*fh
        fh *= (1.0/3.0); f  += fh;      // f  = f  + (h/6)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        F(f1);     f1  = f0;            // fh = F(f1)
        fh *= (0.5*h);   f1 += fh;      // f1 = f0 + (h/2)*fh
        fh *= (2.0/3.0); f  += fh;      // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        F(f1);                          // fh = F(f1)
        fh *= h;          f0 += fh;     // f1 = f0 + h*Yh
        fh *= (1.0/3.0);  f  += fh;     // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
//      Step 4
        F(f0);                          // fh = F(f0)
        fh *= (h/6.0);    f += fh;      // f  = f  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        t += h;
       
        return *this;
    
    }
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//*******************************************************************
//   Definition for the Explicit Integration method for the 
//   energy relaxation of the distribution function
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
    Explicit_EC_FP::Explicit_EC_FP(Stat& Yin, int tout_start)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),         // Initialize time = tout_start * dt_out
        tout(t + Inputdata::IN().cont().dt_out),  // Initialize tout = time + dt_out 
        fc(0.0, Inputdata::IN().list().nump),
        rk4_ee(fc, tout_start),
        Y(Yin) {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        using Inputdata::IN; 
        Nbc = IN().list().RKLevel;
        szx = IN().inp().x.dim();      // size of useful x axis 
        szy = IN().inp().y.dim();      // size of useful y axis

    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
    void Explicit_EC_FP::loop(const double& tnew){    
//-------------------------------------------------------------------
//  This loop scans all the locations in configuration space
//  and calls the RK4 for explicit integration at each location
//-------------------------------------------------------------------

//      Initialization
        rk4_ee.tout()  = tnew;
        // Adam Version of this... rk4_ee.numh()  = static_cast<size_t>(static_cast<int>((tnew-t)/(Inputdata::IN().list().small_dt)-.5))+1;
        //rk4_ee.numh()  = static_cast<size_t>(static_cast<int>((tnew-t)/(Inputdata::IN().list().small_dt)-.5))+1; 
		rk4_ee.numh()  = static_cast<size_t>(static_cast<int>((tnew-t)/(Inputdata::IN().list().small_dt)))+1; 
        rk4_ee.th()    = (tnew-t)/static_cast<double>(rk4_ee.numh());



        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){
				//ScopeTimer("e*e collsion per cell");
				ScopeTimer ___t("e*e collsion per cell");
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                // Copy data for a specific location in space to valarray
                for (size_t ip(0); ip < fc.size(); ++ip){
                    fc[ip] = (Y.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
                }
                // Time loop: Update the valarray
                rk4_ee.time() = t;
                for (size_t h_step(0); h_step < rk4_ee.numh(); ++h_step){
                    rk4_ee.advance(); 
                }
                
                // Return updated data to the harmonic
                for (int ip(0); ip < fc.size(); ++ip){
                    Y.SH(0,0)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            }
        }
        t = tnew;

    }
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//*******************************************************************
//  Definition of scattering with high order harmonics
//*******************************************************************
//*******************************************************************


//*******************************************************************
//--------------------------------------------------------------
    Implicit_SC_step::Implicit_SC_step()
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : 
//       Pre-Calculated Constants
         vr(0.0,  Inputdata::IN().list().nump), 
//       Constants for Integrals
         U4(0.0,  Inputdata::IN().list().nump), 
         U4m1(0.0,Inputdata::IN().list().nump), 
         U2(0.0,  Inputdata::IN().list().nump), 
         U2m1(0.0,Inputdata::IN().list().nump), 
         U1(0.0,  Inputdata::IN().list().nump), 
         U1m1(0.0,Inputdata::IN().list().nump), 
//       Integrals
         J1m(0.0,  Inputdata::IN().list().nump), 
         I2(0.0,  Inputdata::IN().list().nump), 
         I4(0.0,  Inputdata::IN().list().nump), 
         ss1(0.0, Inputdata::IN().list().nump)
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
         double re(2.8179402894e-13);           //classical electron radius
         double kp(sqrt(4.0*M_PI*(Inputdata::IN().list().density_np)*re));
         c_kpre = 4.0*M_PI*re*kp/3.0;

         // Determine vr
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             vr[i] = Inputdata::IN().inp().pr(i);
             vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));         
         }

         // Determine U4, U4m1, U2, U2m1, U1, U1m1 for the integrals
         for (size_t i(1); i < U4.size(); ++i) {
             U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);         
             U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);         
             U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);         
             U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);         
             U1[i]   = 0.5 * vr[i]            * (vr[i]-vr[i-1]);         
             U1m1[i] = 0.5 * vr[i-1]          * (vr[i]-vr[i-1]);         
         }
    }
//--------------------------------------------------------------

//-------------------------------------------------------------------
    void  Implicit_SC_step::reset_coeff(const valarray<double>& fin, const double Dt) {
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------

//     Evaluate the integrals
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

       J1m[J1m.size()-1] = 0;
       for (int n(J1m.size()-2); n > -1; --n) {
           J1m[n]  = U1[n+1]*fin[n+1]+U1m1[n+1]*fin[n]; 
           J1m[n] += J1m[n+1];
       }

//     Electron-electron scattering
       for (int i(0); i < ss1.size(); ++i){
           ss1[i] = 3.0*I2[i] - I4[i]/(vr[i]*vr[i]) + 2.0*vr[i]*J1m[i] ;
       }

//     Multiply with lnL_ee
       ss1 *= LOGee(4.0*M_PI*I2[I2.size()-1],4.0*M_PI*I4[I4.size()-1]);

//     Ion-electron scattering including Z*lnL_ei
       ss1 += ZLOGei(4.0*M_PI*I2[I2.size()-1],4.0*M_PI*I4[I4.size()-1])*3.0*I2[I2.size()-1];

//     multiply by 1/v^3
       for (int i(0); i < ss1.size(); ++i){
           ss1[i] /= pow(vr[i],3);
       }

       ss1 *= (-1.0)* c_kpre;
       ss1 *= (-1.0)*Dt;

    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    void  Implicit_SC_step::advance(valarray< complex<double> >& fin, const int el) {
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------

        double ll1(static_cast<double>(el));
        ll1 *= 0.5*(ll1 + 1.0);

        for (int i(0); i < fin.size(); ++i){
            fin[i] /= 1.0 + ll1 * ss1[i]; 
        }

    }
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//*******************************************************************
//  Definition of collisions with high order harmonics FLM
//*******************************************************************
//*******************************************************************


//*******************************************************************
//--------------------------------------------------------------
    Implicit_FLM_step::Implicit_FLM_step()
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : 
//       Pre-Calculated Constants
         vr(0.0,  Inputdata::IN().list().nump), 
//       Constants for Integrals
         U4(0.0,  Inputdata::IN().list().nump), 
         U4m1(0.0,Inputdata::IN().list().nump), 
         U2(0.0,  Inputdata::IN().list().nump), 
         U2m1(0.0,Inputdata::IN().list().nump), 
         U1(0.0,  Inputdata::IN().list().nump), 
         U1m1(0.0,Inputdata::IN().list().nump), 
//       Integrals
         J1m(0.0,  Inputdata::IN().list().nump), 
         I0(0.0,  Inputdata::IN().list().nump), 
         I2(0.0,  Inputdata::IN().list().nump), 
         Scattering_Term(0.0, Inputdata::IN().list().nump),
         df0(0.0,  Inputdata::IN().list().nump), 
         ddf0(0.0,  Inputdata::IN().list().nump), 
//       Matrices
         Alpha_Tri(Inputdata::IN().list().nump, Inputdata::IN().list().nump),
         if_tridiagonal(Inputdata::IN().list().if_tridiagonal)
//       Make matrix vk/vn
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
        // this is a NOP if the size is the same.
        _implicitCollsionMatrix.resize(Alpha_Tri.dim1(), Alpha_Tri.dim1());
        
         
         double re(2.8179402894e-13);           //classical electron radius
         double kp(sqrt(4.0*M_PI*(Inputdata::IN().list().density_np)*re));
         kpre = re*kp;

         // Determine vr
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             vr[i] = Inputdata::IN().inp().pr(i);
             vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));         
         }

         // Determine U4, U4m1, U2, U2m1, U1, U1m1 for the integrals
         for (size_t i(1); i < U4.size(); ++i) {
             U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);         
             U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);         
             U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);         
             U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);         
             U1[i]   = 0.5 * vr[i]            * (vr[i]-vr[i-1]);         
             U1m1[i] = 0.5 * vr[i-1]          * (vr[i]-vr[i-1]);         
         }
    }
//--------------------------------------------------------------

//-------------------------------------------------------------------
    void  Implicit_FLM_step::reset_coeff(const valarray<double>& fin, const double Delta_t) {
//-------------------------------------------------------------------
//  Reset the coefficients based on f_0^0 
//-------------------------------------------------------------------

//     Calculate Dt
       Dt = Delta_t;

//     INTEGRALS
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     Temperature integral I2 = 4*pi / (v^2) * int_0^v f(u)*u^4du
       I2[0] = 0;
       for (int k(1); k < I2.size(); ++k) {
           I2[k]  = U4[k]*fin[k]+U4m1[k]*fin[k-1]; 
           I2[k] += I2[k-1];
       }
       I2 *= 4.0 * M_PI; 
       I2_temperature = I2[I2.size()-1];
       for (int k(0); k < I2.size(); ++k){
           I2[k] /=  vr[k]*vr[k];
       }

//     Density integral I0 = 4*pi*int_0^v f(u)*u^2du 
       I0[0] = 0;
       for (int k(1); k < I0.size(); ++k) {
           I0[k]  = U2[k]*fin[k]+U2m1[k]*fin[k-1]; 
           I0[k] += I0[k-1];
       }
       I0 *= 4.0 * M_PI;
       I0_density = I0[I0.size()-1];

//     Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
       J1m[J1m.size()-1] = 0;
       for (int k(J1m.size()-2); k > -1; --k) {
           J1m[k]  = U1[k+1]*fin[k+1]+U1m1[k+1]*fin[k]; 
           J1m[k] += J1m[k+1];
       }
       for (int k(0); k < J1m.size(); ++k){
           J1m[k] *= 4.0 * M_PI *vr[k];
       }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


//     COULOMB LOGARITHMS
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       _LOGee  = LOGee(I0_density,I2_temperature);
       _ZLOGei = ZLOGei(I0_density,I2_temperature);
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
       Alpha_Tri = 0.0; _implicitCollsionMatrix.fill(0.0);
       double IvDnDm1, IvDnDp1, Ivsq2Dn;
       Alpha_Tri(0,0) = 8.0 * M_PI * fin[0];                                         //  8*pi*f0[0]
       _implicitCollsionMatrix(0,0) = 8.0 * M_PI * fin[0]; 
       
       
       for (int i(1); i < TriI1.size()-1; ++i){
          IvDnDm1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]  -vr[i-1]));       //  ( v*D_n*D_{n-1/2} )^(-1)
          IvDnDp1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1]-vr[i]  ));       //  ( v*D_n*D_{n+1/2} )^(-1)
          Ivsq2Dn = 1.0 / (vr[i] * vr[i]                 * (vr[i+1]-vr[i-1]));       //  ( v^2 * 2*D_n )^(-1)
          Alpha_Tri(i, i  ) = 8.0 * M_PI * fin[i] - TriI2[i] * (IvDnDm1 + IvDnDp1);
          Alpha_Tri(i, i-1) = TriI2[i] * IvDnDm1 - TriI1[i] * Ivsq2Dn;   
          Alpha_Tri(i, i+1) = TriI2[i] * IvDnDp1 + TriI1[i] * Ivsq2Dn;

          _implicitCollsionMatrix(i, i  ) = 8.0 * M_PI * fin[i] - TriI2[i] * (IvDnDm1 + IvDnDp1);
          _implicitCollsionMatrix(i, i-1) = TriI2[i] * IvDnDm1 - TriI1[i] * Ivsq2Dn;   
          _implicitCollsionMatrix(i, i+1) = TriI2[i] * IvDnDp1 + TriI1[i] * Ivsq2Dn;
	
       }
       Alpha_Tri *=  (-1.0) * _LOGee * kpre * Dt;         // (-1) because the matrix moves to the LHS in the equation
       _implicitCollsionMatrix *= (-1.0) * _LOGee * kpre * Dt; 
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//     CALCULATE DERIVATIVES 
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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

    //   for (int k(0); k < fin.size(); ++k) {
    //       cout <<  vr[k] << "           " << fin[k] <<"           " <<  df0[k] <<"           " <<  ddf0[k] << "\n";
    //   }
    //   exit(1);

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    void  Implicit_FLM_step::advance(valarray< complex<double> >& fin, const int el) {
//-------------------------------------------------------------------
//  Collisions
//-------------------------------------------------------------------
        Matrix2D<double> Alpha(Alpha_Tri);
        valarray< complex<double> > fout(fin);


//      ZEROTH CELL FOR TRIDIAGONAL ARRAY
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        if (el > 1) {
            Alpha(0,0) = 0.0;
			_implicitCollsionMatrix(0,0) = 0.0;
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
               _implicitCollsionMatrix(i,0) += t1 * ( 2.0*M_PI*pow(vr[0]/vr[i],el+2)*vr[0]*vr[0]*(vr[1]-vr[0]) ); 
               _implicitCollsionMatrix(i,0) += t3 * ( 2.0*M_PI*pow(vr[0]/vr[i],el)  *vr[0]*vr[0]*(vr[1]-vr[0]) ); 

               for (int j(1); j < i; ++j){
                   Alpha(i,j) += t1 * ( 2.0*M_PI*pow(vr[j]/vr[i],el+2)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   Alpha(i,j) += t3 * ( 2.0*M_PI*pow(vr[j]/vr[i],el)  *vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   _implicitCollsionMatrix(i,j) += t1 * ( 2.0*M_PI*pow(vr[j]/vr[i],el+2)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   _implicitCollsionMatrix(i,j) += t3 * ( 2.0*M_PI*pow(vr[j]/vr[i],el)  *vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
               }

               Alpha(i,i) += t1 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) ); 
               Alpha(i,i) += t3 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) ); 
               _implicitCollsionMatrix(i,i) += t1 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) ); 
               _implicitCollsionMatrix(i,i) += t3 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i]-vr[i-1]) ); 
               
               Alpha(i,i) += t2 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) ); 
               Alpha(i,i) += t4 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) ); 
               _implicitCollsionMatrix(i,i) += t2 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) ); 
               _implicitCollsionMatrix(i,i) += t4 * ( 2.0*M_PI *vr[i]*vr[i]*(vr[i+1]-vr[i]) ); 

               for (int j(i+1); j < Alpha.dim2()-1; ++j){
                   Alpha(i,j) += t2 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el-1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   Alpha(i,j) += t4 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el+1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   _implicitCollsionMatrix(i,j) += t2 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el-1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
                   _implicitCollsionMatrix(i,j) += t4 * ( 2.0*M_PI*pow(vr[j]/vr[i],-el+1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) ); 
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
			_implicitCollsionMatrix(i,i) += 1.0 - ll1 * Scattering_Term[i];
        }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		cout << "inverting impcit matix" << endl;
		//cout << _implicitCollsionMatrix.inverse();
		cout << "done inverting impcit matix" << endl;

//      SOLVE A * Fout  = Fin
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        //for (int i(0); i < fin.size(); ++i){   // Estimate Fout = Fin /  A(i,i)
        //    fout[i] /= Alpha(i,i); 
        //}

        // cout << "l = " << el << ",        ";       
        if ( if_tridiagonal ) { 
            if ( !(Thomas_Tridiagonal(Alpha, fin, fout)) ) {  // Invert A * fout = fin
                cout << "WARNING: Matrix is not diagonally dominant\n";
            }
        } 
        else {
            if ( !(Gauss_Seidel(Alpha, fin, fout)) ) {  // Invert A * fout = fin
                 cout << "WARNING: Matrix is not diagonally dominant\n";
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
    Anisotropic_Collisions::Anisotropic_Collisions()   
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       : f00(0.0, Inputdata::IN().list().nump),
         fc(0.0, Inputdata::IN().list().nump),
         l0(Inputdata::IN().inp().l0), 
         m0(Inputdata::IN().inp().m0) 
         {
     
        using Inputdata::IN; 
        Nbc = IN().list().RKLevel;
        szx = IN().inp().x.dim(); 
        szy = IN().inp().y.dim();

    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
    Stat& Anisotropic_Collisions::f1_loop1D(Stat& Yin, const double Dt){ 
//-------------------------------------------------------------------
//  This is the calculation for the harmonics f10, f11 
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
               for (size_t ip(0); ip < fc.size(); ++ip){
                   f00[ip] = (Yin.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
               }
               // Reset the integrals and coefficients
               implicit_step.reset_coeff(f00, Dt);

               // Loop over the harmonics for this (x,y)
               for(int m(0); m < 1; ++m){ 

                   // This harmonic --> Valarray
                   for (int ip(0); ip < Yin.SH(1,m).nump(); ++ip) { 
                       fc[ip] = Yin.SH(1,m)(ip,ix+Nbc,iy+Nbc);
                   }

                   // Take an implicit step
                   implicit_step.advance(fc, 1);

                   //  Valarray --> This harmonic
                   for (int ip(0); ip < Yin.SH(1,m).nump(); ++ip) { 
                       Yin.SH(1,m)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                   } 
               }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            }
        }

        return Yin;    
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    Stat& Anisotropic_Collisions::f1_loop(Stat& Yin, const double Dt){ 
//-------------------------------------------------------------------
//  This is the calculation for the harmonics f10, f11 
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
               for (size_t ip(0); ip < fc.size(); ++ip){
                   f00[ip] = (Yin.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
               }
               // Reset the integrals and coefficients
               implicit_step.reset_coeff(f00, Dt);

               // Loop over the harmonics for this (x,y)
               for(int m(0); m < 2; ++m){ 

                   // This harmonic --> Valarray
                   for (int ip(0); ip < Yin.SH(1,m).nump(); ++ip) { 
                       fc[ip] = Yin.SH(1,m)(ip,ix+Nbc,iy+Nbc);
                   }

                   // Take an implicit step
                   implicit_step.advance(fc, 1);

                   //  Valarray --> This harmonic
                   for (int ip(0); ip < Yin.SH(1,m).nump(); ++ip) { 
                       Yin.SH(1,m)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                   } 
               }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            }
        }

        return Yin;    
    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
    Stat& Anisotropic_Collisions::flm_loop1D(Stat& Yin, const double Dt){ 
//-------------------------------------------------------------------
//  This is the calculation for the high order harmonics 
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
               for (size_t ip(0); ip < fc.size(); ++ip){
                   f00[ip] = (Yin.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
               }
               // Reset the integrals and coefficients
               // --> Assume you did that when you called f1_loop
                implicit_step.reset_coeff(f00, Dt);

               // Loop over the harmonics for this (x,y)
               for(int l(2); l < l0+1 ; ++l){ 
                   for(int m(0); m < 1; ++m){ 

                       // This harmonic --> Valarray
                       for (int ip(0); ip < Yin.SH(l,m).nump(); ++ip) { 
                           fc[ip] = Yin.SH(l,m)(ip,ix+Nbc,iy+Nbc);
                       }

                       // Take an implicit step
                       implicit_step.advance(fc, l);

                       //  Valarray --> This harmonic
                       for (int ip(0); ip < Yin.SH(l,m).nump(); ++ip) { 
                           Yin.SH(l,m)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                       } 
                   }
               }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            }
        }

        return Yin;    
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    Stat& Anisotropic_Collisions::flm_loop(Stat& Yin, const double Dt){ 
//-------------------------------------------------------------------
//  This is the calculation for the high order harmonics 
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
               for (size_t ip(0); ip < fc.size(); ++ip){
                   f00[ip] = (Yin.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
               }
               // Reset the integrals and coefficients
               // --> Assume you did that when you called f1_loop
               implicit_step.reset_coeff(f00, Dt);

               // Loop over the harmonics for this (x,y)
               for(int l(2); l < l0+1 ; ++l){ 
                   for(int m(0); m < ((m0 < l)? m0:l)+1; ++m){ 

                       // This harmonic --> Valarray
                       for (int ip(0); ip < Yin.SH(l,m).nump(); ++ip) { 
                           fc[ip] = Yin.SH(l,m)(ip,ix+Nbc,iy+Nbc);
                       }

                       // Take an implicit step
                       implicit_step.advance(fc, l);

                       //  Valarray --> This harmonic
                       for (int ip(0); ip < Yin.SH(l,m).nump(); ++ip) { 
                           Yin.SH(l,m)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                       } 
                   }
               }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            }
        }

        return Yin;    
    }
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//*******************************************************************
//   Definition for the Euler Backward method
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
    Euler_Backward::Euler_Backward(Stat& Yin, int tout_start)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),                    // Initialize time = 0
        if_implicit1D(Inputdata::IN().list().if_implicit1D),       // Initialize time = 0
        Y(Yin){
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    void Euler_Backward::advance_1(const double& tnew){    
//-------------------------------------------------------------------
//  Implicit integration for f10 and f11
//-------------------------------------------------------------------

//      Initialization
        if ( if_implicit1D ) { 
            Flm_Coll.f1_loop1D(Y, tnew - t);
        } 
        else {
            Flm_Coll.f1_loop(Y, tnew - t);
        }

 
        // you do not need to advance to the new t until you call advance_lm
        //t  = tnew;
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    void Euler_Backward::advance_lm(const double& tnew){    
//-------------------------------------------------------------------
//  Implicit integration
//-------------------------------------------------------------------

//      Initialization
        if ( if_implicit1D ) { 
            Flm_Coll.flm_loop1D(Y, tnew - t);
        }
        else {
            Flm_Coll.flm_loop(Y, tnew - t);
        }

        t  = tnew;
    }
//-------------------------------------------------------------------
//*******************************************************************


