/*! \brief Interspecies Collisions - Definitions
 * \date   October 17, 2016
 * \file   interspeciescollisions.h
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

//  Declarations
    #include "input.h"
    #include "state.h"
    #include "formulary.h"
    #include "nmethods.h"
    #include "interspeciescollisions.h"
    // #include "collisions.h"
    

//*******************************************************************
//-------------------------------------------------------------------
    interspecies_f00_explicit_step::interspecies_f00_explicit_step(const DistFunc1D& DF1, const DistFunc1D& DF2, const double& deltat)//, int tout_start) //DistFunc1D& DFin,
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
        : dt(deltat),
//        pgrid_s1(Algorithms::MakeAxis(DF1.pmax()/(2.0*DF1(0,0).nump()-1.0),DF1.pmax(),DF1(0,0).nump())),
//        pgrid_s2(Algorithms::MakeAxis(DF2.pmax()/(2.0*DF2(0,0).nump()-1.0),DF2.pmax(),DF2(0,0).nump())),
        pgrid_s1(Algorithms::MakeCAxis(0.0,DF1.pmax(),DF1(0,0).nump())),
        pgrid_s2(Algorithms::MakeCAxis(0.0,DF2.pmax(),DF2(0,0).nump())),
        // f1(DF1(0,0)),f2(DF2(0,0)),
        fslope(0.0, DF1(0,0).nump()),
        // f0(0.0, DF1(0,0).nump()),
        df0(0.0, DF1(0,0).nump()),

        U4(0.0,  DF2(0,0).nump()), ///< Constants for Integrals
         U4m1(0.0, DF2(0,0).nump()), ///< Constants for Integrals
         U2(0.0, DF2(0,0).nump()), ///< Constants for Integrals
         U2m1(0.0, DF2(0,0).nump()), ///< Constants for Integrals
         U1(0.0, DF2(0,0).nump()), ///< Constants for Integrals
         U1m1(0.0, DF2(0,0).nump()), ///< Constants for Integrals
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         J1_s2(0.0, DF2(0,0).nump()), ///<      Integrals
         I2_s2(0.0, DF2(0,0).nump()), ///<      Integrals
         I4_s2(0.0, DF2(0,0).nump()), ///<      Integrals
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         U3(0.0, DF2(0,0).nump()), ///<      Integrals
         Qn(0.0, DF2(0,0).nump()), ///<      Integrals
         Pn(0.0, DF2(0,0).nump()),  ///<      Integrals

         J1_s1(0.0, DF2(0,0).nump()), ///<      Integrals
         I2_s1(0.0, DF2(0,0).nump()), ///<      Integrals
         I4_s1(0.0, DF2(0,0).nump()) ///<      Integrals
        {
        
        Nbc = Input::List().BoundaryCells;
        szx = Input::List().NxLocal[0];

        m2 = DF2.mass()  ; m1 = DF1.mass();
        z2 = DF2.q();      z1 = DF1.q();

        

        double re(2.8179402894e-13);           //classical electron radius
        double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re));
        kpre = re*kp;

        // f1.reserve(DF1(0,0).nump());
         
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         // Determine U4, U4m1, U2, U2m1, U1, U1m1
         for (size_t i(1); i < DF2(0,0).nump(); ++i) {
             U4[i]   = 0.5 * pow(pgrid_s2[i],4)     * (pgrid_s2[i]-pgrid_s2[i-1]);         
             U4m1[i] = 0.5 * pow(pgrid_s2[i-1],4)   * (pgrid_s2[i]-pgrid_s2[i-1]);         
             U2[i]   = 0.5 * pgrid_s2[i]*pgrid_s2[i]      * (pgrid_s2[i]-pgrid_s2[i-1]);         
             U2m1[i] = 0.5 * pgrid_s2[i-1]*pgrid_s2[i-1]  * (pgrid_s2[i]-pgrid_s2[i-1]);         
             U1[i]   = 0.5 * pgrid_s2[i]            * (pgrid_s2[i]-pgrid_s2[i-1]);         
             U1m1[i] = 0.5 * pgrid_s2[i-1]          * (pgrid_s2[i]-pgrid_s2[i-1]);         
         }
         
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         // Determine U3
         for (size_t i(0); i < U3.size(); ++i) {
             U3[i] = pow(pgrid_s2[i],3);
         }
         // Determine Qn 
         Qn[0] = 1.0 / ((pgrid_s2[0]*pgrid_s2[0]*pgrid_s2[1])/2.0);
         for (size_t i(1); i < Qn.size()-1; ++i) {
             Qn[i] = 1.0 / (pgrid_s2[i]*pgrid_s2[i]*(pgrid_s2[i+1]-pgrid_s2[i-1])/2.0);
         }        
         // Determine Pn 
         for (size_t i(0); i < Pn.size()-1; ++i) {
             Pn[i] = 1.0 / ((pgrid_s2[i+1]-pgrid_s2[i])/2.0*(pgrid_s2[i+1]+pgrid_s2[i]));
         }    

    }
//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief      Collisions between species 1 and 2 in 0,0 harmonic
///
/// @param[in]  f1in      Input distribution function, the change to which is returned
/// @param[in]  f2in      Input distribution function, f1 collides with f2
/// @param[in]  Delta_t   timestep
/// @param[in]  df1       the delta f1
///
    // SHarmonic1D&  interspecies_f00_explicit_step::takestep(const SHarmonic1D& f1in, const SHarmonic1D& f2in) {
    valarray<double>  interspecies_f00_explicit_step::takestep(const valarray<double>& f1in, const valarray<double>& f2in) {

        // valarray<double> temp(f1in.size());
        double temp = 0.0;
        double mu = m2/m1;

        calculateintegrals(f1in,f2in);
        remapintegrals();
    
        for (int n(1); n < f1in.size()-1; ++n) {
           df0[n]  = f1in[n+1]-f1in[n-1];
           df0[n] /= pgrid_s1[n+1]-pgrid_s1[n-1];
        }

        df0[0] = f1in[1]-f1in[0];
        df0[0]/= pgrid_s1[1]-pgrid_s1[0];

        df0[f1in.size()-1] = f1in[f1in.size()-1]-f1in[f1in.size()-2];
        df0[f1in.size()-1]/= pgrid_s1[f1in.size()-1]-pgrid_s1[f1in.size()-2];

        int ipp(1);
        int ipm(0);

        Gamma12 = -1.0*kpre*n2*(z1*z1*z2*z2)*formulary.LOGii(m1,z1,n1,T1,
                                                                  m2,z2,n2,T2);

        temp = (I4_s1[ipp]+J1_s1[ipp])*df0[ipp]*pgrid_s1[ipp]+3.0/mu*I2_s1[ipp]*f1in[ipp];          
        temp -= (I4_s1[ipm]+J1_s1[ipm])*df0[ipm]*pgrid_s1[ipm]+3.0/mu*I2_s1[ipm]*f1in[ipm];          
        temp /= pgrid_s1[ipp]-pgrid_s1[ipm];      
        temp /= 3.0*(0.25*pow(pgrid_s1[ipm]+pgrid_s1[ipp],2.0));
                             
        fslope[0] = temp;
        // std::cout << "\n fslope[0]=" <<  temp*Gamma12*dt << "\n\n";

        for (int ip(1);ip < f1in.size()-2;++ip)
        {
          ipp=ip+1;
          ipm=ip-1;
          temp  = (I4_s1[ipp]+J1_s1[ipp])*df0[ipp]*pgrid_s1[ipp]+3.0/mu*I2_s1[ipp]*f1in[ipp];
          // std::cout << "\n fslope1[ " << ip << "]="<<  temp << "\n\n";
          temp -= (I4_s1[ipm]+J1_s1[ipm])*df0[ipm]*pgrid_s1[ipm]+3.0/mu*I2_s1[ipm]*f1in[ipm];
          // std::cout << "\n fslope2[ " << ip << "]="<<  temp << "\n\n";
          temp /= pgrid_s1[ipp]-pgrid_s1[ipm];
          // std::cout << "\n fslope3[ " << ip << "]="<<  temp << "\n\n";
          temp /= 3.0*(0.25*pow(pgrid_s1[ipm]+pgrid_s1[ipp],2.0));
          // std::cout << "\n fslope4[ " << ip << "]="<<  temp << "\n\n";
          // temp *= Gamma12*dt;
          fslope[ip] = temp;

          // std::cout << "\n fslope[ " << ip << "]="<<  temp*Gamma12*dt << "\n\n";
        }
        // std::cout << "\n fslope5[ " <<  Gamma12 << "\n\n";
        ipp=f1in.size()-1;ipm=f1in.size()-2;

        temp  = (I4_s1[ipp]+J1_s1[ipp])*df0[ipp]*pgrid_s1[ipp]+3.0/mu*I2_s1[ipp]*f1in[ipp];          
        temp -= (I4_s1[ipm]+J1_s1[ipm])*df0[ipm]*pgrid_s1[ipm]+3.0/mu*I2_s1[ipm]*f1in[ipm];          
        temp /= pgrid_s1[ipp]-pgrid_s1[ipm];      
        temp /= 3.0*(0.25*pow(pgrid_s1[ipm]+pgrid_s1[ipp],2.0));
                   
        fslope[f1in.size()-1] = temp;
        // std::cout << "\n fslope[last]=" <<  temp*Gamma12*dt << "\n\n";

        fslope=fslope*Gamma12*dt;
        // fslope[f1in.size()-1] = 0.0;
          // }
        // exit(1);
      return fslope;
  }

//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief      Remap Distribution to momentum grid of colliding particles
///
   void interspecies_f00_explicit_step::calculateintegrals(const valarray<double>& f1in, const valarray<double>& f2in) {
      

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate the integrals in a standard manner
       I4_s2[0] = 0;
       for (int n(1); n < I4_s2.size(); ++n) {
           I4_s2[n]  = U4[n]*f2in[n]+U4m1[n]*f2in[n-1]; 
           I4_s2[n] += I4_s2[n-1];
            // std::cout << "\n I4= " << I4_s2[n] << "\n\n";
       }

       I2_s2[0] = 0;
       for (int n(1); n < I2_s2.size(); ++n) {
           I2_s2[n]  = U2[n]*f2in[n]+U2m1[n]*f2in[n-1]; 
           I2_s2[n] += I2_s2[n-1];
           // std::cout << "\n I2= " << I2_s2[n] << "\n\n";
       }

      /// Pause for density and temperature
       n2 = (4.0*M_PI*I2_s2[I2_s2.size()-1]);
       
       T2 = m2*(I4_s2[I4_s2.size()-1]/3.0/I2_s2[I2_s2.size()-1]);

       // std::cout << "\n n2= " << n2 << "\n\n";
       // std::cout << "\n T2= " << T2 << "\n\n";
       
       J1_s2[J1_s2.size()-1] = 0;
       for (int n(J1_s2.size()-2); n > -1; --n) {
           J1_s2[n]  = U1[n+1]*f2in[n+1]+U1m1[n+1]*f2in[n]; 
           J1_s2[n] += J1_s2[n+1];
       }

       // for (int k(0); k < I2_s2.size(); ++k){
       //     I4_s2[k] /= pgrid_s2[k]*pgrid_s2[k];
       //     J1_s2[k] *= pgrid_s2[k];
       // }

       I4_s2 *= 4.0 * M_PI;
       I2_s2 *= 4.0 * M_PI;
       J1_s2 *= 4.0 * M_PI;

       
       /// The other density and temperature
       // for (size_t ip(0); ip < f1in.size() ; ++ip)
       // {
       //    std::cout << "f1[" << pgrid_s1[ip] << "] = " << f1in[ip] << "\n";
          
       //    // f1[ip]=f1in[ip];
       // }
       n1 = 4.0*M_PI*(Algorithms::moment(f1in,pgrid_s1,2.0));
       T1 = 4.0*M_PI*m1*(Algorithms::moment(f1in,pgrid_s1,4.0))/3.0/n1;

       

       // exit(1);

}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief      Remap Distribution to momentum grid of colliding particles
///
   void interspecies_f00_explicit_step::remapintegrals() {
      
        
        double cell_s2            = 0;
        double dp_s2              = pgrid_s2[2]-pgrid_s2[1];
        double dist               = 0.0;

        /// Initialize 
        I4_s1[0] = 0.0; 
        I2_s1[0] = 0.0;
        J1_s1[J1_s1.size()-1] = 0.0;

        int n(0);

        while ( (n < pgrid_s1.size()) && (pgrid_s1[n] < pgrid_s2[pgrid_s2.size()-1]))
        {            
            dist      = modf(   pgrid_s1[n]/dp_s2 , &cell_s2            ); /// Find nearest cell and distance from it
            // std::cout << "p[" << n << "]=" << pgrid_s1[n] << ", dist=" << dist << ", cell_s2 = " << cell_s2 << "\n";
            I4_s1[n]  = I4_s2[static_cast<size_t>(cell_s2)]*(1.0-dist) + I4_s2[static_cast<size_t>(cell_s2)+1]*dist; /// Interpolate
            I2_s1[n]  = I2_s2[static_cast<size_t>(cell_s2)]*(1.0-dist) + I2_s2[static_cast<size_t>(cell_s2)+1]*dist;
            J1_s1[n]  = J1_s2[static_cast<size_t>(cell_s2)]*(1.0-dist) + J1_s2[static_cast<size_t>(cell_s2)+1]*dist;
            ++n;            
        } 
        while (n < pgrid_s1.size())
        {
            I4_s1[n]  = I4_s2[I4_s2.size()];
            I2_s1[n]  = I2_s2[I2_s2.size()];
            J1_s1[n]  = 0.0;
            ++n;
        }

    }
//-------------------------------------------------------------------
//*******************************************************************
//------------------------------------------------------------------------------
/// @brief      Constructor for RK4 method on f00
///
/// @param      fin         Input distribution
/// @param[in]  tout_start  Hmm...
///
    interspecies_f00_RKfunctor::interspecies_f00_RKfunctor(const DistFunc1D& DF1,const DistFunc1D& DF2,const double& smalldt)
        :collide(DF1,DF2,smalldt){}
//--------------------------------------------------------------
    void interspecies_f00_RKfunctor::operator()(const valarray<double>& fin1, const valarray<double>& fin2, valarray<double>& fslope) { 
        fslope = collide.takestep(fin1,fin2);
    }
    void interspecies_f00_RKfunctor::operator()(const valarray<double>& fin1, valarray<double>& fslope) {}
    void interspecies_f00_RKfunctor::operator()(const valarray<double>& fin1, valarray<double>& fslope, size_t dir) {}

//*******************************************************************
//-------------------------------------------------------------------
    interspecies_f00_explicit_collisions::interspecies_f00_explicit_collisions(const DistFunc1D& DF1,const DistFunc1D& DF2, const double& deltat)//, int tout_start) //DistFunc1D& DFin,
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
      :
        rkf00(DF1, DF2, deltat),
        fin1(0.0,DF1(0,0).nump()), fin2(0.0,DF2(0,0).nump()),
        RK4(fin1),
        num_h(size_t(deltat/Input::List().small_dt)+1)
        // h(Input::List().small_dt)
        {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          h = deltat/static_cast<double>(num_h);
        Nbc = Input::List().BoundaryCells;
        szx = Input::List().NxLocal[0];  //Input::List().numx;      // size of useful x axis
    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
    void interspecies_f00_explicit_collisions::rkloop(SHarmonic1D& SH1, const SHarmonic1D& SH2){    
//-------------------------------------------------------------------
//  This loop scans all the locations in configuration space
//  and calls the RK4 for explicit integration at each location
//-------------------------------------------------------------------
      // std::cout << "\n\n DF1 = " << fin1.size();
      // std::cout << ", DF2 = " << fin2.size() << " \n\n";
      // valarray<double>                        fin1(0.0,SH1.nump()),fin2(0.0,SH2.nump());
      for (size_t ix(0); ix < szx-2*Nbc; ++ix){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          // Copy data for a specific location in space to valarray
          for (size_t ip(0); ip < fin1.size(); ++ip){
              fin1[ip] = (SH1(ip,ix+Nbc)).real();              
          }

          for (size_t ip(0); ip < fin2.size(); ++ip){
              fin2[ip] = (SH2(ip,ix+Nbc)).real();
          }
           
          // // Time loop: Update the valarray
          for (size_t h_step(0); h_step < num_h; ++h_step){
//              fin1 = RK4(fin1,fin2,h,&rkf00);
          }
          // // std::cout << "\n\n 11 \n\n";
          // // Return updated data to the harmonic
          for (int ip(0); ip < fin1.size(); ++ip){
              // std::cout << "\n\nfin[" << ip << "] = "<< fin1[ip];
              SH1(ip,ix+Nbc) = fin1[ip];
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
    interspecies_flm_implicit_step::interspecies_flm_implicit_step(double pmax, size_t nump)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : 
//       Pre-Calculated Constants
//      INCLUDED pmax[0] here only because pmax currently is vector having elements for two species
//         vr(Algorithms::MakeAxis(pmax/(2.0*nump-1.0),pmax,nump)),
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
         if_tridiagonal(Input::List().if_tridiagonal)
//       Make matrix vk/vn
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
         double re(2.8179402894e-13);           //classical electron radius
         double kp(sqrt(4.0*M_PI*(Input::List().density_np)*re));
         kpre = re*kp;

         // Determine vr
         for (size_t i(0); i < vr.size(); ++i) {
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

//------------------------------------------------------------------------------
/// @brief      Resets coefficients and integrals to use in the matrix solve.
///
/// @param[in]  fin      Input distribution function
/// @param[in]  Delta_t  timestep
///
    void  interspecies_flm_implicit_step::reset_coeff(const valarray<double>& fin, const double Delta_t) {
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
       _LOGee  = formulas.LOGee(I0_density,I2_temperature/I0_density);
       _ZLOGei = formulas.Zeta*formulas.LOGei(I0_density,I2_temperature/I0_density,formulas.Zeta);
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
          IvDnDm1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]  -vr[i-1]));       //  ( v*D_n*D_{n-1/2} )^(-1)
          IvDnDp1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1]-vr[i]  ));       //  ( v*D_n*D_{n+1/2} )^(-1)
          Ivsq2Dn = 1.0 / (vr[i] * vr[i]                 * (vr[i+1]-vr[i-1]));       //  ( v^2 * 2*D_n )^(-1)
          Alpha_Tri(i, i  ) = 8.0 * M_PI * fin[i] - TriI2[i] * (IvDnDm1 + IvDnDp1);
          Alpha_Tri(i, i-1) = TriI2[i] * IvDnDm1 - TriI1[i] * Ivsq2Dn;   
          Alpha_Tri(i, i+1) = TriI2[i] * IvDnDp1 + TriI1[i] * Ivsq2Dn;   
       }
       Alpha_Tri *=  (-1.0) * _LOGee * kpre * Dt;         // (-1) because the matrix moves to the LHS in the equation
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


//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    }
//-------------------------------------------------------------------

//------------------------------------------------------------------------------
/// @brief      Perform a matrix solve to calculate effect of collisions on f >= 1
///
/// @param      fin   Input distribution function
/// @param[in]  el    Number of elements in matrix (?)
///
    void  interspecies_flm_implicit_step::advance(valarray<complex<double> >& fin, const int el) {
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


//      SOLVE A * Fout  = Fin
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        //for (int i(0); i < fin.size(); ++i){   // Estimate Fout = Fin /  A(i,i)
        //    fout[i] /= Alpha(i,i); 
        //}

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
    interspecies_flm_implicit_collisions::interspecies_flm_implicit_collisions(const DistFunc1D &DFin, const double& deltat)   
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       : f00(0.0, DFin(0).nump()),
         fc(0.0,DFin(0).nump()),
         l0(DFin.l0()), 
         m0(DFin.m0()),//Yin.SH(0,0).m0()) 
         implicit_step(DFin.pmax(),DFin(0).nump()),
         Dt(deltat)    //
         {
     
        
        Nbc = Input::List().BoundaryCells;
        szx = Input::List().NxLocal[0];

    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
    void interspecies_flm_implicit_collisions::advancef1(DistFunc1D& DF){ 
//-------------------------------------------------------------------
//  This is the calculation for the harmonics f10, f11 
//    To be specific, this routine does just l=1
//    To be specific, this routine differs from 'f1_loop1D' above only in one place: the loop is on m=0 and m=1
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
      
      for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
            for (size_t ip(0); ip < fc.size(); ++ip){
              f00[ip] = (DF(0,0)(ip,ix+Nbc)).real();
            }
               // Reset the integrals and coefficients
            implicit_step.reset_coeff(f00, Dt);

               // Loop over the harmonics for this (x,y)
            for(int m(0); m < 2; ++m){ 

                  // This harmonic --> Valarray
                  for (int ip(0); ip < DF(1,m).nump(); ++ip) { 
                      fc[ip] = DF(1,m)(ip,ix+Nbc);
                  }

        // Take an implicit step
                  implicit_step.advance(fc, 1);

          //  Valarray --> This harmonic
                  for (int ip(0); ip < DF(1,m).nump(); ++ip) { 
                      DF(1,m)(ip,ix+Nbc) = fc[ip];
                  } 
            }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
          }
    }
//-------------------------------------------------------------------



//-------------------------------------------------------------------
    void interspecies_flm_implicit_collisions::advanceflm(DistFunc1D& DF)
    { 
//-------------------------------------------------------------------
//  This is the calculation for the high order harmonics 
//    To be specific, this routine does l=2 to l0
//    To be specific, this routine differs from 'flm_loop1D' aboe only in one place: the loop over m all m from 0 to ((m0 < l)? m0:l)
//-------------------------------------------------------------------

        // For every location in space within the domain of this node
    
        for (size_t ix(0); ix < szx-2*Nbc; ++ix){

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               // "00" harmonic --> Valarray
            for (size_t ip(0); ip < fc.size(); ++ip){
                f00[ip] = (DF(0,0)(ip,ix+Nbc)).real();
            }
            // Reset the integrals and coefficients
            // --> Assume you did that when you called f1_loop
            implicit_step.reset_coeff(f00, Dt);
            // Loop over the harmonics for this (x,y)
            for(int l(2); l < l0+1 ; ++l){ 
                for(int m(0); m < ((m0 < l)? m0:l)+1; ++m){ 

                    // This harmonic --> Valarray
                    for (size_t ip(0); ip < fc.size(); ++ip){
                        fc[ip] = (DF(l,m))(ip,ix+Nbc);
                    }

                    // Take an implicit step
                    implicit_step.advance(fc, l);

                    //  Valarray --> This harmonic
                    for (size_t ip(0); ip < fc.size(); ++ip){
                        DF(l,m)(ip,ix+Nbc) = fc[ip];
                    } 
                }
        }
    }
    
    return;  
    
    }
//-------------------------------------------------------------------


////*******************************************************************
//-------------------------------------------------------------------
    interspecies_collisions::interspecies_collisions(const State1D& Yin, const size_t& sind, const double& deltat)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       // : interspecies_f00_collisions(WS, deltat),
          // interspecies_flm_collisions(DFin, deltat){
    {

        // size_t sind(0);  
        for (size_t s(0); s<Yin.Species(); ++s){
            if (s!=sind) {
              
              interspecies_f00_collisions.push_back(interspecies_f00_explicit_collisions(Yin.DF(sind),Yin.DF(s),deltat));    
              std::cout << "\n\n sind,s = " << sind << "," << s << ", size = " <<interspecies_f00_collisions.size() << "\n\n";
            }
        }

    }  

//-------------------------------------------------------------------

//-------------------------------------------------------------------
//-------------------------------------------------------------------
    // void interspecies_collisions::advanceflm(DistFunc1D& DFin){
//-------------------------------------------------------------------   
      // interspecies_flm_collisions.advanceflm(DFin);
// }


//-------------------------------------------------------------------
    // void interspecies_collisions::advancef1(DistFunc1D& DFin){
//-------------------------------------------------------------------
//       interspecies_flm_collisions.advancef1(DFin);
// }    
//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief      Remap Distribution to momentum grid of colliding particles
///
// interspecies_collisions::interspecies_collisions(const State1D& Yin, const size_t sin, const double& deltat)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------

    
    // }

//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief     Advance f00 collisions for all species
///
void interspecies_collisions::advancef00( State1D& Yin, const size_t& sind){
// //-------------------------------------------------------------------
      // size_t sind(0);

      for (size_t s(0); s<Yin.Species(); ++s){
            if (s!=sind) {
              
              interspecies_f00_collisions[sind].rkloop(Yin.SH(sind,0,0),Yin.SH(s,0,0));
              std::cout << "\n\n sind,s = " << sind << "," << s << "\n\n";
              // ++sind;
          }

      }  
}

//-------------------------------------------------------------------
//------------------------------------------------------------------------------
/// @brief     Advance f1 collisions for all species
///
// void interspecies_collisions::advancef1( State1D& Yin, const size_t sin){
// // //-------------------------------------------------------------------
//       size_t sind(0);
//       for (size_t s(0); s<Yin.Species(); ++s){
//           // sind=0;
//           // for (size_t so(0); so<Yin.Species(); ++s)
//           // {
//             if (s!=sin) {
//               interspecies_flm_collisions[sind].advancef1(Yin.DF(sin,0,0),Yin.SH(s,0,0));
//               ++sind;
//             // }
//           }
//       }  
// }

// //-------------------------------------------------------------------
// //------------------------------------------------------------------------------
// /// @brief     Advance f1 collisions for all species
// ///
// void interspecies_collisions::advanceflm( State1D& Yin, const size_t sin){
// // //-------------------------------------------------------------------
//       size_t sind(0);
//       for (size_t s(0); s<Yin.Species(); ++s){
//           // sind=0;
//           // for (size_t so(0); so<Yin.Species(); ++s)
//           // {
//             if (s!=sin) {
//               interspecies_flm_collisions[sind].advancef1(Yin.SH(sin,0,0),Yin.SH(s,0,0));
//               ++sind;
//             // }
//           }
//       }  
// }
//*******************************************************************

//*******************************************************************