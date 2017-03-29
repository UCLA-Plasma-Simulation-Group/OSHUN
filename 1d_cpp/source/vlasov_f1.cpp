/*! \brief  Vlasov Equation - Definitions
* \author PICKSC
 * \date   October 10, 2016
 * \file   vlasov_f1.cpp
 *
 * Includes spatial advection, electric field advection, and electric field update routines
 * 
 */
//--------------------------------------------------------------
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
#include "nmethods.h"

//  Declerations
#include "state.h"
#include "formulary.h"
#include "input.h"
#include "fluid.h"
#include "vlasov.h"
#include "vlasov_f1.h"


//**************************************************************
//--------------------------------------------------------------
Electric_Field_1D_f1::Electric_Field_1D_f1(size_t Nl, size_t Nm,
                                     double pmin, double pmax, size_t Np,
                                     double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : Hp0(Nl+1),
          H(Np,Nx), G(Np,Nx), TMP(Np,Nx),
          pr(Algorithms::MakeAxis(static_cast<complex<double> >(pmin),
                                  static_cast<complex<double> >(pmax),
                                  Np)),
          invpr(pr)
{

//      - - - - - - - - - - - - - - - - - - - - - - - - - - -
         complex<double> lc, mc;

//       Inverted momentum axis
         for (size_t i(0); i < pr.size(); ++i) { 
             invpr[i] = 1.0/pr[i]; 
         }
         double idp = (-1.0)/ (2.0*(pmax-pmin)/double(Np-1)); // -1/(2dp) 

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        A100 = static_cast<complex<double>>(idp);
        C100 = static_cast<complex<double>>(0.5*idp);
        A210 = static_cast<complex<double>>(1.0/3.0*idp);
        B211 = static_cast<complex<double>>(2.0/3.0);
        A310 = static_cast<complex<double>>(2.0/5.0*idp);
        C311 = static_cast<complex<double>>(-0.5*idp/5.0);

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       -2*Dp*H at the 0 momentum cell
         Hp0[0] = 1.0 / pr[0];
         for (size_t l(1); l < Nl+1; ++l) {
             double ld(l);
             Hp0[l] = Hp0[l-1] * (pr[0]/pr[1]) * (2.0*ld+1.0)/(2.0*ld-1.0);
         }
         Hp0 *= (-2.0)*(pr[1]-pr[0]);

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Electric_Field_1D_f1::operator()(const DistFunc1D& Din,
                                   const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                                   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

    complex<double> ii(0.0,1.0);

    valarray<complex<double> > Ex(FEx.array());
    valarray<complex<double> > Em(FEz.array());
    Em *= (-1.0)*ii;
    Em += FEy.array();
    valarray<complex<double> > Ep(FEz.array());
    Ep *= ii;
    Ep += FEy.array();

    Ex *= Din.q();
    Em *= Din.q();
    Ep *= Din.q();

    size_t l0(Din.l0());
    size_t m0(Din.m0());


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeG00(Din(0,0)); 
    Ex *= A100; TMP = G; Dh(1,0) += G.mxaxis(Ex);
    Em *= C100;            Dh(1,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(1,0),1);
    Ex *= A210 / A100;          Dh(0,0) += H.mxaxis(Ex); 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, 1 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        // Din(2,0)  = Din(0,0);
        // Din(2,0) *= static_cast<complex<double>>(1.0/3.0);

        // MakeGH(Din(2,0),2);
        // Ex *= A310 / A210;      TMP = H;      Dh(1,0) += H.mxaxis(Ex);
        // Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(1,1),1);
        Ep *= B211;              H = H.mxaxis(Ep);     Dh(0,0) += H.Re();

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Electric_Field_1D_f1::Implicit_Ex(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

    valarray<complex<double> > Ex(FEx.array());
    Ex *= Din.q();

//      m = 0, l = 0
        MakeG00(Din(0,0)); 
        Ex *= A100;             Dh(1,0) += G.mxaxis(Ex);     

//      m = 0, l = 1
        MakeGH(Din(1,0),1);
        Ex *= A210 / A100;   Dh(0,0) += H.mxaxis(Ex);
        // Ex *= A1(1,0) / A2(1,0);   Dh(2,0) += G.mxaxis(Ex);

//      m = 0, l = 2
        // Din(2,0)  = Din(0,0);
        // Din(2,0) *= static_cast<complex<double>>(1.0/3.0);

        // MakeGH(Din(2,0),2);
        // Ex *= A310 / A210;      TMP = H;      Dh(1,0) += H.mxaxis(Ex);


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field_1D_f1::Implicit_Ey(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

    valarray<complex<double> > Em(FEy.array());
    valarray<complex<double> > Ep(FEy.array());


    Em *= Din.q();;
    Ep *= Din.q();;

//      m = 0, l = 0
        MakeG00(Din(0,0)); 
        Em *= C100;            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
        // MakeGH(Din(1,0),1);
        // Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
        // Din(2,0)  = Din(0,0);
        // Din(2,0) *= static_cast<complex<double>>(1.0/3.0);

        // MakeGH(Din(2,0),2);
        // Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

//      m = 0,  l = 3 
        // MakeGH(Din(3,0),3);
        // Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),1);
        Ep *= B211;             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
}
//
////--------------------------------------------------------------
void Electric_Field_1D_f1::Implicit_Ez(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------
    complex<double> ii(0.0,1.0);


    valarray<complex<double> > Em(FEz.array());
    Em *= (-1.0)*ii;


    valarray<complex<double> > Ep(FEz.array());
    Ep *= ii;
    //      m = 0, l = 0
        MakeG00(Din(0,0)); 
        Em *= C100;            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
        // MakeGH(Din(1,0),1);
        // Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
        Din(2,0)  = Din(0,0);
        Din(2,0) *= static_cast<complex<double>>(1.0/3.0);

        MakeGH(Din(2,0),2);
        Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

// //      m = 0,  l = 3 
//         MakeGH(Din(3,0),3);
//         Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),1);
        Ep *= B211;             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Make derivatives 2*Dp*(l+1/l)*G and -2*Dp*H for a given f
void Electric_Field_1D_f1::MakeGH(SHarmonic1D& f, size_t el){
//--------------------------------------------------------------
    valarray<complex<double> > invpax(invpr);
    complex<double> ld(el);

    invpax *= (-2.0)*(ld+1.0) * (pr[1]-pr[0]);

    G = f;                   H = f;
    G = G.Dp();
    H  = H.mpaxis(invpax);
    H += G;
    G *= -(2.0*ld+1.0)/ld;
    G += H;

    for (size_t i(0); i < G.numx(); ++i) G(0,i) = 0.0;
    for (size_t i(0); i < H.numx(); ++i) H(0,i) = f(1,i) * Hp0[el];
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Calculation of G00 = -2*Dp* df/dp(p0)
void Electric_Field_1D_f1::MakeG00(SHarmonic1D& f) {
//--------------------------------------------------------------
    G = f; G = G.Dp();

    complex<double> p0p1_sq( pr[0]*pr[0]/(pr[1]*pr[1]) ),
            inv_mp0p1_sq( 1.0/(1.0-p0p1_sq) ),
            g_r = -4.0*(pr[1]-pr[0]) * pr[0]/(pr[1]*pr[1]),
            f00;

    for (size_t i(0); i < f.numx(); ++i) {
        f00    = ( f(0,i) - f(1,i) * p0p1_sq) * inv_mp0p1_sq;
        G(0,i) = ( f(1,i) - f00) * g_r;
    }
}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Magnetic_Field_1D_f1::Magnetic_Field_1D_f1(size_t Nl, size_t Nm,
                                     double pmin, double pmax, size_t Np,
                                     double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : A1(Nm+1), B1(Nl+1), A2(Nl+1,Nm+1), A3(0.5),
          FLM(Np,Nx)
{
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -
    complex<double> lc, mc;
    complex<double> c01(0.0,1.0);

//       Calculate the "A1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t m(0); m < Nm+1; ++m){
        mc = static_cast< complex<double> >(m);
        A1[m] = (-1.0)*c01*mc;
    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // for (size_t l(1); l < Nl+1; ++l)
    //     for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
    //         lc = static_cast< complex<double> >(l);
    //         mc = static_cast< complex<double> >(m);
    //         A2(l,m) = (-0.5)*(lc+1.0-mc)*(lc+mc);
    //     }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    A3 = 0.5;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(0); l < Nl+1; ++l){
        lc = static_cast< complex<double> >(l);
        B1[l] = (-1.0)*lc*(lc+1.0);
    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Magnetic_Field_1D_f1::operator()(const DistFunc1D& Din,
                                   const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                                   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------

    complex<double> ii(0.0,1.0);

    valarray<complex<double> > Bx(FBx.array());
    valarray<complex<double> > Bm(FBy.array());
    Bm *= (-1.0)*ii;
    Bm += FBz.array();
    valarray<complex<double> > Bp(FBy.array());
    Bp *= ii;
    Bp += FBz.array();

    Bx *= Din.q();
    Bm *= Din.q();
    Bp *= Din.q();

    size_t l0(B1.size()-1);
    size_t m0(A1.size()-1);

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, 1 < l < l0+1
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Bp *= A3;
    for (size_t l(1); l < l0+1; ++l){
        FLM = Din(l,0);      Dh(l,1) += FLM.mxaxis(Bp);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FLM = Din(1,1); Bx *= A1[1];                           Dh(1,1) += FLM.mxaxis(Bx);
    FLM = Din(1,1); Bm *= B1[1]; FLM = FLM.mxaxis(Bm);     Dh(1,0) += FLM.Re();

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l > 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(2); l < l0+1; ++l){
//         FLM = Din(l,1);                                                Dh(l,2) += FLM.mxaxis(Bp);
//         FLM = Din(l,1);                                                Dh(l,1) += FLM.mxaxis(Bx);
//         FLM = Din(l,1); Bm *= B1[l]/B1[l-1]; FLM = FLM.mxaxis(Bm); Dh(l,0) += FLM.Re();
//     }
//     Bm *= 1.0/B1[l0];

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m > 1, l = m
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t m(2); m < m0; ++m){
//         FLM = Din(m,m); Bx *= A1[m]/A1[m-1];                   Dh(m,m  ) += FLM.mxaxis(Bx);
//         FLM = Din(m,m); Bm *= A2(m,m);                         Dh(m,m-1) += FLM.mxaxis(Bm);
//         for (size_t l(m+1); l < l0+1; ++l){
//             FLM = Din(l,m);                                    Dh(l,m+1) += FLM.mxaxis(Bp);
//             FLM = Din(l,m);                                    Dh(l,m  ) += FLM.mxaxis(Bx);
//             FLM = Din(l,m); Bm *= A2(l,m)/A2(l-1,m);           Dh(l,m-1) += FLM.mxaxis(Bm);
//         }
//         Bm *= 1.0/A2(l0,m);
//     }

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = m0, l >= m0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     FLM = Din(m0,m0); Bx *= A1[m0]/A1[m0-1];                   Dh(m0,m0)   += FLM.mxaxis(Bx);
//     FLM = Din(m0,m0); Bm *= A2(m0,m0)/*/A2(l0,m0-1)*/;             Dh(m0,m0-1) += FLM.mxaxis(Bm);
//     for (size_t l(m0+1); l < l0+1; ++l){
//         FLM = Din(l,m0);                                     Dh(l,m0  )  += FLM.mxaxis(Bx);
//         FLM = Din(l,m0); Bm *= A2(l,m0)/A2(l-1,m0);          Dh(l,m0-1)  += FLM.mxaxis(Bm);
//     }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Magnetic_Field_1D_f1::implicit(DistFunc1D& Din,
                                 const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                                 double dt) {
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------

    complex<double> ii(0.0,1.0);

    valarray<complex<double> > Bx(FBx.array());
    valarray<complex<double> > Bm(FBy.array());
    Bm *= (-1.0)*ii;
    Bm += FBz.array();
    valarray<complex<double> > Bp(FBy.array());
    Bp *= ii;
    Bp += FBz.array();

    Bx *= Din.q();//Din.mass();
    Bm *= Din.q();//Din.mass();
    Bp *= Din.q();//Din.mass();

    size_t l0(1);
    size_t m0(1);

    int Nm(0);
    size_t ind(0);
    size_t indp(0);

    complex<double> temp(0.0,0.0);

    for (size_t i(0); i<Din(0,0).numx(); ++i)
    {
        for (size_t l(1); l < l0+1; ++l)
        {
            Nm = (l<m0?l:m0);
            valarray<complex<double>> fin((0.0,0.0),2*Nm+1);
            valarray<complex<double>> RHSv((0.0,0.0),2*Nm+1);
            Array2D<complex<double>>  LHS(2*Nm+1,2*Nm+1);
            Array2D<complex<double>>  RHS(2*Nm+1,2*Nm+1);

            //* The matrices corresponding to RHS and LHS are computed
            //  They apply to the vector
            //  [ f_l^m  ]          -> 0 index but applies to Nm th harmonic
            //     ...
            //  [ f_l^2  ]
            //  [ f_l^1  ]
            //  [ f_l^0  ]          -> Nm index but applies to 0th harmonic
            //  [ f_l^-1 ]
            //  [ f_l^-2 ]
            //     ...
            //  [ f_l^-m ]          *//
            RHS(0,0)           = 1.0-ii*double(Nm)*dt/2.0*Bx[i];
            LHS(0,0)           = conj(RHS(0,0));
            RHS(0,1)           =  0.25*dt*Bp[i];
            LHS(0,1)           = -0.25*dt*Bp[i];

            for (int m(Nm-1);m>0;--m)
            {
                ind = Nm-m;

                RHS(ind,ind)   = 1.0-ii*double(m)*dt/2.0*Bx[i];
                LHS(ind,ind)   = conj(RHS(ind,ind));

                RHS(ind,ind+1) =  0.25*dt*Bp[i];
                LHS(ind,ind+1) = -0.25*dt*Bp[i];

                RHS(ind,ind-1) = -0.25*dt*(l-m)*(l+m+1)*Bm[i];
                LHS(ind,ind-1) =  0.25*dt*(l-m)*(l+m+1)*Bm[i];
            }

            RHS(Nm,Nm)         = 1.0;
            LHS(Nm,Nm)         = 1.0;
            RHS(Nm,Nm-1)       = -0.25*dt*(l)*(l+1)*Bm[i];
            LHS(Nm,Nm-1)       =  0.25*dt*(l)*(l+1)*Bm[i];

            //* Fill in bottom half which is a series of reflections and complex conjugates
            RHS(Nm,Nm+1)       = conj(RHS(Nm,Nm-1));
            LHS(Nm,Nm+1)       = conj(LHS(Nm,Nm-1));

            for (int m(-1);m>-Nm;--m)
            {
                ind  = Nm-m;
                indp = Nm-abs(m);

                RHS(ind,ind)   = conj(RHS(indp,indp));
                LHS(ind,ind)   = conj(LHS(indp,indp));

                RHS(ind,ind+1) = conj(RHS(indp,indp-1));
                LHS(ind,ind+1) = conj(LHS(indp,indp-1));

                RHS(ind,ind-1) = conj(RHS(indp,indp+1));
                LHS(ind,ind-1) = conj(LHS(indp,indp+1));
            }

            RHS(2*Nm,2*Nm)     = conj(RHS(0,0));
            LHS(2*Nm,2*Nm)     = conj(LHS(0,0));
            RHS(2*Nm,2*Nm-1)   = conj(RHS(0,1));
            LHS(2*Nm,2*Nm-1)   = conj(LHS(0,1));

            for (size_t k(0); k < Din(0,0).nump(); ++k) {
                /// Distribution function vector. 1 Nmx1 vector per momentum cell
                fin[0] = Din(l, Nm)(k, i);
                for (int m(Nm - 1); m > 1; --m) {
                    ind = Nm - m;
                    fin[ind] = Din(l, m)(k, i);
                }
                fin[Nm] = Din(l, 0)(k, i);
                for (int m(-1); m > -Nm; --m) {
                    ind = Nm - m;
                    fin[ind] = conj(Din(l, abs(m))(k, i));
                }
                fin[2 * Nm] = conj(Din(l, Nm)(k, i));

                /// Multiply Right Side to create right side vector
                RHSv[0] = fin[0] * RHS(0, 0) + fin[1] * RHS(0, 1);
                for (size_t mm(1); mm < 2 * Nm; ++mm) {
                    RHSv[mm] = fin[mm - 1] * RHS(mm, mm - 1) + fin[mm] * RHS(mm, mm) + fin[mm + 1] * RHS(mm, mm + 1);
                }
                RHSv[2 * Nm] = fin[2 * Nm - 1] * RHS(2 * Nm, 2 * Nm - 1) + fin[2 * Nm] * RHS(2 * Nm, 2 * Nm);

//                std::cout << "\n\n LHS = \n";
//                for (size_t i(0); i < LHS.dim1(); ++i) {
//                    std::cout << "i = " << i << " :::: ";
//                    for (size_t j(0); j < LHS.dim2(); ++j) {
//                        std::cout << LHS(i, j) << "   ";
//                    }
//                    std::cout << "\n";
//                }

                Thomas_Tridiagonal(LHS,RHSv,fin);

                /// Unpack
                Din(l,Nm)(k,i)     = fin[0];
                for (int m(Nm-1);m>1;--m)
                {
                    ind = Nm-m;
                    Din(l,m)(k,i)  = fin[ind];
                }
                Din(l,0)(k,i)     = fin[Nm];

            }

        }
    }

}
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Spatial_Advection_1D_f1::Spatial_Advection_1D_f1(size_t Nl, size_t Nm,
                                           double pmin, double pmax, size_t Np,
                                           double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : fd1(Np,Nx),
          vr(Algorithms::MakeAxis(static_cast<complex<double> >(pmin),
                                  static_cast<complex<double> >(pmax),
                                  Np)) {
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -

    complex<double> lc, mc;

    for (size_t i(0); i < vr.size(); ++i) {
        vr[i] = vr[i]/(sqrt(1.0+vr[i]*vr[i]));
    }

         double idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx)); // -1/(2dx)

//       Calculate the "A1, A2" parameters
        A00 = static_cast<complex<double>>(idx);
        A10 = static_cast<complex<double>>(idx/3.0);
        A20 = static_cast<complex<double>>(idx*2.0/5.0);

}
//--------------------------------------------------------------

//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection_1D_f1::operator()(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vt(vr); vt *= 1.0/Din.mass();

    fd1 = Din(0,0);     
    fd1 = fd1.Dx();     vt *= A00;      
    Dh(1,0)+=(fd1.mpaxis(vt));
    

    fd1 = Din(1,0);     
    fd1 = fd1.Dx();     vt *= A10/A00; 
    Dh(0,0)+=(fd1.mpaxis(vt));

    fd1 = Din(0,0);     fd1 *= static_cast<complex<double>>(1.0/3.0);
    fd1 = fd1.Dx();     vt *= A20/A10; 
    Dh(1,0)+=(fd1.mpaxis(vt));          


}
//--------------------------------------------------------------


//**************************************************************
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with explicit 
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_f1_explicitE::VlasovFunctor1D_f1_explicitE(vector<size_t> Nl,vector<size_t> Nm,vector<double> pmax, vector<size_t> Np,
                                                     double xmin, double xmax, size_t Nx) {
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );

        SA.push_back( Spatial_Advection_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        EF.push_back( Electric_Field_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        JX.push_back( Current_1D(pmin, pmax[s], Np[s], Nx) );

        BF.push_back( Magnetic_Field_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

//        HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        AM.push_back( Ampere_1D(xmin, xmax, Nx) );

        FA.push_back( Faraday_1D(xmin, xmax, Nx) );

    }
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor1D_f1_explicitE::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        SA[s](Yin.DF(s),Yslope.DF(s));

        EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

//        if (Input::List().hydromotion) HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));

        BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

        AM[s](Yin.EMF(),Yslope.EMF());

        FA[s](Yin.EMF(),Yslope.EMF());

        Yslope.DF(s) = Yslope.DF(s).Filterp();
    }

}

void VlasovFunctor1D_f1_explicitE::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_f1_explicitE::operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope){}

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_f1_explicitEB::VlasovFunctor1D_f1_explicitEB(vector<size_t> Nl,vector<size_t> Nm,vector<double> pmax, vector<size_t> Np,
                                                       double xmin, double xmax, size_t Nx) {
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );

        SA.push_back( Spatial_Advection_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        EF.push_back( Electric_Field_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        JX.push_back( Current_1D(pmin, pmax[s], Np[s], Nx) );

//        BF.push_back( Magnetic_Field_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

//        HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        AM.push_back( Ampere_1D(xmin, xmax, Nx) );

        FA.push_back( Faraday_1D(xmin, xmax, Nx) );

    }
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor1D_f1_explicitEB::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        SA[s](Yin.DF(s),Yslope.DF(s));

        EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

//        if (Input::List().hydromotion) HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));

//        BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

        AM[s](Yin.EMF(),Yslope.EMF());

        FA[s](Yin.EMF(),Yslope.EMF());

        // std::cout << "\n";

        // for (int p(0);p < Yslope.SH(0,1,0).nump(); ++p)
        // {
        //     std::cout << Yslope.SH(0,1,0)(p,43).real() << "\n";
        // }

        // exit(1);
        Yslope.DF(s) = Yslope.DF(s).Filterp();
    }

}

void VlasovFunctor1D_f1_explicitEB::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_f1_explicitEB::operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope){}
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with implicit
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_f1_implicitE_p1::VlasovFunctor1D_f1_implicitE_p1(vector<size_t> Nl,vector<size_t> Nm,vector<double> pmax, vector<size_t> Np,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );

        SA.push_back( Spatial_Advection_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        // EF.push_back( Electric_Field_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        // JX.push_back( Current_1D(pmin, pmax[s], Np[s], Nx) );

        // AM.push_back( Ampere_1D(xmin, xmax, Nx) );

//        HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        BF.push_back( Magnetic_Field_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

//        FA.push_back( Faraday_1D(xmin, xmax, Nx) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor1D_f1_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        // Yslope.DF(s).checknan();std::cout << "Vlasov 1 \n";
        SA[s](Yin.DF(s),Yslope.DF(s));
        // Yslope.DF(s).checknan();std::cout << "Vlasov 2 \n";

        // Yslope.DF(s).checknan();std::cout << "Vlasov 3 \n";

        // EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

        // AM[s](Yin.EMF(),Yslope.EMF());

//       if (Input::List().hydromotion) HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));

        BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
        // Yslope.DF(s).checknan();std::cout << "Vlasov 4 \n";
//       FA[s](Yin.EMF(),Yslope.EMF());

        // Yslope.DF(s).checknan();std::cout << "Vlasov 5 \n";
        Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor1D_f1_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_f1_implicitE_p1::operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope){}
//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_f1_implicitEB_p1::VlasovFunctor1D_f1_implicitEB_p1(vector<size_t> Nl,vector<size_t> Nm,vector<double> pmax, vector<size_t> Np,
                                                             double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );

        SA.push_back( Spatial_Advection_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        // EF.push_back( Electric_Field_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

        // JX.push_back( Current_1D(pmin, pmax[s], Np[s], Nx) );

        // AM.push_back( Ampere_1D(xmin, xmax, Nx) );

//        HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

//        BF.push_back( Magnetic_Field_1D(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

//        FA.push_back( Faraday_1D(xmin, xmax, Nx) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor1D_f1_implicitEB_p1::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        // Yslope.DF(s).checknan();std::cout << "Vlasov 1 \n";
        SA[s](Yin.DF(s),Yslope.DF(s));
        // Yslope.DF(s).checknan();std::cout << "Vlasov 2 \n";

        // Yslope.DF(s).checknan();std::cout << "Vlasov 3 \n";

        // EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

        // AM[s](Yin.EMF(),Yslope.EMF());

//       if (Input::List().hydromotion) HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));

//        BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
        // Yslope.DF(s).checknan();std::cout << "Vlasov 4 \n";
//        FA[s](Yin.EMF(),Yslope.EMF());

        // Yslope.DF(s).checknan();std::cout << "Vlasov 5 \n";
        Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor1D_f1_implicitEB_p1::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
void VlasovFunctor1D_f1_implicitEB_p1::operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope){}
//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_f1_implicitE_p2::VlasovFunctor1D_f1_implicitE_p2(vector<size_t> Nl,vector<size_t> Nm,vector<double> pmax, vector<size_t> Np,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );
        FA.push_back( Faraday_1D(xmin, xmax, Nx) );
        EF.push_back( Electric_Field_1D_f1(Nl[s], Nm[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );

    }

}


//------------------------------------------------------------------------------------------------------
void VlasovFunctor1D_f1_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope){
// //---------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        FA[s](Yin.EMF(),Yslope.EMF());
        EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));
        Yslope.DF(s) = Yslope.DF(s).Filterp();

    }

}

// //--------------------------------------------------------------------------------------------------
void VlasovFunctor1D_f1_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope, size_t direction){
// //--------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    if (direction == 1)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            EF[s].Implicit_Ex(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else if (direction == 2)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            EF[s].Implicit_Ey(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else
    {
        // complex<double> ii(0.0,1.0);
        for (size_t s(0); s < Yin.Species(); ++s) {
            EF[s].Implicit_Ez(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }

}
void VlasovFunctor1D_f1_implicitE_p2::operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope){}
// //**************************************************************
