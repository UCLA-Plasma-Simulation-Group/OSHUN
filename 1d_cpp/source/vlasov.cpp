/*! \brief  Vlasov Equation - Definitions
 * \author PICKSC
 * \file   vlasov.cpp
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



//**************************************************************
//--------------------------------------------------------------
//  Current

Current_1D::Current_1D( double pmin, double pmax, size_t Np,
                        size_t Nx )
        : Jx(Nx), Jy(Nx), Jz(Nx), small(0.5*pmax/Np) {
//          pr(Algorithms::MakeAxis(complex<double>(pmin),
//                                  complex<double>(pmax),
//                                  Np)),
//          invg(pr) {
    small *= small; small *= small; small *= 0.5*pmax/Np; 
    small *= 0.2; small *= 1.0/(1.5*pmax/Np); 
//    for (size_t i(0); i < pr.size(); ++i) {
//        invg[i] = 1.0 / (sqrt(1.0+pr[i]*pr[i]));
//    }
};
//--------------------------------------------------------------

void Current_1D::operator()(const DistFunc1D& Din, Field1D& Exh, Field1D& Eyh, Field1D& Ezh) {

    Array2D<double> temp(3,Din(0).numx());

    temp = Din.getcurrent();

    for (size_t i(0); i < Jx.numx(); ++i) {
        Jx(i) = complex<double >(temp(0,i));
        Jy(i) = complex<double >(temp(1,i));
        Jz(i) = complex<double >(temp(2,i));
    }

    Exh += Jx;
    Eyh += Jy;
    Ezh += Jz;

}

void Current_1D::es1d(const DistFunc1D& Din, Field1D& Exh) {

    valarray<double> temp(Din(0).numx());

    temp = Din.getcurrent(0);

    for (size_t i(0); i < Jx.numx(); ++i) {
        Jx(i) = complex<double >(temp[i]);//+4.0/3.0*3.1415926*small*Din(1,0)(1,i);
    }

    Exh += Jx;

}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Electric_Field_1D::Electric_Field_1D(size_t Nl, size_t Nm,
                                     double pmin, double pmax, size_t Np,
                                     double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : A1(Nl+1,Nm+1), A2(Nl+1,Nm+1),
          B1(Nl+1), B2(Nl+1),
          C1(Nl+1), C2(Nl+1,Nm+1), C3(Nl+1), C4(Nl+1,Nm+1),
          Hp0(Nl+1),
          H(Np,Nx), G(Np,Nx), TMP(Np,Nx),
         // pr(Algorithms::MakeAxis(complex<double>(complex<double>(pmax/2.0/Np)),
         //                         complex<double>(pmax),
         //                         Np)),

          pr(Algorithms::MakeCAxis(complex<double>(0.0),
                                  complex<double>(pmax),
                                  Np)),
          invpr(pr)
{
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -
    complex<double> lc, mc;

//       Inverted momentum axis
    for (size_t i(0); i < pr.size(); ++i) {
        invpr[i] = 1.0/pr[i];
    }
    double idp = (-1.0)/ (2.0*(pmax-pmin)/double(Np)); // -1/(2dp)
    // double idp = (-1.0)/ (2.0*(pmax-pmax/2.0/double(Np))/double(Np-1));

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Calculate the A1 * l/(l+1) * 1/(2Dp), A2 * (-1)/(2Dp) parameters
    for (size_t l(1); l < Nl+1; ++l){
        for (size_t m(0); m<((Nm<l)?Nm:l)+1; ++m){
            lc = complex<double>(l);
            mc = complex<double>(m);
            A1(l,m) = idp *(-1.0) *  (lc+1.0-mc)/(2.0*lc+1.0)  * lc/(lc+1.0);
            A2(l,m) = idp *              (lc+mc)/(2.0*lc+1.0);
        }
    }
    A1(0,0) = idp;
    // A2[0] is not used

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1, B2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(1); l<Nl+1; ++l){
        lc = complex<double>(l);
        B1[l] = idp * lc* lc     /(2.0*lc+1.0);
        B2[l] = idp * lc*(lc+1.0)/(2.0*lc+1.0);
    }
    B2[0] = 1.0;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C1, C3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(1); l<Nl+1; ++l){
        lc = complex<double>(l);
        C1[l] = (-0.5) * idp * lc /((2.0*lc+1.0)*(lc+1.0));
        C3[l] = (-0.5) * idp /(2.0*lc+1.0);
    }
    C1[0] = 0.5 * idp;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C2, C4" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(2); l<Nl+1; ++l){
        for (size_t m(2); m<((Nm<l)?Nm:l)+1; ++m){
            lc = complex<double>(l);
            mc = complex<double>(m);
            C2(l,m) = (0.5) * idp * lc * (lc-mc+2.0)*(lc-mc+1.0)/((2.0*lc+1.0)*(lc+1.0));
            C4(l,m) = (0.5) * idp * (lc+mc-1.0)*(lc+mc)/(2.0*lc+1.0);
        }
    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       -2*Dp*H at the 0 momentum cell
    Hp0[0] = 1.0 / pr[0];
    for (size_t l(1); l < Nl+1; ++l) {
        double ld(l);
        Hp0[l] = Hp0[l-1] * (pr[0]/pr[1]) * (2.0*ld+1.0)/(2.0*ld-1.0);
    }
    Hp0 *= (-2.0)*(pr[1]-pr[0]);

    A100 = complex<double>(idp);
    C100 = complex<double>(0.5*idp);
    A210 = complex<double>(1.0/3.0*idp);
    B211 = complex<double>(2.0/3.0);
    A310 = complex<double>(2.0/5.0*idp);
    C311 = complex<double>(-0.5*idp/5.0);

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Electric_Field_1D::operator()(const DistFunc1D& Din,
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
    Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxaxis(Ex);
    Em *= C1[0];            Dh(1,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(1,0),1);
    Ex *= A2(1,0) / A1(0,0);          Dh(0,0) += H.mxaxis(Ex);
    Ex *= A1(1,0) / A2(1,0); TMP = G; Dh(2,0) += G.mxaxis(Ex);
    Em *= C1[1]   / C1[0]  ;          Dh(2,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, 1 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(2); l < l0; ++l){
        MakeGH(Din(l,0),l);
        Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Dh(l-1,0) += H.mxaxis(Ex);
        Em *= C3[l]   / C1[l-1];             Dh(l-1,1) += TMP.mxaxis(Em);
        Ex *= A1(l,0) / A2(l,0);   TMP = G;  Dh(l+1,0) += G.mxaxis(Ex);
        Em *= C1[l]   / C3[l];               Dh(l+1,1) += TMP.mxaxis(Em);
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0,  l = l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(l0,0),l0);
    Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Dh(l0-1,0) += H.mxaxis(Ex);
    Em *= C3[l0]   / C1[l0-1];             Dh(l0-1,1) += TMP.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(1,1),1);
    Ep *= B2[1];              H = H.mxaxis(Ep);     Dh(0,0) += H.Re();
    Ex *= A1(1,1) / A2(l0,0); TMP = G;              Dh(2,1) += G.mxaxis(Ex);
    Em *= C1[1] / C3[l0];     G = TMP;              Dh(2,2) += TMP.mxaxis(Em);
    Ep *= B1[1] / B2[1];      G = G.mxaxis(Ep);     Dh(2,0) += G.Re();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 2
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(2,1),2);
    Ex *= A2(2,1) / A1(1,1); TMP = H;               Dh(1,1) += TMP.mxaxis(Ex);
    Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep);      Dh(1,0) += H.Re();
    Ex *= A1(2,1) / A2(2,1); TMP = G;               Dh(3,1) += G.mxaxis(Ex);
    Em *= C1[2]   / C1[1];   G = TMP;               Dh(3,2) += TMP.mxaxis(Em);
    Ep *= B1[2]   / B2[2];   G = G.mxaxis(Ep);      Dh(3,0) += G.Re();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, 1 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(3); l < l0; ++l){
        MakeGH(Din(l,1),l);
        Ex *= A2(l,1) / A1(l-1,1); TMP = H;             Dh(l-1,1) += H.mxaxis(Ex);
        Em *= C3[l]   / C1[l-1];   H = TMP;             Dh(l-1,2) += TMP.mxaxis(Em);
        Ep *= B2[l]   / B1[l-1];   H = H.mxaxis(Ep);    Dh(l-1,0) += H.Re();
        Ex *= A1(l,1) / A2(l,1);   TMP = G;             Dh(l+1,1) += G.mxaxis(Ex);
        Em *= C1[l]   / C3[l];     G = TMP;             Dh(l+1,2) += TMP.mxaxis(Em);
        Ep *= B1[l]   / B2[l];     G = G.mxaxis(Ep);    Dh(l+1,0) += G.Re();
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       m = 1,  l = l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(l0,1),l0);
    Ex *= A2(l0,1) / A1(l0-1,1); TMP = H;              Dh(l0-1,1) += H.mxaxis(Ex);
    Em *= C3[l0]   / C1[l0-1];   H = TMP;              Dh(l0-1,2) += TMP.mxaxis(Em);
    Ep *= B2[l0]   / B1[l0-1];   H = H.mxaxis(Ep);     Dh(l0-1,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    C4(l0,1) = B2[l0];
    for (size_t m(2); m < m0; ++m){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//          m > 1 , l = m
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(m,m),m);
        Ep *= C4(m,m) / C4(l0,m-1);              Dh(m-1,m-1) += H.mxaxis(Ep);
        Ex *= A1(m,m) / A2(l0,m-1); TMP = G;     Dh(m+1,m  ) += G.mxaxis(Ex);
        Em *= C1[m]   / C3[l0];     G = TMP;     Dh(m+1,m+1) += TMP.mxaxis(Em);
        Ep *= C2(m,m) / C4(m,m);                 Dh(m+1,m-1) += G.mxaxis(Ep);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//          m > 1 , l = m+1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(m+1,m),m+1);
        Ex *= A2(m+1,m) / A1(m,m);   TMP = H;     Dh(m  ,m  ) += TMP.mxaxis(Ex);
        Ep *= C4(m+1,m) / C2(m,m);                Dh(m  ,m-1) += H.mxaxis(Ep);
        if ( m+1 < l0) { //always true except when m = m0-1 = l0-1
            Ex *= A1(m+1,m) / A2(m+1,m); TMP = G;     Dh(m+2,m  ) += G.mxaxis(Ex);
            Em *= C1[m+1]   / C1[m];     G = TMP;     Dh(m+2,m+1) += TMP.mxaxis(Em);
            Ep *= C2(m+1,m) / C4(m+1,m);              Dh(m+2,m-1) += G.mxaxis(Ep);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//              m > 1, 1 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t l(m+2); l < l0; ++l){
                MakeGH(Din(l,m),l);
                Ex *= A2(l,m) / A1(l-1,m); TMP = H;   Dh(l-1,m  ) += H.mxaxis(Ex);
                Em *= C3[l]   / C1[l-1];   H = TMP;   Dh(l-1,m+1) += TMP.mxaxis(Em);
                Ep *= C4(l,m) / C2(l-1,m);            Dh(l-1,m-1) += H.mxaxis(Ep);
                Ex *= A1(l,m) / A2(l,m);   TMP = G;   Dh(l+1,m  ) += G.mxaxis(Ex);
                Em *= C1[l]   / C3[l];     G = TMP;   Dh(l+1,m+1) += TMP.mxaxis(Em);
                Ep *= C2(l,m) / C4(l,m);              Dh(l+1,m-1) += G.mxaxis(Ep);
            }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//               m > 1,  l = l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeGH(Din(l0,m),l0);
            Ex *= A2(l0,m) / A1(l0-1,m); TMP = H;    Dh(l0-1,m  ) += H.mxaxis(Ex);
            Em *= C3[l0]   / C1[l0-1];   H = TMP;    Dh(l0-1,m+1) += TMP.mxaxis(Em);
            Ep *= C4(l0,m) / C2(l0-1,m);             Dh(l0-1,m-1) += H.mxaxis(Ep);
        }
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(m0,m0),m0);
    Ep *= C4(m0,m0) / C4(l0,m0-1);              Dh(m0-1,m0-1) += H.mxaxis(Ep);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       m = m0, l0 > l > m0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ( m0 < l0) {
        Ex *= A1(m0,m0) / A2(l0,m0-1); TMP = G;  Dh(m0+1,m0  ) += TMP.mxaxis(Ex);
        Ep *= C2(m0,m0) / C4(m0,m0);             Dh(m0+1,m0-1) += G.mxaxis(Ep);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//          m = m0 , l = m0+1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(m0+1,m0),m0+1);
        Ex *= A2(m0+1,m0) / A1(m0,m0); TMP = H;        Dh(m0,m0  )  += TMP.mxaxis(Ex);
        Ep *= C4(m0+1,m0) / C2(m0,m0);                 Dh(m0,m0-1)  += H.mxaxis(Ep);
        if ( m0+1 < l0) {
            Ex *= A1(m0+1,m0) / A2(m0+1,m0); TMP = G;  Dh(m0+2,m0  )+= TMP.mxaxis(Ex);
            Ep *= C2(m0+1,m0) / C4(m0+1,m0);           Dh(m0+2,m0-1)+= G.mxaxis(Ep);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//              m = m0, m0+2 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t l(m0+2); l < l0; ++l){
                MakeGH(Din(l,m0),l);
                Ex *= A2(l,m0) / A1(l-1,m0); TMP = H;  Dh(l-1,m0)   += TMP.mxaxis(Ex);
                Ep *= C4(l,m0) / C2(l-1,m0);           Dh(l-1,m0-1) += H.mxaxis(Ep);
                Ex *= A1(l,m0) / A2(l,  m0); TMP = G;  Dh(l+1,m0  ) += TMP.mxaxis(Ex);
                Ep *= C2(l,m0) / C4(l,  m0);           Dh(l+1,m0-1) += G.mxaxis(Ep);
            }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//               m > 1,  l = l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeGH(Din(l0,m0),l0);
            Ex *= A2(l0,m0) / A1(l0-1,m0); TMP = H;    Dh(l0-1,m0)   += TMP.mxaxis(Ex);
            Ep *= C4(l0,m0) / C2(l0-1,m0);             Dh(l0-1,m0-1) += H.mxaxis(Ep);
        }
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field_1D::es1d(const DistFunc1D& Din,
                                   const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                                   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

    complex<double> ii(0.0,1.0);

    valarray<complex<double> > Ex(FEx.array());
    // valarray<complex<double> > Em(FEz.array());
    // Em *= (-1.0)*ii;
    // Em += FEy.array();
    // valarray<complex<double> > Ep(FEz.array());
    // Ep *= ii;
    // Ep += FEy.array();

    Ex *= Din.q();
    // Em *= Din.q();
    // Ep *= Din.q();

    size_t l0(Din.l0());
    // size_t m0(Din.m0());


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeG00(Din(0,0));
    Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxaxis(Ex);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(1,0),1);
    Ex *= A2(1,0) / A1(0,0);          Dh(0,0) += H.mxaxis(Ex);
    Ex *= A1(1,0) / A2(1,0); TMP = G; Dh(2,0) += G.mxaxis(Ex);


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, 1 < l < l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(2); l < l0; ++l){
        MakeGH(Din(l,0),l);
        Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Dh(l-1,0) += H.mxaxis(Ex);

        Ex *= A1(l,0) / A2(l,0);   TMP = G;  Dh(l+1,0) += G.mxaxis(Ex);

    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0,  l = l0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(l0,0),l0);
    Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Dh(l0-1,0) += H.mxaxis(Ex);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


}
//--------------------------------------------------------------



//--------------------------------------------------------------
void Electric_Field_1D::Implicit_Ex(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

    valarray<complex<double> > Ex(FEx.array());
    Ex *= Din.q();

//      m = 0, l = 0
    MakeG00(Din(0,0));
    Ex *= A1(0,0);             Dh(1,0) += G.mxaxis(Ex);

//      m = 0, l = 1
    MakeGH(Din(1,0),1);
    Ex *= A2(1,0) / A1(0,0);   Dh(0,0) += H.mxaxis(Ex);
    Ex *= A1(1,0) / A2(1,0);   Dh(2,0) += G.mxaxis(Ex);

//      m = 0, l = 2
    MakeGH(Din(2,0),2);
    Ex *= A2(2,0) / A1(1,0);   Dh(1,0) += H.mxaxis(Ex);

//      m = 0, l = 3
    MakeGH(Din(3,0),3);
    Ex *= A2(3,0) / A2(2,0);   Dh(2,0) += H.mxaxis(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
    MakeGH(Din(1,1),1);
    Ex *= A1(1,1) / A2(3,0);   Dh(2,1) += G.mxaxis(Ex);

//      m = 1, l = 2
    MakeGH(Din(2,1),2);
    Ex *= A2(2,1) / A1(1,1);   Dh(1,1) += H.mxaxis(Ex);

//      m = 1, l = 3
    MakeGH(Din(3,1),3);
    Ex *= A2(3,1) / A2(2,1);   Dh(2,1) += H.mxaxis(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field_1D::Implicit_Ey(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh) {
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
    Em *= C1[0];            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
    MakeGH(Din(1,0),1);
    Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
    MakeGH(Din(2,0),2);
    Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxaxis(Em);

//      m = 0,  l = 3
    MakeGH(Din(3,0),3);
    Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
    MakeGH(Din(1,1),1);
    Ep *= B2[1];             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
    Em *= C1[1] / C3[3];     TMP = G;          Dh(2,2) += TMP.mxaxis(Em);
    Ep *= B1[1] / B2[1];     G = G.mxaxis(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
    MakeGH(Din(2,1),2);
    Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
    MakeGH(Din(3,1),3);
    Em *= C3[3]   / C1[1];   TMP = H;          Dh(2,2) += TMP.mxaxis(Em);
    Ep *= B2[3]   / B2[2];   H = H.mxaxis(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
    MakeGH(Din(2,2),2);
    Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxaxis(Ep);

//      m = 2, l = 3
    MakeGH(Din(3,2),3);
    Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxaxis(Ep);

    size_t m0(Din.m0());
    if ( m0 > 2) {
//          m = 3, l = 3
        MakeGH(Din(3,3),3);
        Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxaxis(Ep);
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
//
////--------------------------------------------------------------
void Electric_Field_1D::Implicit_Ez(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------
    complex<double> ii(0.0,1.0);


    valarray<complex<double> > Em(FEz.array());
    Em *= (-1.0)*ii;


    valarray<complex<double> > Ep(FEz.array());
    Ep *= ii;
    // Ep += FEy.array();
    //
    // Ex *= Din.q();;
    Em *= Din.q();
    Ep *= Din.q();



//      m = 0, l = 0
    MakeG00(Din(0,0));
    Em *= C1[0];            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
    MakeGH(Din(1,0),1);
    Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
    MakeGH(Din(2,0),2);
    Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxaxis(Em);

//      m = 0,  l = 3
    MakeGH(Din(3,0),3);
    Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
    MakeGH(Din(1,1),1);
    Ep *= B2[1];             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
    Em *= C1[1] / C3[3];     TMP = G;              Dh(2,2) += TMP.mxaxis(Em);
    Ep *= B1[1] / B2[1];     G = G.mxaxis(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
    MakeGH(Din(2,1),2);
    Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
    MakeGH(Din(3,1),3);
    Em *= C3[3]   / C1[1];   TMP = H;              Dh(2,2) += TMP.mxaxis(Em);
    Ep *= B2[3]   / B2[2];   H = H.mxaxis(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
    MakeGH(Din(2,2),2);
    Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxaxis(Ep);

//      m = 2, l = 3
    MakeGH(Din(3,2),3);
    Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxaxis(Ep);

    size_t m0(Din.m0());
    if ( m0 > 2) {
//          m = 3, l = 3
        MakeGH(Din(3,3),3);
        Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxaxis(Ep);
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field_1D::f1only(const DistFunc1D& Din,
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

    SHarmonic1D f20(Din(0,0));
    f20 *= complex<double>(1.0/3.0);

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

//     MakeGH(f20,2);
//     Ex *= A310 / A210;      TMP = H;      Dh(1,0) += H.mxaxis(Ex);
//     Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    MakeGH(Din(1,1),1);
    Ep *= B211;              H = H.mxaxis(Ep);     Dh(0,0) += H.Re();

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Electric_Field_1D::Implicit_Ex_f1only(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

    valarray<complex<double> > Ex(FEx.array());
    Ex *= Din.q();

    SHarmonic1D f20(Din(0,0));
    f20 *= complex<double>(1.0/3.0);

//      m = 0, l = 0
    MakeG00(Din(0,0));
    Ex *= A100;             Dh(1,0) += G.mxaxis(Ex);

//      m = 0, l = 1
    MakeGH(Din(1,0),1);
    Ex *= A210 / A100;   Dh(0,0) += H.mxaxis(Ex);
    // Ex *= A1(1,0) / A2(1,0);   Dh(2,0) += G.mxaxis(Ex);

//      m = 0, l = 2
//     MakeGH(f20,2);
//     Ex *= A310 / A210;      TMP = H;      Dh(1,0) += H.mxaxis(Ex);


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field_1D::Implicit_Ey_f1only(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

    valarray<complex<double> > Em(FEy.array());
    valarray<complex<double> > Ep(FEy.array());

    SHarmonic1D f20(Din(0,0));
    f20 *= complex<double>(1.0/3.0);

    Em *= Din.q();
    Ep *= Din.q();

//      m = 0, l = 0
    MakeG00(Din(0,0));
    Em *= C100;            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
    // MakeGH(Din(1,0),1);
    // Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
//     MakeGH(f20,2);
//     Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

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
void Electric_Field_1D::Implicit_Ez_f1only(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------
    complex<double> ii(0.0,1.0);
    SHarmonic1D f20(Din(0,0));
    f20 *= complex<double>(1.0/3.0);

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
//    Din(2,0)  = Din(0,0);
//    Din(2,0) *= complex<double>(1.0/3.0);
//
//    MakeGH(f20,2);
//    Em *= C311 / C100;                    Dh(1,1) += TMP.mxaxis(Em);

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
void Electric_Field_1D::MakeGH(SHarmonic1D& f, size_t el){
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
void Electric_Field_1D::MakeG00(SHarmonic1D& f) {
//--------------------------------------------------------------
    G = f; G = G.Dp();

    complex<double> p0p1_sq( pr[0]*pr[0]/(pr[1]*pr[1]) ),
            inv_mp0p1_sq( 1.0/(1.0-p0p1_sq) ),
            g_r = -4.0*(pr[1]-pr[0]) * pr[0]/(pr[1]*pr[1]),
            // g_r = 2.0*(pr[0]/pr[1]/pr[1]),
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
Magnetic_Field_1D::Magnetic_Field_1D(size_t Nl, size_t Nm,
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
        mc = complex<double>(m);
        A1[m] = (-1.0)*c01*mc;
    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(1); l < Nl+1; ++l)
        for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
            lc = complex<double>(l);
            mc = complex<double>(m);
            A2(l,m) = (-0.5)*(lc+1.0-mc)*(lc+mc);
        }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    A3 = 0.5;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(0); l < Nl+1; ++l){
        lc = complex<double>(l);
        B1[l] = (-1.0)*lc*(lc+1.0);
    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Magnetic_Field_1D::operator()(const DistFunc1D& Din,
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
    FLM = Din(1,1); Bm *= B1[1]; FLM = FLM.mxaxis(Bm); Dh(1,0) += FLM.Re();

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l > 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t l(2); l < l0+1; ++l){
        FLM = Din(l,1);                                                Dh(l,2) += FLM.mxaxis(Bp);
        FLM = Din(l,1);                                                Dh(l,1) += FLM.mxaxis(Bx);
        FLM = Din(l,1); Bm *= B1[l]/B1[l-1]; FLM = FLM.mxaxis(Bm); Dh(l,0) += FLM.Re();
    }
    Bm *= 1.0/B1[l0];

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m > 1, l = m
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (size_t m(2); m < m0; ++m){
        FLM = Din(m,m); Bx *= A1[m]/A1[m-1];                   Dh(m,m  ) += FLM.mxaxis(Bx);
        FLM = Din(m,m); Bm *= A2(m,m);                         Dh(m,m-1) += FLM.mxaxis(Bm);
        for (size_t l(m+1); l < l0+1; ++l){
            FLM = Din(l,m);                                    Dh(l,m+1) += FLM.mxaxis(Bp);
            FLM = Din(l,m);                                    Dh(l,m  ) += FLM.mxaxis(Bx);
            FLM = Din(l,m); Bm *= A2(l,m)/A2(l-1,m);           Dh(l,m-1) += FLM.mxaxis(Bm);
        }
        Bm *= 1.0/A2(l0,m);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = m0, l >= m0
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FLM = Din(m0,m0); Bx *= A1[m0]/A1[m0-1];                   Dh(m0,m0)   += FLM.mxaxis(Bx);
    FLM = Din(m0,m0); Bm *= A2(m0,m0)/*/A2(l0,m0-1)*/;             Dh(m0,m0-1) += FLM.mxaxis(Bm);
    for (size_t l(m0+1); l < l0+1; ++l){
        FLM = Din(l,m0);                                     Dh(l,m0  )  += FLM.mxaxis(Bx);
        FLM = Din(l,m0); Bm *= A2(l,m0)/A2(l-1,m0);          Dh(l,m0-1)  += FLM.mxaxis(Bm);
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Magnetic_Field_1D::implicit(DistFunc1D& Din,
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
//--------------------------------------------------------------
void Magnetic_Field_1D::f1only(const DistFunc1D& Din,
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
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Spatial_Advection_1D::Spatial_Advection_1D(size_t Nl, size_t Nm,
                                           double pmin, double pmax, size_t Np,
                                           double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : A1(Nl+1,Nm+1), A2(Nl+1,Nm+1),
          fd1(Np,Nx), fd2(Np,Nx),
         // vr(Algorithms::MakeAxis(complex<double>(pmax/2.0/Np),
         //                         complex<double>(pmax),
         //                         Np)) {

          vr(Algorithms::MakeCAxis(complex<double>(0.0),
                                  complex<double>(pmax),
                                  Np)){
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -

    complex<double> lc, mc;

    for (size_t i(0); i < vr.size(); ++i) {
        vr[i] = vr[i]/(sqrt(1.0+vr[i]*vr[i]));
    }

    double idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx)); // -1/(2dx)


//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Calculate the "A1, A2" parameters
    for (size_t l(0); l < Nl+1; ++l){
        for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
            lc = complex<double>(l);
            mc = complex<double>(m);
            A1(l,m) = idx *(-1.0) * (lc-mc+1.0) / (2.0*lc+1.0);
            A2(l,m) = idx *(-1.0) * (lc+mc)     / (2.0*lc+1.0);
        }
    }
    A2(0,0) = 1.0;

    //       Calculate the "A1, A2" parameters
    A00 = complex<double>(idx);
    A10 = complex<double>(idx/3.0);
    A20 = complex<double>(idx*2.0/5.0);



}
//--------------------------------------------------------------

//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection_1D::operator()(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vt(vr); vt *= 1.0/Din.mass();

    size_t l0(Din.l0());
    size_t m0(Din.m0());

    for (size_t m(0); m < ((m0<l0)?(m0+1):m0); ++m){

        fd1 = Din(m,m);     fd1 = fd1.Dx();

        vt *= A1(m,m);                          Dh(m+1,m) += fd1.mpaxis(vt);

        for (size_t l(m+1); l < l0; ++l) {
            fd1 = Din(l,m);            fd1 = fd1.Dx();

            vt *= A2(l,m)/A1(l-1,m);    fd2 = fd1;  Dh(l-1,m) += fd1.mpaxis(vt);
            vt *= A1(l,m)/A2(l  ,m);                Dh(l+1,m) += fd2.mpaxis(vt);
        }

        fd1 = Din(l0,m);            fd1 = fd1.Dx();
        vt *= A2(l0,m)/A1(l0-1,m);                 Dh(l0-1,m) += fd1.mpaxis(vt);
        vt *= 1.0/A2(l0,m);
    }


}
//--------------------------------------------------------------
//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection_1D::es1d(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vt(vr); vt *= 1.0/Din.mass();

    size_t l0(Din.l0());

    fd1 = Din(0,0);     fd1 = fd1.Dx();

    vt *= A1(0,0);                          Dh(1,0) += fd1.mpaxis(vt);

    for (size_t l(1); l < l0; ++l) {
        fd1 = Din(l,0);  //std::cout << "\n \n before dx, l = " << l << " \n";          
        fd1 = fd1.Dx(); //std::cout << " \n after dx\n";

        vt *= A2(l,0)/A1(l-1,0);    fd2 = fd1;  Dh(l-1,0) += fd1.mpaxis(vt);
        vt *= A1(l,0)/A2(l  ,0);                Dh(l+1,0) += fd2.mpaxis(vt);
    }
    
    fd1 = Din(l0,0);            fd1 = fd1.Dx();
    vt *= A2(l0,0)/A1(l0-1,0);                 Dh(l0-1,0) += fd1.mpaxis(vt);
    // vt *= 1.0/A2(l0,0);
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection_1D::f1only(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vt(vr); vt *= 1.0/Din.mass();

    fd1 = Din(0,0);
    fd1 = fd1.Dx();     vt *= A00;
    Dh(1,0) += (fd1.mpaxis(vt));

    fd1 = Din(1,0);
    fd1 = fd1.Dx();     vt *= A10/A00;
    Dh(0,0) += (fd1.mpaxis(vt));

//    fd1 = Din(0,0);     fd1 *= complex<double>(1.0/3.0);
//    fd1 = fd1.Dx();     vt *= A20/A10;
//    Dh(1,0) += (fd1.mpaxis(vt));

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
//  Update B with E term from Faraday's Law
//-------------------------------------------------------------- 
Faraday_1D::Faraday_1D( double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
// Constructor
//--------------------------------------------------------------
        : tmpE(Nx), numx(Nx) {

    idx = complex<double>((-1.0)/ (2.0*(xmax-xmin)/double(Nx))); // -1/(2dx)

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Faraday_1D::operator()(EMF1D& EMFin, EMF1D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Faraday's Law 
//--------------------------------------------------------------

//      dBx/dt += - dEz/dy  
//         tmpE          = Fin.Ez(); 
//         tmpE         *= (-1.0) * idy;
//         Yh.EMF().Bx() += tmpE.Dy();

//      dBy/dt +=   dEz/dx       
    tmpE                 = EMFin.Ez();
    tmpE                *= idx;
    EMFh.By()           += tmpE.Dx();
    // EMFh.By()(numx-1)    = 0.0;

//      dBz/dt +=   dEx/dy       
//         tmpE          = Ein.Ex(); 
//         tmpE         *= idy;
//         Yh.EMF().Bz() += tmpE.Dy();    

//      dBz/dt += - dEy/dx       
    tmpE                 = EMFin.Ey();
    tmpE                *= (-1.0) * idx;
    EMFh.Bz()           += tmpE.Dx();

    // EMFh.Bz()(numx-1)    = 0.0;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Update E with B term from Faraday's Law
//--------------------------------------------------------------
Ampere_1D::Ampere_1D( double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
// Constructor
//--------------------------------------------------------------
        : tmpB(Nx), numx(Nx) {

    idx = complex<double>((-1.0)/ (2.0*(xmax-xmin)/double(Nx))); // -1/(2dx)

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Ampere_1D::operator()(EMF1D& EMFin, EMF1D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Ampere's Law 
//--------------------------------------------------------------

//      dEx/dt +=   dBz/dy       
//         tmpB          = Fin.Bz(); 
//         tmpB         *= idy;
//         Fh.Ex() += tmpB.Dy();

//      dEy/dt +=  - dBz/dx       
    tmpB                 = EMFin.Bz();
    tmpB                *= (-1.0) * idx;
    EMFh.Ey()           += tmpB.Dx();

    //.Dx();
    // EMFh.Ey()(numx-1)    = 0.0;

//      dEz/dt +=  - dBx/dy       
//         tmpB          = Fin.Bx(); 
//         tmpB         *= (-1.0) * idy;
//         Fh.Ez() += tmpB.Dy();    

//      dEz/dt += dBy/dx       
    tmpB                 = EMFin.By();
    tmpB                *=  idx;
    EMFh.Ez()           += tmpB.Dx();
    // EMFh.Ez()(numx-1)    = 0.0;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//**************************************************************
