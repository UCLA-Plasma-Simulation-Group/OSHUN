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
#include <omp.h>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
// #include "nmethods.h"

//  Declerations
#include "state.h"
#include "formulary.h"
#include "input.h"
#include "fluid.h"
#include "vlasov.h"

//**************************************************************
//--------------------------------------------------------------
//  Current

Current::Current(){};
//--------------------------------------------------------------

void Current::operator()(const DistFunc1D& Din, Field1D& Exh, Field1D& Eyh, Field1D& Ezh) {

    Array2D<double> temp(3,Din(0).numx());

    if (Input::List().relativity)   temp = Din.getrelativisticcurrent();
    else                            temp = Din.getcurrent();

    for (size_t i(0); i < Exh.numx(); ++i) {
        Exh(i) += static_cast<complex<double> >(temp(0,i));
        Eyh(i) += static_cast<complex<double> >(temp(1,i));
        Ezh(i) += static_cast<complex<double> >(temp(2,i));
    }

}
void Current::operator()(const DistFunc2D& Din, Field2D& Exh, Field2D& Eyh, Field2D& Ezh) {

    Array3D<double> temp(3,Din(0).numx(),Din(0).numy());

    if (Input::List().relativity)   temp = Din.getrelativisticcurrent();
    else                            temp = Din.getcurrent();

    for (size_t ix(0); ix < Exh.numx(); ++ix) 
    {
        for (size_t iy(0); iy < Exh.numy(); ++iy) 
        {
            Exh(ix,iy) += static_cast<complex<double> >(temp(0,ix,iy));
            Eyh(ix,iy) += static_cast<complex<double> >(temp(1,ix,iy));
            Ezh(ix,iy) += static_cast<complex<double> >(temp(2,ix,iy));
        }
    }

}
void Current::es1d(const DistFunc1D& Din, Field1D& Exh) {

    valarray<double> temp(Din(0).numx());

    if (Input::List().relativity)   temp = Din.getrelativisticcurrent(0);
    else                            temp = Din.getcurrent(0);

    for (size_t i(0); i < Exh.numx(); ++i) {
        Exh(i) += static_cast<complex<double> >(temp[i]);
    }

}
//--------------------------------------------------------------
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Electric_Field::Electric_Field(size_t Nl, size_t Nm, valarray<double> dp)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
:   A1(Nl+1,Nm+1), A2(Nl+1,Nm+1),
    B1(Nl+1), B2(Nl+1),
    C1(Nl+1), C3(Nl+1), C2(Nl+1,Nm+1), C4(Nl+1,Nm+1),
    Hp0(Nl+1),

    pr(Algorithms::MakeCAxis(complex<double>(0.0),
      complex<double>(1.0),
      dp.size())),
    invdp(Algorithms::MakeCAxis(complex<double>(0.0),
      complex<double>(1.0),
      dp.size())),
    invpr(pr),
    f_start(Input::List().ompthreads),f_end(Input::List().ompthreads),
    dist_il((Nm+1)*(2*Nl-Nm+2)/2),dist_im((Nm+1)*(2*Nl-Nm+2)/2),
    nwsediag_il((Nm+1)*(2*Nl-Nm+2)/2),nwsediag_im((Nm+1)*(2*Nl-Nm+2)/2),
    neswdiag_il((Nm+1)*(2*Nl-Nm+2)/2),neswdiag_im((Nm+1)*(2*Nl-Nm+2)/2)
{
    //      - - - - - - - - - - - - - - - - - - - - - - - - - - -
    complex<double> lc, mc;

    // ------------------------------------------------------------------------ // 
    // Non-uniform velocity grids
    // ------------------------------------------------------------------------ //  
    pr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        pr[ip] = pr[ip-1] + dp[ip];
            // std::cout << "\n Epr[" << ip << "] = " << pr[ip] << std::endl;
    }
    // ------------------------------------------------------------------------ // 


    // ------------------------------------------------------------------------ // 
    // Since .Dp() calculates df/dx[i]  =  -(f[i+1] - f[i-1])
    // We still need to divide by delta x
    // for a nonuniform, cell-centered, center difference scheme
    // this is -(1/(x[i]-x[i-1]))
    // ------------------------------------------------------------------------ // 

    invdp[0] =  -1.0/((dp[1]+dp[0]));
    for (size_t i(1); i < pr.size()-1; ++i)
    {
        invdp[i] = -1.0/(dp[i]+dp[i+1]);
    }
    invdp[pr.size()-1] = -1.0/((dp[pr.size()-1]+dp[pr.size()-2]));
    // ------------------------------------------------------------------------ // 
    
    
    
    // ------------------------------------------------------------------------ // 
    //       Inverted momentum axis
    for (size_t i(0); i < pr.size(); ++i) invpr[i] = 1.0/pr[i];
    // ------------------------------------------------------------------------ // 

    // ------------------------------------------------------------------------ // 
    //       Calculate the A1 * -l/(l+1), A2 parameters
    // ------------------------------------------------------------------------ // 
        for (size_t l(1); l < Nl+1; ++l){
            for (size_t m(0); m<((Nm<l)?Nm:l)+1; ++m){
                lc = complex<double>(l);
                mc = complex<double>(m);
                A1(l,m) = (-1.0) *  (lc+1.0-mc)/(2.0*lc+1.0)  * lc/(lc+1.0);
                A2(l,m) =               (lc+mc)/(2.0*lc+1.0);

                A1(l,m) = (-1.0) *  (lc+1.0-mc)/(2.0*lc+1.0)  * lc/(lc+1.0);
                A2(l,m) =               (lc+mc)/(2.0*lc+1.0);
            }
        }
        A1(0,0) = 1.0;
    // A2[0] is not used

    // ------------------------------------------------------------------------ // 
    //       Calculate the "B1, B2" parameters
    // ------------------------------------------------------------------------ // 
        for (size_t l(1); l<Nl+1; ++l){
            lc = complex<double>(l);
            B1[l] =  lc* lc     /(2.0*lc+1.0);
            B2[l] =  lc*(lc+1.0)/(2.0*lc+1.0);
        }
        B1[0] = 1.0;   // Need dummy for division later
        B2[0] = 1.0;

    // ------------------------------------------------------------------------ // 
    //       Calculate the "C1, C3" parameters
    // ------------------------------------------------------------------------ // 
        for (size_t l(1); l<Nl+1; ++l){
            lc = complex<double>(l);
            C1[l] = (-0.5) *  lc /((2.0*lc+1.0)*(lc+1.0));
            C3[l] = (-0.5) /(2.0*lc+1.0);
        }
        C1[0] = 0.5;

    // ------------------------------------------------------------------------ // 
    //       Calculate the "C2, C4" parameters
    // ------------------------------------------------------------------------ // 
        for (size_t l(0); l<Nl+1; ++l)
        {
            for (size_t m(0); m<((Nm<l)?Nm:l)+1; ++m)
            {
                if (l < 3 && m < 2)
                {
                    C2(l,m) = static_cast<complex<double> > (1.0);
                }
                else
                {
                    lc = complex<double>(l);
                    mc = complex<double>(m);
                    C2(l,m) = (0.5) *  lc * (lc-mc+2.0)*(lc-mc+1.0)/((2.0*lc+1.0)*(lc+1.0));
                    C4(l,m) = (0.5) *  (lc+mc-1.0)*(lc+mc)/(2.0*lc+1.0);    
                }
                
            }
        }
    // ------------------------------------------------------------------------ // 

    // ------------------------------------------------------------------------ // 
    //       H at the 0 momentum cell
        Hp0[0] = 1.0 / pr[0];
        for (size_t l(1); l < Nl+1; ++l) {
            double ld(l);
            Hp0[l] = Hp0[l-1] * (pr[0]/pr[1]) * (2.0*ld+1.0)/(2.0*ld-1.0);
        }

        A100 = complex<double>(1.0);
        C100 = complex<double>(0.5);
        A210 = complex<double>(1.0/3.0);
        B211 = complex<double>(2.0/3.0);
        A310 = complex<double>(2.0/5.0);
        C311 = complex<double>(-0.5/5.0);

// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
    // OpenMP stuff
        size_t num_threads = Input::List().ompthreads;
        size_t num_dists = (((Nm+1)*(2*Nl-Nm+2))/2);
        size_t fsperthread = static_cast<size_t>(num_dists/num_threads - 2);

        if (num_threads > 272 || num_threads < 1) num_threads = 1;

    // // std::cout << "\n numthreads = " << num_threads << "\n\n";    
    //     size_t lsperthread = static_cast<size_t>(Nl/num_threads - 2);
        
        if (num_threads > 1) 
        {    
            f_start[0] = 1;
            f_end[0]   = fsperthread - 2;           // Remember that f_end isn't actually processed
                                                // until the boundaries so that a gap of 2 is 
                                                // actually a gap of 3

        // if (num_threads > 1) std::cout << "\nls[" << 0 << "]=" << f_start[0] << "\n\n";
        // if (num_threads > 1) std::cout << "le[" << 0 << "]=" << f_end[0] << "\n\n";

        for (size_t i(1); i<num_threads-1; ++i)
        {        
            f_start[i] = f_end[i-1] + 2;
            f_end[i]   = f_start[i] + fsperthread;
            // if (num_threads > 1)    std::cout << "\nls[" << i << "]=" << f_start[i] << "\n\n";
            // if (num_threads > 1)    std::cout << "le[" << i << "]=" << f_end[i] << "\n\n";
        }

        f_start[num_threads-1] = f_end[num_threads-2] + 2;
        f_end[num_threads-1]   = num_dists; 
    }
    else
    {
        f_start[0] = 1;
        f_end[0]   = num_dists; 
    }

    

    size_t il(0), im(0);
    for (size_t id(0); id < num_dists; ++id)
    {
        dist_il[id] = il;
        dist_im[id] = im;

        if (il < Nl)
        {
            ++il;
        }
        else
        {
            ++im;
            il = im;
        }
    }
    

    il = 0; im = 0;
    size_t diag(0);
    for (size_t id(0); id < num_dists; ++id)
    {
        // std::cout << "\n1(id,l,m) =  " << id << "," << il << "," << im << "," << "\n";
        nwsediag_il[id] = il;
        nwsediag_im[id] = im;

        if (il < Nl && im < Nm && im <= il)
        {
            ++il; ++im;
        }
        else
        {
            ++diag;   
            im = 0;
            il = diag;                
        }
    }

    


    il = 0; im = 0; diag = 0;

    for (size_t id(0); id < num_dists; ++id)
    {
        // std::cout << "\n2(id,l,m) =  " << id << "," << il << "," << im << "," << "\n";
        neswdiag_il[id] = il;
        neswdiag_im[id] = im;

        if (il < Nl && im > 0)
        {
            ++il; --im;
            diag++;
        }
        else
        {
            if (neswdiag_il[id-diag] == neswdiag_im[id-diag])
            {
                im = neswdiag_im[id-diag];
                il = neswdiag_il[id-diag] + 1;
                diag = 0;
            }
            else if (neswdiag_il[id-diag] == neswdiag_im[id-diag] + 1)
            {
                if (neswdiag_im[id-diag] < Nm)
                {
                    im = neswdiag_im[id-diag] + 1;
                    il = neswdiag_il[id-diag];
                }
                else 
                {
                    im = Nm;
                    il = neswdiag_il[id-diag] + 1;
                }
                diag = 0;
            }
            else
            {
                im = Nm;
                il = neswdiag_il[id-diag] + 1;
                diag = 0;
            }
        }
    }
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // -----     
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Electric_Field::operator()(const DistFunc1D& Din,
   const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

//     complex<double> ii(0.0,1.0);

//     valarray<complex<double> > Ex(FEx.array());
//     valarray<complex<double> > Em(FEz.array());
//     Em *= (-1.0)*ii;
//     Em += FEy.array();
//     valarray<complex<double> > Ep(FEz.array());
//     Ep *= ii;
//     Ep += FEy.array();

//     Ex *= Din.q();
//     Em *= Din.q();
//     Ep *= Din.q();

//     size_t l0(Din.l0());
//     size_t m0(Din.m0());

//     SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx()), TMP(pr.size(),FEx.numx());

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, l = 0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeG00(Din(0,0),G);
//     Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxaxis(Ex);
//     Em *= C1[0];            Dh(1,1) += TMP.mxaxis(Em);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, l = 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(1,0),G,H,1);
//     Ex *= A2(1,0) / A1(0,0);          Dh(0,0) += H.mxaxis(Ex);
//     Ex *= A1(1,0) / A2(1,0); TMP = G; Dh(2,0) += G.mxaxis(Ex);
//     Em *= C1[1]   / C1[0]  ;          Dh(2,1) += TMP.mxaxis(Em);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(2); l < l0; ++l){
//         MakeGH(Din(l,0),G,H,l);
//         Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Dh(l-1,0) += H.mxaxis(Ex);
//         Em *= C3[l]   / C1[l-1];             Dh(l-1,1) += TMP.mxaxis(Em);
//         Ex *= A1(l,0) / A2(l,0);   TMP = G;  Dh(l+1,0) += G.mxaxis(Ex);
//         Em *= C1[l]   / C3[l];               Dh(l+1,1) += TMP.mxaxis(Em);
//     }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(l0,0),G,H,l0);
//     Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Dh(l0-1,0) += H.mxaxis(Ex);
//     Em *= C3[l0]   / C1[l0-1];             Dh(l0-1,1) += TMP.mxaxis(Em);
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l = 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(1,1),G,H,1);
//     Ep *= B2[1];              H = H.mxaxis(Ep);     Dh(0,0) += H.Re();
//     Ex *= A1(1,1) / A2(l0,0); TMP = G;              Dh(2,1) += G.mxaxis(Ex);
//     Em *= C1[1] / C3[l0];     G = TMP;              Dh(2,2) += TMP.mxaxis(Em);
//     Ep *= B1[1] / B2[1];      G = G.mxaxis(Ep);     Dh(2,0) += G.Re();

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l = 2
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(2,1),G,H,2);
//     Ex *= A2(2,1) / A1(1,1); TMP = H;               Dh(1,1) += TMP.mxaxis(Ex);
//     Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep);      Dh(1,0) += H.Re();
//     Ex *= A1(2,1) / A2(2,1); TMP = G;               Dh(3,1) += G.mxaxis(Ex);
//     Em *= C1[2]   / C1[1];   G = TMP;               Dh(3,2) += TMP.mxaxis(Em);
//     Ep *= B1[2]   / B2[2];   G = G.mxaxis(Ep);      Dh(3,0) += G.Re();

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(3); l < l0; ++l){
//         MakeGH(Din(l,1),G,H,l);
//         Ex *= A2(l,1) / A1(l-1,1); TMP = H;             Dh(l-1,1) += H.mxaxis(Ex);
//         Em *= C3[l]   / C1[l-1];   H = TMP;             Dh(l-1,2) += TMP.mxaxis(Em);
//         Ep *= B2[l]   / B1[l-1];   H = H.mxaxis(Ep);    Dh(l-1,0) += H.Re();
//         Ex *= A1(l,1) / A2(l,1);   TMP = G;             Dh(l+1,1) += G.mxaxis(Ex);
//         Em *= C1[l]   / C3[l];     G = TMP;             Dh(l+1,2) += TMP.mxaxis(Em);
//         Ep *= B1[l]   / B2[l];     G = G.mxaxis(Ep);    Dh(l+1,0) += G.Re();
//     }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //       m = 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(l0,1),G,H,l0);
//     Ex *= A2(l0,1) / A1(l0-1,1); TMP = H;              Dh(l0-1,1) += H.mxaxis(Ex);
//     Em *= C3[l0]   / C1[l0-1];   H = TMP;              Dh(l0-1,2) += TMP.mxaxis(Em);
//     Ep *= B2[l0]   / B1[l0-1];   H = H.mxaxis(Ep);     Dh(l0-1,0) += H.Re();
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     C4(l0,1) = B2[l0];
//     for (size_t m(2); m < m0; ++m){
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m > 1 , l = m
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m,m),G,H,m);
//         Ep *= C4(m,m) / C4(l0,m-1);              Dh(m-1,m-1) += H.mxaxis(Ep);
//         Ex *= A1(m,m) / A2(l0,m-1); TMP = G;     Dh(m+1,m  ) += G.mxaxis(Ex);
//         Em *= C1[m]   / C3[l0];     G = TMP;     Dh(m+1,m+1) += TMP.mxaxis(Em);
//         Ep *= C2(m,m) / C4(m,m);                 Dh(m+1,m-1) += G.mxaxis(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m > 1 , l = m+1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m+1,m),G,H,m+1);
//         Ex *= A2(m+1,m) / A1(m,m);   TMP = H;     Dh(m  ,m  ) += TMP.mxaxis(Ex);
//         Ep *= C4(m+1,m) / C2(m,m);                Dh(m  ,m-1) += H.mxaxis(Ep);
//         if ( m+1 < l0) { //always true except when m = m0-1 = l0-1
//             Ex *= A1(m+1,m) / A2(m+1,m); TMP = G;     Dh(m+2,m  ) += G.mxaxis(Ex);
//             Em *= C1[m+1]   / C1[m];     G = TMP;     Dh(m+2,m+1) += TMP.mxaxis(Em);
//             Ep *= C2(m+1,m) / C4(m+1,m);              Dh(m+2,m-1) += G.mxaxis(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //              m > 1, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//             for (size_t l(m+2); l < l0; ++l){
//                 MakeGH(Din(l,m),G,H,l);
//                 Ex *= A2(l,m) / A1(l-1,m); TMP = H;   Dh(l-1,m  ) += H.mxaxis(Ex);
//                 Em *= C3[l]   / C1[l-1];   H = TMP;   Dh(l-1,m+1) += TMP.mxaxis(Em);
//                 Ep *= C4(l,m) / C2(l-1,m);            Dh(l-1,m-1) += H.mxaxis(Ep);
//                 Ex *= A1(l,m) / A2(l,m);   TMP = G;   Dh(l+1,m  ) += G.mxaxis(Ex);
//                 Em *= C1[l]   / C3[l];     G = TMP;   Dh(l+1,m+1) += TMP.mxaxis(Em);
//                 Ep *= C2(l,m) / C4(l,m);              Dh(l+1,m-1) += G.mxaxis(Ep);
//             }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //               m > 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//             MakeGH(Din(l0,m),G,H,l0);
//             Ex *= A2(l0,m) / A1(l0-1,m); TMP = H;    Dh(l0-1,m  ) += H.mxaxis(Ex);
//             Em *= C3[l0]   / C1[l0-1];   H = TMP;    Dh(l0-1,m+1) += TMP.mxaxis(Em);
//             Ep *= C4(l0,m) / C2(l0-1,m);             Dh(l0-1,m-1) += H.mxaxis(Ep);
//         }
//     }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(m0,m0),G,H,m0);
//     Ep *= C4(m0,m0) / C4(l0,m0-1);              Dh(m0-1,m0-1) += H.mxaxis(Ep);
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //       m = m0, l0 > l > m0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     if ( m0 < l0) {
//         Ex *= A1(m0,m0) / A2(l0,m0-1); TMP = G;  Dh(m0+1,m0  ) += TMP.mxaxis(Ex);
//         Ep *= C2(m0,m0) / C4(m0,m0);             Dh(m0+1,m0-1) += G.mxaxis(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m = m0 , l = m0+1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m0+1,m0),G,H,m0+1);
//         Ex *= A2(m0+1,m0) / A1(m0,m0); TMP = H;        Dh(m0,m0  )  += TMP.mxaxis(Ex);
//         Ep *= C4(m0+1,m0) / C2(m0,m0);                 Dh(m0,m0-1)  += H.mxaxis(Ep);
//         if ( m0+1 < l0) {
//             Ex *= A1(m0+1,m0) / A2(m0+1,m0); TMP = G;  Dh(m0+2,m0  )+= TMP.mxaxis(Ex);
//             Ep *= C2(m0+1,m0) / C4(m0+1,m0);           Dh(m0+2,m0-1)+= G.mxaxis(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //              m = m0, m0+2 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//             for (size_t l(m0+2); l < l0; ++l){
//                 MakeGH(Din(l,m0),G,H,l);
//                 Ex *= A2(l,m0) / A1(l-1,m0); TMP = H;  Dh(l-1,m0)   += TMP.mxaxis(Ex);
//                 Ep *= C4(l,m0) / C2(l-1,m0);           Dh(l-1,m0-1) += H.mxaxis(Ep);
//                 Ex *= A1(l,m0) / A2(l,  m0); TMP = G;  Dh(l+1,m0  ) += TMP.mxaxis(Ex);
//                 Ep *= C2(l,m0) / C4(l,  m0);           Dh(l+1,m0-1) += G.mxaxis(Ep);
//             }

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //               m > 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//             MakeGH(Din(l0,m0),G,H,l0);
//             Ex *= A2(l0,m0) / A1(l0-1,m0); TMP = H;    Dh(l0-1,m0)   += TMP.mxaxis(Ex);
//             Ep *= C4(l0,m0) / C2(l0-1,m0);             Dh(l0-1,m0-1) += H.mxaxis(Ep);
//         }
//     }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    /////  Vertical Iteration
//  -------------------------------------------------------- //
    //   Because each iteration in the loop modifies + and - 1
    //   The parallelization is performed in chunks and boundaries
    //   are taken care of later
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(Input::List().ompthreads)
    {   
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();
        size_t f_start_thread(f_start[this_thread]); ///< Chunk starts here
        size_t f_end_thread(f_end[this_thread]);     ///< Chunk ends here

        /// Local variables for each thread
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

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx()), TMP(pr.size(),FEx.numx());

        size_t l(0),m(0);


        if (this_thread == 0)
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 0, l = 0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeG00(Din(0,0),G);
            Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxaxis(Ex);    Ex *= 1.0/A1(0,0); 
            Em *= C1[0];            Dh(1,1) += TMP.mxaxis(Em);  Em *= 1.0/C1[0];            
            
            // std::cout << "\n Checkpoint #0";   Dh.checknan();
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 1 loop, 1 <= l < l0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Ep *= B1[0];
            for (size_t il = 1; il < l0; ++il)
            {
                MakeGH(Din(il,1),G,H,il);
                Ep *= B2[il]/B1[il-1];      H = H.mxaxis(Ep);       Dh(il-1,0) += H.Re();
                Ep *= B1[il]/B2[il];        G = G.mxaxis(Ep);       Dh(il+1,0) += G.Re();
                
            }

            MakeGH(Din(l0,1),G,H,l0);     Ep *= B2[l0]/B1[l0-1];          H = H.mxaxis(Ep);       Dh(l0-1,0) += H.Re();   Ep *= 1.0/B2[l0];

            
        }
        // std::cout << "\n Checkpoint #1 \n";   Dh.checknan();
        

        //  -------------------------------------------------------- //
        //  Do the chunks
        //  -------------------------------------------------------- //       
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Vertical loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {

            l = dist_il[id];
            m = dist_im[id];

            // std::cout << "\n (l,m) = " << l << ", " << m << " \n";

            MakeGH(Din(l,m),G,H,l);

            if (l == m)         // Diagonal, no l - 1
            {
                if (l < l0) {   Ex *= A1(l,m);   Dh(l+1,m) += G.mxaxis(Ex);     Ex *= 1.0 / A1(l,m);}
            }
            else if (l == l0)   // Last l, no l + 1
            {
                                Ex *= A2(l0,m);  Dh(l0-1,0) += H.mxaxis(Ex);      Ex *= 1.0 / A2(l0,m);
            }
            else
            {
                                Ex *= A2(l,m);              Dh(l-1,m) += H.mxaxis(Ex);
                                Ex *= A1(l,m) / A2(l,m);    Dh(l+1,m) += G.mxaxis(Ex);  Ex *= 1.0 / A1(l,m);
            }
        }

        // std::cout << "\n Checkpoint #2 \n";   Dh.checknan();
        
        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = nwsediag_il[id];
            m = nwsediag_im[id];

            MakeGH(Din(l,m),G,H,l);

            if (m == 0)         // Top or Left, no l - 1, m - 1
            {
                // std::cout << "1: ( " << l << "," << m << ")";
                if (l < l0) {   Em *= C1[m];    Dh(l+1,m+1) += G.mxaxis(Em);        Em *= 1.0/C1[m];}
            }
            else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
            {
                // std::cout << "2: ( " << l << "," << m << ")";
                                Ep *= C4(l,m);   Dh(l-1,m-1) += H.mxaxis(Ep);       Ep *= 1.0/C4(l,m);
            }
            else
            {
                // std::cout << "3: ( " << l << "," << m << ")";
                                Em *= C1[m];    Dh(l+1,m+1) += G.mxaxis(Em);        Em *= 1.0/C1[m];
                if (m > 1)  {   Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxaxis(Ep);        Ep *= 1.0/C4(l,m);}
            }

            // std::cout << "\n Checkpoint ( " << l << "," << m << ") \n\n";   Dh.checknan();
        }

        // std::cout << "\n Checkpoint #3 \n";   Dh.checknan();        
        
        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Anti-Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = neswdiag_il[id];
            m = neswdiag_im[id];

            // std::cout << "Before ( " << l << "," << m << ") \n";

            MakeGH(Din(l,m),G,H,l);

            

            if (m == 0)         // Left wall, no l + 1, m - 1
            {
                // std::cout << "1: ( " << l << "," << m << ")" << std::endl;
                if (l > 1)                  {    Em *= C3[l];    Dh(l-1,m+1) += H.mxaxis(Em);        Em *= 1.0/C3[l];}
            }
            else if (m == m0)   // Right boundary, no l - 1, m + 1
            {
                // std::cout << "2: ( " << l << "," << m << ")" << std::endl;
                if (l < l0)                 {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxaxis(Ep);        Ep *= 1.0/C2(l,m);}
            }
            else
            {
                // std::cout << "3: ( " << l << "," << m << ")" << std::endl;

                if (m > 1 && l < l0)        {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxaxis(Ep);        Ep *= 1.0/C2(l,m);}
                if (l - 1 != m && l != m)   {   Em *= C3[l];    Dh(l-1,m+1) += H.mxaxis(Em);        Em *= 1.0/C3[l];}
            }
            // std::cout << "\n Checkpoint ( " << l << "," << m << ") ... ";   Dh.checknan();  std::cout << "passed \n" << std::endl;
        }
        // std::cout << "Checkpoint #4";   Dh.checknan();
    }

    

    //  -------------------------------------------------------- //
    //  Do the boundaries between the chunks
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(f_start.size()-1)
    // for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {  
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();

        /// Local variables for each thread
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

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx()), TMP(pr.size(),FEx.numx());

        size_t l(0),m(0);

        if (this_thread < f_start.size() - 1) 
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      Vertical loop, boundaries between threads
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {

                l = dist_il[id];
                m = dist_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (l == m)         // Diagonal, no l - 1
                {
                    if (l < l0) {   Ex *= A1(l,0);   Dh(l+1,0) += G.mxaxis(Ex);     Ex *= 1.0 / A1(l,0);}
                }
                else if (l == l0)   // Last l, no l + 1
                {
                                    Ex *= A2(l,0);  Dh(l-1,0) += H.mxaxis(Ex);      Ex *= 1.0 / A2(l,0);
                }
                else
                {
                                    Ex *= A2(l,0);              Dh(l-1,0) += H.mxaxis(Ex);
                                    Ex *= A1(l,0) / A2(l,0);    Dh(l+1,0) += G.mxaxis(Ex);  Ex *= 1.0 / A1(l,0);
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = nwsediag_il[id];
                m = nwsediag_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (m == 0)         // Top or Left, no l - 1, m - 1
                {
                                    Em *= C1[m];    Dh(l+1,m+1) += G.mxaxis(Em);        Em *= 1.0/C1[m];
                }
                else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
                {
                                    Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxaxis(Ep);        Ep *= 1.0/C4(l,m);
                }
                else
                {
                                    Em *= C1[m];    Dh(l+1,m+1) += G.mxaxis(Em);        Em *= 1.0/C1[m];
                    if (m > 1)  {   Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxaxis(Ep);        Ep *= 1.0/C4(l,m);}
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Anti-Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = neswdiag_il[id];
                m = neswdiag_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (m == 0)         // Left wall, no l + 1, m - 1
                {
                    if (l > 1)                  {   Em *= C3[l];    Dh(l-1,1) += H.mxaxis(Em);        Em *= 1.0/C3[l];}
                }
                else if (m == m0)   // Right boundary, no l - 1, m + 1
                {
                    if (l < l0)                 {   Ep *= C2(l,m0);  Dh(l+1,m0-1) += G.mxaxis(Ep);        Ep *= 1.0/C2(l,m0);}
                }
                else
                {
                    if (m > 1 && l < l0)        {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxaxis(Ep);        Ep *= 1.0/C2(l,m);}
                    if (l - 1 != m && l != m)   {   Em *= C3[l];    Dh(l-1,m+1) += H.mxaxis(Em);        Em *= 1.0/C3[l];}
                }
            }
        }
    }

}
//--------------------------------------------------------------
//
void Electric_Field::es1d(const DistFunc1D& Din,
   const Field1D& FEx,
   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

    //  -------------------------------------------------------- //
    //   Because each iteration in the loop modifies + and - 1
    //   The parallelization is performed in chunks and boundaries
    //   are taken care of later
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(Input::List().ompthreads)
    {   
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();
        size_t f_start_thread(f_start[this_thread]); ///< Chunk starts here
        size_t f_end_thread(f_end[this_thread]);     ///< Chunk ends here

        /// Local variables for each thread
        size_t l0(Din.l0());
        valarray<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx());

        //  -------------------------------------------------------- //
        //   First thread takes the boundary conditions (l = 0, 1)
        //   Last thread takes the boundary condition (l = l0)
        //   Rest of the threads proceed to chunks
        //  -------------------------------------------------------- //
        if (this_thread==0)
        {       
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 0, l = 0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeG00(Din(0,0),G);
            Ex *= A1(0,0);  Dh(1,0) += G.mxaxis(Ex);
        }

        if (this_thread==Input::List().ompthreads - 1)
        {                       
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      m = 0,  l = l0
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeGH(Din(l0,0),G,H,l0);
            Ex *= A2(l0,0);  Dh(l0-1,0) += H.mxaxis(Ex);
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Ex /= A2(l0,0);             // Reset Ex
                                        // 
            f_end_thread -= 1;
        }

        //  -------------------------------------------------------- //
        //  Do the chunks
        //  Initialize Ex so that it its ready for loop iteration l
        //  -------------------------------------------------------- //        
        Ex *= A1(f_start_thread-1,0);

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      m = 0, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t l = f_start_thread; l < f_end_thread; ++l)
        {
            MakeGH(Din(l,0),G,H,l);

            Ex *= A2(l,0) / A1(l-1,0);  Dh(l-1,0) += H.mxaxis(Ex);
            Ex *= A1(l,0) / A2(l,0);   Dh(l+1,0) += G.mxaxis(Ex);
        }
    }

    //  -------------------------------------------------------- //
    //  Do the boundaries between the chunks
    //  -------------------------------------------------------- //
    #pragma omp parallel for num_threads(f_start.size()-1)
    for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {    
        SHarmonic1D G(pr.size(),FEx.numx()),H(pr.size(),FEx.numx());
        valarray<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

    //  Initialize Ex so that it its ready for loop iteration l
        Ex *= A1(f_end[threadboundaries]-1,0);

        for (size_t l = f_end[threadboundaries]; l < f_start[threadboundaries+1]; ++l)
        {
            MakeGH(Din(l,0),G,H,l);

            Ex *= A2(l,0) / A1(l-1,0);     Dh(l-1,0) += H.mxaxis(Ex);
            Ex *= A1(l,0) / A2(l,0);       Dh(l+1,0) += G.mxaxis(Ex);
        }
    }


    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    void Electric_Field::Implicit_Ex(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        valarray<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx()), TMP(pr.size(),FEx.numx());

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Ex *= A1(0,0);             Dh(1,0) += G.mxaxis(Ex);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Ex *= A2(1,0) / A1(0,0);   Dh(0,0) += H.mxaxis(Ex);
        Ex *= A1(1,0) / A2(1,0);   Dh(2,0) += G.mxaxis(Ex);

//      m = 0, l = 2
        MakeGH(Din(2,0),G,H,2);
        Ex *= A2(2,0) / A1(1,0);   Dh(1,0) += H.mxaxis(Ex);

//      m = 0, l = 3
        MakeGH(Din(3,0),G,H,3);
        Ex *= A2(3,0) / A2(2,0);   Dh(2,0) += H.mxaxis(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ex *= A1(1,1) / A2(3,0);   Dh(2,1) += G.mxaxis(Ex);

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ex *= A2(2,1) / A1(1,1);   Dh(1,1) += H.mxaxis(Ex);

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Ex *= A2(3,1) / A2(2,1);   Dh(2,1) += H.mxaxis(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::Implicit_Ey(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        valarray<complex<double> > Em(FEy.array());
        valarray<complex<double> > Ep(FEy.array());

        SHarmonic1D G(pr.size(),FEy.numx()), H(pr.size(),FEy.numx()), TMP(pr.size(),FEy.numx());

        Em *= Din.q();
        Ep *= Din.q();

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Em *= C1[0];            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
        MakeGH(Din(2,0),G,H,2);
        Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxaxis(Em);

//      m = 0,  l = 3
        MakeGH(Din(3,0),G,H,3);
        Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ep *= B2[1];             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;          Dh(2,2) += TMP.mxaxis(Em);
        Ep *= B1[1] / B2[1];     G = G.mxaxis(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Em *= C3[3]   / C1[1];   TMP = H;          Dh(2,2) += TMP.mxaxis(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxaxis(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(Din(2,2),G,H,2);
        Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxaxis(Ep);

//      m = 2, l = 3
        MakeGH(Din(3,2),G,H,3);
        Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxaxis(Ep);

        size_t m0(Din.m0());
        if ( m0 > 2) {
//          m = 3, l = 3
            MakeGH(Din(3,3),G,H,3);
            Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxaxis(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
//
////--------------------------------------------------------------
    void Electric_Field::Implicit_Ez(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh) {
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
        SHarmonic1D G(pr.size(),FEz.numx()), H(pr.size(),FEz.numx()), TMP(pr.size(),FEz.numx());
    // Ex *= Din.q();;
        Em *= Din.q();
        Ep *= Din.q();



//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Em *= C1[0];            Dh(1,1) += G.mxaxis(Em);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxaxis(Em);

//      m = 0, 1 < l < l0
        MakeGH(Din(2,0),G,H,2);
        Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxaxis(Em);

//      m = 0,  l = 3
        MakeGH(Din(3,0),G,H,3);
        Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxaxis(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ep *= B2[1];             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;              Dh(2,2) += TMP.mxaxis(Em);
        Ep *= B1[1] / B2[1];     G = G.mxaxis(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ep *= B2[2]   / B1[1];   H = H.mxaxis(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Em *= C3[3]   / C1[1];   TMP = H;              Dh(2,2) += TMP.mxaxis(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxaxis(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(Din(2,2),G,H,2);
        Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxaxis(Ep);

//      m = 2, l = 3
        MakeGH(Din(3,2),G,H,3);
        Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxaxis(Ep);

        size_t m0(Din.m0());
        if ( m0 > 2) {
//          m = 3, l = 3
            MakeGH(Din(3,3),G,H,3);
            Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxaxis(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::f1only(const DistFunc1D& Din,
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

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx()), TMP(pr.size(),FEx.numx());

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeG00(Din(0,0),G);
        Ex *= A100; TMP = G; Dh(1,0) += G.mxaxis(Ex);
        Em *= C100;            Dh(1,1) += TMP.mxaxis(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(1,0),G,H,1);
        Ex *= A210 / A100;          Dh(0,0) += H.mxaxis(Ex);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;              H = H.mxaxis(Ep);     Dh(0,0) += H.Re();

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field::Implicit_Ex_f1only(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        valarray<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

        SHarmonic1D G(pr.size(),FEx.numx()), H(pr.size(),FEx.numx());

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Ex *= A100;             Dh(1,0) += G.mxaxis(Ex);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Ex *= A210 / A100;   Dh(0,0) += H.mxaxis(Ex);


    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::Implicit_Ey_f1only(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh) {
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

        SHarmonic1D G(pr.size(),FEy.numx()), H(pr.size(),FEy.numx());

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
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
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
    }
//
////--------------------------------------------------------------
    void Electric_Field::Implicit_Ez_f1only(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh) {
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


        SHarmonic1D G(pr.size(),FEz.numx()), H(pr.size(),FEz.numx());

        Em *= Din.q();
        Ep *= Din.q();

    //      m = 0, l = 0
        MakeG00(Din(0,0),G);
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
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;             H = H.mxaxis(Ep); Dh(0,0) += H.Re();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Make derivatives -(l+1/l)*G and H for a given f , used in openMP routine
    void Electric_Field::MakeGH(const SHarmonic1D& f, SHarmonic1D& G, SHarmonic1D& H, size_t el)
    {
//--------------------------------------------------------------
        valarray<complex<double> > invpax(invpr);
        valarray<complex<double> > invdp_local(invdp);
        complex<double> ld(el);

    // invpax *= (-2.0)*(ld+1.0) * (pr[1]-pr[0]);
    invpax *= (ld+1.0);   // Non-uniform grid

    G = f;                   H = f;
    G = G.Dp();              G = G.mpaxis(invdp_local);  // Non-uniform grid
    H  = H.mpaxis(invpax);
    H += G;
    G *= -(2.0*ld+1.0)/ld;
    G += H;

    for (size_t i(0); i < G.numx(); ++i) G(0,i) = 0.0;
        for (size_t i(0); i < H.numx(); ++i) H(0,i) = f(1,i) * Hp0[el];
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Make derivatives -(l+1/l)*G and H for a given f , used in openMP routine
void Electric_Field::MakeGH(const SHarmonic2D& f, SHarmonic2D& G, SHarmonic2D& H, size_t el)
{
//--------------------------------------------------------------
    valarray<complex<double> > invpax(invpr);
    valarray<complex<double> > invdp_local(invdp);
    complex<double> ld(el);

    // invpax *= (-2.0)*(ld+1.0) * (pr[1]-pr[0]);
    invpax *= (ld+1.0);   // Non-uniform grid

    G = f;                   H = f;
    G = G.Dp();              G = G.mpaxis(invdp_local);  // Non-uniform grid
    H  = H.mpaxis(invpax);
    H += G;
    G *= -(2.0*ld+1.0)/ld;
    G += H;

    for (size_t ix(0); ix < G.numx(); ++ix)
    {
        for (size_t iy(0); iy < G.numy(); ++iy)
        {
            G(0,ix,iy) = 0.0;
        }
    } 
    
    for (size_t ix(0); ix < G.numx(); ++ix)
    {
        for (size_t iy(0); iy < G.numy(); ++iy)
        {
            H(0,ix,iy) = f(1,ix,iy) * Hp0[el];
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Calculation of G00 = df/dp(p0)
void Electric_Field::MakeG00(const SHarmonic1D& f, SHarmonic1D& G) {
//--------------------------------------------------------------
    G = f; G = G.Dp();  G = G.mpaxis(invdp);

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
//  Calculation of G00 = df/dp(p0)
void Electric_Field::MakeG00(const SHarmonic2D& f, SHarmonic2D& G) {
//--------------------------------------------------------------
    G = f; G = G.Dp();  G = G.mpaxis(invdp);

    complex<double> p0p1_sq( pr[0]*pr[0]/(pr[1]*pr[1]) ),
    inv_mp0p1_sq( 1.0/(1.0-p0p1_sq) ),
    g_r = -4.0*(pr[1]-pr[0]) * pr[0]/(pr[1]*pr[1]),
    f00;

    for (size_t ix(0); ix < G.numx(); ++ix)
    {
        for (size_t iy(0); iy < G.numy(); ++iy)
        {
            f00    = ( f(0,ix,iy) - f(1,ix,iy) * p0p1_sq) * inv_mp0p1_sq;
            G(0,ix,iy) = ( f(1,ix,iy) - f00) * g_r;
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Electric_Field::operator()(const DistFunc2D& Din, 
   const Field2D& FEx, const Field2D& FEy, const Field2D& FEz, 
   DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

//     complex<double> ii(0.0,1.0);

//     Array2D<complex<double> > Ex(FEx.array());
//     Array2D<complex<double> > Em(FEz.array());
//     Em *= (-1.0)*ii;
//     Em += FEy.array();
//     Array2D<complex<double> > Ep(FEz.array());
//     Ep *= ii;
//     Ep += FEy.array();

//     Ex *= Din.q();;
//     Em *= Din.q();;
//     Ep *= Din.q();;
    
//     size_t l0(Din.l0());
//     size_t m0(Din.m0());


//     SHarmonic2D G(pr.size(),FEx.numx(),FEx.numy()), H(pr.size(),FEx.numx(),FEx.numy()),
//     TMP(pr.size(),FEx.numx(),FEx.numy());

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, l = 0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeG00(Din(0,0),G); 
//     Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxy_matrix(Ex);
//     Em *= C1[0];            Dh(1,1) += TMP.mxy_matrix(Em);



// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, l = 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(1,0),G,H,1);  H=0.0;
//     Ex *= A2(1,0) / A1(0,0);          Dh(0,0) += H.mxy_matrix(Ex);
//     Ex *= A1(1,0) / A2(1,0); TMP = G; Dh(2,0) += G.mxy_matrix(Ex);
//     Em *= C1[1]   / C1[0]  ;          Dh(2,1) += TMP.mxy_matrix(Em);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(2); l < l0; ++l){
//         MakeGH(Din(l,0),G,H,l);
//         Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Dh(l-1,0) += H.mxy_matrix(Ex);
//         Em *= C3[l]   / C1[l-1];             Dh(l-1,1) += TMP.mxy_matrix(Em);
//         Ex *= A1(l,0) / A2(l,0);   TMP = G;  Dh(l+1,0) += G.mxy_matrix(Ex);
//         Em *= C1[l]   / C3[l];               Dh(l+1,1) += TMP.mxy_matrix(Em);
//     }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(l0,0),G,H,l0);
//     Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Dh(l0-1,0) += H.mxy_matrix(Ex);
//     Em *= C3[l0]   / C1[l0-1];             Dh(l0-1,1) += TMP.mxy_matrix(Em);
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l = 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(1,1),G,H,1);
//     Ep *= B2[1];              H = H.mxy_matrix(Ep);     Dh(0,0) += H.Re();
//     Ex *= A1(1,1) / A2(l0,0); TMP = G;              Dh(2,1) += G.mxy_matrix(Ex); 
//     Em *= C1[1] / C3[l0];     G = TMP;              Dh(2,2) += TMP.mxy_matrix(Em); 
//     Ep *= B1[1] / B2[1];      G = G.mxy_matrix(Ep);     Dh(2,0) += G.Re();

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l = 2
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(2,1),G,H,2);
//     Ex *= A2(2,1) / A1(1,1); TMP = H;               Dh(1,1) += TMP.mxy_matrix(Ex);
//     Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep);      Dh(1,0) += H.Re();
//     Ex *= A1(2,1) / A2(2,1); TMP = G;               Dh(3,1) += G.mxy_matrix(Ex);
//     Em *= C1[2]   / C1[1];   G = TMP;               Dh(3,2) += TMP.mxy_matrix(Em);
//     Ep *= B1[2]   / B2[2];   G = G.mxy_matrix(Ep);      Dh(3,0) += G.Re();

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(3); l < l0; ++l){
//         MakeGH(Din(l,1),G,H,l);
//         Ex *= A2(l,1) / A1(l-1,1); TMP = H;             Dh(l-1,1) += H.mxy_matrix(Ex);
//         Em *= C3[l]   / C1[l-1];   H = TMP;             Dh(l-1,2) += TMP.mxy_matrix(Em);
//         Ep *= B2[l]   / B1[l-1];   H = H.mxy_matrix(Ep);    Dh(l-1,0) += H.Re();
//         Ex *= A1(l,1) / A2(l,1);   TMP = G;             Dh(l+1,1) += G.mxy_matrix(Ex);
//         Em *= C1[l]   / C3[l];     G = TMP;             Dh(l+1,2) += TMP.mxy_matrix(Em);
//         Ep *= B1[l]   / B2[l];     G = G.mxy_matrix(Ep);    Dh(l+1,0) += G.Re();
//     }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //       m = 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     MakeGH(Din(l0,1),G,H,l0);
//     Ex *= A2(l0,1) / A1(l0-1,1); TMP = H;              Dh(l0-1,1) += H.mxy_matrix(Ex);
//     Em *= C3[l0]   / C1[l0-1];   H = TMP;              Dh(l0-1,2) += TMP.mxy_matrix(Em);
//     Ep *= B2[l0]   / B1[l0-1];   H = H.mxy_matrix(Ep);     Dh(l0-1,0) += H.Re();
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     C4(l0,1) = B2[l0];
//     for (size_t m(2); m < m0; ++m){
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m > 1 , l = m
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m,m),G,H,m);
//         Ep *= C4(m,m) / C4(l0,m-1);              Dh(m-1,m-1) += H.mxy_matrix(Ep); 
//         Ex *= A1(m,m) / A2(l0,m-1); TMP = G;     Dh(m+1,m  ) += G.mxy_matrix(Ex); 
//         Em *= C1[m]   / C3[l0];     G = TMP;     Dh(m+1,m+1) += TMP.mxy_matrix(Em); 
//         Ep *= C2(m,m) / C4(m,m);                 Dh(m+1,m-1) += G.mxy_matrix(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m > 1 , l = m+1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m+1,m),G,H,m+1);
//         Ex *= A2(m+1,m) / A1(m,m);   TMP = H;     Dh(m  ,m  ) += TMP.mxy_matrix(Ex);
//         Ep *= C4(m+1,m) / C2(m,m);                Dh(m  ,m-1) += H.mxy_matrix(Ep);
//             if ( m+1 < l0) { //always true except when m = m0-1 = l0-1 
//                 Ex *= A1(m+1,m) / A2(m+1,m); TMP = G;     Dh(m+2,m  ) += G.mxy_matrix(Ex);
//                 Em *= C1[m+1]   / C1[m];     G = TMP;     Dh(m+2,m+1) += TMP.mxy_matrix(Em);
//                 Ep *= C2(m+1,m) / C4(m+1,m);              Dh(m+2,m-1) += G.mxy_matrix(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //              m > 1, 1 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                 for (size_t l(m+2); l < l0; ++l){
//                     MakeGH(Din(l,m),G,H,l);
//                     Ex *= A2(l,m) / A1(l-1,m); TMP = H;   Dh(l-1,m  ) += H.mxy_matrix(Ex);
//                     Em *= C3[l]   / C1[l-1];   H = TMP;   Dh(l-1,m+1) += TMP.mxy_matrix(Em);
//                     Ep *= C4(l,m) / C2(l-1,m);            Dh(l-1,m-1) += H.mxy_matrix(Ep);
//                     Ex *= A1(l,m) / A2(l,m);   TMP = G;   Dh(l+1,m  ) += G.mxy_matrix(Ex);
//                     Em *= C1[l]   / C3[l];     G = TMP;   Dh(l+1,m+1) += TMP.mxy_matrix(Em);
//                     Ep *= C2(l,m) / C4(l,m);              Dh(l+1,m-1) += G.mxy_matrix(Ep); 
//                 }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //               m > 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                 MakeGH(Din(l0,m),G,H,l0);
//                 Ex *= A2(l0,m) / A1(l0-1,m); TMP = H;    Dh(l0-1,m  ) += H.mxy_matrix(Ex);
//                 Em *= C3[l0]   / C1[l0-1];   H = TMP;    Dh(l0-1,m+1) += TMP.mxy_matrix(Em);
//                 Ep *= C4(l0,m) / C2(l0-1,m);             Dh(l0-1,m-1) += H.mxy_matrix(Ep);
//             }
//         }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         MakeGH(Din(m0,m0),G,H,m0);
//         Ep *= C4(m0,m0) / C4(l0,m0-1);              Dh(m0-1,m0-1) += H.mxy_matrix(Ep); 
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //       m = m0, l0 > l > m0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         if ( m0 < l0) { 
//             Ex *= A1(m0,m0) / A2(l0,m0-1); TMP = G;  Dh(m0+1,m0  ) += TMP.mxy_matrix(Ex); 
//             Ep *= C2(m0,m0) / C4(m0,m0);             Dh(m0+1,m0-1) += G.mxy_matrix(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //          m = m0 , l = m0+1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//             MakeGH(Din(m0+1,m0),G,H,m0+1);
//             Ex *= A2(m0+1,m0) / A1(m0,m0); TMP = H;        Dh(m0,m0  )  += TMP.mxy_matrix(Ex);
//             Ep *= C4(m0+1,m0) / C2(m0,m0);                 Dh(m0,m0-1)  += H.mxy_matrix(Ep);
//             if ( m0+1 < l0) { 
//                 Ex *= A1(m0+1,m0) / A2(m0+1,m0); TMP = G;  Dh(m0+2,m0  )+= TMP.mxy_matrix(Ex);
//                 Ep *= C2(m0+1,m0) / C4(m0+1,m0);           Dh(m0+2,m0-1)+= G.mxy_matrix(Ep);

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //              m = m0, m0+2 < l < l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                 for (size_t l(m0+2); l < l0; ++l){
//                     MakeGH(Din(l,m0),G,H,l);
//                     Ex *= A2(l,m0) / A1(l-1,m0); TMP = H;  Dh(l-1,m0)   += TMP.mxy_matrix(Ex);
//                     Ep *= C4(l,m0) / C2(l-1,m0);           Dh(l-1,m0-1) += H.mxy_matrix(Ep);
//                     Ex *= A1(l,m0) / A2(l,  m0); TMP = G;  Dh(l+1,m0  ) += TMP.mxy_matrix(Ex);
//                     Ep *= C2(l,m0) / C4(l,  m0);           Dh(l+1,m0-1) += G.mxy_matrix(Ep); 
//                 }

// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //               m > 1,  l = l0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                 MakeGH(Din(l0,m0),G,H,l0);
//                 Ex *= A2(l0,m0) / A1(l0-1,m0); TMP = H;    Dh(l0-1,m0)   += TMP.mxy_matrix(Ex);
//                 Ep *= C4(l0,m0) / C2(l0-1,m0);             Dh(l0-1,m0-1) += H.mxy_matrix(Ep);
//             } 
//         }
// // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    /////  Vertical Iteration
//  -------------------------------------------------------- //
    //   Because each iteration in the loop modifies + and - 1
    //   The parallelization is performed in chunks and boundaries
    //   are taken care of later
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(Input::List().ompthreads)
    {   
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();
        size_t f_start_thread(f_start[this_thread]); ///< Chunk starts here
        size_t f_end_thread(f_end[this_thread]);     ///< Chunk ends here

        /// Local variables for each thread
        complex<double> ii(0.0,1.0);

        Array2D<complex<double> > Ex(FEx.array());
        Array2D<complex<double> > Em(FEz.array());
        Em *= (-1.0)*ii;
        Em += FEy.array();
        Array2D<complex<double> > Ep(FEz.array());
        Ep *= ii;
        Ep += FEy.array();

        Ex *= Din.q();
        Em *= Din.q();
        Ep *= Din.q();

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        SHarmonic2D G(Din(0,0)), H(Din(0,0)), TMP(Din(0,0));

        size_t l(0),m(0);


        if (this_thread == 0)
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 0, l = 0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            MakeG00(Din(0,0),G);
            Ex *= A1(0,0); TMP = G; Dh(1,0) += G.mxy_matrix(Ex);    Ex *= 1.0/A1(0,0); 
            Em *= C1[0];            Dh(1,1) += TMP.mxy_matrix(Em);  Em *= 1.0/C1[0];            
            
            // std::cout << "\n Checkpoint #0";   Dh.checknan();
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 1 loop, 1 <= l < l0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Ep *= B1[0];
            for (size_t il = 1; il < l0; ++il)
            {
                MakeGH(Din(il,1),G,H,il);
                Ep *= B2[il]/B1[il-1];      H = H.mxy_matrix(Ep);       Dh(il-1,0) += H.Re();
                Ep *= B1[il]/B2[il];        G = G.mxy_matrix(Ep);       Dh(il+1,0) += G.Re();
                
            }

            MakeGH(Din(l0,1),G,H,l0);     Ep *= B2[l0]/B1[l0-1];          H = H.mxy_matrix(Ep);       Dh(l0-1,0) += H.Re();   Ep *= 1.0/B2[l0];

            
        }
        // std::cout << "\n Checkpoint #1 \n";   Dh.checknan();
        

        //  -------------------------------------------------------- //
        //  Do the chunks
        //  -------------------------------------------------------- //       
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Vertical loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {

            l = dist_il[id];
            m = dist_im[id];

            // std::cout << "\n (l,m) = " << l << ", " << m << " \n";

            MakeGH(Din(l,m),G,H,l);

            if (l == m)         // Diagonal, no l - 1
            {
                if (l < l0) {   Ex *= A1(l,m);   Dh(l+1,m) += G.mxy_matrix(Ex);     Ex *= 1.0 / A1(l,m);}
            }
            else if (l == l0)   // Last l, no l + 1
            {
                                Ex *= A2(l0,m);  Dh(l0-1,0) += H.mxy_matrix(Ex);      Ex *= 1.0 / A2(l0,m);
            }
            else
            {
                                Ex *= A2(l,m);              Dh(l-1,m) += H.mxy_matrix(Ex);
                                Ex *= A1(l,m) / A2(l,m);    Dh(l+1,m) += G.mxy_matrix(Ex);  Ex *= 1.0 / A1(l,m);
            }
        }

        // std::cout << "\n Checkpoint #2 \n";   Dh.checknan();
        
        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = nwsediag_il[id];
            m = nwsediag_im[id];

            MakeGH(Din(l,m),G,H,l);

            if (m == 0)         // Top or Left, no l - 1, m - 1
            {
                // std::cout << "1: ( " << l << "," << m << ")";
                if (l < l0) {   Em *= C1[m];    Dh(l+1,m+1) += G.mxy_matrix(Em);        Em *= 1.0/C1[m];}
            }
            else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
            {
                // std::cout << "2: ( " << l << "," << m << ")";
                                Ep *= C4(l,m);   Dh(l-1,m-1) += H.mxy_matrix(Ep);       Ep *= 1.0/C4(l,m);
            }
            else
            {
                // std::cout << "3: ( " << l << "," << m << ")";
                                Em *= C1[m];    Dh(l+1,m+1) += G.mxy_matrix(Em);        Em *= 1.0/C1[m];
                if (m > 1)  {   Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxy_matrix(Ep);        Ep *= 1.0/C4(l,m);}
            }

            // std::cout << "\n Checkpoint ( " << l << "," << m << ") \n\n";   Dh.checknan();
        }

        // std::cout << "\n Checkpoint #3 \n";   Dh.checknan();        
        
        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Anti-Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = neswdiag_il[id];
            m = neswdiag_im[id];

            // std::cout << "Before ( " << l << "," << m << ") \n";

            MakeGH(Din(l,m),G,H,l);

            if (m == 0)         // Left wall, no l + 1, m - 1
            {
                // std::cout << "1: ( " << l << "," << m << ")" << std::endl;
                if (l > 1)                  {    Em *= C3[l];    Dh(l-1,m+1) += H.mxy_matrix(Em);        Em *= 1.0/C3[l];}
            }
            else if (m == m0)   // Right boundary, no l - 1, m + 1
            {
                // std::cout << "2: ( " << l << "," << m << ")" << std::endl;
                if (l < l0)                 {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxy_matrix(Ep);        Ep *= 1.0/C2(l,m);}
            }
            else
            {
                // std::cout << "3: ( " << l << "," << m << ")" << std::endl;

                if (m > 1 && l < l0)        {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxy_matrix(Ep);        Ep *= 1.0/C2(l,m);}
                if (l - 1 != m && l != m)   {   Em *= C3[l];    Dh(l-1,m+1) += H.mxy_matrix(Em);        Em *= 1.0/C3[l];}
            }
            // std::cout << "\n Checkpoint ( " << l << "," << m << ") ... ";   Dh.checknan();  std::cout << "passed \n" << std::endl;
        }
        // std::cout << "Checkpoint #4";   Dh.checknan();
    }

    

    //  -------------------------------------------------------- //
    //  Do the boundaries between the chunks
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(f_start.size()-1)
    // for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {  
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();

        /// Local variables for each thread
        complex<double> ii(0.0,1.0);

        Array2D<complex<double> > Ex(FEx.array());
        Array2D<complex<double> > Em(FEz.array());
        Em *= (-1.0)*ii;
        Em += FEy.array();
        Array2D<complex<double> > Ep(FEz.array());
        Ep *= ii;
        Ep += FEy.array();

        Ex *= Din.q();
        Em *= Din.q();
        Ep *= Din.q();

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        SHarmonic2D G(Din(0,0)), H(Din(0,0)), TMP(Din(0,0));

        size_t l(0),m(0);

        if (this_thread < f_start.size() - 1) 
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      Vertical loop, boundaries between threads
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {

                l = dist_il[id];
                m = dist_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (l == m)         // Diagonal, no l - 1
                {
                    if (l < l0) {   Ex *= A1(l,0);   Dh(l+1,0) += G.mxy_matrix(Ex);     Ex *= 1.0 / A1(l,0);}
                }
                else if (l == l0)   // Last l, no l + 1
                {
                                    Ex *= A2(l,0);  Dh(l-1,0) += H.mxy_matrix(Ex);      Ex *= 1.0 / A2(l,0);
                }
                else
                {
                                    Ex *= A2(l,0);              Dh(l-1,0) += H.mxy_matrix(Ex);
                                    Ex *= A1(l,0) / A2(l,0);    Dh(l+1,0) += G.mxy_matrix(Ex);  Ex *= 1.0 / A1(l,0);
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = nwsediag_il[id];
                m = nwsediag_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (m == 0)         // Top or Left, no l - 1, m - 1
                {
                                    Em *= C1[m];    Dh(l+1,m+1) += G.mxy_matrix(Em);        Em *= 1.0/C1[m];
                }
                else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
                {
                                    Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxy_matrix(Ep);        Ep *= 1.0/C4(l,m);
                }
                else
                {
                                    Em *= C1[m];    Dh(l+1,m+1) += G.mxy_matrix(Em);        Em *= 1.0/C1[m];
                    if (m > 1)  {   Ep *= C4(l,m);  Dh(l-1,m-1) += H.mxy_matrix(Ep);        Ep *= 1.0/C4(l,m);}
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Anti-Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = neswdiag_il[id];
                m = neswdiag_im[id];

                MakeGH(Din(l,m),G,H,l);

                if (m == 0)         // Left wall, no l + 1, m - 1
                {
                    if (l > 1)                  {   Em *= C3[l];    Dh(l-1,1) += H.mxy_matrix(Em);        Em *= 1.0/C3[l];}
                }
                else if (m == m0)   // Right boundary, no l - 1, m + 1
                {
                    if (l < l0)                 {   Ep *= C2(l,m0);  Dh(l+1,m0-1) += G.mxy_matrix(Ep);        Ep *= 1.0/C2(l,m0);}
                }
                else
                {
                    if (m > 1 && l < l0)        {   Ep *= C2(l,m);  Dh(l+1,m-1) += G.mxy_matrix(Ep);        Ep *= 1.0/C2(l,m);}
                    if (l - 1 != m && l != m)   {   Em *= C3[l];    Dh(l-1,m+1) += H.mxy_matrix(Em);        Em *= 1.0/C3[l];}
                }
            }
        }
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::f1only(const DistFunc2D& Din,
      const Field2D& FEx, const Field2D& FEy, const Field2D& FEz,
      DistFunc2D& Dh) 
    {
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------

        complex<double> ii(0.0,1.0);

        Array2D<complex<double> > Ex(FEx.array());
        Array2D<complex<double> > Em(FEz.array());
        Em *= (-1.0)*ii;
        Em += FEy.array();
        Array2D<complex<double> > Ep(FEz.array());
        Ep *= ii;
        Ep += FEy.array();

        Ex *= Din.q();
        Em *= Din.q();
        Ep *= Din.q();

        SHarmonic2D G(Din(0,0)), H(Din(0,0)), TMP(Din(0,0));

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeG00(Din(0,0),G);
        Ex *= A100; TMP = G;    Dh(1,0) += G.mxy_matrix(Ex);
        Em *= C100;             Dh(1,1) += TMP.mxy_matrix(Em);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(1,0),G,H,1);
        Ex *= A210 / A100;      Dh(0,0) += H.mxy_matrix(Ex);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;              H = H.mxy_matrix(Ep);     Dh(0,0) += H.Re();

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field::Implicit_Ex(const DistFunc2D& Din, const Field2D& FEx, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        Array2D<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

        SHarmonic2D G(pr.size(),FEx.numx(),FEx.numy()), H(pr.size(),FEx.numx(),FEx.numy());

//      m = 0, l = 0
        MakeG00(Din(0,0),G); 
        Ex *= A1(0,0);             Dh(1,0) += G.mxy_matrix(Ex);     

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Ex *= A2(1,0) / A1(0,0);   Dh(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0);   Dh(2,0) += G.mxy_matrix(Ex);

//      m = 0, l = 2
        MakeGH(Din(2,0),G,H,2);
        Ex *= A2(2,0) / A1(1,0);   Dh(1,0) += H.mxy_matrix(Ex);

//      m = 0, l = 3
        MakeGH(Din(3,0),G,H,3);
        Ex *= A2(3,0) / A2(2,0);   Dh(2,0) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ex *= A1(1,1) / A2(3,0);   Dh(2,1) += G.mxy_matrix(Ex); 

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ex *= A2(2,1) / A1(1,1);   Dh(1,1) += H.mxy_matrix(Ex);

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Ex *= A2(3,1) / A2(2,1);   Dh(2,1) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::Implicit_Ey(const DistFunc2D& Din, const Field2D& FEy, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        Array2D<complex<double> > Em(FEy.array());
        Array2D<complex<double> > Ep(FEy.array());


        Em *= Din.q();
        Ep *= Din.q();


        SHarmonic2D G(pr.size(),FEy.numx(),FEy.numy()), H(pr.size(),FEy.numx(),FEy.numy()),
        TMP(pr.size(),FEy.numx(),FEy.numy());

//      m = 0, l = 0
        MakeG00(Din(0,0),G); 
        Em *= C1[0];            Dh(1,1) += G.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        MakeGH(Din(2,0),G,H,2);
        Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxy_matrix(Em);

//      m = 0,  l = 3 
        MakeGH(Din(3,0),G,H,3);
        Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ep *= B2[1];             H = H.mxy_matrix(Ep); Dh(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;          Dh(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];     G = G.mxy_matrix(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Em *= C3[3]   / C1[1];   TMP = H;          Dh(2,2) += TMP.mxy_matrix(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxy_matrix(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(Din(2,2),G,H,2);
        Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxy_matrix(Ep); 

//      m = 2, l = 3
        MakeGH(Din(3,2),G,H,3);
        Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxy_matrix(Ep);

        size_t m0(Din.m0());
        if ( m0 > 2) { 
//          m = 3, l = 3
            MakeGH(Din(3,3),G,H,3);
            Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxy_matrix(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
//
////--------------------------------------------------------------
    void Electric_Field::Implicit_Ez(const DistFunc2D& Din, const Field2D& FEz, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);

        
        Array2D<complex<double> > Em(FEz.array());
        Em *= (-1.0)*ii;
        
        
        Array2D<complex<double> > Ep(FEz.array());
        Ep *= ii;
        // Ep += FEy.array();
        // 
        // Ex *= Din.q();;
        Em *= Din.q();;
        Ep *= Din.q();;
        
        SHarmonic2D G(pr.size(),FEz.numx(),FEz.numy()), H(pr.size(),FEz.numx(),FEz.numy()),
        TMP(pr.size(),FEz.numx(),FEz.numy());

//      m = 0, l = 0
        MakeG00(Din(0,0),G); 
        Em *= C1[0];            Dh(1,1) += G.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Em *= C1[1]   / C1[0];  Dh(2,1) += G.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        MakeGH(Din(2,0),G,H,2);
        Em *= C3[2]   / C1[1];  Dh(1,1) += H.mxy_matrix(Em);

//      m = 0,  l = 3 
        MakeGH(Din(3,0),G,H,3);
        Em *= C3[3]   / C3[2];  Dh(2,1) += H.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(Din(1,1),G,H,1);
        Ep *= B2[1];             H = H.mxy_matrix(Ep); Dh(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;              Dh(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];     G = G.mxy_matrix(Ep); Dh(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(Din(2,1),G,H,2);
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Dh(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(Din(3,1),G,H,3);
        Em *= C3[3]   / C1[1];   TMP = H;              Dh(2,2) += TMP.mxy_matrix(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxy_matrix(Ep); Dh(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(Din(2,2),G,H,2);
        Ep *= C4(2,2) / B2[3];        Dh(1,1) += H.mxy_matrix(Ep); 

//      m = 2, l = 3
        MakeGH(Din(3,2),G,H,3);
        Ep *= C4(3,2) / C4(2,2);      Dh(2,1)  += H.mxy_matrix(Ep);

        size_t m0(Din.m0());
        if ( m0 > 2) { 
//          m = 3, l = 3
            MakeGH(Din(3,3),G,H,3);
            Ep *= C4(3,3) / C4(3,2);  Dh(2,2)  += H.mxy_matrix(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::Implicit_Ex_f1only(const DistFunc2D& Din, const Field2D& FEx, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        Array2D<complex<double> > Ex(FEx.array());
        Ex *= Din.q();

        SHarmonic2D G(Din(0,0)), H(Din(0,0));

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Ex *= A100;             Dh(1,0) += G.mxy_matrix(Ex);

//      m = 0, l = 1
        MakeGH(Din(1,0),G,H,1);
        Ex *= A210 / A100;   Dh(0,0) += H.mxy_matrix(Ex);


    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field::Implicit_Ey_f1only(const DistFunc2D& Din, const Field2D& FEy, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

        Array2D<complex<double> > Em(FEy.array());
        Array2D<complex<double> > Ep(FEy.array());

        Em *= Din.q();
        Ep *= Din.q();

        SHarmonic2D G(Din(0,0)), H(Din(0,0));

//      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Em *= C100;            Dh(1,1) += G.mxy_matrix(Em);

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
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;             H = H.mxy_matrix(Ep); Dh(0,0) += H.Re();
    }
//
////--------------------------------------------------------------
    void Electric_Field::Implicit_Ez_f1only(const DistFunc2D& Din, const Field2D& FEz, DistFunc2D& Dh) {
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);

        Array2D<complex<double> > Em(FEz.array());
        Em *= (-1.0)*ii;


        Array2D<complex<double> > Ep(FEz.array());
        Ep *= ii;


        SHarmonic2D G(Din(0,0)), H(Din(0,0));

        Em *= Din.q();
        Ep *= Din.q();

    //      m = 0, l = 0
        MakeG00(Din(0,0),G);
        Em *= C100;            Dh(1,1) += G.mxy_matrix(Em);

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
        MakeGH(Din(1,1),G,H,1);
        Ep *= B211;             H = H.mxy_matrix(Ep); Dh(0,0) += H.Re();
    }
//--------------------------------------------------------------


//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Magnetic_Field::Magnetic_Field(size_t Nl, size_t Nm,
        valarray<double> dp)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
    : A1(Nm+1), B1(Nl+1), A2(Nl+1,Nm+1), A3(0.5),
        f_start(Input::List().ompthreads),f_end(Input::List().ompthreads),//, sigma(Nx)//, killedbyPML((1.0,0.0),Nx)
        dist_il((Nm+1)*(2*Nl-Nm+2)/2),dist_im((Nm+1)*(2*Nl-Nm+2)/2)
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
        {
            for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m)
            {
                lc = complex<double>(l);
                mc = complex<double>(m);
                A2(l,m) = (-0.5)*(lc+1.0-mc)*(lc+mc);
            }
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

        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        // Prepare chunk indices for OpenMP 
        size_t num_threads = Input::List().ompthreads;
        size_t num_dists = (Nm+1)*(2*Nl-Nm+2)/2;
        size_t fsperthread = static_cast<size_t>(num_dists/num_threads - 2);
        
        if (num_threads > 1)
        { 
            // std::cout << "\n\n numthreads = " << num_threads;

            f_start[0] = 3;
            f_end[0]   = fsperthread-2;

            // std::cout << "\nls[" << 0 << "]=" << f_start[0] << "\n\n";
            // std::cout << "le[" << 0 << "]=" << f_end[0] << "\n\n";

            for (size_t i(1); i<num_threads-1; ++i)
            {        
                f_start[i] = f_end[i-1] + 2;
                f_end[i]   = f_start[i] + fsperthread;
                
                // std::cout << "\nls[" << i << "]=" << f_start[i] << "\n\n";
                // std::cout << "le[" << i << "]=" << f_end[i] << "\n\n";
            }

            f_start[num_threads-1] = f_end[num_threads-2] + 2;
            f_end[num_threads-1]   = num_dists;     

            // std::cout << "\nls[" << num_threads-1 << "]=" << f_start[num_threads-1] << "\n\n";
            // std::cout << "le[" << num_threads-1 << "]=" << f_end[num_threads-1] << "\n\n";
        }
        else
        {
            f_start[0] = 3;
            f_end[0] = num_dists;
        }

        size_t il(0), im(0);
        for (size_t id(0); id < num_dists; ++id)
        {
            // std::cout << "\n0(id,l,m) =  " << id << "," << il << "," << im << "," << "\n";

            dist_il[id] = il;
            dist_im[id] = im;

            if (im < il && im < Nm)
            {
                ++im;
            }
            else
            {
                ++il;
                im = 0;
            }
        }


    }
//--------------------------------------------------------------


//--------------------------------------------------------------
  void Magnetic_Field::operator()(const DistFunc1D& Din,
   const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
   DistFunc1D& Dh) {
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------

    complex<double> ii(0.0,1.0);

    #pragma omp parallel num_threads(Input::List().ompthreads)
    {
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

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        size_t this_thread  = omp_get_thread_num();

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);

        SHarmonic1D FLM(Din(0,0)), TMP(Din(0,0));

        Bp *= A3;

        size_t l(0),m(0);

        if (this_thread == 0)
        {
            FLM = Din(1,0);                 Dh(1,1) += FLM.mxaxis(Bp);
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      l = 1, m = 1
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            FLM = Din(1,1); Bx *= A1[1];    Dh(1,1) += FLM.mxaxis(Bx);  Bx /= A1[1];
        
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 1, l = 1
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            FLM = Din(1,1); Bm *= B1[1];    FLM = FLM.mxaxis(Bm);     Dh(1,0) += FLM.Re();  Bm /= B1[1];
        }

        // ----------------------------------------- //
        //              Do the chunks
        // ----------------------------------------- //       
        
        for (size_t id = f_start_thread; id < f_end_thread ; ++id)
        {   

            l = dist_il[id];
            m = dist_im[id];

            FLM = Din(l,m);
            TMP = FLM;

            if (l == m || m == m0)         // Diagonal, no m + 1
            {
                            Bx *= A1[m];        Dh(l,m  ) += TMP.mxaxis(Bx);        Bx /= A1[m];
                            Bm *= A2(l,m);      Dh(l,m-1) += FLM.mxaxis(Bm);        Bm /= A2(l,m);                   
            }
            else if (m == 0)    // m = 0, no m - 1
            {   
                                                Dh(l,1) += TMP.mxaxis(Bp);
            }
            else if (m == 1)
            {
                            Bm *= B1[l];        TMP = TMP.mxaxis(Bm);               Dh(l,0) += TMP.Re();        Bm /= B1[l];
                TMP = FLM;                      Dh(l,2) += TMP.mxaxis(Bp);
                TMP = FLM;  Bx *= A1[1];        Dh(l,1) += TMP.mxaxis(Bx);          Bx /= A1[1];
            }
            else
            {
                                                Dh(l,m+1) += TMP.mxaxis(Bp);
                TMP = FLM; Bx *= A1[m];         Dh(l,m  ) += TMP.mxaxis(Bx);        Bx /= A1[m];
                TMP = FLM; Bm *= A2(l,m);       Dh(l,m-1) += TMP.mxaxis(Bm);        Bm /= A2(l,m);                   
            } 
        }
    }

    // ----------------------------------------- //
    //          Boundaries between chunks
    // ----------------------------------------- //

    #pragma omp parallel for num_threads(f_start.size()-1)
    for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {
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

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        size_t this_thread  = omp_get_thread_num();

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);

        SHarmonic1D FLM(Din(0,0)), TMP(Din(0,0));

        Bp *= A3;

        size_t l(0),m(0);

        for (size_t id = f_end[threadboundaries]; id < f_start[threadboundaries+1]; ++id)
        {   
            l = dist_il[id];
            m = dist_im[id];

            FLM = Din(l,m);
            TMP = FLM;

            if (l == m || m == m0)         // Diagonal, no m + 1
            {
                            Bx *= A1[m];        Dh(l,m  ) += TMP.mxaxis(Bx);        Bx /= A1[m];
                            Bm *= A2(l,m);      Dh(l,m-1) += FLM.mxaxis(Bm);        Bm /= A2(l,m);                   
            }
            else if (m == 0)    // m = 0, no m - 1
            {
                                                Dh(l,1) += TMP.mxaxis(Bp);
            }
            else if (m == 1)
            {
                            Bm *= B1[l];        TMP = TMP.mxaxis(Bm);               Dh(l,0) += TMP.Re();        Bm /= B1[l];           
                TMP = FLM;                      Dh(l,2) += TMP.mxaxis(Bp);
                TMP = FLM;  Bx *= A1[1];        Dh(l,1) += TMP.mxaxis(Bx);          Bx /= A1[1];        
            }
            else
            {
                                                Dh(l,m+1) += TMP.mxaxis(Bp);
                TMP = FLM; Bx *= A1[m];         Dh(l,m  ) += TMP.mxaxis(Bx);        Bx /= A1[m];
                TMP = FLM; Bm *= A2(l,m);       Dh(l,m-1) += TMP.mxaxis(Bm);        Bm /= A2(l,m);                   
            }
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Magnetic_Field::operator()(const DistFunc2D& Din, 
       const Field2D& FBx, const Field2D& FBy, const Field2D& FBz, 
       DistFunc2D& Dh) 
{
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------

//     complex<double> ii(0.0,1.0);

//     Array2D<complex<double> > Bx(FBx.array());
//     Array2D<complex<double> > Bm(FBy.array());
//     Bm *= (-1.0)*ii;
//     Bm += FBz.array();
//     Array2D<complex<double> > Bp(FBy.array());
//     Bp *= ii;
//     Bp += FBz.array();

//     Bx *= Din.q();;
//     Bm *= Din.q();;
//     Bp *= Din.q();;

//     size_t l0(Din.l0());
//     size_t m0(Din.m0());

//     SHarmonic2D FLM(Din(0,0));

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 0, 1 < l < l0+1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     Bp *= A3;
//     for (size_t l(1); l < l0+1; ++l){
//         FLM = Din(l,0);      Dh(l,1) += FLM.mxy_matrix(Bp);
//     }

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l = 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     FLM = Din(1,1); Bx *= A1[1];                           Dh(1,1) += FLM.mxy_matrix(Bx); 
//     FLM = Din(1,1); Bm *= B1[1]; FLM = FLM.mxy_matrix(Bm); Dh(1,0) += FLM.Re(); 

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = 1, l > 1
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t l(2); l < l0+1; ++l){
//         FLM = Din(l,1);                                                Dh(l,2) += FLM.mxy_matrix(Bp); 
//         FLM = Din(l,1);                                                Dh(l,1) += FLM.mxy_matrix(Bx); 
//         FLM = Din(l,1); Bm *= B1[l]/B1[l-1]; FLM = FLM.mxy_matrix(Bm); Dh(l,0) += FLM.Re(); 
//     }        
//     Bm *= 1.0/B1[l0];        

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m > 1, l = m
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t m(2); m < m0; ++m){
//         FLM = Din(m,m); Bx *= A1[m]/A1[m-1];                   Dh(m,m  ) += FLM.mxy_matrix(Bx); 
//         FLM = Din(m,m); Bm *= A2(m,m);                         Dh(m,m-1) += FLM.mxy_matrix(Bm); 
//         for (size_t l(m+1); l < l0+1; ++l){
//             FLM = Din(l,m);                                    Dh(l,m+1) += FLM.mxy_matrix(Bp); 
//             FLM = Din(l,m);                                    Dh(l,m  ) += FLM.mxy_matrix(Bx); 
//             FLM = Din(l,m); Bm *= A2(l,m)/A2(l-1,m);           Dh(l,m-1) += FLM.mxy_matrix(Bm); 
//         }
//         Bm *= 1.0/A2(l0,m);
//     }

// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
// //      m = m0, l >= m0
// // - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     FLM = Din(m0,m0); Bx *= A1[m0]/A1[m0-1];                   Dh(m0,m0)   += FLM.mxy_matrix(Bx); 
//     FLM = Din(m0,m0); Bm *= A2(m0,m0)/*/A2(l0,m0-1)*/;             Dh(m0,m0-1) += FLM.mxy_matrix(Bm);

//     for (size_t l(m0+1); l < l0+1; ++l)
//     {
//         FLM = Din(l,m0);                                     Dh(l,m0  )  += FLM.mxy_matrix(Bx); 
//         FLM = Din(l,m0); Bm *= A2(l,m0)/A2(l-1,m0);          Dh(l,m0-1)  += FLM.mxy_matrix(Bm); 
//     }


    

    #pragma omp parallel num_threads(Input::List().ompthreads)
    {
        complex<double> ii(0.0,1.0);

        Array2D<complex<double> > Bx(FBx.array());
        Array2D<complex<double> > Bm(FBy.array());
        Bm *= (-1.0)*ii;
        Bm += FBz.array();
        Array2D<complex<double> > Bp(FBy.array());
        Bp *= ii;
        Bp += FBz.array();

        Bx *= Din.q();
        Bm *= Din.q();
        Bp *= Din.q();

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        size_t this_thread  = omp_get_thread_num();

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);

        SHarmonic2D FLM(Din(0,0)), TMP(Din(0,0));

        Bp *= A3;

        size_t l(0),m(0);

        if (this_thread == 0)
        {
            FLM = Din(1,0);                 Dh(1,1) += FLM.mxy_matrix(Bp);
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      l = 1, m = 1
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            FLM = Din(1,1); Bx *= A1[1];    Dh(1,1) += FLM.mxy_matrix(Bx);  Bx *= 1.0/A1[1];

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 1, l = 1
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            FLM = Din(1,1); Bm *= B1[1];    FLM = FLM.mxy_matrix(Bm);     Dh(1,0) += FLM.Re();  Bm *= 1.0/B1[1];
        }

        // ----------------------------------------- //
        //              Do the chunks
        // ----------------------------------------- //       
        
        for (size_t id = f_start_thread; id < f_end_thread ; ++id)
        {   


            l = dist_il[id];
            m = dist_im[id];

            FLM = Din(l,m);
            TMP = FLM;

            if (l == m || m == m0)         // Diagonal or last m, no m + 1
            {
                            Bx *= A1[m];        Dh(l,m  ) += TMP.mxy_matrix(Bx);        Bx *= 1.0/A1[m];
                            Bm *= A2(l,m);      Dh(l,m-1) += FLM.mxy_matrix(Bm);        Bm *= 1.0/A2(l,m);                   
            }
            else if (m == 0)    // m = 0, no m - 1
            {   
                                                Dh(l,1) += TMP.mxy_matrix(Bp);
            }
            else if (m == 1)
            {
                            Bm *= B1[l];        TMP = TMP.mxy_matrix(Bm);               Dh(l,0) += TMP.Re();        Bm *= 1.0/B1[l];           
                TMP = FLM;                      Dh(l,2) += TMP.mxy_matrix(Bp);
                TMP = FLM;  Bx *= A1[1];        Dh(l,1) += TMP.mxy_matrix(Bx);          Bx *= 1.0/A1[1];
            }
            else
            {
                                                Dh(l,m+1) += TMP.mxy_matrix(Bp);
                TMP = FLM; Bx *= A1[m];         Dh(l,m  ) += TMP.mxy_matrix(Bx);        Bx *= 1.0/A1[m];
                TMP = FLM; Bm *= A2(l,m);       Dh(l,m-1) += TMP.mxy_matrix(Bm);        Bm *= 1.0/A2(l,m);                   
            } 
        }
    }

    // ----------------------------------------- //
    //          Boundaries between chunks
    // ----------------------------------------- //

    #pragma omp parallel for num_threads(f_start.size()-1)
    for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {
        /// Local variables for each thread
        complex<double> ii(0.0,1.0);
        
        Array2D<complex<double> > Bx(FBx.array());
        Array2D<complex<double> > Bm(FBy.array());
        Bm *= (-1.0)*ii;
        Bm += FBz.array();
        Array2D<complex<double> > Bp(FBy.array());
        Bp *= ii;
        Bp += FBz.array();

        Bx *= Din.q();
        Bm *= Din.q();
        Bp *= Din.q();

        size_t l0(Din.l0());
        size_t m0(Din.m0());

        size_t this_thread  = omp_get_thread_num();

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);

        SHarmonic2D FLM(Din(0,0)), TMP(Din(0,0));

        Bp *= A3;

        size_t l(0),m(0);

        for (size_t id = f_end[threadboundaries]; id < f_start[threadboundaries+1]; ++id)
        {   
            l = dist_il[id];
            m = dist_im[id];

            FLM = Din(l,m);
            TMP = FLM;

            if (l == m || m == m0)         // Diagonal, no m + 1
            {
                            Bx *= A1[m];        Dh(l,m  ) += TMP.mxy_matrix(Bx);        Bx *= 1.0/A1[m];
                            Bm *= A2(l,m);      Dh(l,m-1) += FLM.mxy_matrix(Bm);        Bm *= 1.0/A2(l,m);                   
            }
            else if (m == 0)    // m = 0, no m - 1
            {
                                                Dh(l,1) += TMP.mxy_matrix(Bp);
            }
            else if (m == 1)
            {
                            Bm *= B1[l];        TMP = TMP.mxy_matrix(Bm);               Dh(l,0) += TMP.Re();        Bm *= 1.0/B1[l];           
                TMP = FLM;                      Dh(l,2) += TMP.mxy_matrix(Bp);
                TMP = FLM;  Bx *= A1[1];        Dh(l,1) += TMP.mxy_matrix(Bx);          Bx *= 1.0/A1[1];        
            }
            else
            {
                                                Dh(l,m+1) += TMP.mxy_matrix(Bp);
                TMP = FLM; Bx *= A1[m];         Dh(l,m  ) += TMP.mxy_matrix(Bx);        Bx *= 1.0/A1[m];
                TMP = FLM; Bm *= A2(l,m);       Dh(l,m-1) += TMP.mxy_matrix(Bm);        Bm *= 1.0/A2(l,m);                   
            }
        }
    }

}
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
void Magnetic_Field::f1only(const DistFunc1D& Din,
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

    SHarmonic1D FLM(Din(0,0));

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


}
//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
void Magnetic_Field::f1only(const DistFunc2D& Din, 
    const Field2D& FBx, const Field2D& FBy, const Field2D& FBz, 
    DistFunc2D& Dh) 
{
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------
    complex<double> ii(0.0,1.0);

    Array2D<complex<double> > Bx(FBx.array());
    Array2D<complex<double> > Bm(FBy.array());
    Bm *= (-1.0)*ii;
    Bm += FBz.array();
    Array2D<complex<double> > Bp(FBy.array());
    Bp *= ii;
    Bp += FBz.array();

    Bx *= Din.q();;
    Bm *= Din.q();;
    Bp *= Din.q();;

    size_t l0(B1.size()-1);
    size_t m0(A1.size()-1);

    SHarmonic2D FLM(Din(0,0));

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, 1 < l < l0+1
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Bp *= A3;
    for (size_t l(1); l < l0+1; ++l)
    {
        FLM = Din(l,0);      Dh(l,1) += FLM.mxy_matrix(Bp);
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
// - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FLM = Din(1,1); Bx *= A1[1];                           Dh(1,1) += FLM.mxy_matrix(Bx); 
    FLM = Din(1,1); Bm *= B1[1]; FLM = FLM.mxy_matrix(Bm); Dh(1,0) += FLM.Re(); 
}

//**************************************************************
//--------------------------------------------------------------
Spatial_Advection::Spatial_Advection(size_t Nl, size_t Nm,
    valarray<double> dp,
    double xmin, double xmax, size_t Nx,
    double ymin, double ymax, size_t Ny)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        :   A1(Nl+1,Nm+1), A2(Nl+1,Nm+1),
            C2(Nl+1,Nm+1), C4(Nl+1,Nm+1),
            B1(Nl+1), B2(Nl+1), C1(Nl+1), C3(Nl+1), 
            vr(Algorithms::MakeCAxis(complex<double>(0.0),
              complex<double>(1.0),
              dp.size())),
            f_start(Input::List().ompthreads),f_end(Input::List().ompthreads),//, sigma(Nx)//, killedbyPML((1.0,0.0),Nx)
            dist_il((Nm+1)*(2*Nl-Nm+2)/2),dist_im((Nm+1)*(2*Nl-Nm+2)/2),
            nwsediag_il((Nm+1)*(2*Nl-Nm+2)/2),nwsediag_im((Nm+1)*(2*Nl-Nm+2)/2),
            neswdiag_il((Nm+1)*(2*Nl-Nm+2)/2),neswdiag_im((Nm+1)*(2*Nl-Nm+2)/2)
    {
// ------------------------------------------------------------------------ // 

    // ------------------------------------------------------------------------ // 
    // Make velocity grid
    // ------------------------------------------------------------------------ // 
            vr[0] = 0.5*dp[0];
            for (size_t ip(1); ip < dp.size(); ++ip)
            {
                vr[ip]  = dp[ip];
                vr[ip] += vr[ip-1];
            // std::cout << "\n Spr[" << ip << "] = " << vr[ip] << std::endl;
            }

            if (Input::List().relativity)
            {
                for (size_t i(0); i < vr.size(); ++i) {
                    vr[i] = vr[i]/(sqrt(1.0+vr[i]*vr[i]));   // Turn pr into vr
                }
            }

            double idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx)); // -1/(2dx)
            
            complex<double> lc, mc;

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
            A00 = complex<double>(-idx);
            A10 = complex<double>(-idx/3.0);
            A20 = complex<double>(-idx*2.0/5.0);

        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        //// 2D begins here
        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
        // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
            
            double idy = (-1.0) / (2.0*(ymax-ymin)/double(Ny)); // -1/(2dy) 
        
        //       - - - - - - - - - - - - - - - - - - - - - - - - - - -        
        //       Calculate the "B1, B2" parameters
        //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t l(0); l<Nl+1; ++l){
               lc = static_cast< complex<double> >(l);
               B1[l] = idy * (lc + 1.0) * lc / (2.0*lc + 1.0);
               B2[l] = (-1.0)*B1[l];
           }
           B1[0] = 1.0;

    //       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //       Calculate the "C1, C3" parameters
    //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
           for (size_t l(0); l<Nl+1; ++l){
               lc = static_cast< complex<double> >(l);
               C1[l] = (-0.5) * idy / (2.0*lc + 1.0);
               C3[l] = (-1.0) * C1[l];
           }
    //       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //       Calculate the "C2, C4" parameters
    //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
           for (size_t l(0); l<Nl+1; ++l){
               for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
                   lc = static_cast< complex<double> >(l);
                   mc = static_cast< complex<double> >(m);
                   C2(l,m) = idy * 0.5 * (lc + 2.0 - mc)*(lc - mc + 1.0) / (2.0*lc + 1.0);
                   C4(l,m) = idy * (-0.5) * (lc + mc - 1.0)*(lc + mc) / (2.0*lc + 1.0);
               }
           }
    // ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
    // OpenMP stuff
        size_t num_threads = Input::List().ompthreads;
        size_t num_dists = (((Nm+1)*(2*Nl-Nm+2))/2);
        size_t fsperthread = static_cast<size_t>(num_dists/num_threads - 2);

        if (num_threads > 272 || num_threads < 1) num_threads = 1;

    // // std::cout << "\n numthreads = " << num_threads << "\n\n";    
    //     size_t lsperthread = static_cast<size_t>(Nl/num_threads - 2);
    
        if (num_threads > 1) 
        {       
            f_start[0] = 1;
            f_end[0]   = fsperthread - 2;           // Remember that f_end isn't actually processed
                                                    // until the boundaries so that a gap of 2 is 
                                                    // actually a gap of 3

            // if (num_threads > 1) std::cout << "\nls[" << 0 << "]=" << f_start[0] << "\n\n";
            // if (num_threads > 1) std::cout << "le[" << 0 << "]=" << f_end[0] << "\n\n";

            for (size_t i(1); i<num_threads-1; ++i)
            {        
                f_start[i] = f_end[i-1] + 2;
                f_end[i]   = f_start[i] + fsperthread;
                // if (num_threads > 1)    std::cout << "\nls[" << i << "]=" << f_start[i] << "\n\n";
                // if (num_threads > 1)    std::cout << "le[" << i << "]=" << f_end[i] << "\n\n";
            }

            f_start[num_threads-1] = f_end[num_threads-2] + 2;
            f_end[num_threads-1]   = num_dists; 
        }
        else
        {
            f_start[0] = 1;
            f_end[0]   = num_dists; 
        }

        

        size_t il(0), im(0);
        for (size_t id(0); id < num_dists; ++id)
        {
            dist_il[id] = il;
            dist_im[id] = im;

            if (il < Nl)
            {
                ++il;
            }
            else
            {
                ++im;
                il = im;
            }
        }
        

        il = 0; im = 0;
        size_t diag(0);
        for (size_t id(0); id < num_dists; ++id)
        {
            // std::cout << "\n1(id,l,m) =  " << id << "," << il << "," << im << "," << "\n";
            nwsediag_il[id] = il;
            nwsediag_im[id] = im;

            if (il < Nl && im < Nm && im <= il)
            {
                ++il; ++im;
            }
            else
            {
                ++diag;   
                im = 0;
                il = diag;                
            }
        }


        il = 0; im = 0; diag = 0;

        for (size_t id(0); id < num_dists; ++id)
        {
            // std::cout << "\n2(id,l,m) =  " << id << "," << il << "," << im << "," << "\n";
            neswdiag_il[id] = il;
            neswdiag_im[id] = im;

            if (il < Nl && im > 0)
            {
                ++il; --im;
                diag++;
            }
            else
            {
                if (neswdiag_il[id-diag] == neswdiag_im[id-diag])
                {
                    im = neswdiag_im[id-diag];
                    il = neswdiag_il[id-diag] + 1;
                    diag = 0;
                }
                else if (neswdiag_il[id-diag] == neswdiag_im[id-diag] + 1)
                {
                    if (neswdiag_im[id-diag] < Nm)
                    {
                        im = neswdiag_im[id-diag] + 1;
                        il = neswdiag_il[id-diag];
                    }
                    else 
                    {
                        im = Nm;
                        il = neswdiag_il[id-diag] + 1;
                    }
                    diag = 0;
                }
                else
                {
                    im = Nm;
                    il = neswdiag_il[id-diag] + 1;
                    diag = 0;
                }
            }
        }
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // ----- 
// ----- // ----- // ----- // ----- // ----- // ----- // ----- // -----     
                
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//   Advection in x
   void Spatial_Advection::operator()(const DistFunc2D& Din, DistFunc2D& Dh) {
//--------------------------------------------------------------

//     valarray<complex<double> > vt(vr); vt *= 1.0/Din.mass(); 
//         // vt += fluidvelocity;       
//     size_t l0(Din.l0());
//     size_t m0(Din.m0());

//     SHarmonic2D fd1(Din(0,0)), fd2(Din(0,0));

//     for (size_t m(0); m < ((m0<l0)?(m0+1):m0); ++m){

//         fd1 = Din(m,m);     fd1 = fd1.Dx();
//         vt *= A1(m,m);                          Dh(m+1,m) += fd1.mpaxis(vt);

//         for (size_t l(m+1); l < l0; ++l) {
//             fd1 = Din(l,m);            fd1 = fd1.Dx();  
//             vt *= A2(l,m)/A1(l-1,m);    fd2 = fd1;  Dh(l-1,m) += fd1.mpaxis(vt);  
//             vt *= A1(l,m)/A2(l  ,m);                Dh(l+1,m) += fd2.mpaxis(vt);  
//         }

//         fd1 = Din(l0,m);            fd1 = fd1.Dx();
//         vt *= A2(l0,m)/A1(l0-1,m);                 Dh(l0-1,m) += fd1.mpaxis(vt);  
//         vt *= 1.0/A2(l0,m); 
//     }

//         //       m = 0, advection in y
// //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     fd1 = Din(0,0);     fd1 = fd1.Dy();
//     vt *= C1[0];                            Dh(1,1) += fd1.mpaxis(vt);
//     fd1 = Din(1,0);     fd1 = fd1.Dy();
//     vt *= C1[1]/C1[0];                      Dh(2,1) += fd1.mpaxis(vt);
//     for (size_t l(2); l < l0; ++l) 
//     {
//         fd1 = Din(l,0);     fd1 = fd1.Dy();
//         vt *= C3[l]/C1[l-1];   fd2 = fd1;   Dh(l-1,1) += fd1.mpaxis(vt);
//         vt *= C1[l]/C3[l];                  Dh(l+1,1) += fd2.mpaxis(vt);
//     }
//     fd1 = Din(l0,0);    fd1 = fd1.Dy();
//     vt *= C3[l0]/C1[l0-1];              Dh(l0-1,1) += fd1.mpaxis(vt);
//     vt *= 1.0   /C3[l0];
// //       - - - - - - - - - - - - - - - - - - - - - - - - - - -


// //       m = 1, advection in y
// //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     fd1 = Din(1,1);        fd1 = fd1.Dy();
//     vt *= C1[1];           fd2 = fd1;                                 Dh(2,2) += fd1.mpaxis(vt);
//     vt *= B2[1]/C1[1];     fd1 = fd2;  fd2 = fd2.mpaxis(vt);  Dh(0,0) += fd2.Re();
//     vt *= B1[1]/B2[1];                 fd1 = fd1.mpaxis(vt);  Dh(2,0) += fd1.Re();

//     fd1 = Din(2,1);        fd1 = fd1.Dy();
//     vt *= C1[2]/B1[1];     fd2 = fd1;                                 Dh(3,2) += fd1.mpaxis(vt);
//     vt *= B2[2]/C1[2];     fd1 = fd2;  fd2 = fd2.mpaxis(vt);  Dh(1,0) += fd2.Re();
//     vt *= B1[2]/B2[2];                 fd1 = fd1.mpaxis(vt);  Dh(3,0) += fd1.Re();

//     for (size_t l(3); l < l0; ++l)
//     {
//         fd1 = Din(l,1);        fd1 = fd1.Dy();                           
//         vt *= C3[l]/B1[l-1];   fd2 = fd1;                                Dh(l-1,2) += fd1.mpaxis(vt);
//         vt *= C1[l]/C3[l];     fd1 = fd2;                                Dh(l+1,2) += fd2.mpaxis(vt);
//         vt *= B2[l]/C1[l];     fd2 = fd1;  fd1 = fd1.mpaxis(vt); Dh(l-1,0) += fd1.Re();
//         vt *= B1[l]/B2[l];                 fd2 = fd2.mpaxis(vt); Dh(l+1,0) += fd2.Re();
//     }

//     fd1 = Din(l0,1);     fd1 = fd1.Dy();                           
//     vt *= C3[l0]/B1[l0-1];     fd2 = fd1;                                Dh(l0-1,2) += fd1.mpaxis(vt);
//     vt *= B2[l0]/C3[l0];       fd2 = fd2.mpaxis(vt);             Dh(l0-1,0) += fd2.Re();
//     vt *= 1.0   /B2[l0];
//     //       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//     //       m > 1, advection in y
//     //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     for (size_t m(2); m < m0; ++m)
//     {
//         //           m > 1, l = m
//         fd1 = Din(m,m);       fd1 = fd1.Dy();
//         vt *= C4(m,m);           fd2 = fd1;     Dh(m-1,m-1) += fd1.mpaxis(vt);
//         vt *= C2(m,m)/C4(m,m);   fd1 = fd2;     Dh(m+1,m-1) += fd2.mpaxis(vt);
//         vt *= C1[m]  /C2(m,m);                  Dh(m+1,m+1) += fd1.mpaxis(vt);

//     //           m > 1, l = m+1
//         fd1 = Din(m+1,m);       fd1 = fd1.Dy();
//         vt *= C4(m+1,m)/C1[m];     fd2 = fd1;     Dh(m  ,m-1) += fd1.mpaxis(vt);
        
//         if (m+1 < l0) {
//             vt *= C2(m+1,m)/C4(m+1,m);   fd1 = fd2;   Dh(m+2,m-1) += fd2.mpaxis(vt);
//             vt *= C1[m+1]  /C2(m+1,m);                Dh(m+2,m+1) += fd1.mpaxis(vt);

// //               m > 1, 3 < l < l0
//             for (size_t l(m+2); l < l0; ++l){
//                 fd1 = Din(l,m);       fd1 = fd1.Dy();
//                 vt *= C4(l,m)/C1[m+1];   fd2 = fd1;     Dh(l-1,m-1) += fd1.mpaxis(vt);
//                 vt *= C2(l,m)/C4(l,m);   fd1 = fd2;     Dh(l+1,m-1) += fd2.mpaxis(vt);
//                 vt *= C3[l]  /C2(l,m);   fd2 = fd1;     Dh(l-1,m+1) += fd1.mpaxis(vt);
//                 vt *= C1[l]  /C3[l];                    Dh(l+1,m+1) += fd2.mpaxis(vt);
//             }

//             fd1 = Din(l0,m);       fd1 = fd1.Dy();
//             vt *= C4(l0,m)/C1[l0-1];  fd2 = fd1;       Dh(l0-1,m-1) += fd1.mpaxis(vt);
//             vt *= C3[l0]/C4(l0,m);                     Dh(l0-1,m+1) += fd2.mpaxis(vt);
//             vt *= 1.0/C3[l0];
//         }
//         else {
//            vt *= 1.0/C4(m+1,m);
//         }
//    }
// //       - - - - - - - - - - - - - - - - - - - - - - - - - - -

// //       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     fd1 = Din(m0,m0);          fd1 = fd1.Dy();
//     vt *= C4(m0,m0);           fd2 = fd1;     Dh(m0-1,m0-1) += fd1.mpaxis(vt);
//     if (m0 < l0) {
//         vt *= C2(m0,m0)/C4(m0,m0);            Dh(m0+1,m0-1) += fd2.mpaxis(vt);
//         for (size_t l(m0+1); l < l0; ++l){
//             fd1 = Din(l,m0);          fd1 = fd1.Dy();
//             vt *= C4(l,m0)/C2(l-1,m0);   fd2 = fd1;     Dh(l-1,m0-1) += fd1.mpaxis(vt);
//             vt *= C2(l,m0)/C4(l,  m0);                  Dh(l+1,m0-1) += fd2.mpaxis(vt);
//         }
//         fd1 = Din(l0,m0);           fd1 = fd1.Dy();
//         vt *= C4(l0,m0)/C2(l0-1,m0);   Dh(l0-1,m0-1) += fd1.mpaxis(vt);
//     }
    
    //  ------------------------------------------------------- //
    //   Because each iteration in the loop modifies + and - 1
    //   The parallelization is performed in chunks and boundaries
    //   are taken care of later
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(Input::List().ompthreads)
    {   
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();
        size_t f_start_thread(f_start[this_thread]); ///< Chunk starts here
        size_t f_end_thread(f_end[this_thread]);     ///< Chunk ends here

        /// Local variables for each thread
        valarray<complex<double> > vtemp(vr); 
        vtemp /= Din.mass();
        size_t l(0),m(0);
        size_t l0(Din.l0());
        size_t m0(Din.m0());

        SHarmonic2D fd1(Din(0,0)),fd2(Din(0,0));

        if (this_thread == 0)
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 0, l = 0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            fd1 = Din(0,0);                         fd1 = fd1.Dx(Input::List().dbydx_order);
            vtemp *= A1(0,0);                       Dh(1,0) += fd1.mpaxis(vtemp);
            vtemp /= A1(0,0);
        
            // std::cout << "\n Checkpoint #0 \n";   Dh.checknan(); 
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      m = 1 loop, 1 <= l < l0
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            vtemp *= B1[0];
            for (size_t il = 1; il < l0; ++il)
            {
                fd1 = Din(il,1);                fd1 = fd1.Dy(Input::List().dbydy_order);;
                vtemp *= B2[il]/B1[il-1];       fd2 = fd1;      fd1 = fd1.mpaxis(vtemp);    Dh(il-1,0) += fd1.Re();
                vtemp *= B1[il]/B2[il];                         fd2 = fd2.mpaxis(vtemp);    Dh(il+1,0) += fd2.Re();

                // std::cout << "\n Checkpoint (" << il  << ")\n";   Dh.checknan(); std::cout << ".. passed \n";
            }
            fd1 = Din(l0,1);                            fd1 = fd1.Dy(Input::List().dbydy_order);;
            vtemp *= B2[l0]/B1[l0-1];                   fd1 = fd1.mpaxis(vtemp);    Dh(l0-1,0) += fd1.Re();

            vtemp *= 1.0/B2[l0];
        }

        // std::cout << "\n Checkpoint #1 \n";   Dh.checknan(); std::cout << ".. passed \n";
        
        //  -------------------------------------------------------- //
        //                      Do the chunks
        //  -------------------------------------------------------- //       
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Vertical loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {

            l = dist_il[id];
            m = dist_im[id];

            // std::cout << "\n (l,m) = " << l << ", " << m << " \n";

            fd1 = Din(l,m);     fd1 = fd1.Dx(Input::List().dbydx_order);  

            if (l == m)         // Diagonal, no l - 1
            {
                if (l < l0){    vtemp *= A1(m,m);           Dh(m+1,m) += fd1.mpaxis(vtemp);     vtemp /= A1(m,m);}
            }
            else if (l == l0)   // Last l, no l + 1
            {                
                                vtemp *= A2(l0,m);          Dh(l0-1,m) += fd1.mpaxis(vtemp);    vtemp /= A2(l0,m);
            }
            else
            {   
                fd2 = fd1;      vtemp *= A2(l,m);           Dh(l-1,m) += fd1.mpaxis(vtemp);
                                vtemp *= A1(l,m)/A2(l  ,m); Dh(l+1,m) += fd2.mpaxis(vtemp);     vtemp /= A1(l,m);
            }
        }
        // std::cout << "\n Checkpoint #2 \n";   Dh.checknan();    std::cout << ".. passed \n";

        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = nwsediag_il[id];
            m = nwsediag_im[id];

            fd1 = Din(l,m);     fd1 = fd1.Dy(Input::List().dbydy_order);;

            if (m == 0)         // Top or Left, no l - 1, m - 1
            {
                if (l < l0) {   vtemp *= C1[l];             Dh(l+1,m+1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C1[l];}
            }
            else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
            {
                                vtemp *= C4(l,m);           Dh(l-1,m-1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C4(l,m);
            }
            else
            {       
                fd2 = fd1;      vtemp *= C1[l];             Dh(l+1,m+1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C1[l]; 
                if (m>1)    {   vtemp *= C4(l,m);           Dh(l-1,m-1) += fd2.mpaxis(vtemp);   vtemp *= 1.0/C4(l,m);}            
            }
        }
    
        // std::cout << "\n Checkpoint #3 \n";   Dh.checknan();    std::cout << ".. passed \n";        
        #pragma omp barrier
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //      Anti-Diagonal loop, f_start < l < f_end
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {
            l = neswdiag_il[id];
            m = neswdiag_im[id];

            fd1 = Din(l,m);     fd1 = fd1.Dy(Input::List().dbydy_order);;

            if (m == 0)         // Left wall, no l + 1, m - 1
            {
                if (l > 1)  {               vtemp *= C3[l];         Dh(l-1,m+1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C3[l];}
            }
            else if (m == m0)   // Right boundary, no l - 1, m + 1
            {
                if (l < l0) {               vtemp *= C2(l,m);       Dh(l+1,m-1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C2(l,m);}
            }
            else
            {
                if (m > 1 && l < l0)        
                {   
                    fd2 = fd1;              vtemp *= C2(l,m);       Dh(l+1,m-1) += fd2.mpaxis(vtemp);    vtemp *= 1.0/C2(l,m);                    
                }
                if (l - 1 != m && l != m){  vtemp *= C3[l];         Dh(l-1,m+1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C3[l];}
            }
        }
    }

    

    //  -------------------------------------------------------- //
    //  Do the boundaries between the chunks
    //  -------------------------------------------------------- //
    #pragma omp parallel num_threads(f_start.size()-1)
    // for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {  
        /// Determine which chunk to do
        size_t this_thread  = omp_get_thread_num();

        /// Local variables for each thread
        valarray<complex<double> > vtemp(vr); 
        vtemp /= Din.mass();
        size_t l(0),m(0);
        size_t l0(Din.l0());
        size_t m0(Din.m0());

        SHarmonic2D fd1(Din(0,0)),fd2(Din(0,0));

        if (this_thread < f_start.size() - 1) 
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //      Vertical loop, boundaries between threads
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {

                l = dist_il[id];
                m = dist_im[id];

                // std::cout << "\n (l,m) = " << l << ", " << m << " \n";

                fd1 = Din(l,m);     fd1 = fd1.Dx(Input::List().dbydx_order);  

                if (l == m)         // Diagonal, no l - 1
                {
                    if (l < l0){    vtemp *= A1(m,m);                           Dh(m+1,m) += fd1.mpaxis(vtemp);     vtemp /= A1(m,m);}
                }
                else if (l == l0)   // Last l, no l + 1
                {                
                                    vtemp *= A2(l0,m);                          Dh(l0-1,m) += fd1.mpaxis(vtemp);    vtemp /= A2(l0,m);
                }
                else
                {
                                    vtemp *= A2(l,m);               fd2 = fd1;  Dh(l-1,m) += fd1.mpaxis(vtemp);
                                    vtemp *= A1(l,m)/A2(l  ,m);                 Dh(l+1,m) += fd2.mpaxis(vtemp);     vtemp /= A1(l,m);                    
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = nwsediag_il[id];
                m = nwsediag_im[id];

                fd1 = Din(l,m);     fd1 = fd1.Dy(Input::List().dbydy_order);;

                if (m == 0)         // Top or Left, no l - 1, m - 1
                {
                    if (l < l0) {   vtemp *= C1[l];             Dh(l+1,m+1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C1[l];}
                }
                else if (m == m0 || l == l0)   // Bottom or right, no l + 1, m + 1
                {
                                    vtemp *= C4(l,m);           Dh(l-1,m-1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C4(l,m);
                }
                else
                {       
                    fd2 = fd1;      vtemp *= C1[l];             Dh(l+1,m+1) += fd1.mpaxis(vtemp);   vtemp *= 1.0/C1[l]; 
                    if (m>1)    {   vtemp *= C4(l,m);           Dh(l-1,m-1) += fd2.mpaxis(vtemp);   vtemp *= 1.0/C4(l,m);}
                }
            }

            #pragma omp barrier
            /// -------------------------------------------------------------------------------------- ///
            /// Anti-Diagonal
            /// -------------------------------------------------------------------------------------- ///
            for (size_t id = f_end[this_thread]; id < f_start[this_thread+1]; ++id)
            {
                l = neswdiag_il[id];
                m = neswdiag_im[id];

                fd1 = Din(l,m);     fd1 = fd1.Dy(Input::List().dbydy_order);;

                if (m == 0)         // Left wall, no l + 1, m - 1
                {
                    if (l > 1)  {               vtemp *= C3[l];         Dh(l-1,m+1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C3[l];}
                }
                else if (m == m0)   // Right boundary, no l - 1, m + 1
                {
                    if (l < l0) {               vtemp *= C2(l,m);       Dh(l+1,m-1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C2(l,m);}
                }
                else
                {          
                    if (m > 1 && l < l0)        
                    {   
                        fd2 = fd1;              vtemp *= C2(l,m);       Dh(l+1,m-1) += fd2.mpaxis(vtemp);    vtemp *= 1.0/C2(l,m);                  
                    }
                    if (l - 1 != m && l != m){  vtemp *= C3[l];         Dh(l-1,m+1) += fd1.mpaxis(vtemp);    vtemp *= 1.0/C3[l];}
                }
            }
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------

//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection::operator()(const DistFunc1D& Din, DistFunc1D& Dh) 
{
//--------------------------------------------------------------
    size_t l0(Din.l0());
    size_t m0(Din.m0());

    #pragma omp parallel num_threads(Input::List().ompthreads)
    {
        size_t this_thread  = omp_get_thread_num();

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);

        valarray<complex<double> > vtemp(vr); 
        vtemp /= Din.mass();
        size_t l(0),m(0);

        SHarmonic1D fd1(vr.size(),Din(0,0).numx()),fd2(vr.size(),Din(0,0).numx());
        
        if (this_thread == 0)
        {
            fd1 = Din(0,0);                         fd1 = fd1.Dx(Input::List().dbydx_order);
            vtemp *= A1(0,0);                       Dh(1,0) += fd1.mpaxis(vtemp);
            vtemp /= A1(0,0);
        }

        // ----------------------------------------- //
        //              Do the chunks
        // ----------------------------------------- //        
        for (size_t id = f_start_thread; id < f_end_thread; ++id)
        {   
            l = dist_il[id];    m = dist_im[id];
            
            fd1 = Din(l,m);     fd1 = fd1.Dx(Input::List().dbydx_order);

            if (l == m)         // Diagonal, no l - 1
            {
                if (l < l0) {   vtemp *= A1(m,m);           Dh(m+1,m) += fd1.mpaxis(vtemp);     vtemp /= A1(m,m);}
            }
            else if (l == l0)   // Last l, no l + 1
            {
                                vtemp *= A2(l0,m);          Dh(l0-1,m) += fd1.mpaxis(vtemp);    vtemp /= A2(l0,m);
            }
            else
            {
                fd2 = fd1;      vtemp *= A2(l,m);           Dh(l-1,m) += fd1.mpaxis(vtemp);
                                vtemp *= A1(l,m)/A2(l,m);   Dh(l+1,m) += fd2.mpaxis(vtemp);     vtemp /= A1(l,m);             
            }
        }
    }


    // ----------------------------------------- //
    //          Boundaries between chunks
    // ----------------------------------------- //

    #pragma omp parallel for num_threads(f_start.size()-1)
    for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {
        valarray<complex<double> > vtemp(vr);
        vtemp /= Din.mass();        
        size_t l(0),m(0);

        SHarmonic1D fd1(vr.size(),Din(0,0).numx()),fd2(vr.size(),Din(0,0).numx());

        for (size_t id = f_end[threadboundaries]; id < f_start[threadboundaries+1]; ++id)
        {   
            l = dist_il[id];    m = dist_im[id];

            fd1 = Din(l,m);     fd1 = fd1.Dx(Input::List().dbydx_order);

            if (l == m)         // Diagonal, no l - 1
            {
                if (l < l0) {   vtemp *= A1(m,m);           Dh(m+1,m) += fd1.mpaxis(vtemp);     vtemp /= A1(m,m);}
            }
            else if (l == l0)   // Last l, no l + 1
            {
                                vtemp *= A2(l0,m);          Dh(l0-1,m) += fd1.mpaxis(vtemp);    vtemp /= A2(l0,m);
            }
            else
            {
                fd2 = fd1;      vtemp *= A2(l,m);           Dh(l-1,m) += fd1.mpaxis(vtemp);
                                vtemp *= A1(l,m)/A2(l,m);   Dh(l+1,m) += fd2.mpaxis(vtemp);     vtemp /= A1(l,m);             
            }
        }
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection::es1d(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    // valarray<complex<double> > vtemp(vr);

    size_t l0(Din.l0());

    #pragma omp parallel num_threads(Input::List().ompthreads)
    {   
        size_t this_thread  = omp_get_thread_num();
        // std::cout << "\n hi i'm " << this_thread << "\n";

        size_t f_start_thread(f_start[this_thread]);
        size_t f_end_thread(f_end[this_thread]);


        // std::cout << "\n els[ " << this_thread << "] = " << f_start_thread << "....\n";
        // std::cout << "\n ele[ " << this_thread << "] = " << f_end_thread << "....\n";
        //  Initialize work variables
        SHarmonic1D fd1(vr.size(),Din(0,0).numx()),fd2(vr.size(),Din(0,0).numx());
        valarray<complex<double> > vtemp(vr);
        vtemp /= Din.mass();
        

        //  -------------------------------------------------------- //
        //   First thread takes the boundary conditions (l = 0)
        //   Last thread takes the boundary condition (l = l0)
        //  Rest proceed to chunks
        //  -------------------------------------------------------- //
        if (this_thread == 0)
        {
            fd1 = Din(0,0);                         fd1 = fd1.Dx(Input::List().dbydx_order);
            vtemp *= A1(0,0);                       Dh(1,0) += fd1.mpaxis(vtemp);
            vtemp /= A1(0,0);

            f_start_thread = 1;
        }

        if (this_thread == Input::List().ompthreads - 1)    
        {    
            fd1 = Din(l0,0);                        fd1 = fd1.Dx(Input::List().dbydx_order);
            vtemp *= A2(l0,0);                      Dh(l0-1,0) += fd1.mpaxis(vtemp);
            vtemp /= A2(l0,0);

            f_end_thread -= 1;
        }

        //  -------------------------------------------------------- //
        //  Do the chunks
        //  Initialize vtemp so that it starts correctly
        //  -------------------------------------------------------- //
        vtemp *= A1(f_start_thread-1,0);

        for (size_t l = f_start_thread; l < f_end_thread; ++l)
        {
            fd1 = Din(l,0);  //std::cout << "\n \n before dx, l = " << l << " \n";          
            fd1 = fd1.Dx(Input::List().dbydx_order);  //std::cout << " \n after dx\n";

            vtemp *= A2(l,0)/A1(l-1,0);    fd2 = fd1;  Dh(l-1,0) += fd1.mpaxis(vtemp);
            vtemp *= A1(l,0)/A2(l  ,0);                Dh(l+1,0) += fd2.mpaxis(vtemp);
        }    
    }

    //  -------------------------------------------------------- //
    //  Do the boundaries between the chunks
    //  -------------------------------------------------------- //
    #pragma omp parallel for num_threads(f_start.size()-1)
    for (size_t threadboundaries = 0; threadboundaries < f_start.size()-1; ++threadboundaries)
    {
        SHarmonic1D fd1(vr.size(),Din(0,0).numx()),fd2(vr.size(),Din(0,0).numx());
        valarray<complex<double> > vtemp(vr);
        vtemp /= Din.mass();        
        
        vtemp *= A1(f_end[threadboundaries]-1,0);

        for (size_t l = f_end[threadboundaries]; l < f_start[threadboundaries+1]; ++l)
        {   
            fd1 = Din(l,0);  //std::cout << "\n \n before dx, l = " << l << " \n";          
            fd1 = fd1.Dx(Input::List().dbydx_order); //std::cout << " \n after dx\n";

            vtemp *= A2(l,0)/A1(l-1,0);    fd2 = fd1;  Dh(l-1,0) += fd1.mpaxis(vtemp);
            vtemp *= A1(l,0)/A2(l  ,0);                Dh(l+1,0) += fd2.mpaxis(vtemp);
        }
    }         

    // }
    
    
    // vt *= 1.0/A2(l0,0);
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection::f1only(const DistFunc1D& Din, DistFunc1D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vtemp(vr);
    vtemp /= Din.mass();    

    SHarmonic1D fd1(vr.size(),Din(0,0).numx()),fd2(vr.size(),Din(0,0).numx());

    fd1 = Din(0,0);
    fd1 = fd1.Dx(Input::List().dbydx_order);     vtemp *= A00;
    Dh(1,0) += (fd1.mpaxis(vtemp));

    fd1 = Din(1,0);
    fd1 = fd1.Dx(Input::List().dbydx_order);     vtemp *= A10/A00;
    Dh(0,0) += (fd1.mpaxis(vtemp));

}

//--------------------------------------------------------------
//   Advection in x
void Spatial_Advection::f1only(const DistFunc2D& Din, DistFunc2D& Dh) {
//--------------------------------------------------------------

    valarray<complex<double> > vtemp(vr);
    vtemp /= Din.mass();    

    SHarmonic2D fd1(Din(0,0));

    fd1 = Din(0,0);
    fd1 = fd1.Dx(Input::List().dbydx_order);     vtemp *= A00;
    Dh(1,0) += (fd1.mpaxis(vtemp));

    fd1 = Din(1,0);
    fd1 = fd1.Dx(Input::List().dbydx_order);     vtemp *= A10/A00;
    Dh(0,0) += (fd1.mpaxis(vtemp));

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //       m = 0, advection in y
    //  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fd1 = Din(0,0);                 fd1 = fd1.Dy(Input::List().dbydy_order);;
    vtemp *= C1[0]/A10;             Dh(1,1) += fd1.mpaxis(vtemp);

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //       m = 1, advection in y
    //  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fd1 = Din(1,1);        fd1 = fd1.Dy(Input::List().dbydy_order);;
    
    vtemp *= B2[1]/C1[0];  fd1 = fd1.mpaxis(vtemp);  Dh(0,0) += fd1.Re();

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
//  Update B with E term from Faraday's Law
//-------------------------------------------------------------- 
Faraday::Faraday( double xmin, double xmax, size_t Nx,
    double ymin, double ymax, size_t Ny)
//--------------------------------------------------------------
// Constructor
//--------------------------------------------------------------

{
    idx = complex<double>((-1.0)/ (2.0*(xmax-xmin)/double(Nx))); // -1/(2dx)
    idy = complex<double>((-1.0)/ (2.0*(ymax-ymin)/double(Ny))); // -1/(2dy)
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Faraday::operator()(EMF1D& EMFin, EMF1D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Faraday's Law 
//--------------------------------------------------------------
    Field1D tmpE(EMFin.Ez());
//      dBy/dt +=   dEz/dx       
    // tmpE                 = EMFin.Ez();
    tmpE                *= idx;
    EMFh.By()           += tmpE.Dx(Input::List().dbydx_order);
    
//      dBz/dt += - dEy/dx       
    tmpE                 = EMFin.Ey();
    tmpE                *= (-1.0) * idx;
    EMFh.Bz()           += tmpE.Dx(Input::List().dbydx_order);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Faraday::operator()(EMF2D& EMFin, EMF2D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Faraday's Law 
//--------------------------------------------------------------

    Field2D tmpE(EMFin.Ez());
//      dBx/dt += - dEz/dy  
    // tmpE          = EMFin.Ez(); 
    tmpE         *= (-1.0) * idy;
    EMFh.Bx()    += tmpE.Dy(Input::List().dbydy_order);;

//      dBy/dt +=   dEz/dx       
    tmpE                 = EMFin.Ez(); 
    tmpE                *= idx;
    EMFh.By()           += tmpE.Dx(Input::List().dbydx_order);
        // EMFh.By()(numx-1)    = 0.0;        

//      dBz/dt +=   dEx/dy       
    tmpE          = EMFin.Ex(); 
    tmpE         *= idy;
    EMFh.Bz()    += tmpE.Dy(Input::List().dbydy_order);;    

//      dBz/dt += - dEy/dx       
    tmpE                 = EMFin.Ey(); 
    tmpE                *= (-1.0) * idx;
    EMFh.Bz()           += tmpE.Dx(Input::List().dbydx_order);  

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Update E with B term from Faraday's Law
//--------------------------------------------------------------
Ampere::Ampere(double xmin, double xmax, size_t Nx,
    double ymin, double ymax, size_t Ny)
//--------------------------------------------------------------
// Constructor
//--------------------------------------------------------------
{
    idx = complex<double>((-1.0)/ (2.0*(xmax-xmin)/double(Nx))); // -1/(2dx)
    idy = complex<double>((-1.0)/ (2.0*(ymax-ymin)/double(Ny))); // -1/(2dx)
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Ampere::operator()(EMF1D& EMFin, EMF1D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Ampere's Law 
//--------------------------------------------------------------
    Field1D tmpB(EMFin.Bz());
//      dEy/dt +=  - dBz/dx       
    // tmpB                 = EMFin.Bz();
    tmpB                *= (-1.0) * idx;
    EMFh.Ey()           += tmpB.Dx(Input::List().dbydx_order);

//      dEz/dt += dBy/dx       
    tmpB                 = EMFin.By();
    tmpB                *=  idx;
    EMFh.Ez()           += tmpB.Dx(Input::List().dbydx_order);    

}
//--------------------------------------------------------------
void Ampere::operator()(EMF2D& EMFin, EMF2D& EMFh) {
//--------------------------------------------------------------
//  This is the core calculation for Ampere's Law 
//--------------------------------------------------------------

    Field2D tmpB(EMFin.Bz());
//      dEx/dt +=   dBz/dy       
    // tmpB                 = EMFin.Bz(); 
    tmpB                *= idy;
    EMFh.Ex()           += tmpB.Dy(Input::List().dbydy_order);;

//      dEy/dt +=  - dBz/dx       
    tmpB                 = EMFin.Bz(); 
    tmpB                *= (-1.0) * idx;
    EMFh.Ey()           += tmpB.Dx(Input::List().dbydx_order);

//      dEz/dt +=  - dBx/dy       
    tmpB                 = EMFin.Bx(); 
    tmpB                *= (-1.0) * idy;
    EMFh.Ez()           += tmpB.Dy(Input::List().dbydy_order);;    

//      dEz/dt += dBy/dx       
    tmpB                 = EMFin.By(); 
    tmpB                *=  idx;
    EMFh.Ez()           += tmpB.Dx(Input::List().dbydx_order);   
        // EMFh.Ez()(numx-1)    = 0.0;

}
//--------------------------------------------------------------
//**************************************************************
