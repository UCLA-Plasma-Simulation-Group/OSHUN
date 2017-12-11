/*! \brief Functors for various methods
 * \author PICKSC
 * \file   functors.cpp
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

//  Declerations
#include "state.h"
// #include "input.h"
#include "fluid.h"
#include "vlasov.h"
#include "functors.h"

//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with explicit
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_explicitE::VlasovFunctor1D_explicitE(vector<size_t> Nl,vector<size_t> Nm,
                                                    // vector<double> pmax, vector<size_t> Np,
                                                    vector<valarray<double> > dp,
                                                     double xmin, double xmax, size_t Nx) {
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, 0., 1., 1) );

        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );        

        JX.push_back( Current() );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

        AM.push_back( Ampere(xmin, xmax, Nx, 0., 1., 1) );

        FA.push_back( Faraday(xmin, xmax, Nx, 0., 1., 1) );

    }
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor1D_explicitE::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------
    bool debug(0);

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) {

            SA[s].f1only(Yin.DF(s),Yslope.DF(s));

            EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            AM[s](Yin.EMF(),Yslope.EMF());

            FA[s](Yin.EMF(),Yslope.EMF());

        }
        else if (Yin.DF(s).m0() == 0) {

            // GA[s].es1d(Yin.DF(s),Yslope.EMF().Ex());
            
            // if (debug) 
            // {
            //     std::cout << "\n\n f at start:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yin.SH(0,1,0)(ip,4);
            //     }

            //     std::cout << "\n\n E at start:";
            //     for (size_t ix(0); ix < Yin.SH(0,0,0).numx(); ++ix){
            //         std::cout << "\nEx(" << ix << ") = " << Yin.EMF().Ex()(ix);
            //     }            
            // }
            EF[s].es1d(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));


            // if (debug) 
            // {
            //     std::cout << "\n\nf after E:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yslope.SH(0,1,0)(ip,4);
            //     }
            // }
            JX[s].es1d(Yin.DF(s),Yslope.EMF().Ex());



            // if (debug) 
            // {
            //     std::cout << "\n\n after J:";
            //     for (size_t ix(0); ix < Yin.SH(0,0,0).numx(); ++ix){
            //         std::cout << "\nEx(" << ix << ") = " << Yslope.EMF().Ex()(ix);
            //     }            
            // }
            
            SA[s].es1d(Yin.DF(s),Yslope.DF(s));

            // if (debug) 
            // {
            //     std::cout << "\n\n after SA:";
            //     for (size_t ip(0); ip < Yin.SH(0,0,0).nump(); ++ip){
            //         std::cout << "\nf(" << ip << ") = " << Yslope.SH(0,1,0)(ip,4);
            //     }
            // }
        }

        else {

            SA[s](Yin.DF(s),Yslope.DF(s));

            EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            
        }

    }

}

void VlasovFunctor1D_explicitE::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}
//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_explicitE::VlasovFunctor2D_explicitE(vector<size_t> Nl,vector<size_t> Nm,
                                                    // vector<double> pmax, vector<size_t> Np,
                                                    vector<valarray<double> > dp,
                                                     double xmin, double xmax, size_t Nx,
                                                     double ymin, double ymax, size_t Ny) {
//--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, ymin, ymax, Ny) );

        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );        

        JX.push_back( Current() );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

        AM.push_back( Ampere(xmin, xmax, Nx, ymin, ymax, Ny) );

        FA.push_back( Faraday(xmin, xmax, Nx, ymin, ymax, Ny) );

    }
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Collect all of the terms
void VlasovFunctor2D_explicitE::operator()(const State2D& Yin, State2D& Yslope){
//--------------------------------------------------------------
    bool debug(0);

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) {

            SA[s].f1only(Yin.DF(s),Yslope.DF(s));

            EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            AM[s](Yin.EMF(),Yslope.EMF());

            FA[s](Yin.EMF(),Yslope.EMF());

        }
        
        else {

            SA[s](Yin.DF(s),Yslope.DF(s));

            EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));

            JX[s](Yin.DF(s),Yslope.EMF().Ex(),Yslope.EMF().Ey(),Yslope.EMF().Ez());

            
        }

    }

}

void VlasovFunctor2D_explicitE::operator()(const State2D& Yin, State2D& Yslope, size_t direction){}
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods with implicit
//  E-field solver

//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_implicitE_p1::VlasovFunctor1D_implicitE_p1(vector<size_t> Nl,vector<size_t> Nm,
                                                            // vector<double> pmax, vector<size_t> Np,
                                                    vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, 0., 1., 1) );

        // HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], dp[s], xmin, xmax, Nx) );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor1D_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) 
        {
            SA[s].f1only(Yin.DF(s),Yslope.DF(s));
            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }
        else {
            SA[s](Yin.DF(s),Yslope.DF(s));
            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }

        // if (Input::List().filterdistribution)  Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor1D_implicitE_p1::operator()(const State1D& Yin, State1D& Yslope, size_t direction){}

//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_implicitE_p1::VlasovFunctor2D_implicitE_p1(vector<size_t> Nl,vector<size_t> Nm,
                                            
                                                    vector<valarray<double> > dp,
                                                    double xmin, double xmax, size_t Nx,
                                                    double ymin, double ymax, size_t Ny) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        SA.push_back( Spatial_Advection(Nl[s], Nm[s], dp[s], xmin, xmax, Nx, ymin, ymax, Ny) );

        // HA.push_back( Hydro_Advection_1D(Nl[s], Nm[s], dp[s], xmin, xmax, Nx) );

        BF.push_back( Magnetic_Field(Nl[s], Nm[s], dp[s]) );

    }

}

//--------------------------------------------------------------
//
//
void VlasovFunctor2D_implicitE_p1::operator()(const State2D& Yin, State2D& Yslope){
//--------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {

        if (Yin.DF(s).l0() == 1) 
        {
            SA[s].f1only(Yin.DF(s),Yslope.DF(s));
            BF[s].f1only(Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }
        else {
            SA[s](Yin.DF(s),Yslope.DF(s));
            BF[s](Yin.DF(s),Yin.EMF().Bx(),Yin.EMF().By(),Yin.EMF().Bz(),Yslope.DF(s));
            // HA[s](Yin.DF(s),Yin.HYDRO(),Yslope.DF(s));
        }

        // if (Input::List().filterdistribution)  Yslope.DF(s) = Yslope.DF(s).Filterp();


    }

}
void VlasovFunctor2D_implicitE_p1::operator()(const State2D& Yin, State2D& Yslope, size_t direction){}


//--------------------------------------------------------------
//  Constructor
VlasovFunctor1D_implicitE_p2::VlasovFunctor1D_implicitE_p2(vector<size_t> Nl,vector<size_t> Nm,
                                                                vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        FA.push_back( Faraday(xmin, xmax, Nx, 0., 1., 1) );
        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );

    }

}
//------------------------------------------------------------------------------------------------------
void VlasovFunctor1D_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope){
// //---------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        FA[s](Yin.EMF(),Yslope.EMF());

        if (Yin.DF(s).l0() == 1) EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));
        else                     EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();

    }

}

//--------------------------------------------------------------------------------------------------
void VlasovFunctor1D_implicitE_p2::operator()(const State1D& Yin, State1D& Yslope, size_t direction){
//--------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    if (direction == 1)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ex_f1only(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            else                     EF[s].Implicit_Ex(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else if (direction == 2)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ey_f1only(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            else                     EF[s].Implicit_Ey(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else
    {
        // complex<double> ii(0.0,1.0);
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ez_f1only(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            else                     EF[s].Implicit_Ez(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }

}

//--------------------------------------------------------------
//  Constructor
VlasovFunctor2D_implicitE_p2::VlasovFunctor2D_implicitE_p2(vector<size_t> Nl,vector<size_t> Nm,
                                                                vector<valarray<double> > dp,
                                                           double xmin, double xmax, size_t Nx,
                                                           double ymin, double ymax, size_t Ny) {
// //--------------------------------------------------------------

    for (size_t s(0); s < Nl.size(); ++s){

        FA.push_back( Faraday(xmin, xmax, Nx, ymin, ymax, Ny) );
        EF.push_back( Electric_Field(Nl[s], Nm[s], dp[s]) );

    }

}
//------------------------------------------------------------------------------------------------------
void VlasovFunctor2D_implicitE_p2::operator()(const State2D& Yin, State2D& Yslope){
// //---------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    for (size_t s(0); s < Yin.Species(); ++s) {
        FA[s](Yin.EMF(),Yslope.EMF());

        if (Yin.DF(s).l0() == 1) EF[s].f1only(Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));
        else                     EF[s](Yin.DF(s),Yin.EMF().Ex(),Yin.EMF().Ey(),Yin.EMF().Ez(),Yslope.DF(s));

        // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();

    }

}

//--------------------------------------------------------------------------------------------------
void VlasovFunctor2D_implicitE_p2::operator()(const State2D& Yin, State2D& Yslope, size_t direction){
//--------------------------------------------------------------------------------------------------

    Yslope = 0.0;

    if (direction == 1)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ex_f1only(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            else                     EF[s].Implicit_Ex(Yin.DF(s),Yin.EMF().Ex(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else if (direction == 2)
    {
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ey_f1only(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            else                     EF[s].Implicit_Ey(Yin.DF(s),Yin.EMF().Ey(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }
    else
    {
        // complex<double> ii(0.0,1.0);
        for (size_t s(0); s < Yin.Species(); ++s) {
            if (Yin.DF(s).l0() == 1) EF[s].Implicit_Ez_f1only(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            else                     EF[s].Implicit_Ez(Yin.DF(s),Yin.EMF().Ez(),Yslope.DF(s));
            // if (Input::List().filterdistribution) Yslope.DF(s) = Yslope.DF(s).Filterp();
        }
    }

}



