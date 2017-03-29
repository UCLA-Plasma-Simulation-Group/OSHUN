/*! \brief Input reader - Definitions
 * \author PICKSC
 * \date   March 2017
 * \file   input.cpp
 * 
 * Contains:
 * 1) input reader
 * 2) default values for input variables.
 */
//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <cfloat>
#include <fstream>
#include <string>

#include <math.h>
#include <map>

//  Declarations
#include "input.h"

//**************************************************************
//--------------------------------------------------------------
Input::Input_List::Input_List(): 
    isthisarestart(0),
    NnodesX(4),
    numsp(1),
    l0(6),
    m0(4), 
    nump(96), 
    Nx(32),
    xmin(-1000.0), 
    xmax(1000.0),
    clf_dp(1.0),
    if_tridiagonal(1),
    implicit_E(1),
    implicit_B(0),
    collisions(1),
    f00_implicitorexplicit(2),
    flm_collisions(0),
    BoundaryCells(3),
    bndX(0),
    n_outsteps(1000), 
    n_distoutsteps(100),
    t_stop(8000),
    restart_time(10000.0),
    n_restarts(100),

//          Output
    o_EHist(0),
    o_Ex(0), o_Ey(0), o_Ez(0), o_Bx(0), o_By(0), o_Bz(0), o_x1x2(0), o_pth(0), o_p2p1x1(0), o_p1p2p3(0),
    o_G(0), o_Px(0), o_PxPx(0), o_Py(0), o_PxPy(0), o_PyPy(0), o_Pz(0), o_PxPz(0), o_PyPz(0), o_PzPz(0),
    o_Vx(0), o_VxVx(0), o_Vy(0), o_VxVy(0), o_VyVy(0), o_VxVz(0), o_VyVz(0), o_VzVz(0),
    o_Vsq(0), o_Qx(0), o_Qy(0), o_Qz(0),
    o_vNx(0), o_vNy(0), o_vNz(0),
    o_Jx(0), o_Jy(0), o_Jz(0),
    o_Pressure(0), o_Temperature(0), o_ND(0), o_Nu(0), o_p1x1(0), o_f0x1(0), o_f10x1(0), o_f11x1(0), o_f20x1(0), o_fl0x1(0),
    o_Ux(0), o_Uy(0), o_Uz(0), o_Z(0), o_ni(0), o_Ti(0),
    
    numpx(96), nump1(96), nump2(96), nump3(96),

//          Electron-ion collisions
    density_np(1.0e21),

//          Electron-electron collisions
    RB_D_itmax(100),
    RB_D_tolerance(1e-12),
    small_dt(1e-1),
    smaller_dt(1e-5),
    NB_algorithms(4),


//          Hydro parameters
    hydromotion(0),
    hydromass(100), hydrocharge(79),
    polarization_direction(0),
    init_f1(0), init_f2(0),
    MX_cooling(0),
    super_gaussian_m(2.0),

    pth_ref(0.025),
        

    hydro_dens_profile_str("cst{0.0}"),
    hydro_temp_profile_str("cst{0.0}"),
    hydro_vel_profile_str("cst{0.0}"),
    hydro_Z_profile_str("cst{0.0}"),


// ----------------------------------------------------------------------
        /// External fields
//
    ext_fields(0),trav_wave(0),num_waves(0),
    IB_heating(0),
    I_0(0.0), lambda_0(0.351),

    intensity_profile_str("cst{0.0}"),
    intensity_time_profile_str("cst{0.0}"),
    ex_time_profile_str("cst{0.0}"),
    ey_time_profile_str("cst{0.0}"),
    ez_time_profile_str("cst{0.0}"),
    bx_time_profile_str("cst{0.0}"),
    by_time_profile_str("cst{0.0}"),
    bz_time_profile_str("cst{0.0}"),
    ex_profile_str("cst{0.0}"),
    ey_profile_str("cst{0.0}"),
    ez_profile_str("cst{0.0}"),
    bx_profile_str("cst{0.0}"),
    by_profile_str("cst{0.0}"),
    bz_profile_str("cst{0.0}")
    // ex_wave_profile_str("cst{0.0}"),
    // ey_wave_profile_str("cst{0.0}"),
    // ez_wave_profile_str("cst{0.0}"),
    // bx_wave_profile_str("cst{0.0}"),
    // by_wave_profile_str("cst{0.0}"),
    // bz_wave_profile_str("cst{0.0}"),
    // wave_time_envelope_str("cst{0.0}")

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
        // std::vector< double > pmax;
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
        /// Initialization

        // std::vector<std::string> dens_profile_str;
        // std::vector<std::string> temp_profile_str;
        // std::vector<std::string> f10x_profile_str;
        // std::vector<std::string> f20x_profile_str;
        // std::vector<std::string> f_pedestal;

        // std::vector<double> qs;
        // std::vector<double> mass;
        // std::vector<size_t> ls;
        // std::vector<size_t> ms;
        // std::vector<size_t> ps;
        // std::vector<double> pth;

        // std::vector<size_t> Npx;
        // std::vector<size_t> Npy;
        // std::vector<size_t> Npz;


        // std::vector<size_t> NxGlobal;
        // std::vector<size_t> NxLocalnobnd;
        // std::vector<size_t> NxLocal;

        // std::vector<double> xminGlobal,
        // std::vector<double> xmaxGlobal,
        // std::vector<double> xminLocal,
        // std::vector<double> xmaxLocal,
        // std::vector<double> globdx
         {
//--------------------------------------------------------------
//  The constructor for the input_list structure
//--------------------------------------------------------------

    std::ifstream deckfile("inputdeck");
    std::string deckstring, deckequalssign, deckstringbool;
    double deckreal;
    size_t tempint;


    if (deckfile.is_open()) {

        while (deckfile >> deckstring) {

            if (deckstring == "if_restart") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                isthisarestart = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }

            if (deckstring == "n_restarts") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> n_restarts;
            }

            if (deckstring == "restart_time") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1); 
                }
                deckfile >> restart_time;
            }

            if (deckstring == "NnodesX") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> NnodesX;
            }

            if (deckstring == "numsp") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> numsp;
            }
            if (deckstring == "l0") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> tempint;
                    ls.push_back(tempint);
                    // std::cout<< "l0 = " << tempint << "\n";
                }
            }
            if (deckstring == "m0") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> tempint;
                    ms.push_back(tempint);
                    // std::cout<< "m0 = " << tempint << "\n";
                }
            }
            if (deckstring == "nump") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> tempint;
                    ps.push_back(tempint);
                    Npx.push_back(2*tempint+1);
                    Npy.push_back(2*tempint+1);
                    Npz.push_back(2*tempint+1);
                    // std::cout<< "nump = " << tempint << "\n";
                }
            }
            if (deckstring == "pmax") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                // deckfile >> deckreal;
                // TWICE HERE FOR TWO SPECIES -- SHOULD BE CHANGED
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckreal;
                    pmax.push_back(deckreal);
                    // std::cout<< "pmax = " << deckreal << "\n";
                }
            }
            if (deckstring == "charge") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckreal;
                    qs.push_back(-1.0*deckreal);
                    // std::cout<< "q = " << deckreal << "\n";
                }
            }
            if (deckstring == "mass") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckreal;
                    mass.push_back(deckreal);
                    // std::cout<< "mass = " << deckreal << "\n";
                }
            }
            if (deckstring == "Nx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> Nx;
                NxGlobal.push_back(Nx);
            }

            if (deckstring == "xmin") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> xmin;
                xminGlobal.push_back(xmin);
            }
            if (deckstring == "xmax") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> xmax;
                xmaxGlobal.push_back(xmax);
            }
            if (deckstring == "clf_dp") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> clf_dp;
            }
            if (deckstring == "hydro") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                hydromotion = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');

            }
            if (deckstring == "hydroatomicmass") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> hydromass;
                hydromass *= 1836;
            }
            if (deckstring == "hydrocharge") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> hydrocharge;
            }
            if (deckstring == "ext_fields") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                ext_fields = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');

            }
            if (deckstring == "traveling_wave") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                trav_wave = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "num_waves") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> num_waves;
            }
            if (deckstring == "small_dt") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> small_dt;
            }
            if (deckstring == "smaller_dt") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> smaller_dt;
            }
            if (deckstring == "Rosenbluth_D_tolerance") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> RB_D_tolerance;
            }
            if (deckstring == "Rosenbluth_D_maximum_iterations") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> RB_D_itmax;
            }
            if (deckstring == "NB_algorithms") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> NB_algorithms;
            }
            if (deckstring == "assume_tridiagonal_flm_collisions") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                if_tridiagonal = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "implicit_E") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                implicit_E = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "implicit_B_push") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                implicit_B = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "collisions_switch") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                collisions = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "f00_collisions") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                if (deckstringbool[0] == 'i' || deckstringbool[0] == 'I'){
                    f00_implicitorexplicit = 2;
                }
                else if (deckstringbool[0] == 'e' || deckstringbool[0] == 'E'){
                    f00_implicitorexplicit = 1;
                }
                if (deckstringbool[0] == 'f' || deckstringbool[0] == 'F'){
                    f00_implicitorexplicit = 0;
                }
            }
            if (deckstring == "flm_collisions") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                flm_collisions = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "MX_cooling") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                MX_cooling = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }



            if (deckstring == "t_stop") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> t_stop;
            }
            if (deckstring == "n_outsteps") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> n_outsteps;
                // if (fmod(t_stop,n_outsteps) != 0){
                //     std::cout << "\n Need n_outsteps to be a n integer divisor of t_stop i.e. remainder (t_stop / n_outsteps) = 0 \n";
                //     exit(1);
                // }
            }

            if (deckstring == "n_distoutsteps") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> n_distoutsteps;
                // if (fmod(t_stop,n_distoutsteps) != 0){
                //     std::cout << "\n Need n_distoutsteps to be a n integer divisor of t_stop i.e. remainder (t_stop / n_distoutsteps) = 0 \n";
                //     exit(1);
                // }
            }

            if (deckstring == "bndX") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> bndX;
            }
            if (deckstring == "o_EHist") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_EHist = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Ex") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Ex = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Ey") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Ey = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Ez") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Ez = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Bx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Bx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_By") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_By = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Bz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Bz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_x1x2") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_x1x2 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_pth") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_pth = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_G") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_G = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Px") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Px = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PxPx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PxPx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Py") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Py = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PxPy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PxPy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PyPy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PyPy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Pz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Pz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PxPz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PxPz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PyPz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PyPz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_PzPz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_PzPz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Jx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Jx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Jy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Jy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Jz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Jz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_vNx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_vNx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_vNy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_vNy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_vNz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_vNz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Vx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Vx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VxVx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VxVx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Vy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Vy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VxVy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VxVy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VyVy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VyVy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VxVz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VxVz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VyVz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VyVz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_VzVz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_VzVz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Vsq") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Vsq = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Qx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Qx = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Qy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Qy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Qz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Qz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Temperature") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Temperature = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Pressure") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Pressure = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_ND") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_ND = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Nu") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Nu = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_p1x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_p1x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_f0x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_f0x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_f10x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_f10x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_f11x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_f11x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_f20x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_f20x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_fl0x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_fl0x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_p2p1x1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_p2p1x1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_p1p2p3") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_p1p2p3 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Ux") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Ux = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Uy") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Uy = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Uz") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Uz = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Z") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Z = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_ni") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_ni = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "o_Ti") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                o_Ti = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "nump1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> nump1;
                // Npx.push_back(nump1);

            }
            if (deckstring == "nump2") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> nump2;
                // Npx.push_back(nump2);
                // Npx.push_back(2*ps[0]);
            }
            if (deckstring == "nump3") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> nump3;
            }
            if (deckstring == "numpx") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> numpx;
            }
            if (deckstring == "only_output") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                only_output = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "lnLambda") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> lnLambda;
            }
            // if (deckstring == "Zeta") {
            //     deckfile >> deckequalssign;
            //     if(deckequalssign != "=") {
            //         std::cout << "Error reading " << deckstring << std::endl;
            //         exit(1);
            //     }
            //     deckfile >> Zeta;
            // }
            if (deckstring == "density_np") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> density_np;
            }
            if (deckstring == "super_gaussian_distribution") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> super_gaussian_m;
            }
            if (deckstring == "init_f1") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                init_f1 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "init_f2") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                init_f2 = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "polarization_direction") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> polarization_direction;
            }
            if (deckstring == "inverse_bremsstrahlung") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstringbool;
                IB_heating = (deckstringbool[0] == 't' || deckstringbool[0] == 'T');
            }
            if (deckstring == "I_0") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> I_0;
            }
            if (deckstring == "lambda_0") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> lambda_0;
            }
            if (deckstring == "pth_ref") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> pth_ref;
            }
            if (deckstring == "n(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckstring;
                    dens_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "T(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckstring;
                    temp_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "f_pedestal") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckstring;
                    f_pedestal.push_back(deckstring);
                }
            }
            if (deckstring == "multiplier-f10(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckstring;
                    f10x_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "multiplier-f20(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t s(0);s<numsp;++s)
                {
                    deckfile >> deckstring;
                    f20x_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "ni(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                hydro_dens_profile_str = deckstring;
            }
            if (deckstring == "Ti(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                hydro_temp_profile_str = deckstring;
            }
            if (deckstring == "Ux(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                hydro_vel_profile_str = deckstring;
            }
            if (deckstring == "Z(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                hydro_Z_profile_str = deckstring;
            }
            if (deckstring == "I(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                intensity_profile_str = deckstring;
            }
            if (deckstring == "I(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                intensity_time_profile_str = deckstring;
            }
            if (deckstring == "Ex(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ex_profile_str = deckstring;
            }
            if (deckstring == "Ey(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ey_profile_str = deckstring;
            }
            if (deckstring == "Ez(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ez_profile_str = deckstring;
            }
            if (deckstring == "Bx(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                bx_profile_str = deckstring;
            }
            if (deckstring == "By(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                by_profile_str = deckstring;
            }
            if (deckstring == "Bz(x)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                bz_profile_str = deckstring;
            }
            if (deckstring == "Ex(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ex_time_profile_str = deckstring;
            }
            if (deckstring == "Ey(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ey_time_profile_str = deckstring;
            }
            if (deckstring == "Ez(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                ez_time_profile_str = deckstring;
            }
            if (deckstring == "Bx(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                bx_time_profile_str = deckstring;
            }
            if (deckstring == "By(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                by_time_profile_str = deckstring;
            }
            if (deckstring == "Bz(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                deckfile >> deckstring;
                bz_time_profile_str = deckstring;
            }
            if (deckstring == "Ex(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    ex_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "Ey(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    ey_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "Ez(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    ez_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "Bx(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {    
                    deckfile >> deckstring;
                    bx_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "By(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    by_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "Bz(x,t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    bz_wave_profile_str.push_back(deckstring);
                }
            }
            if (deckstring == "envelope(t)") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckstring;
                    wave_time_envelope_str.push_back(deckstring);
                }
            }

            if (deckstring == "rise_flat_fall_center") {
                deckfile >> deckequalssign;
                if(deckequalssign != "=") {
                    std::cout << "Error reading " << deckstring << std::endl;
                    exit(1);
                }
                for (size_t n(0);n<num_waves;++n)
                {
                    deckfile >> deckreal;
                    trav_wave_rise.push_back(deckreal);
                    

                    deckfile >> deckreal;
                    trav_wave_flat.push_back(deckreal);
                    

                    deckfile >> deckreal;
                    trav_wave_fall.push_back(deckreal);
                    
                    deckfile >> deckreal;
                    trav_wave_center.push_back(deckreal);                    
                }
            }



        }

        deckfile.close();

        if (((NnodesX%2)!=0) && (NnodesX !=1))
            std::cout << "The number of nodes " <<NnodesX<< " is not even" << std::endl;

//            numx      = numx_glob / NnodesX;
//            numx     += 2*BoundaryCells;
//            numx_glob = NnodesX * (numx-2*BoundaryCells);

        //  INPUT PARAMETERS

        oTags.push_back("Time");
        oTags.push_back("Space");
        oTags.push_back("px-x");
        // oTags.push_back("px-x_1");
        oTags.push_back("f0-x");
        oTags.push_back("f10-x");
        oTags.push_back("f11-x");
        oTags.push_back("f20-x");
        oTags.push_back("fl0-x");
        oTags.push_back("pxpy-x");
        oTags.push_back("Ex");
        oTags.push_back("Ey");
        oTags.push_back("Ez");
        oTags.push_back("Bx");
        oTags.push_back("By");
        oTags.push_back("Bz");
        oTags.push_back("n");
        oTags.push_back("T_eV");
        oTags.push_back("T");
        oTags.push_back("Jx");
        oTags.push_back("Jy");
        oTags.push_back("Jz");
        oTags.push_back("Qx");
        oTags.push_back("Qy");
        oTags.push_back("Qz");
        oTags.push_back("vNx");
        oTags.push_back("vNy");
        oTags.push_back("vNz");
        oTags.push_back("Ux");
        oTags.push_back("Uy");
        oTags.push_back("Uz");
        oTags.push_back("Z");
        oTags.push_back("ni");
        oTags.push_back("Ti");



        // Determination of the local computational domain (i.e. the x-axis and the y-axis)

        BoundaryCells = 4;

        for (size_t i(0); i < xminGlobal.size(); ++i){
            NxLocalnobnd.push_back(NxGlobal[i]/ NnodesX) ;
            NxLocal.push_back(NxLocalnobnd[i] + 2 * BoundaryCells);
            xminLocal.push_back(0.0);
            xmaxLocal.push_back(0.0);
            globdx.push_back((xmaxGlobal[i]-xminGlobal[i])/(double (NxGlobal[i]) ));
        }
    } 
    else {
        std::cout << "Unable to open inputdeck" << std::endl;
        exit(1);
    }

}
//--------------------------------------------------------------

//--------------------------------------------------------------

//--------------------------------------------------------------

Input::Input_List& Input::List(){
    static Input::Input_List incoming;
    return incoming;
}

//--------------------------------------------------------------

//**************************************************************

