/*! \brief Input reader - Declarations
 * \author PICKSC
 * \date   March 2017
 * \file   input.h
 * 
 * Contains:
 * 1) input reader
 * 2) default values for input variables.
 */

#ifndef DECL_IMPORT_H
#define DECL_IMPORT_H



//**************************************************************
//--------------------------------------------------------------
namespace Input{
//**************************************************************
//--------------------------------------------------------------
    class Input_List {
    public:
        Input_List();

//--------------------------------------------------------------
//  declaration of input list variables (input deck)
//--------------------------------------------------------------
        bool isthisarestart;

        size_t NnodesX;

        size_t numsp;

        size_t l0, m0, nump, Nx;
        double xmin, xmax;


        std::vector< double > pmax;
        double clf_dp;


//          Algorithms
        bool if_tridiagonal;
        bool implicit_E;
        bool implicit_B;
        bool collisions;
        int f00_implicitorexplicit;
        bool flm_collisions;

        int BoundaryCells;
        int bndX;
        size_t n_outsteps, n_distoutsteps;
        double t_stop;
        int restart_time;  int n_restarts;

//          Output
        bool o_EHist;
        bool o_Ex, o_Ey, o_Ez, o_Bx, o_By, o_Bz, o_x1x2, o_pth, o_p2p1x1, o_p1p2p3;
        bool o_G, o_Px, o_PxPx, o_Py, o_PxPy, o_PyPy, o_Pz, o_PxPz, o_PyPz, o_PzPz;
        bool o_Vx, o_VxVx, o_Vy, o_VxVy, o_VyVy, o_VxVz, o_VyVz, o_VzVz;
        bool o_Vsq, o_Qx, o_Qy, o_Qz;
        bool o_vNx, o_vNy, o_vNz;
        bool o_Jx, o_Jy, o_Jz;
        bool o_Pressure, o_Temperature, o_ND, o_Nu, o_p1x1, o_f0x1, o_f10x1, o_f11x1,
                o_f20x1, o_fl0x1;
        bool o_Ux, o_Uy, o_Uz, o_Z, o_ni, o_Ti;
        bool only_output;
        size_t numpx, nump1, nump2, nump3;

//          Electron-ion collisions
        double lnLambda, density_np;

//          Electron-electron collisions
        size_t RB_D_itmax;
        double RB_D_tolerance;
        double small_dt;
        double smaller_dt;
        int NB_algorithms;


//          Hydro parameters
        bool hydromotion;
        double hydromass, hydrocharge;

        int polarization_direction;

        bool init_f1, init_f2;

        bool MX_cooling;

        double super_gaussian_m;


        

        double pth_ref;

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
        /// Initialization

        std::vector<std::string> dens_profile_str;
        std::vector<std::string> temp_profile_str;
        std::vector<std::string> f10x_profile_str;
        std::vector<std::string> f20x_profile_str;
        std::vector<std::string> f_pedestal;

        std::string hydro_dens_profile_str;
        std::string hydro_temp_profile_str;
        std::string hydro_vel_profile_str;
        std::string hydro_Z_profile_str;


// ----------------------------------------------------------------------
        /// External fields
//
        bool ext_fields,trav_wave;
        int num_waves;
        bool IB_heating;
        double I_0, lambda_0;

        std::string intensity_profile_str;
        std::string intensity_time_profile_str;
        std::string ex_time_profile_str;
        std::string ey_time_profile_str;
        std::string ez_time_profile_str;
        std::string bx_time_profile_str;
        std::string by_time_profile_str;
        std::string bz_time_profile_str;

        std::string ex_profile_str;
        std::string ey_profile_str;
        std::string ez_profile_str;
        std::string bx_profile_str;
        std::string by_profile_str;
        std::string bz_profile_str;

        std::vector< std::string > ex_wave_profile_str;
        std::vector< std::string > ey_wave_profile_str;
        std::vector< std::string > ez_wave_profile_str;
        std::vector< std::string > bx_wave_profile_str;
        std::vector< std::string > by_wave_profile_str;
        std::vector< std::string > bz_wave_profile_str;
        std::vector< std::string > wave_time_envelope_str;
        std::vector< double > trav_wave_rise;
        std::vector< double > trav_wave_flat;
        std::vector< double > trav_wave_fall;
        std::vector< double > trav_wave_center;

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

        std::vector< std::string > oTags;
        std::vector<double> qs;
        std::vector<double> mass;
        std::vector<size_t> ls;
        std::vector<size_t> ms;
        std::vector<size_t> ps;
        std::vector<double> pth;

        std::vector<size_t> Npx;
        std::vector<size_t> Npy;
        std::vector<size_t> Npz;


        std::vector<size_t> NxGlobal;
        std::vector<size_t> NxLocalnobnd;
        std::vector<size_t> NxLocal;

        std::vector<double> xminGlobal;
        std::vector<double> xmaxGlobal;
        std::vector<double> xminLocal;
        std::vector<double> xmaxLocal;
        std::vector<double> globdx;

    };

//--------------------------------------------------------------
    Input_List& List();

}
//--------------------------------------------------------------
//**************************************************************

#endif
