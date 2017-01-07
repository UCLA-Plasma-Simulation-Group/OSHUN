///////////////////////////////////////////////////////////
//   Contributing authors : Michail Tzoufras, Benjamin Winjum
//
//  Last Modified:  September 1, 2016
///////////////////////////////////////////////////////////

//   
//   This header file contains the declaration for the 
//   structures that are required for importing the input
//   parameters:
//
//   1. Import:
//     Input data
//
///////////////////////////////////////////////////////////

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
            size_t NnodesX, NnodesY;

            size_t numsp;

            size_t l0, m0, nump, numx_glob, numy_glob, numx, numy;
            std::vector< double > xmin, xmax;


            double ymin, ymax;
            // double pmin;
            std::vector< double > pmax;
            double clf_dp;
            double small_dt;
            double smaller_dt;
            int NB_algorithms;

//          Algorithms
            bool if_tridiagonal;
            bool implicit_E;
            bool implicit_B;

            int RKLevel;
            int bndX, bndY;
            size_t n_outsteps; double t_stop;
            int restart_time;  int n_restarts;

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
            double sigma_p, center_p;
//            double sigma_x, center_x;
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//          Output
            bool o_EHist;
            bool o_Ex, o_Ey, o_Ez, o_Bx, o_By, o_Bz, o_x1x2, o_pth, o_p2p1x1, o_p1p2p3;
            bool o_G, o_Px, o_PxPx, o_Py, o_PxPy, o_PyPy, o_Pz, o_PxPz, o_PyPz, o_PzPz;
            bool o_Vx, o_VxVx, o_Vy, o_VxVy, o_VyVy, o_VxVz, o_VyVz, o_VzVz;
            bool o_Vsq, o_Qx, o_Qy, o_Qz;
            bool o_vNx;
            bool o_Jx, o_Jy, o_Jz;
            bool o_Pressure, o_Temperature, o_ND, o_Nu, o_p1x1, o_f0x1, o_f10x1, o_f11x1;
            bool o_Ux, o_ni, o_Ti;
            bool only_output;
            size_t numpx, nump1, nump2, nump3;

//          Electron-ion collisions
            double lnLambda, density_np;


//          Hydro parameters
            bool hydromotion;
            double hydromass, hydrocharge;

//          Phenomenological laser parameters
            int polarization_direction;


            int inverse_bremsstrahlung;
            double I_0, lambda_0;
            double rise_time, fall_time, flat_time;
            
            bool ext_fields;

//          density initialization

//          temperature initialization
            double pth_ref;

            // Strings

            std::vector<std::string> dens_profile_str;
            std::vector<std::string> temp_profile_str;

            std::string hydro_dens_profile_str;
            std::string hydro_temp_profile_str;
            std::string hydro_vel_profile_str;


            std::string intensity_profile_str;
            
            std::string ex_time_profile_str;
            std::string ey_time_profile_str;
            std::string ez_time_profile_str;
            std::string bx_time_profile_str;
            std::string by_time_profile_str;
            std::string bz_time_profile_str;

            std::string Ex_profile_str;
            std::string Ey_profile_str;
            std::string Ez_profile_str;
            std::string Bx_profile_str;
            std::string By_profile_str;
            std::string Bz_profile_str;

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
        std::vector< size_t > Nx; 

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
