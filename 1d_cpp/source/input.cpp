///////////////////////////////////////////////////////////
//   Contributing authors : Michail Tzoufras, Benjamin Winjum
//
//  Last Modified:  September 1 2016
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the input
//
///////////////////////////////////////////////////////////
//
// 
//   class Export_Formatted_Data::
//
//   This class receives the output matrices and saves the 
//   data in txt files with recognizable name after it attaches
//   a small header with information necessary for plotting. 
//
// 
//   class Restart_Facility::
//
//   This class writes restart files from each node, and 
//   reads restart files for each node.  
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

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
    Input::Input_List::Input_List() {
//--------------------------------------------------------------
//  The constructor for the input_list structure
//--------------------------------------------------------------

        std::ifstream deckfile("inputdeck");
        std::string deckstring, deckequalssign, deckstringbool;
        double deckreal;
        size_t tempint;


        if (deckfile.is_open()) {

            while (deckfile >> deckstring) {

                if (deckstring == "NnodesX") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> NnodesX;
                }
                if (deckstring == "NnodesY") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> NnodesY;
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
                        qs.push_back(deckreal);
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
                if (deckstring == "numx_glob") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> numx_glob;
                    Nx.push_back(numx_glob);
                }
                if (deckstring == "numy_glob") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> numy_glob;
                }
                numy_glob = 1;
                if (deckstring == "xmin") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckreal;
                    xmin.push_back(deckreal);
                }
                if (deckstring == "xmax") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckreal;
                    xmax.push_back(deckreal);
                }
                if (deckstring == "ymin") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> ymin;
                }
                if (deckstring == "ymax") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> ymax;
                }

     //            if (deckstring == "pmax") {
     //                deckfile >> deckequalssign;
     //                if(deckequalssign != "=") {
     //                    std::cout << "Error reading " << deckstring << std::endl;
     //                    exit(1);
     //                }
     //                deckfile >> deckreal;
     //                // TWICE HERE FOR TWO SPECIES -- SHOULD BE CHANGED
                    // pmax.push_back(deckreal); 
                    // // pmax.push_back(deckreal);
     //            }

                // if (deckstring == "pmin") {
                //     deckfile >> deckequalssign;
                //     if(deckequalssign != "=") {
                //         std::cout << "Error reading " << deckstring << std::endl;
                //         exit(1);
                //     }
                //     deckfile >> pmin;
                // }
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
                    // if (hydromotion){
                    //     if (ls( < 4){
                    //         std::cout << "\n ERROR: Need l0 > 4 for hydrosolver. l0 = " << l0 << " \n";
                    //         break;
                    //     }
                    // }
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
                    // if (hydromotion){
                    //     if (ls( < 4){
                    //         std::cout << "\n ERROR: Need l0 > 4 for hydrosolver. l0 = " << l0 << " \n";
                    //         break;
                    //     }
                    // }
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
                if (deckstring == "NB_algorithms") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> NB_algorithms;
                }
                if (deckstring == "if_tridiagonal") {
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
                // if (deckstring == "RKLevel") {
                //     deckfile >> deckequalssign;
                //     if(deckequalssign != "=") {
                //         std::cout << "Error reading " << deckstring << std::endl;
                //         exit(1);
                //     }
                //     deckfile >> RKLevel;
                // }
                RKLevel = 3;
                if (deckstring == "n_outsteps") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> n_outsteps;
                }
                if (deckstring == "t_stop") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> t_stop;
                }
                if (deckstring == "restart_time") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> restart_time;
                }
                if (deckstring == "n_restarts") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> n_restarts;
                }
                if (deckstring == "bndX") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> bndX;
                }
                if (deckstring == "bndY") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> bndY;
                }
                if (deckstring == "sigma_p") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> sigma_p;
                }
                if (deckstring == "center_p") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> center_p;
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
                    deckfile >> inverse_bremsstrahlung;
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
                if (deckstring == "U(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    hydro_vel_profile_str = deckstring;
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
                if (deckstring == "Ex(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    Ex_profile_str = deckstring;
                }
                if (deckstring == "Ey(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    Ey_profile_str = deckstring;
                }
                if (deckstring == "Ez(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    Ez_profile_str = deckstring;
                }
                if (deckstring == "Bx(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    Bx_profile_str = deckstring;
                }
                if (deckstring == "By(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    By_profile_str = deckstring;
                }
                if (deckstring == "Bz(x)") {
                    deckfile >> deckequalssign;
                    if(deckequalssign != "=") {
                        std::cout << "Error reading " << deckstring << std::endl;
                        exit(1);
                    }
                    deckfile >> deckstring;
                    Bz_profile_str = deckstring;
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



            }

            deckfile.close();

            if (((NnodesX%2)!=0) && (NnodesX !=1)) 
                std::cout << "The number of nodes " <<NnodesX<< " is not even" << std::endl;

            numx      = numx_glob / NnodesX;
            numx     += 2*RKLevel;        
            numx_glob = NnodesX * (numx-2*RKLevel);

            if (((NnodesY%2)!=0) && (NnodesY !=1)) 
                std::cout << "The number of nodes " <<NnodesY<< " is not even" << std::endl;
            numy      = numy_glob / NnodesY;
            numy     += 2*RKLevel;        
            numy_glob = NnodesY * (numy-2*RKLevel);

            //  INPUT PARAMETERS

                oTags.push_back("Time");
                oTags.push_back("Space");
                oTags.push_back("px-x");
                // oTags.push_back("px-x_1");
                oTags.push_back("f0-x");
                oTags.push_back("f10-x");
                oTags.push_back("f11-x");
                oTags.push_back("pxpy-x");
                oTags.push_back("Ex");
                oTags.push_back("Ey");
                oTags.push_back("Ez");
                oTags.push_back("Bx");
                oTags.push_back("By");
                oTags.push_back("Bz");
                oTags.push_back("n");
                oTags.push_back("T_eV");
                oTags.push_back("Jx");
                oTags.push_back("Jy");
                oTags.push_back("Jz");
                oTags.push_back("Qx");
                oTags.push_back("vNx");
                oTags.push_back("Ux");
                oTags.push_back("ni");
                oTags.push_back("Ti");

            

                // Determination of the local computational domain (i.e. the x-axis and the y-axis) 

                for(size_t i(0); i < xmin.size(); ++i) {
                    NxGlobal.push_back( Nx[0] );
                    NxLocalnobnd.push_back( NxGlobal[i] / NnodesX );
                    NxLocal.push_back( NxLocalnobnd[i]  + 2 * RKLevel );

                    xminGlobal.push_back( xmin[0] );
                    xmaxGlobal.push_back( xmax[0] );
                    xminLocal.push_back( 0.0 );
                    xmaxLocal.push_back( 0.0 );
                    globdx.push_back( (xmaxGlobal[i]-xminGlobal[i])/(static_cast<double>(NxGlobal[i])) );
                }


        } else {
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

