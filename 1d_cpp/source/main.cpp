/*! \brief Main Loop
 * \author PICKSC
 * \date   September 1, 2016
 * \file   main.cpp
 * 
 * This file includes
 * - the main loop.
 * - the definition for the clock class
 * - the definition of periodic boundaries
 *
 */  

// Standard libraries 
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdarg.h>
 

#include <map>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <ctime>


// My libraries 
#include "lib-array.h"
#include "lib-algorithms.h"
#include "H5Cpp.h"
#include "exprtk.hpp"

// Misc Declerations
#include "input.h"
#include "state.h"
#include "formulary.h"
#include "setup.h"
#include "fluid.h"
#include "vlasov.h"
#include "collisions.h"
#include "functors.h"
#include "parallel.h"
#include "implicitE.h"

#include "export.h"



//--------------------------------------------------------------
//**************************************************************

//----------------------------------------------------------------------------------------------------------------------------
    double startmessages(State1D& Y){

        Formulary formulas;
        // double Tref = Input::List().pth_ref*Input::List().pth_ref;
        double ND = 1.72e9*sqrt(pow(formulas.T0,3.0)/formulas.n);
        double nuei_wp = sqrt(2.0/Formulary::pi)/9.0/ND*
        formulas.LOGei(1.0,pow(Input::List().pth_ref,2.0),formulas.Zeta)*formulas.Zeta;

        double nu_ei = nuei_wp * formulas.wp;
        double dt_out(Input::List().t_stop / Input::List().n_outsteps);
        double CLF(Input::List().clf_dp);
        size_t n_outsteps(Input::List().n_outsteps);
        double numh  = size_t(static_cast<int>(dt_out/CLF))+1; 
        double h     = dt_out/static_cast<double>(numh);


        std::cout << "\n\n"; 
        std::cout<< "------- ---------------- */*/*/*/*/* ---------------- -------\n";                          
        std::cout<< "------- ,-----.  ,---.  ,--.  ,--.,--. ,--.,--.  ,--. -------\n";
        std::cout<< "-------'  .-.  ''   .-' |  '--'  ||  | |  ||  ,'.|  | -------\n";
        std::cout<< "-------|  | |  |`.  `-. |  .--.  ||  | |  ||  |' '  | -------\n";
        std::cout<< "-------'  '-'  '.-'    ||  |  |  |'  '-'  '|  | `   | -------\n";
        std::cout<< "------- `-----' `-----' `--'  `--' `-----' `--'  `--' -------\n";
        std::cout<< "------- ---------------- */*/*/*/*/* ---------------- -------\n";                          

        

        std::cout << "\n\n";
        std::cout << "---------------- OSHUN Beta - 1D 3P ---------------- \n";
        std::cout << "    Particle-in-Cell and Kinetic Simulation Center   \n";
        std::cout << "------------------- UCLA - 2017 -------------------- \n";
        
        std::cout << "----------------- Reference Units ------------------ \n";
        std::cout << "normalizing density                   = " << formulas.n << " / cc \n";
        std::cout << "plasma frequency                      = " << formulas.wp << " Hz \n";
        std::cout << "background ions - Z                   = " << formulas.Zeta << " \n";
        std::cout << "normalizing magnetic field            = " << formulas.B0 << " T \n";
        std::cout << "Reference Temperature                 = " << formulas.T0 << " eV \n";
        std::cout << "Corresponding thermal velocity        = " << Input::List().pth_ref*Formulary::cL  << " m/s \n";

        std::cout << "\n";
        std::cout << "Corresponding e-i log Lambda          = " << formulas.LOGei(1.0,pow(Input::List().pth_ref,2.0),formulas.Zeta) << " \n";
        std::cout << "Corresponding e-i collision frequency = " << nu_ei << " Hz \n";
        std::cout << "Corresponding tau_e-i / tau_p         = " << 1.0/nuei_wp << " \n";
        std::cout << "Corresponding e-i mean free path      = " << Input::List().pth_ref*Formulary::cL/nu_ei *1e6<< " microns \n";

        std::cout << "\n";
        std::cout << "Corresponding e-e log Lambda          = " << formulas.LOGee(1.0,pow(Input::List().pth_ref,2.0)) << " \n";
        std::cout << "Corresponding e-e collision frequency = " << 1.0/formulas.Tau_e(1.0,pow(Input::List().pth_ref,2.0)) << " Hz \n";
            // std::cout << "Corresponding thermal mean free path =" << Input::List().pth_ref/cL / nu_ei  << "\n";
            // std::cout << "Corresponding thermal velocity  c =" << pow(Input::List().pth_ref*cL,2.0)*0.5*me << "\n";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "--------------  Simulation parameters -------------- \n";
        std::cout << "skin depth (normalizing distance)     = " << formulas.skindepth*1e6 << " microns \n";
        std::cout << "plasma period (normalizing time)      = " << 1.0/formulas.wp*1e15 << " fs \n";
        std::cout << "Time step     (normalizing time)      = " << h << " plasma periods \n";
        std::cout << "\n";
        std::cout << "\n";
        


        return 1.0/formulas.wp*1e15;
    }

//**************************************************************
int main(int argc, char** argv) {

    MPI_Init(&argc,&argv);

///  Initiate the Parallel Environment and decompose the Computational Domain
    Parallel_Environment_1D PE;

    time_t tstart, tend;
    tstart = time(0);


///  Set up the grid
///	Moves all of the relevant data from the input deck into a single container
    Grid_Info grid(Input::List().ls, Input::List().ms,
    				Input::List().mass, Input::List().qs,
					Input::List().xminLocal, Input::List().xmaxLocal, Input::List().NxLocal,
					Input::List().xminGlobal, Input::List().xmaxGlobal, Input::List().NxGlobal,
                    Input::List().pmax, Input::List().ps, Input::List().Npx, Input::List().Npy, Input::List().Npz);


///  Clock
    int tout_start;
    if (Input::List().isthisarestart) tout_start = Input::List().restart_time;
    else  tout_start = 0;

    double dt_out(Input::List().t_stop / Input::List().n_outsteps);
    double CLF(Input::List().clf_dp);
    size_t n_outsteps(Input::List().n_outsteps);
    double numh  = size_t(static_cast<int>(dt_out/CLF))+1;
    double h     = dt_out/static_cast<double>(numh);

///  INITIALIZATION
	if (PE.RANK() == 0) {
		Export_Files::Folders();
	}

    Export_Files::Restart_Facility Re(PE.RANK());

    State1D Y( grid.axis.Nx(0), Input::List().ls, Input::List().ms, grid.Np, Input::List().pmax, grid.charge, grid.mass, Input::List().hydromass, Input::List().hydrocharge);
    Y = 0.0;
    Setup_Y::initialize(Y, grid);


    collisions collide(Y,h);
//    InverseBremsstrahlung IB(Y.DF(0).pmax(),Y.SH(0,0,0).nump(),Y.SH(0,0,0).numx(),tout_start,grid.axis.x(0));
    Hydro_Functor         HydroFunc(grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

    Algorithms::RK3<State1D> RK(Y);

    // Algorithms::LEAPs<State1D> LEAP(Y);
    // Algorithms::PEFRL<State1D> MAGIC(Y);


   if (Input::List().isthisarestart) Re.Read(PE.RANK(),tout_start,Y);

    Output_Data::Output_Preprocessor_1D  output( grid, Input::List().oTags);
    output( Y, grid, tout_start, PE );
    output.distdump(Y, grid, tout_start, PE );


    double plasmaperiod;
    if (!(PE.RANK())){
        plasmaperiod = startmessages(Y);
    }

// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
//  ITERATION LOOP 
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
    if (Input::List().implicit_B) {
        if (Input::List().implicit_E) {
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // IMPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_implicitE_implicitB_p1 impE_p1_Functor(Input::List().ls, Input::List().ms, Input::List().pmax,
                                                         Input::List().ps,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

            VlasovFunctor1D_implicitE_p2 impE_p2_Functor(Input::List().ls, Input::List().ms, Input::List().pmax,
                                                         Input::List().ps,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

//            Magnetic_Field_1D Bfield(Input::List().ls[0], Input::List().ms[0], Input::List().pmax[0]/(2.0*Input::List().ps[0]-1), Input::List().pmax[0],
//                                     Input::List().ps[0],
//                                     grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));


            Magnetic_Field_1D Bfield(Input::List().ls[0], Input::List().ms[0],
                                     0.0, Input::List().pmax[0], Input::List().ps[0],
                                     grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));



            // --------------------------------------------------------------------------------------------------------------------------------
            using Electric_Field_Methods::Efield_Method;
            Electric_Field_Methods::Implicit_E_Field eim(Y.EMF(), h,
                                                         (grid.axis.xmax(0) - grid.axis.xmin(0)) / grid.axis.Nx(0));

            if (!Input::List().collisions) {
                if (!PE.RANK())
                    std::cout << "\n Need collisions for implicit E field solver. \n Exiting. \n";
                exit(0);
            }

            for (size_t t_out(tout_start + 1); t_out < n_outsteps + 1; ++t_out) {
                // --------------------------------------------------------------------------------------------------------------------------------
                for (Clock W(t_out - 1, dt_out, CLF); W.tick() < W.numh(); ++W) {
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // if (!(PE.RANK())) {
                    //     // cout << "Time = " <<  W.time()<< "\n";
                    //     printf("\r Time = %4.4f tau_p = %4.4f fs", W.time(), W.time() * plasmaperiod);
                    //     fflush(stdout);
                    // }
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    //                                   the guts
                    // --------------------------------------------------------------------------------------------------------------------------------
                    if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, W.time());

                    Y = RK(Y, W.h(),&impE_p1_Functor);                                                              /// Vlasov - Updates the distribution function: Spatial Advection and B Field "action".
                    PE.Neighbor_ImplicitE_Communications(Y);                                                        /// Boundaries
                    eim.advance(&RK, Y, collide,&impE_p2_Functor);                                                  /// Finds new electric field
                    Y = RK(Y, W.h(),&impE_p2_Functor);

                    /// Uses new electric field to push distribution functions
                    for (size_t s(0); s < Y.Species(); ++s)
                    {
                        Bfield.implicit(Y.DF(s),Y.EMF().Bx(),Y.EMF().By(),Y.EMF().Bz(),W.h());
                        Y.DF(s) = Y.DF(s).Filterp();
                    }

                    collide.advance(Y,W);                                                                             /// Fokker-Planck
                                                                                                                      
                    if (Input::List().hydromotion)
                        Y = RK(Y, W.h(), &HydroFunc);                                                               /// Hydro Motion
                                                                                                    /// Error Checking
                    PE.Neighbor_Communications(Y);                                                                  /// Boundaries
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                }

                if (!(PE.RANK())) cout << " \n Output #" << t_out << "\n";
                output(Y, grid, t_out, PE);
                Y.checknan();

                if (!(t_out%Input::List().n_distoutsteps))
                    output.distdump(Y, grid, t_out, PE);
            }
        } else {
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // EXPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_explicitE_implicitB rkF(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps,
                                          grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

            Magnetic_Field_1D Bfield(Input::List().ls[0], Input::List().ms[0],
                                     0.0, Input::List().pmax[0], Input::List().ps[0],
                                     grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

            // --------------------------------------------------------------------------------------------------------------------------------
            for (size_t t_out(tout_start + 1); t_out < n_outsteps + 1; ++t_out) {
                for (Clock W(t_out - 1, dt_out, CLF); W.tick() < W.numh(); ++W) {
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // if (!(PE.RANK())) {

                    //     printf("\r Time = %4.4f tau_p = %4.4f fs", W.time(), W.time() * plasmaperiod);
                    //     fflush(stdout);
                    // }
                    // --------------------------------------------------------------------------------------------------------------------------------
                    //                                   the guts
                    // --------------------------------------------------------------------------------------------------------------------------------

                    if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, W.time());

                    Y = RK(Y, W.h(), &rkF);                                                                         ///  Vlasov          //

                    for (size_t s(0); s < Y.Species(); ++s)
                    {
                        Bfield.implicit(Y.DF(s),Y.EMF().Bx(),Y.EMF().By(),Y.EMF().Bz(),W.h());
                        Y.DF(s) = Y.DF(s).Filterp();
                    }

                    if (Input::List().collisions)
                        collide.advance(Y,W);                                               ///  Fokker-Planck   //

                    if (Input::List().hydromotion)
                        Y = RK(Y, W.h(), &HydroFunc);                                      ///  Hydro           //

                    PE.Neighbor_Communications(Y);                                                                  ///  Boundaries      //
                }

                if (!(PE.RANK())) cout << " \n Output #" << t_out << "\n";
                output(Y, grid, t_out, PE);
                Y.checknan();

                if (!(t_out%Input::List().n_distoutsteps))
                    output.distdump(Y, grid, t_out, PE);


            }
        }
    }
    else {
        if (Input::List().implicit_E) {
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // IMPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_implicitE_p1 impE_p1_Functor(Input::List().ls, Input::List().ms, Input::List().pmax,
                                                         Input::List().ps,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            VlasovFunctor1D_implicitE_p2 impE_p2_Functor(Input::List().ls, Input::List().ms, Input::List().pmax,
                                                         Input::List().ps,
                                                         grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            // --------------------------------------------------------------------------------------------------------------------------------
            using Electric_Field_Methods::Efield_Method;
            Electric_Field_Methods::Implicit_E_Field eim(Y.EMF(), h,
                                                         (grid.axis.xmax(0) - grid.axis.xmin(0)) / grid.axis.Nx(0));

            if (!Input::List().collisions) {
                if (!PE.RANK())
                    std::cout << "\n Need collisions for implicit E field solver. \n Exiting. \n";
                exit(0);
            }

            for (size_t t_out(tout_start + 1); t_out < n_outsteps + 1; ++t_out) {
                // --------------------------------------------------------------------------------------------------------------------------------
                for (Clock W(t_out - 1, dt_out, CLF); W.tick() < W.numh(); ++W) {
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // if (!(PE.RANK())) {
                    //     // cout << "Time = " <<  W.time()<< "\n";
                    //     printf("\r Time = %4.4f tau_p = %4.4f fs", W.time(), W.time() * plasmaperiod);
                    //     fflush(stdout);
                    // }
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    //                                   the guts
                    // --------------------------------------------------------------------------------------------------------------------------------
                    if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, W.time());

                    Y = RK(Y, W.h(),
                           &impE_p1_Functor);                                                             /// Vlasov - Updates the distribution function: Spatial Advection and B Field "action".
                    PE.Neighbor_ImplicitE_Communications(Y);                                            /// Boundaries
                    eim.advance(&RK, Y, collide,
                                &impE_p2_Functor);                                                 /// Finds new electric field
                    Y = RK(Y, W.h(),
                           &impE_p2_Functor);                                                             /// Uses new electric field to push distribution functions
                    collide.advance(Y,W);                                                                             /// Fokker-Planck
//                    if (Input::List().inverse_bremsstrahlung)
//                        IB.loop(Y.SH(0, 0, 0), Y.HYDRO().Zarray(), W.time());     /// Explicit IB heating for f00
                    if (Input::List().hydromotion)
                        Y = RK(Y, W.h(), &HydroFunc);                                      /// Hydro Motion
                    // Y.checknan();                                                                                /// Error Checking
                    PE.Neighbor_Communications(
                            Y);                                                                  /// Boundaries
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // --------------------------------------------------------------------------------------------------------------------------------
                }

                if (!(PE.RANK())) cout << " \n Output #" << t_out << "\n";
                output(Y, grid, t_out, PE);
                Y.checknan();

                if (!(t_out%Input::List().n_distoutsteps))
                    output.distdump(Y, grid, t_out, PE);
            }
        } else {
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------------------------------
            // EXPLICIT E-FIELD
            // --------------------------------------------------------------------------------------------------------------------------------
            VlasovFunctor1D_explicitE rkF(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps,
                                          grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

            // VlasovFunctor1D_spatialpush SA(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps,
            //                               grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));

            // VlasovFunctor1D_momentumpush PA(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps,
            //                               grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
            // --------------------------------------------------------------------------------------------------------------------------------
            for (size_t t_out(tout_start + 1); t_out < n_outsteps + 1; ++t_out) {
                for (Clock W(t_out - 1, dt_out, CLF); W.tick() < W.numh(); ++W) {
            
                    // --------------------------------------------------------------------------------------------------------------------------------
                    // if (!(PE.RANK())) {

                    //     printf("\r Time = %4.4f tau_p = %4.4f fs", W.time(), W.time() * plasmaperiod);
                    

                    // }
                    // --------------------------------------------------------------------------------------------------------------------------------
                    //                                   the guts
                    // --------------------------------------------------------------------------------------------------------------------------------
                
                    if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, W.time());
                    // if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid, Y, t);
                    if (Input::List().trav_wave) Setup_Y::applytravelingwave(grid, Y, W.time());

                    Y = RK(Y, W.h(), &rkF);                                                 ///  Vlasov          //
                    // Y = LEAP(Y, W.h(), &SA, &PA);                                                 ///  Vlasov          //
                    // Y = MAGIC(Y, W.h(), &SA, &PA);                                                 ///  Vlasov          //
                    if (Input::List().collisions)
                        collide.advance(Y,W);                                               ///  Fokker-Planck   //

                    if (Input::List().hydromotion)
                        Y = RK(Y, W.h(), &HydroFunc);                                      ///  Hydro           //

                    PE.Neighbor_Communications(Y);                                         ///  Boundaries      //
                }
                

                if (!(PE.RANK())) cout << " \n Output #" << t_out << "\n";

                output(Y, grid, t_out, PE);
                Y.checknan();

                if (!(t_out%Input::List().n_distoutsteps))
                    output.distdump(Y, grid, t_out, PE);

                if (!(t_out%Input::List().n_restarts))
                    Re.Write(PE.RANK(), t_out, Y);
            }
        }
    }
    tend = time(0); 
    if (!(PE.RANK())){
        cout << "Simulation took "<< difftime(tend, tstart) <<" second(s)."<< endl;
    }
    MPI_Finalize();
	return 0;
}
//**************************************************************
//--------------------------------------------------------------------------- 
