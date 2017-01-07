/*! \brief Main Loop
 * \author Michail Tzoufras, Archis Joglekar, Benjamin Winjum
 * \date   September 1, 2016
 * \file   main.cpp
 * 
 * This file includes
 * - the main loop.
 * - the definition for the clock class
 * - the definition of periodic boundaries
 * 
 * \todo Change clock to have fixed timestep from input deck and all corresponding changes. Rename CLF to timestep. Implicit solver timestep is FUBAR
 * \todo Cleanup home directory - new makefile, source directory, version output
 * \todo Choose output files by strings rather than toggles
 * \todo Implement various useful moments i.e. Resistivity, Nernst velocity, Anisotropic Pressure
 * \todo 2D - Field Solver
 * \todo Hydro ion motion
 * \todo Kinetic ions
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
#include "interspeciescollisions.h"
#include "parallel.h"
#include "implicitE.h"
#include "laser.h"
#include "export.h"




//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can 
//  be incremented; 

    class Clock {
    public:
//      Constructor
        Clock(int tout_start, double dt_out, double CFL) {
            _hstep = 0;
            _t_start = double(tout_start)*dt_out;
            _numh  = size_t(static_cast<int>(dt_out/CFL))+1;
            _h     = dt_out/static_cast<double>(_numh);
        }

//      Clock readings
        double h()     const {return _h;}                    // Step size
        size_t numh()  const {return _numh;}                 // # steps 
        size_t tick()  const {return _hstep;}                // Current step
        double time()  const {return tick()*h() + _t_start;} // Current time
             
//      Increment time
        Clock& operator++() { ++_hstep; return *this;}  

    private:
        double _t_start;   // Initial output timestep
        size_t _numh;      // # of steps derived from this CFL
        double _h;         // resulting time-step from CFL
        size_t _hstep;
    };
//--------------------------------------------------------------
//**************************************************************

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
        std::cout << "------------------- UCLA - 2016 -------------------- \n";
        
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

//  Initiate the Paralell Environment / Decompose the Computational Domain  
    Parallel_Environment_1D PE; 

    time_t tstart, tend; 
    tstart = time(0);

    // std::cout << "xminLocal = " << Input::List().xminLocal[0] << " , xmaxLocal = " << Input::List().xmaxLocal[0] << "\n";    
//  Set up the grid
//	Moves all of the relevant data from the input deck into a single container 
    Grid_Info grid(Input::List().ls, Input::List().ms, 
    				Input::List().mass, Input::List().qs, 
					Input::List().xminLocal, Input::List().xmaxLocal, Input::List().NxLocal,
					Input::List().xminGlobal, Input::List().xmaxGlobal, Input::List().NxGlobal,
                    Input::List().pmax, Input::List().ps, Input::List().Npx, Input::List().Npy, Input::List().Npz);


//  Clock
    int tout_start(0);

    double dt_out(Input::List().t_stop / Input::List().n_outsteps);
    double CLF(Input::List().clf_dp);
    size_t n_outsteps(Input::List().n_outsteps);
    double numh  = size_t(static_cast<int>(dt_out/CLF))+1; 
    double h     = dt_out/static_cast<double>(numh);


// WARNING::::
// XMIN AND XMAX AND NX CURRENTLY HAVE TWO ELEMENTS BUT SECOND IS MEANINGLESS SO BEWARE THEIR ACCESS
// ALSO, XMIN AND XMAX AND NX FROM THE INPUT DECK ARE EXPANDED TO INCLUDE GUARD CELLS...
// ... SO ACCESSING THEM IN THE GRID OBJECT WILL GIVE THE ENTIRE SPACE + GUARD CELLS.
//.... USE XGMIN AND XGMAX AND NXG TO ACCESS ONLY THE REAL PHYSICAL SPACE PARAMETERS

//  INITIALIZATION
	if (PE.RANK() == 0) {
		Export_Files::Folders();
	    Export_Files::Restart_Facility Re;
	}

    State1D Y( grid.axis.Nx(0), Input::List().ls, Input::List().ms, grid.Np, Input::List().pmax, grid.charge, grid.mass, Input::List().hydromass, Input::List().hydrocharge);
    Y = 0.0;    
    Setup_Y::initialize(Y, grid);
    collisions collide(Y,h);  
    InverseBremsstrahlung IB(Y.DF(0).pmax(),Y.SH(0,0,0).nump(),Y.SH(0,0,0).numx(),tout_start,grid.axis.x(0));
    Hydro_Functor         HydroFunc(grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
 
    Algorithms::RK2<State1D> RK(Y);
    Output_Data::Output_Preprocessor_1D  output( grid, Input::List().oTags); 
    output( Y, grid, tout_start, PE );
    

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
    if (Input::List().implicit_E)
    {   
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// IMPLICIT E-FIELD
// --------------------------------------------------------------------------------------------------------------------------------
        VlasovFunctor1D_implicitE_p1 impE_p1_Functor(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps, 
                        grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
        VlasovFunctor1D_implicitE_p2 impE_p2_Functor(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps, 
                        grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));
// --------------------------------------------------------------------------------------------------------------------------------
        using Electric_Field_Methods::Efield_Method;
        Electric_Field_Methods::Implicit_E_Field eim(Y.EMF(),h, (grid.axis.xmax(0)-grid.axis.xmin(0))/grid.axis.Nx(0));
        
        for (size_t  t_out(tout_start+1); t_out < n_outsteps+1; ++t_out) {
// --------------------------------------------------------------------------------------------------------------------------------
            for (Clock W(t_out-1,dt_out,CLF); W.tick() < W.numh(); ++W) {
// --------------------------------------------------------------------------------------------------------------------------------                
                if (!(PE.RANK())){
                    // cout << "Time = " <<  W.time()<< "\n";    
                    printf("\r Time = %4.4f tau_p = %4.4f fs",W.time(),W.time()*plasmaperiod);
                    fflush(stdout); 
                }
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
//                                   the guts
// --------------------------------------------------------------------------------------------------------------------------------                
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid,Y,W.time());

                Y = RK(Y, W.h(), &impE_p1_Functor);                                             // Vlasov - Updates the distribution function: Spatial Advection and B Field "action".  
                                                        PE.Neighbor_ImplicitE_Communications(Y);
                collide.advancef0(Y);                                                           // Fokker-Planck
                eim.advance(&RK, Y, collide, &impE_p2_Functor);                                 // Finds new electric field
                Y = RK(Y, W.h(), &impE_p2_Functor);                                             // Uses new electric field to push distribution functions
                collide.advancef1(Y);                                                           // Rest of F-P
                collide.advanceflm(Y);
                
                if (Input::List().inverse_bremsstrahlung) IB.loop(Y.SH(0,0,0),W.time());        // Explicit IB heating for f00
                if (Input::List().hydromotion) Y = RK(Y,W.h(),&HydroFunc);                      // Hydro Motion

                // Y.checknan();                                                                // Error Checking
                
                // Boundaries and communications between processors
			    PE.Neighbor_Communications(Y);
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
            }
    		
            if (!(PE.RANK())){
                cout << " \n Output...\n";
            }
	        output( Y, grid, t_out, PE );
        }
    }
    else 
    {   
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------
// EXPLICIT E-FIELD
// --------------------------------------------------------------------------------------------------------------------------------        
        VlasovFunctor1D_explicitE rkF(Input::List().ls, Input::List().ms, Input::List().pmax, Input::List().ps, 
                        grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0)); 
// --------------------------------------------------------------------------------------------------------------------------------
        for (size_t  t_out(tout_start+1); t_out < n_outsteps+1; ++t_out) {
            for (Clock W(t_out-1,dt_out,CLF); W.tick() < W.numh(); ++W) {
// --------------------------------------------------------------------------------------------------------------------------------            
                if (!(PE.RANK())){
                    // cout << "Time = " <<  W.time()<< "\n";    
                    printf("\r Time = %4.4f tau_p = %4.4f fs",W.time(),W.time()*plasmaperiod);
                    fflush(stdout); 
                }
// --------------------------------------------------------------------------------------------------------------------------------
//                                   the guts
// --------------------------------------------------------------------------------------------------------------------------------
                
                if (Input::List().ext_fields) Setup_Y::applyexternalfields(grid,Y,W.time());
                // Setup_Y::applyexternalfields(Y,W.time());

                Y = RK(Y, W.h(), &rkF);                                                     ///  Vlasov          //
                
                collide.advancef0(Y);
                collide.advancef1(Y);                                                       ///  Fokker-Planck   //
                collide.advanceflm(Y);
            
                if (Input::List().hydromotion) Y = RK(Y,W.h(),&HydroFunc);                  ///  Hydro           //
                
                if (Input::List().inverse_bremsstrahlung) IB.loop(Y.SH(0,0,0),W.time());    ///  Inverse-Brem    //

                // Y.checknan();                                                               ///  Error-checking  //
                
                PE.Neighbor_Communications(Y);                                              ///  Boundaries      //
            }
            
            if (!(PE.RANK())){
                // std::cout << "\nBz = " << Y.EMF().Bz()(4);
                cout << " \n Output...\n";
            }
            output( Y, grid, t_out, PE );
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
