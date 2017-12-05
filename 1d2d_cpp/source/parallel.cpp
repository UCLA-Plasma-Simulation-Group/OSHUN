/*!\brief  Parallelization routines - Definitions
* \author  PICKSC
 * \date   2017
 * \file   parallel.cpp
 *
 * In here are the structures that enable parallelization
 * 
 * Periodic and reflecting have been implemented.
 * 
 * 
 */


//  Standard libraries
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <valarray>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <float.h>
#include <string>
#include <iomanip>
#include <fstream>

//  My libraries
#include "lib-array.h"
#include <map>

//  Declarations
#include "input.h"
#include "state.h"
#include "parallel.h"


//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Node_ImplicitE_Communications_1D:: Node_ImplicitE_Communications_1D() :
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Input::List().BoundaryCells),      // # of boundary cells
        bndX(Input::List().bndX) {        // Type of boundary in X

    // 3 components for Bx, By, Bz
    msg_sizeX = 3;

    msg_sizeX *= Nbc; //(IN().inp().y.dim()*Nbc);
    msg_bufX = new complex<double>[msg_sizeX];

}
//--------------------------------------------------------------

//--------------------------------------------------------------
Node_ImplicitE_Communications_1D:: ~Node_ImplicitE_Communications_1D(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] msg_bufX;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
int Node_ImplicitE_Communications_1D:: BNDX()  const {return bndX;}
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::Send_right_X(State1D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

    static size_t step_f(Nbc);
    size_t bufind(0);

    // Fields:   x0 "Right-Bound --> "
    for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.FLD(i)(Y.FLD(i).numx()-2*Nbc+e);
        }
        bufind += step_f;
    }

    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::Recv_from_left_X(State1D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

    static size_t step_f(Nbc);
    size_t bufind(0);
    MPI_Status status;

    // Receive Data
    MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

    // Fields:   x0-"---> Left-Guard"
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        for(int e(0); e < Nbc; e++) {
            Y.FLD(i)(e) = msg_bufX[bufind + e];
        }
        bufind += step_f;
    }
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::Send_left_X(State1D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

    static size_t step_f(Nbc);
    size_t bufind(0);

    // Fields:   x0 " <--- Left-Bound "
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.FLD(i)(Nbc+e);
        }
        bufind += step_f;
    }

    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::Recv_from_right_X(State1D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

    static size_t step_f(Nbc);
    size_t bufind(0);
    MPI_Status status;

    // Receive Data
    MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

    // Fields:   x0-"Right-Guard <--- "
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        for(int e(0); e < Nbc; e++) {
            Y.FLD(i)(Y.FLD(i).numx()-Nbc+e) = msg_bufX[bufind + e];
        }
        bufind += step_f;
    }
}
//--------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//*************************************************************************************
//-------------------------------------------------------------------------------------
void Node_ImplicitE_Communications_1D::mirror_bound_Xleft(State1D& Y) {
//-------------------------------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//-------------------------------------------------------------------------------------

    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the fields
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(c) = Y.FLD(i)(2*Nbc-c-1);
    }

    //Bx
    for(int c(0); c < Nbc; ++c) {
        Y.EMF().Bx()(c) *= -1.0; // left  boundary
    }

}
//-------------------------------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::mirror_bound_Xright(State1D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------

    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the fields
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(Nx-c-1) = Y.FLD(i)(Nx-2*Nbc+c);
    }

    //Bx
    for(int c(0); c < Nbc; ++c) {
        Y.EMF().Bx()(Nx-1-c) *= -1.0; // right  boundary
    }

}
//--------------------------------------------------------------


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::sameNode_bound_X(State1D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
    switch (BNDX()) {
        case 0:                   // periodic
            sameNode_periodic_X(Y);

            break;
        case 1:                   // mirror boundary
            sameNode_mirror_X(Y);
            break;
        default:
            cout<<"Not a valid boundary condition." << endl;
            break;
    }
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::sameNode_periodic_X(State1D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------

    // Fields:   x0 "Right-Bound ---> Left-Guard"
    for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
        for(int c(0); c < Nbc; c++) {
            Y.FLD(i)(c) = Y.FLD(i)(Y.EMF().Ex().numx()-2*Nbc+c);
        }

    // Fields:   x0 "Left-Bound ---> Right-Guard"
    for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
        for(int c(0); c < Nbc; c++) {
            Y.FLD(i)(Y.EMF().Ex().numx()-Nbc+c) = Y.FLD(i)(Nbc+c);
        }

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_ImplicitE_Communications_1D::sameNode_mirror_X(State1D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------

    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the fields
    for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
        // right boundary
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(Nx-c-1) = Y.FLD(i)(Nx-2*Nbc+c);
        // left boundary
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(c) = Y.FLD(i)(2*Nbc-c-1);
    }

    //Bx
    for (int c(0); c < Nbc; ++c) {
        Y.EMF().Bx()(c) *= -1.0; // left  boundary
        Y.EMF().Bx()(Nx-Nbc+c) *= -1.0; // right  boundary
    }

}
//--------------------------------------------------------------


//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Node_Communications_1D:: Node_Communications_1D() :
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Input::List().BoundaryCells),      // # of boundary cells
        bndX(Input::List().bndX) {        // Type of boundary in X

    numspec = Input::List().ls.size();
    // numpmax = Input::List().ps;

    double temp;

    // # of harmonics
    msg_sizeX = 0; temp = 0;
    for(int s(0); s < numspec; ++s) {
        temp = ((Input::List().ms[s]+1)*(2*Input::List().ls[s]-Input::List().ms[s]+2))/2;
        msg_sizeX += temp*(Input::List().dp[s]).size();
    }
    // (# of harmonics) * (# cells in p)
    // msg_sizeX *= numpmax
    // 6 fields: Ex, Ey, Ez, Bx, By, Bz
    msg_sizeX += 6;

    // 6 Hydro Quantities, density, vx, vy, vz, temp, Z
    if (Input::List().hydromotion)
    {
        msg_sizeX += 6;
    }

    msg_sizeX *= Nbc;  //(IN().inp().y.dim()*Nbc);
                       //
    if (Input::List().particlepusher)
    {
        msg_sizeX += Input::List().numparticles; // Going Left or Right?
        msg_sizeX += Input::List().numparticles; // Position
        msg_sizeX += Input::List().numparticles; // X-momentum
        msg_sizeX += Input::List().numparticles; // Y-momentum
        msg_sizeX += Input::List().numparticles; // Z-momentum
    }
    
    msg_bufX = new complex<double>[msg_sizeX];

}
//--------------------------------------------------------------

//--------------------------------------------------------------
Node_Communications_1D:: ~Node_Communications_1D(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] msg_bufX;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
int Node_Communications_1D:: BNDX()  const {return bndX;}
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_Communications_1D::Send_right_X(State1D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

    // static size_t step_h(numpmax*Nbc);
    static size_t step_f(Nbc);
    size_t bufind(0);

    // Harmonics:x0 "Right-Bound ---> "
    for(int s(0); s < Y.Species(); ++s) {
        for(size_t i = 0; i < Y.DF(s).dim(); ++i){
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int e(0); e < Nbc; e++) {
                    msg_bufX[bufind + e] = (Y.DF(s)(i))(p, Y.FLD(0).numx()-2*Nbc+e);
                }
                bufind += step_f;
            }
        }
    }
    // Fields:   x0 "Right-Bound --> "
    for(size_t i = 0; i < Y.EMF().dim(); ++i){
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.FLD(i)(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;
    }

    if (Input::List().hydromotion)
    {
        // Hydro-velocity:   x0 "Right-Bound --> "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().density(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;
        // Hydro-velocity:   x0 "Right-Bound --> "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vx(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;

        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vy(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;

        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vz(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;

        // Hydro-temperature:   x0 "Right-Bound --> "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().temperature(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;

        // Hydro-chargefraction:   x0 "Right-Bound --> "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().Z(Y.FLD(0).numx()-2*Nbc+e);
        }
        bufind += step_f;
    }

    if (Input::List().particlepusher)
    {
        for (int ip(0); ip < Y.particles().numpar(); ++ip)
        {          
            // if (Y.particles().x(ip) > Input::List().xmaxGlobal[0]){
            //     std::cout << " global max = " << Input::List().xmaxGlobal[0] << "\n";
            //     Y.particles().x(ip) -= Input::List().xmaxGlobal[0]-Input::List().xminGlobal[0];
            // }
            // if (Y.particles().goingright(ip) == 1)
            // {
                msg_bufX[bufind]   = Y.particles().goingright(ip);
                msg_bufX[bufind+1] = Y.particles().x(ip);
                msg_bufX[bufind+2] = Y.particles().px(ip);
                msg_bufX[bufind+3] = Y.particles().py(ip);
                msg_bufX[bufind+4] = Y.particles().pz(ip);    

                // Y.particles().goingright(ip) = 0;
            // }
            

            bufind += 4;
        }
    }

        
    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);


}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::Recv_from_left_X(State1D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

    // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
    static size_t step_f(Nbc);
    size_t bufind(0);
    
    MPI_Status status;

    // Receive Data
    MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

    // Harmonics:x0-"---> Left-Guard"
    for(int s(0); s < Y.Species(); ++s) {
        for(int i = 0; i < Y.DF(s).dim(); ++i){
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int e(0); e < Nbc; e++) {
                    (Y.DF(s)(i))(p, e) = msg_bufX[bufind + e];
                }
                bufind += step_f;
            }
        }
    }


    // Fields:   x0-"---> Left-Guard"
    for(int i = 0; i < Y.EMF().dim(); ++i){
        for(int e(0); e < Nbc; e++) {
            Y.FLD(i)(e) = msg_bufX[bufind + e];
        }
        bufind += step_f;
    }

    if (Input::List().hydromotion)
    {
        // Hydro-velocity:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().density(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vx(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vy(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vz(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-temperature:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().temperature(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-charge fraction:   x0-"---> Left-Guard"
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().Z(e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;
    }


    if (Input::List().particlepusher)
    {   
        
        // if (Input::List().xminLocalnobnd[0] == 0.0)
        // {
        //     std::cout << "\n is here = " << Y.particles().ishere(2) << "\n";
        // }

        for (int ip(0); ip < Y.particles().numpar(); ++ip){
            


            // See if particle is now in local box, but wasn't before. 
            // Only overwrite if so.
            // If not, step over buffer and go to next particle.
            
            // if (((msg_bufX[bufind]).real() < Input::List().xmaxLocalnobnd[0] && (msg_bufX[bufind]).real() >= Input::List().xminLocalnobnd[0])
            //     // && (Y.particles().x(ip) >= Input::List().xmaxLocal || Y.particles().x(ip) < Input::List().xminLocal[0]))
            //     && !(Y.particles().ishere(ip)) && (Y.particles().goingleft(ip)))
            //     
            if (int (msg_bufX[bufind].real()) == 1)
            {
                
                // if (ip == 2){
                // std::cout << "\n\n\n origin = " << origin << "\n";
                // std::cout << "in cell = " << Input::List().xminLocalnobnd[0] << "\n" ;
                // std::cout << "is it here? "<< Y.particles().ishere(ip) << "\n\n";  

                if (msg_bufX[bufind+1].real() >= Input::List().xmaxGlobal[0]){
                    Y.particles().x(ip) = (msg_bufX[bufind+1]).real() - Input::List().xmaxGlobal[0] + Input::List().xminGlobal[0];
                }
                else Y.particles().x(ip)  = (msg_bufX[bufind+1]).real();


                Y.particles().px(ip) = (msg_bufX[bufind+2]).real();
                Y.particles().py(ip) = (msg_bufX[bufind+3]).real();
                Y.particles().pz(ip) = (msg_bufX[bufind+4]).real();
                Y.particles().ishere(ip) = 1;
            }
            
            bufind+=4;  // Go to next particle in buffer
        }
    }



}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::Send_left_X(State1D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

    // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
    static size_t step_f(Nbc);
    size_t bufind(0);

    // Harmonics:x0 " <--- Left-Bound "
    for(int s(0); s < Y.Species(); ++s) {
        for(int i = 0; i < Y.DF(s).dim(); ++i){
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int e(0); e < Nbc; e++) {
                    msg_bufX[bufind + e] = (Y.DF(s)(i))(p, Nbc+e);
                }
                bufind += step_f;
            }
        }
    }
    // Fields:   x0 " <--- Left-Bound "
    for(int i = 0; i < Y.EMF().dim(); ++i){
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.FLD(i)(Nbc+e);
        }
        bufind += step_f;
    }

    if (Input::List().hydromotion)
    {
        // Hydro-velocity:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().density(Nbc+e);
        }
        bufind += step_f;
        // Hydro-velocity:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vx(Nbc+e);
        }
        bufind += step_f;

        // Hydro-velocity:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vy(Nbc+e);
        }
        bufind += step_f;

        // Hydro-velocity:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().vz(Nbc+e);
        }
        bufind += step_f;

        // Hydro-temperature:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().temperature(Nbc+e);
        }
        bufind += step_f;

        // Hydro-charge-fraction:   x0 " <--- Left-Bound "
        for(int e(0); e < Nbc; e++) {
            msg_bufX[bufind + e] = Y.HYDRO().Z(Nbc+e);
        }
        bufind += step_f;
    }

    

    if (Input::List().particlepusher)
    {   

        // Send all
        // Receive routine is the discriminator
        for (int ip(0); ip < Y.particles().numpar(); ++ip)
        {          
            // if (Y.particles().goingright(ip) == -1){
            //     Y.particles().x(ip) += Input::List().xmaxGlobal[0]-Input::List().xminGlobal[0];

            msg_bufX[bufind] = Y.particles().goingright(ip);
            msg_bufX[bufind+1] = Y.particles().x(ip);
            msg_bufX[bufind+2] = Y.particles().px(ip);
            msg_bufX[bufind+3] = Y.particles().py(ip);
            msg_bufX[bufind+4] = Y.particles().pz(ip);

            bufind += 4;
        }
    }

    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::Recv_from_right_X(State1D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

    // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
    static size_t step_f(Nbc);
    size_t bufind(0);
    MPI_Status status;

    // Receive Data
    MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

    // Harmonics:x0-"Right-Guard <--- "
    for(int s(0); s < Y.Species(); ++s) {
        for(int i = 0; i < Y.DF(s).dim(); ++i){
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int e(0); e < Nbc; e++) {
                    (Y.DF(s)(i))(p, Y.FLD(0).numx()-Nbc+e) = msg_bufX[bufind + e];
                }
                bufind += step_f;
            }
        }
    }
    // Fields:   x0-"Right-Guard <--- "
    for(int i = 0; i < Y.EMF().dim(); ++i){
        for(int e(0); e < Nbc; e++) {
            Y.FLD(i)(Y.FLD(0).numx()-Nbc+e) = msg_bufX[bufind + e];
        }
        bufind += step_f;
    }

    if (Input::List().hydromotion)
    {
        // Hydro-density:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().density(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vx(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vy(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-velocity:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().vz(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-temperature:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().temperature(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;

        // Hydro-charge fraction:   x0-"Right-Guard <--- "
        for(int e(0); e < Nbc; e++) {
            Y.HYDRO().Z(Y.FLD(0).numx()-Nbc+e) = (msg_bufX[bufind + e]).real();
        }
        bufind += step_f;
    }


    
    if (Input::List().particlepusher)
    {   
        for (int ip(0); ip < Y.particles().numpar(); ++ip){
            

            // See if particle is now in local box, but wasn't before. 
            // Only overwrite if so.
            // If not, step over buffer and go to next particle.
            
            // if (((msg_bufX[bufind]).real() < Input::List().xmaxLocalnobnd[0] && (msg_bufX[bufind]).real() >= Input::List().xminLocalnobnd[0])
            //     // && (Y.particles().x(ip) >= Input::List().xmaxLocal || Y.particles().x(ip) < Input::List().xminLocal[0]))
            //     && !(Y.particles().ishere(ip)))
            // {

                // if (ip == 2){
                //     std::cout << "\n\n\n origin = " << origin << "\n";
                //     std::cout << "in cell = " << Input::List().xminLocalnobnd[0] << "\n" ;
                //     std::cout << "is it here? "<< Y.particles().ishere(ip) << "\n\n";  
                // } 
            if (msg_bufX[bufind].real() == -1)
            {

                if (msg_bufX[bufind+1].real() < Input::List().xminGlobal[0]){
                    Y.particles().x(ip) = (msg_bufX[bufind+1]).real() + Input::List().xmaxGlobal[0] - Input::List().xminGlobal[0];
                }
                else Y.particles().x(ip)  = (msg_bufX[bufind+1]).real();

                
                Y.particles().px(ip) = (msg_bufX[bufind+2]).real();
                Y.particles().py(ip) = (msg_bufX[bufind+3]).real();
                Y.particles().pz(ip) = (msg_bufX[bufind+4]).real();
                Y.particles().ishere(ip) = 1;
            }
            
            bufind+=4;  // Go to next particle in buffer
        }
    }

}
//--------------------------------------------------------------


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//**************************************************************
//--------------------------------------------------------------
void Node_Communications_1D::mirror_bound_Xleft(State1D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//--------------------------------------------------------------
    int sign(1);
    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the harmonics
    for(int s(0); s < Y.Species(); ++s) {
        for(int l(0); l < Y.DF(s).l0(); ++l){
            for(int m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                for (int c(0); c < Nbc; ++c) {
                    for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        Y.SH(s,l,m)(p, c) = Y.SH(s,l,m)(p, 2*Nbc-c-1);
                        Y.SH(s,l,m)(p, c) *= sign;
                    }
                }

            }
        }
    }

    // Mirror the fields
    for(int i(0); i < Y.EMF().dim(); ++i){
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(c) = Y.FLD(i)(2*Nbc-c-1);
    }

    for (int c(0); c < Nbc; ++c) {
// 			Ey
        Y.EMF().Ey()(c) *= -1.0; // left  boundary
// 			Ez
        Y.EMF().Ez()(c) *= -1.0; // left  boundary
// 			Bx
        Y.EMF().Bx()(c) *= -1.0; // left  boundary

//        Y.EMF().By()(c) = 1e-3; // left  boundary
    }

    if (Input::List().hydromotion)
    {
        // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().density(c) =  Y.HYDRO().density(2*Nbc-c-1);
            Y.HYDRO().temperature(c) =  Y.HYDRO().temperature(2*Nbc-c-1);
            Y.HYDRO().Z(c) =  Y.HYDRO().Z(2*Nbc-c-1);

            Y.HYDRO().vx(c) *= -1.0;
        }
    }

    if (Input::List().particlepusher)
    {
        for (int ip(0); ip < Y.particles().numpar(); ++ip){
            if (Y.particles().x(ip) < Input::List().xminLocalnobnd[0])  
            {
                Y.particles().x(ip) = Input::List().xminLocalnobnd[0] + (Input::List().xminLocalnobnd[0] - Y.particles().x(ip));
                Y.particles().px(ip) *= -1.0;
                Y.particles().ishere(ip) = 1;
            }
        }
    }


}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::mirror_bound_Xright(State1D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------
    int sign(1);
    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the harmonics
    for(int s(0); s < Y.Species(); ++s) {
        for(int l(0); l < Y.DF(s).l0(); ++l){
            for(int m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                    for (int c(0); c < Nbc; ++c) {
                        Y.SH(s,l,m)(p, Nx-c-1) = Y.SH(s,l,m)(p, Nx-2*Nbc+c);
                        Y.SH(s,l,m)(p, Nx-c-1) *= sign;
                    }
                }

            }
        }
    }

    // Mirror the fields
    for(int i(0); i < Y.EMF().dim(); ++i){
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(Nx-c-1) = Y.FLD(i)(Nx-2*Nbc+c);
    }

    for (int c(0); c < Nbc; ++c) {
        //Ey
        Y.EMF().Ey()(Nx-c-1) *= -1.0; // right boundary
        //Ez
        Y.EMF().Ez()(Nx-c-1) *= -1.0; // right boundary
        //Bx
        Y.EMF().Bx()(Nx-c-1) *= -1.0; // right boundary

//        Y.EMF().By()(Nx-c-1) = 1e-3; // right  boundary
    }

    if (Input::List().hydromotion)
    {
        // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().density(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().density(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().temperature(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().temperature(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().Z(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().Z(Y.EMF().Ex().numx()-2*Nbc+c);
            

            Y.HYDRO().vx(Y.EMF().Ex().numx()-Nbc+c) *=  -1.0;
        }
    }

    if (Input::List().particlepusher)
    {
        for (int ip(0); ip < Y.particles().numpar(); ++ip){
            if (Y.particles().x(ip) > Input::List().xmaxLocalnobnd[0])  
            {
                Y.particles().x(ip) = Input::List().xmaxLocalnobnd[0] - (Y.particles().x(ip) - Input::List().xmaxLocalnobnd[0]);
                Y.particles().px(ip) *= -1.0;
                Y.particles().ishere(ip) = 1;
            }
        }
    }

}
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_Communications_1D::sameNode_bound_X(State1D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------

    switch (BNDX()) {
        case 0:                   // periodic
            sameNode_periodic_X(Y);
            break;
        case 1:                   // mirror boundary
            sameNode_mirror_X(Y);
            break;
        // case 2:                   // mirror boundary
            // sameNode_PML_X(Y);
            // break;
        default:
            cout<<"Not a valid boundary condition." << endl;
            break;
    }
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::sameNode_periodic_X(State1D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------

    // Harmonics:x0 "Right-Bound ---> Left-Guard"
    for(int s(0); s < Y.Species(); ++s) {
        for(int i(0); i < Y.DF(s).dim(); ++i) {
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int c(0); c < Nbc; c++) {
                    (Y.DF(s)(i))(p, c) = (Y.DF(s)(i))(p, Y.EMF().Ex().numx()-2*Nbc+c);
                    // std::cout << "\n 1: DF(" << i << "," << p << "," << c << ") = " << (Y.DF(s)(i))(p, Y.EMF().Ex().numx()-2*Nbc+c) << "\n";
                }
            }
        }
    }
    // Fields:   x0 "Right-Bound ---> Left-Guard"
    for(int i(0); i < Y.EMF().dim(); ++i) {
        for(int c(0); c < Nbc; c++) {
            Y.FLD(i)(c) = Y.FLD(i)(Y.EMF().Ex().numx()-2*Nbc+c);
            // std::cout << "\n 2: FLD(" << i << "," << c << ") = " << Y.FLD(i)(Y.EMF().Ex().numx()-2*Nbc+c) << "\n";
        }
    }



    // Harmonics:x0 "Left-Bound ---> Right-Guard"
    for(int s(0); s < Y.Species(); ++s) {
        for(int i(0); i < Y.DF(s).dim(); ++i) {
            for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                for(int c(0); c < Nbc; c++) {
                    // std::cout << "\n 3: DF(" << i << "," << p << "," << Y.EMF().Ex().numx()-Nbc+c << ") = " << (Y.DF(s)(i))(p, Y.EMF().Ex().numx()-Nbc+c) << "\n";

                    (Y.DF(s)(i))(p, Y.EMF().Ex().numx()-Nbc+c) = (Y.DF(s)(i))(p, Nbc+c);

                    // std::cout << "\n 3n: DF(" << i << "," << p << "," << Y.EMF().Ex().numx()-Nbc+c << ") = " << (Y.DF(s)(i))(p, Y.EMF().Ex().numx()-Nbc+c-100.0) << "\n";
                    // if (c==2)


                }
            }
        }
    }

    // Fields:   x0 "Left-Bound ---> Right-Guard"
    for(int i(0); i < Y.EMF().dim(); ++i) {
        for(int c(0); c < Nbc; c++) {
            // if (i == 0)                 std::cout << "\n 4: FLD(" << i << "," << Y.EMF().Ex().numx()-Nbc+c << ") = " << Y.FLD(i)(Y.EMF().Ex().numx()-Nbc+c) << "\n";
            Y.FLD(i)(Y.EMF().Ex().numx()-Nbc+c) = Y.FLD(i)(Nbc+c);
            // if (i == 0)                 std::cout << "\n 4n: FLD(" << i << "," << Y.EMF().Ex().numx()-Nbc+c << ") = " << Y.FLD(i)(Y.EMF().Ex().numx()-Nbc+c) << "\n";


        }
        // Y.FLD(i)(Y.EMF().Ex().numx()-1) = 0.0;
    }

    if (Input::List().hydromotion)
    {
        // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().density(c) = Y.HYDRO().density(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().vx(c) = Y.HYDRO().vx(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().vy(c) = Y.HYDRO().vy(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().vz(c) = Y.HYDRO().vz(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().temperature(c) = Y.HYDRO().temperature(Y.EMF().Ex().numx()-2*Nbc+c);
            Y.HYDRO().Z(c) = Y.HYDRO().Z(Y.EMF().Ex().numx()-2*Nbc+c);
        }

        // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().density(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().density(Nbc+c);
            Y.HYDRO().vx(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vx(Nbc+c);
            Y.HYDRO().vy(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vy(Nbc+c);
            Y.HYDRO().vz(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vz(Nbc+c);
            Y.HYDRO().temperature(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().temperature(Nbc+c);
            Y.HYDRO().Z(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().Z(Nbc+c);
        }
    }

    if (Input::List().particlepusher)
    {
        for (int ip(0); ip < Y.particles().numpar(); ++ip)
        {
            if (Y.particles().x(ip) > Input::List().xmaxLocalnobnd[0])
            {
                Y.particles().x(ip) = Input::List().xminLocalnobnd[0] + (Y.particles().x(ip) - Input::List().xmaxLocalnobnd[0]); 
                Y.particles().ishere(ip) = 1;  
            }
            else if (Y.particles().x(ip) < Input::List().xminLocalnobnd[0])
            {
                Y.particles().x(ip) = Input::List().xmaxLocalnobnd[0] - (Input::List().xminLocalnobnd[0] - Y.particles().x(ip));   
                Y.particles().ishere(ip) = 1;
            }
        }
    }

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications_1D::sameNode_mirror_X(State1D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
    int sign(1);
    size_t Nx(Y.SH(0,0,0).numx());

    // Mirror the harmonics
    for(int s(0); s < Y.Species(); ++s) {
        for(int l(0); l < Y.DF(s).l0(); ++l){
            for(int m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                // right boundary
                for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                    for (int c(0); c < Nbc; ++c) {
                        Y.SH(s,l,m)(p, Nx-c-1) = Y.SH(s,l,m)(p, Nx-2*Nbc+c);
                        Y.SH(s,l,m)(p, Nx-c-1) *= sign;
                    }
                }

                // left boundary
                for(int p(0); p < Y.SH(s,0,0).nump(); ++p) {
                    for (int c(0); c < Nbc; ++c) {
                        Y.SH(s,l,m)(p, c) = Y.SH(s,l,m)(p, 2*Nbc-c-1);
                        Y.SH(s,l,m)(p, c) *= sign;
                    }
                }

            }
        }
    }

    // Mirror the fields
    for(int i(0); i < Y.EMF().dim(); ++i){
        // right boundary
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(Nx-c-1) = Y.FLD(i)(Nx-2*Nbc+c);
        // left boundary
        for (int c(0); c < Nbc; ++c)
            Y.FLD(i)(c) = Y.FLD(i)(2*Nbc-c-1);
    }

    for(int c(0); c < Nbc; c++) {
        //Ey
        Y.EMF().Ey()(Nx-Nbc+c) *= -1.0; // right boundary
        Y.EMF().Ey()(c) *= -1.0;  // left  boundary

        //Ez
        Y.EMF().Ez()(Nx-Nbc+c) *= -1.0; // right boundary
        Y.EMF().Ez()(c) *= -1.0;  // left  boundary

        //Bx
        Y.EMF().Bx()(Nx-Nbc+c) *= -1.0; // right boundary
        Y.EMF().Bx()(c) *= -1.0;  // left  boundary
    }


    if (Input::List().hydromotion)
    {
        // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().vx(c) *= -1.0;
        }

        // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().vx(Y.EMF().Ex().numx()-Nbc+c) *= -1.0;


        }
    }

    if (Input::List().particlepusher)
    {
        for (int ip(0); ip < Y.particles().numpar(); ++ip)
        {
            if (Y.particles().x(ip) > Input::List().xmaxLocalnobnd[0])
            {
                Y.particles().x(ip) = Input::List().xmaxLocalnobnd[0] - (Y.particles().x(ip) - Input::List().xmaxLocalnobnd[0]);
                Y.particles().px(ip) *= -1.0;
            }
            else if (Y.particles().x(ip) < Input::List().xminLocalnobnd[0])
            {
                Y.particles().x(ip) = Input::List().xminLocalnobnd[0] + (Input::List().xminLocalnobnd[0] - Y.particles().x(ip));
                Y.particles().px(ip) *= -1.0;
            }

        }
    }

}
//--------------------------------------------------------------


//**************************************************************
//**************************************************************
//   Definition of the Parallel Environment
//**************************************************************
//**************************************************************


//--------------------------------------------------------------
Parallel_Environment_1D:: Parallel_Environment_1D() :
//--------------------------------------------------------------
//  Constructor, domain decomposition
//--------------------------------------------------------------
        bndX(Input::List().bndX),           // Type of boundary
        MPI_Procs(Input::List().MPI_X[0])   // Number of nodes in X-direction
{
    // Determination of the rank and size of the run
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_Procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (error_check()) 
    {
        std::cout << "PE error check failed" << std::endl;
        MPI_Finalize(); exit(1);
    }

    // Determination of the local computational domain (i.e. the x-axis and the y-axis)
    for(size_t i(0); i < Input::List().xminLocal.size(); ++i) {

        Input::List().xminLocal[i] = Input::List().xminGlobal[i]
                                     + rank * Input::List().NxLocalnobnd[i] * Input::List().globdx[i]
                                     - Input::List().BoundaryCells * Input::List().globdx[i];
        Input::List().xmaxLocal[i] = Input::List().xminLocal[i]
                                     + (Input::List().NxLocal[i]) * Input::List().globdx[i];

        Input::List().xminLocalnobnd[i] = Input::List().xminLocal[i] + Input::List().BoundaryCells * Input::List().globdx[i];
        Input::List().xmaxLocalnobnd[i] = Input::List().xmaxLocal[i] - Input::List().BoundaryCells * Input::List().globdx[i];

    }

    double  xval_lastcell = Input::List().xmaxLocal[0] - 0.5*Input::List().globdx[0];
    double xval_firstcell = Input::List().xminLocal[0] + 0.5*Input::List().globdx[0];

    
    // Restart files will be generated when output_step % restart_step == 0
//          if ( ( (Input::List().n_outsteps + 1) >  Input::List().n_restarts ) &&
//                (Input::List().n_restarts > 0)  ) { 
//               restart_step = Input::List().n_outsteps / Input::List().n_restarts ;
//          } 
//          else restart_step = -1;
//          if ( Input::List().restart_time > Input::List().n_restarts ) {
//              restart_time = 0; 
//          }
//          else restart_time = Input::List().restart_time;

}
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
Parallel_Environment_1D:: ~Parallel_Environment_1D(){ }
//--------------------------------------------------------------

//--------------------------------------------------------------
bool Parallel_Environment_1D:: error_check() {
//--------------------------------------------------------------
//  Temporary error check
//--------------------------------------------------------------

    // Test the number of cells
    if (rank == 0){
        if (Input::List().NxLocal[0] < 2*Input::List().BoundaryCells+2){
            std::cout << "Not enough cells per processor" << endl;
            return true;
        }
        if ( MPI_Procs != (Input::List().MPI_X[0]) ) {
            std:: cout << "the number of nodes in the input deck is "
                       << (Input::List().MPI_X[0]) << ", terminating ..." << endl;
            return true;
        }
    }
    return false;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Basic information
//--------------------------------------------------------------
int Parallel_Environment_1D:: RANK()  const {return rank;}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int Parallel_Environment_1D:: MPI_Processes() const {return MPI_Procs;}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int Parallel_Environment_1D:: BNDX()  const {return bndX;}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Parallel_Environment_1D::Neighbor_ImplicitE_Communications(State1D& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------

    int moduloX(RANK()%2);
    int RNx((RANK()+1)%MPI_Processes()),         // This is the right neighbor
            LNx((RANK()-1+MPI_Processes())%MPI_Processes()); // This is the left  neighbor

    if (MPI_Processes() > 1) {
        //even nodes
        if (moduloX==0){
            Bfield_Data.Send_right_X(Y,RNx);                  //   (Send) 0 --> 1
            if ((RANK() != 0) || (BNDX()==0)){
                Bfield_Data.Recv_from_left_X(Y,LNx);          //          1 --> 0 (Receive)
                Bfield_Data.Send_left_X(Y,LNx);               //          1 <-- 0 (Send)
            }
            else {
                if (BNDX()==1) {
                    Bfield_Data.mirror_bound_Xleft(Y);        // Update node "0" in the x direction
                }
                else {
                    cout<<"Invalid Boundary." << endl;
                }
            }
            Bfield_Data.Recv_from_right_X(Y,RNx);               // (Receive) 0 <-- 1
        }
            //odd nodes 
        else {
            Bfield_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)
            if ((RANK()!=(MPI_Processes()-1)) || (BNDX()==0)){
                Bfield_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0
                Bfield_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
            }
            else {
                if (BNDX()==1) {
                    Bfield_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
                }
                else {
                    cout<<"Invalid Boundary." << endl;;
                }
            }
            Bfield_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)
        }
    }
    else { Bfield_Data.sameNode_bound_X(Y); }

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Parallel_Environment_1D::Neighbor_Communications(State1D& Y) {
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------

    int moduloX(RANK() % 2);
    int RNx((RANK() + 1) % MPI_Processes()),         // This is the right neighbor
            LNx((RANK() - 1 + MPI_Processes()) % MPI_Processes()); // This is the left  neighbor

    if (MPI_Processes() > 1) {
        //even nodes
//        if (moduloX==0){
        // std::cout << "\n before boundary on " << Input::List().xminLocalnobnd[0] << "? " << Y.particles().x(2) << "\n\n";  
        if (BNDX() == 0) {
            if (((RANK() != 0) && (RANK() != MPI_Processes() - 1))) {
                X_Data.Send_right_X(Y, RNx);                  //   (Send) 0 --> 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
                X_Data.Recv_from_right_X(Y, RNx);               // (Receive) 0 <-- 1
            } else if (RANK() == 0) {                         /// Update node "0" in the x direction
                X_Data.Recv_from_left_X(Y, MPI_Processes() - 1);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, MPI_Processes() - 1);               //          1 <-- 0 (Send)
                X_Data.Send_right_X(Y, RNx);                ///   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, RNx);           /// (Receive) 0 <-- 1
            } else if (RANK() == MPI_Processes() - 1) {               ///        // Update node "MPI_Processes()" in the x direction
                X_Data.Send_right_X(Y, 0);                  //   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, 0);           /// (Receive) 0 <-- 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
            }

            Y.particles().par_goingright_array() = 0.0;

        } 
        else if (BNDX() == 1) {
            if (((RANK() != 0) && (RANK() != MPI_Processes() - 1))) {
                X_Data.Send_right_X(Y, RNx);                  //   (Send) 0 --> 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
                X_Data.Recv_from_right_X(Y, RNx);               // (Receive) 0 <-- 1
            } else if (RANK() == 0) {                         /// Update node "0" in the x direction
                X_Data.mirror_bound_Xleft(Y);
                X_Data.Send_right_X(Y, RNx);                ///   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, RNx);           /// (Receive) 0 <-- 1
            } else if (RANK() == MPI_Processes() - 1) {               ///        // Update node "MPI_Processes()" in the x direction
                X_Data.mirror_bound_Xright(Y);
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
            }

        }
        else if (BNDX() == 2) {
            if (((RANK() != 0) && (RANK() != MPI_Processes() - 1))) {
                X_Data.Send_right_X(Y, RNx);                  //   (Send) 0 --> 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
                X_Data.Recv_from_right_X(Y, RNx);               // (Receive) 0 <-- 1
            } else if (RANK() == 0) {                         /// Update node "0" in the x direction
                X_Data.Recv_from_left_X(Y, MPI_Processes() - 1);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, MPI_Processes() - 1);               //          1 <-- 0 (Send)
                X_Data.Send_right_X(Y, RNx);                ///   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, RNx);           /// (Receive) 0 <-- 1
            } else if (RANK() == MPI_Processes() - 1) {               ///        // Update node "MPI_Processes()" in the x direction
                X_Data.Send_right_X(Y, 0);                  //   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, 0);           /// (Receive) 0 <-- 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
            }

        }
        else {
            cout << "Invalid Boundary." << endl;
        }

//            //odd nodes
//        else {
//            X_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)
//            if ((RANK()!=(MPI_Processes()-1)) || (BNDX()==0)){
//                X_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0
//                X_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
//            }
//            else {
//                if (BNDX()==1) {
//                    X_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
//                }
//                else {
//                    cout<<"Invalid Boundary." << endl;
//                }
//            }
//            X_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)
//        }
    } else { X_Data.sameNode_bound_X(Y); }

}
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
    Node_ImplicitE_Communications_2D:: Node_ImplicitE_Communications_2D() : 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Input::List().BoundaryCells),      // # of boundary cells
        bndX(Input::List().bndX),        // Type of boundary in X
        bndY(Input::List().bndY),
        Nx_local(Input::List().NxLocal[0]),
        Ny_local(Input::List().NxLocal[1]) {       // Type of boundary in X
     
        // 3 components for Bx, By, Bz
        msg_sizeX = 3;  
        msg_sizeY = 3;  
        
        msg_sizeX *= Nbc * Ny_local;   
        msg_sizeY *= Nbc * Nx_local;   
        msg_bufX = new complex<double>[msg_sizeX];
        msg_bufY = new complex<double>[msg_sizeY];

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Node_ImplicitE_Communications_2D:: ~Node_ImplicitE_Communications_2D(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
        delete[] msg_bufX;
        delete[] msg_bufY;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int Node_ImplicitE_Communications_2D:: BNDX()  const {return bndX;} 
//--------------------------------------------------------------
//--------------------------------------------------------------
    int Node_ImplicitE_Communications_2D:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Send_right_X(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);

        // Fields:   x0 "Right-Bound --> "
        
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                    msg_bufX[bufind + e] = Y.FLD(i)(Nx_local-2*Nbc+e,iy);
                } 
                bufind += step_f;
            }
        } 

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Recv_from_left_X(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Fields:   x0-"---> Left-Guard"
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(e,iy) = msg_bufX[bufind + e];
                }     
                bufind += step_f;
            }
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Send_left_X(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0); 

        // Fields:   x0 " <--- Left-Bound "
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   msg_bufX[bufind + e] = Y.FLD(i)(Nbc+e,iy);
                } 
                
                bufind += step_f;
            }
        }

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Recv_from_right_X(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Fields:   x0-"Right-Guard <--- "
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(Nx_local-Nbc+e,iy) = msg_bufX[bufind + e];
                } 
                bufind += step_f;
            }
        }
    }
//--------------------------------------------------------------
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the Y direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Send_right_Y(State2D& Y, int dest) {
//--------------------------------------------------------------
//  Y-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);

        // Fields:   x0 "Right-Bound --> "
        
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t e(0); e < Nbc; e++) {
                    msg_bufY[bufind + e] = Y.FLD(i)(ix,Ny_local-2*Nbc+e);
                } 
                bufind += step_f;
            }
        } 

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Recv_from_left_Y(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Fields:   x0-"---> Left-Guard"
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(ix,e) = msg_bufY[bufind + e];
                }     
                bufind += step_f;
            }
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Send_left_Y(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0); 

        // Fields:   x0 " <--- Left-Bound "
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t e(0); e < Nbc; e++) {
                   msg_bufY[bufind + e] = Y.FLD(i)(ix,Nbc+e);
                } 
                
                bufind += step_f;
            }
        }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::Recv_from_right_Y(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Fields:   x0-"Right-Guard <--- "
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(ix,Ny_local-Nbc+e) = msg_bufY[bufind + e];
                } 
                bufind += step_f;
            }
        }
    }
//--------------------------------------------------------------


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//*************************************************************************************
//-------------------------------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::mirror_bound_Xleft(State2D& Y) {
//-------------------------------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//-------------------------------------------------------------------------------------



        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the x cells                
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(c,iy) = Y.FLD(i)(2*Nbc-c-1,iy);
                }
            }
        }

        //Bx
        for(size_t iy(0); iy < Ny_local; ++iy){  // All the x cells                
            for (size_t c(0); c < Nbc; ++c) {
               Y.EMF().Bx()(c,iy) *= -1.0; // left  boundary
            }
        }

    }
//-------------------------------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::mirror_bound_Xright(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------



        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the x cells
                for (size_t c(0); c < Nbc; ++c){ 
                    Y.FLD(i)(Nx_local-c-1,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
                }       
            }
        }

        //Bx
        for(size_t iy(0); iy < Ny_local; ++iy){  // All the x cells                
            for (size_t c(0); c < Nbc; ++c) {
                Y.EMF().Bx()(Nx_local-1-c,iy) *= -1.0; // right  boundary
            }
        }

    }
//--------------------------------------------------------------

//-------------------------------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::mirror_bound_Yleft(State2D& Y) {
//-------------------------------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//-------------------------------------------------------------------------------------

        

        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(ix,c) = Y.FLD(i)(ix,2*Nbc-c-1);
                }
            }
        }

        //Bx
        for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
            for (size_t c(0); c < Nbc; ++c) {
                Y.EMF().By()(ix,c) *= -1.0; // left  boundary
            }
        }

    }
//-------------------------------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::mirror_bound_Yright(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------



        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(ix,Ny_local-c-1) = Y.FLD(i)(ix,Ny_local-2*Nbc+c);
                }
            }
        }

        //Bx
        for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
            for (size_t c(0); c < Nbc; ++c) {
                Y.EMF().By()(ix,Ny_local-1-c) *= -1.0; // right  boundary
            }
        }

    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_bound_X(State2D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
        switch (BNDX()) {
            case 0:                   // periodic
                sameNode_periodic_X(Y);

                break;
            case 1:                   // mirror boundary
                sameNode_mirror_X(Y);
                break;
            default:
                cout<<"Not a valid boundary condition." << endl;
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_periodic_X(State2D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------

        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(c,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
                }
            }
        }

        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(Nx_local-Nbc+c,iy) = Y.FLD(i)(Nbc+c,iy);
                }
            }
        }
       
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_mirror_X(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
 


        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
               // right boundary
               for (size_t c(0); c < Nbc; ++c) 
                   Y.FLD(i)(Nx_local-c-1,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
               // left boundary
               for (size_t c(0); c < Nbc; ++c) 
                  Y.FLD(i)(c,iy) = Y.FLD(i)(2*Nbc-c-1,iy);
          }
        }

        //Bx
        for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
            for (size_t c(0); c < Nbc; ++c) {
                Y.EMF().Bx()(c,iy) *= -1.0; // left  boundary
                Y.EMF().Bx()(Nx_local-Nbc+c,iy) *= -1.0; // right  boundary
            }
        }

    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_bound_Y(State2D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
        switch (BNDY()) {
            case 0:                   // periodic
                sameNode_periodic_Y(Y);

                break;
            case 1:                   // mirror boundary
                sameNode_mirror_Y(Y);
                break;
            default:
                cout<<"Not a valid boundary condition." << endl;
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_periodic_Y(State2D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------

        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for (size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(ix,c) = Y.FLD(i)(ix,Ny_local-2*Nbc+c);
                }
            }
        }

        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for (size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(ix,Ny_local-Nbc+c) = Y.FLD(i)(ix,Nbc+c);
                }
            }
        }
       
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications_2D::sameNode_mirror_Y(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
 


        // Mirror the fields 
        for (size_t i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
               // right boundary
               for (size_t c(0); c < Nbc; ++c) 
                   Y.FLD(i)(ix,Ny_local-c-1) = Y.FLD(i)(ix,Ny_local-2*Nbc+c);
               // left boundary
               for (size_t c(0); c < Nbc; ++c) 
                  Y.FLD(i)(ix,c) = Y.FLD(i)(ix,2*Nbc-c-1);
          }
        }

        //Bx
        for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
            for (size_t c(0); c < Nbc; ++c) {
                Y.EMF().By()(ix,c) *= -1.0; // left  boundary
                Y.EMF().By()(ix,Ny_local-Nbc+c) *= -1.0; // right  boundary
            }
        }

    }
//--------------------------------------------------------------

//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Node_Communications_2D:: Node_Communications_2D() : 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Input::List().BoundaryCells),      // # of boundary cells
        bndX(Input::List().bndX),
        bndY(Input::List().bndY),
        Nx_local(Input::List().NxLocal[0]),
        Ny_local(Input::List().NxLocal[1]){        // Type of boundary in X
         
        numspec = Input::List().ls.size();
        // numpmax = Input::List().ps;

        double temp;

        // # of harmonics
        msg_sizeX = 0; msg_sizeY = 0; temp = 0;
        for (size_t s(0); s < numspec; ++s) {
            temp = ((Input::List().ms[s]+1)*(2*Input::List().ls[s]-Input::List().ms[s]+2))/2;
            msg_sizeX += temp*(Input::List().dp[s]).size();
        }
        // (# of harmonics) * (# cells in p)
        // msg_sizeX *= numpmax
        // 6 fields: Ex, Ey, Ez, Bx, By, Bz
        msg_sizeX += 6;  
        msg_sizeY  = msg_sizeX;  

        // 5 Hydro Quantities, density, vel, temp
        // if (Input::List().hydromotion)
        // {
        //     msg_sizeX += 3;
        // }
        
        msg_sizeX *= Nbc*Input::List().NxLocal[1];   
        msg_bufX = new complex<double>[msg_sizeX];

        msg_sizeY *= Nbc*Input::List().NxLocal[0];   
        msg_bufY = new complex<double>[msg_sizeY];

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Node_Communications_2D:: ~Node_Communications_2D(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
        delete[] msg_bufX;
        delete[] msg_bufY;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int Node_Communications_2D:: BNDX()  const {return bndX;} 
//--------------------------------------------------------------
//--------------------------------------------------------------
    int Node_Communications_2D:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications_2D::Send_right_X(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

        // static size_t step_h(numpmax*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);

        // Harmonics:x0 "Right-Bound ---> " 
        for (size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                    for (size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for (size_t e(0); e < Nbc; e++) {
                          msg_bufX[bufind + e] = (Y.DF(s)(i))(p, Nx_local-2*Nbc+e,iy);
                        }
                        bufind += step_f;
                    }
                } 
            }
        }
        // Fields:   x0 "Right-Bound --> "
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t e(0); e < Nbc; e++) {
                    msg_bufX[bufind + e] = Y.FLD(i)(Nx_local-2*Nbc+e,iy);
                }
               bufind += step_f;
            } 
         } 
        
      //   if (Input::List().hydromotion)
      //   {
      //       // Hydro-velocity:   x0 "Right-Bound --> "
      //       for (size_t e(0); e < Nbc; e++) {
      //           msg_bufX[bufind + e] = Y.HYDRO().density(Nx_local-2*Nbc+e);
      //       } 
      //       bufind += step_f;
      //       // Hydro-velocity:   x0 "Right-Bound --> "
            // for (size_t e(0); e < Nbc; e++) {
            //  msg_bufX[bufind + e] = Y.HYDRO().velocity(Nx_local-2*Nbc+e);
            // } 
            // bufind += step_f;

      //       // Hydro-temperature:   x0 "Right-Bound --> "
            // for (size_t e(0); e < Nbc; e++) {
            //  msg_bufX[bufind + e] = Y.HYDRO().temperature(Nx_local-2*Nbc+e);
            // } 
            // bufind += step_f;
        // }
        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Recv_from_left_X(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        for (size_t s(0); s < Y.Species(); ++s) {
            for (size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                    for (size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for (size_t e(0); e < Nbc; e++) {
                            (Y.DF(s)(i))(p, e, iy) = msg_bufX[bufind + e];
                        }
                        bufind += step_f;
                    }
                }
            } 
        }


        // Fields:   x0-"---> Left-Guard"
        for (size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(e,iy) = msg_bufX[bufind + e];
                } 
                bufind += step_f;
            }            
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro-velocity:   x0-"---> Left-Guard"
        //     for (size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().density(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;        

        //     // Hydro-velocity:   x0-"---> Left-Guard"
        //     for (size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().velocity(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;        
            
        //     // Hydro-temperature:   x0-"---> Left-Guard"        
        //     for (size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().temperature(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;   
        // }     
        
        
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Send_left_X(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        for (size_t s(0); s < Y.Species(); ++s) {
            for (size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                    for (size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for (size_t e(0); e < Nbc; e++) {
                          msg_bufX[bufind + e] = (Y.DF(s)(i))(p, Nbc+e,iy);
                        }
                        bufind += step_f;
                    }
                }
            } 
        } 
        // Fields:   x0 " <--- Left-Bound "
        for (size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t e(0); e < Nbc; e++) {
                    msg_bufX[bufind + e] = Y.FLD(i)(Nbc+e,iy);
                } 
                bufind += step_f;
            }
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro-velocity:   x0 " <--- Left-Bound "
        //     for (size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().density(Nbc+e);
        //     } 
        //     bufind += step_f;
        //     // Hydro-velocity:   x0 " <--- Left-Bound "
        //     for (size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().velocity(Nbc+e);
        //     } 
        //     bufind += step_f;

        //     // Hydro-temperature:   x0 " <--- Left-Bound "
        //     for (size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().temperature(Nbc+e);
        //     } 
        //     bufind += step_f;
        // }

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Recv_from_right_X(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t e(0); e < Nbc; e++) {
                            (Y.DF(s)(i))(p, Nx_local-Nbc+e,iy) = msg_bufX[bufind + e];
                        }
                        bufind += step_f;
                    }
                }
            } 
        }
        // Fields:   x0-"Right-Guard <--- "
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(Nx_local-Nbc+e,iy) = msg_bufX[bufind + e];
                } 
                bufind += step_f;
            }
        }

      //   if (Input::List().hydromotion)
      //   {
      //       // Hydro-density:   x0-"Right-Guard <--- "
      //       for(size_t e(0); e < Nbc; e++) {
      //           Y.HYDRO().density(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
      //       } 
      //       bufind += step_f;

      //       // Hydro-velocity:   x0-"Right-Guard <--- "
            // for(size_t e(0); e < Nbc; e++) {
            //  Y.HYDRO().velocity(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
            // } 
            // bufind += step_f;
      //       // Hydro-temperature:   x0-"Right-Guard <--- "
            // for(size_t e(0); e < Nbc; e++) {
            //  Y.HYDRO().temperature(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
            // } 
            // bufind += step_f;
      //   }
        
    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the Y direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications_2D::Send_right_Y(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------

        // static size_t step_h(numpmax*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);

        // Harmonics:x0 "Right-Bound ---> " 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t e(0); e < Nbc; e++) {
                          msg_bufY[bufind + e] = (Y.DF(s)(i))(p, ix, Ny_local-2*Nbc+e);
                        }
                        bufind += step_f;
                    }
                } 
            }
        }
        // Fields:   x0 "Right-Bound --> "
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                    msg_bufY[bufind + e] = Y.FLD(i)(ix,Ny_local-2*Nbc+e);
                }
               bufind += step_f;
            } 
         } 
        
      //   if (Input::List().hydromotion)
      //   {
      //       // Hydro-velocity:   x0 "Right-Bound --> "
      //       for(size_t e(0); e < Nbc; e++) {
      //           msg_bufX[bufind + e] = Y.HYDRO().density(Nx_local-2*Nbc+e);
      //       } 
      //       bufind += step_f;
      //       // Hydro-velocity:   x0 "Right-Bound --> "
            // for(size_t e(0); e < Nbc; e++) {
            //  msg_bufX[bufind + e] = Y.HYDRO().velocity(Nx_local-2*Nbc+e);
            // } 
            // bufind += step_f;

      //       // Hydro-temperature:   x0 "Right-Bound --> "
            // for(size_t e(0); e < Nbc; e++) {
            //  msg_bufX[bufind + e] = Y.HYDRO().temperature(Nx_local-2*Nbc+e);
            // } 
            // bufind += step_f;
        // }
        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Recv_from_left_Y(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t e(0); e < Nbc; e++) {
                            (Y.DF(s)(i))(p, ix, e) = msg_bufY[bufind + e];
                        }
                        bufind += step_f;
                    }
                }
            } 
        }


        // Fields:   x0-"---> Left-Guard"
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(ix, e) = msg_bufY[bufind + e];
                } 
                bufind += step_f;
            }            
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro-velocity:   x0-"---> Left-Guard"
        //     for(size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().density(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;        

        //     // Hydro-velocity:   x0-"---> Left-Guard"
        //     for(size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().velocity(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;        
            
        //     // Hydro-temperature:   x0-"---> Left-Guard"        
        //     for(size_t e(0); e < Nbc; e++) {
        //         Y.HYDRO().temperature(e) = (msg_bufX[bufind + e]).real();
        //     } 
        //     bufind += step_f;   
        // }     
        
        
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Send_left_Y(State2D& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t e(0); e < Nbc; e++) {
                          msg_bufY[bufind + e] = (Y.DF(s)(i))(p, ix, Nbc+e);
                        }
                        bufind += step_f;
                    }
                }
            } 
        } 
        // Fields:   x0 " <--- Left-Bound "
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                    msg_bufY[bufind + e] = Y.FLD(i)(ix,Nbc+e);
                } 
                bufind += step_f;
            }
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro-velocity:   x0 " <--- Left-Bound "
        //     for(size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().density(Nbc+e);
        //     } 
        //     bufind += step_f;
        //     // Hydro-velocity:   x0 " <--- Left-Bound "
        //     for(size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().velocity(Nbc+e);
        //     } 
        //     bufind += step_f;

        //     // Hydro-temperature:   x0 " <--- Left-Bound "
        //     for(size_t e(0); e < Nbc; e++) {
        //         msg_bufX[bufind + e] = Y.HYDRO().temperature(Nbc+e);
        //     } 
        //     bufind += step_f;
        // }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::Recv_from_right_Y(State2D& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------

        // static size_t step_h(Y.SH(s,0,0).nump()*Nbc);
        static size_t step_f(Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i){
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t e(0); e < Nbc; e++) {
                            (Y.DF(s)(i))(p, ix, Ny_local-Nbc+e) = msg_bufY[bufind + e];
                        }
                        bufind += step_f;
                    }
                }
            } 
        }
        // Fields:   x0-"Right-Guard <--- "
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the y cells
                for(size_t e(0); e < Nbc; e++) {
                   Y.FLD(i)(ix,Ny_local-Nbc+e) = msg_bufY[bufind + e];
                } 
                bufind += step_f;
            }
        }

      //   if (Input::List().hydromotion)
      //   {
      //       // Hydro-density:   x0-"Right-Guard <--- "
      //       for(size_t e(0); e < Nbc; e++) {
      //           Y.HYDRO().density(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
      //       } 
      //       bufind += step_f;

      //       // Hydro-velocity:   x0-"Right-Guard <--- "
            // for(size_t e(0); e < Nbc; e++) {
            //  Y.HYDRO().velocity(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
            // } 
            // bufind += step_f;
      //       // Hydro-temperature:   x0-"Right-Guard <--- "
            // for(size_t e(0); e < Nbc; e++) {
            //  Y.HYDRO().temperature(Nx_local-Nbc+e) = (msg_bufX[bufind + e]).real();
            // } 
            // bufind += step_f;
      //   }
        
    }
//--------------------------------------------------------------

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//**************************************************************
//--------------------------------------------------------------
    void Node_Communications_2D::mirror_bound_Xleft(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//--------------------------------------------------------------
        int sign(1);
        

        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)
                
                        for (size_t c(0); c < Nbc; ++c) {
                            for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                                Y.SH(s,l,m)(p, c, iy) = Y.SH(s,l,m)(p, 2*Nbc-c-1,iy);
                                Y.SH(s,l,m)(p, c, iy) *= sign;
                            }
                        }

                    } 
                }
            }
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(c,iy) = Y.FLD(i)(2*Nbc-c-1,iy);
                }
            }
        }

        for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
            for (size_t c(0); c < Nbc; ++c) {
//          Ey
                 Y.EMF().Ey()(c,iy) *= -1.0; // left  boundary
//          Ez
                 Y.EMF().Ez()(c,iy) *= -1.0; // left  boundary
//          Bx
                 Y.EMF().Bx()(c,iy) *= -1.0; // left  boundary
            }
        }

        // if (Input::List().hydromotion)
        // {
        //   // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(c) *= -1.0;
        //     }    
        // }


    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::mirror_bound_Xright(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------
        int sign(1);
        

        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, Nx_local-c-1, iy) = Y.SH(s,l,m)(p, Nx_local-2*Nbc+c, iy);
                                Y.SH(s,l,m)(p, Nx_local-c-1, iy) *= sign;
                            }
                        }
     
                    }
                }
            }
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
           for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(Nx_local-c-1,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
                }
            }
        }
        for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells   
            for (size_t c(0); c < Nbc; ++c) {
                //Ey
                Y.EMF().Ey()(Nx_local-c-1,iy) *= -1.0; // right boundary 
                //Ez
                Y.EMF().Ez()(Nx_local-c-1,iy) *= -1.0; // right boundary 
                //Bx
                Y.EMF().Bx()(Nx_local-c-1,iy) *= -1.0; // right boundary 
            }
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(Nx_local-Nbc+c) *=  -1.0;
        //     }
        // }

    }
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
    void Node_Communications_2D::mirror_bound_Yleft(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the y direction on the left
//--------------------------------------------------------------
        int sign(1);


        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)
                
                        for (size_t c(0); c < Nbc; ++c) {
                            for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                                Y.SH(s,l,m)(p, ix, c) = Y.SH(s,l,m)(p, ix, 2*Nbc-c-1);
                                Y.SH(s,l,m)(p, ix, c) *= sign;
                            }
                        }

                    } 
                }
            }
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(ix,c) = Y.FLD(i)(ix, 2*Nbc-c-1);
                }
            }
        }

        for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
            for (size_t c(0); c < Nbc; ++c) {
//          Ey
                 Y.EMF().Ey()(ix,c) *= -1.0; // left  boundary
//          Ez
                 Y.EMF().Ez()(ix,c) *= -1.0; // left  boundary
//          Bx
                 Y.EMF().Bx()(ix,c) *= -1.0; // left  boundary
            }
        }

        // if (Input::List().hydromotion)
        // {
        //   // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(c) *= -1.0;
        //     }    
        // }


    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::mirror_bound_Yright(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------
        int sign(1);


        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, ix, Ny_local-c-1) = Y.SH(s,l,m)(p, ix, Ny_local-2*Nbc+c);
                                Y.SH(s,l,m)(p, ix, Ny_local-c-1) *= sign;
                            }
                        }
     
                    }
                }
            }
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
           for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for (size_t c(0); c < Nbc; ++c) {
                    Y.FLD(i)(ix,Ny_local-c-1) = Y.FLD(i)(ix,Ny_local-2*Nbc+c);
                }
            }
        }

        for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
            for (size_t c(0); c < Nbc; ++c) {
                //Ey
                Y.EMF().Ey()(ix,Ny_local-c-1) *= -1.0; // right boundary 
                //Ez
                Y.EMF().Ez()(ix,Ny_local-c-1) *= -1.0; // right boundary 
                //Bx
                Y.EMF().Bx()(ix,Ny_local-c-1) *= -1.0; // right boundary 
            }
        }

        // if (Input::List().hydromotion)
        // {
        //     // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        //     for (size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(Nx_local-Nbc+c) *=  -1.0;
        //     }
        // }

    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_bound_X(State2D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
    
        switch (BNDX()) {
            case 0:                   // periodic

                sameNode_periodic_X(Y);

                break;
            case 1:                   // mirror boundary
                sameNode_mirror_X(Y);
                break;
            default:
                cout<<"Not a valid boundary condition." << endl;
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_periodic_X(State2D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------
        // valarray<double> temp(0.0,Y.SH(0,0,0).numx());
  
       

       // temp = Y.DF(0).getcurrent(0);
       // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    //     {        
    //         std::cout<<"\n original[" << ix  << "] = " << temp[ix] << "\n";      
            
            
    //     }



        // Harmonics:x0 "Right-Bound ---> Left-Guard" 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i) {
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                    for(size_t c(0); c < Nbc; c++) {
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            (Y.DF(s)(i))(p, c, iy) = (Y.DF(s)(i))(p, Nx_local-2*Nbc+c, iy);
                            // std::cout << "\n 1: DF(" << i << "," << p << "," << c << ") = " << (Y.DF(s)(i))(p, Nx_local-2*Nbc+c) << "\n";
                        }
                    }
                }
            }
        }
        
        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for(size_t i(0); i < Y.EMF().dim(); ++i) {
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(c,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
                    // std::cout << "\n 2: FLD(" << i << "," << c << ") = " << Y.FLD(i)(Nx_local-2*Nbc+c) << "\n";
                }
            }
        }
        


        // Harmonics:x0 "Left-Bound ---> Right-Guard" 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i) {
                for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells                    
                    for(size_t c(0); c < Nbc; c++) {
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            // std::cout << "\n 3: DF(" << i << "," << p << "," << Nx_local-Nbc+c << ") = " << (Y.DF(s)(i))(p, Nx_local-Nbc+c) << "\n";

                            (Y.DF(s)(i))(p, Nx_local-Nbc+c,iy) = (Y.DF(s)(i))(p, Nbc+c,iy);

                            // std::cout << "\n 3n: DF(" << i << "," << p << "," << Nx_local-Nbc+c << ") = " << (Y.DF(s)(i))(p, Nx_local-Nbc+c-100.0) << "\n";
                            // if (c==2)

                            
                        }
                    }
                }
            }
        }
        
        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for(size_t i(0); i < Y.EMF().dim(); ++i) {
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells
                for(size_t c(0); c < Nbc; c++) {
                    // if (i == 0)                 std::cout << "\n 4: FLD(" << i << "," << Nx_local-Nbc+c << ") = " << Y.FLD(i)(Nx_local-Nbc+c) << "\n";
                    Y.FLD(i)(Nx_local-Nbc+c,iy) = Y.FLD(i)(Nbc+c,iy);
                    // if (i == 0)                 std::cout << "\n 4n: FLD(" << i << "," << Nx_local-Nbc+c << ") = " << Y.FLD(i)(Nx_local-Nbc+c) << "\n";
                }
            }
            // Y.FLD(i)(Nx_local-1) = 0.0;
        }
        

        // temp = Y.DF(0).getcurrent(0);
        // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
        // {        
        //     std::cout<<"\n      new[" << ix  << "] = " << temp[ix] << "\n";      
            
        // }


   //      if (Input::List().hydromotion)
   //      {
   //          // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
            // for(size_t c(0); c < Nbc; c++) {
            //  Y.HYDRO().density(c) = Y.HYDRO().density(Nx_local-2*Nbc+c);
   //              Y.HYDRO().velocity(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            //  Y.HYDRO().temperature(c) = Y.HYDRO().temperature(Nx_local-2*Nbc+c);
            //  // Y.HYDRO().kpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            //  // Y.HYDRO().mpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            // }    
    
   //          // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
            // for(size_t c(0); c < Nbc; c++) {
            //  Y.HYDRO().density(Nx_local-Nbc+c) =  Y.HYDRO().density(Nbc+c);
   //              Y.HYDRO().velocity(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
            //  Y.HYDRO().temperature(Nx_local-Nbc+c) =  Y.HYDRO().temperature(Nbc+c);
            //  // Y.HYDRO().kpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
            //  // Y.HYDRO().mpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);

            // }
   //      }

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_mirror_X(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
        int sign(1);


        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells                
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                        // right boundary
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, Nx_local-c-1, iy) = Y.SH(s,l,m)(p, Nx_local-2*Nbc+c, iy);
                                Y.SH(s,l,m)(p, Nx_local-c-1, iy) *= sign;
                            }
                        }

                        // left boundary
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, c, iy) = Y.SH(s,l,m)(p, 2*Nbc-c-1, iy);
                                Y.SH(s,l,m)(p, c, iy) *= sign;
                            }
                        }

                    }
                }
            }           
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells    
               // right boundary
               for (size_t c(0); c < Nbc; ++c) 
                   Y.FLD(i)(Nx_local-c-1,iy) = Y.FLD(i)(Nx_local-2*Nbc+c,iy);
               // left boundary
               for (size_t c(0); c < Nbc; ++c) 
                  Y.FLD(i)(c,iy) = Y.FLD(i)(2*Nbc-c-1,iy);
            }
        }

        for(size_t iy(0); iy < Ny_local; ++iy){  // All the y cells    
            for(size_t c(0); c < Nbc; c++) {
                //Ey
                Y.EMF().Ey()(Nx_local-Nbc+c,iy) *= -1.0; // right boundary 
                Y.EMF().Ey()(c,iy) *= -1.0;  // left  boundary

                //Ez
                Y.EMF().Ez()(Nx_local-Nbc+c,iy) *= -1.0; // right boundary 
                Y.EMF().Ez()(c,iy) *= -1.0;  // left  boundary
           
                //Bx
                Y.EMF().Bx()(Nx_local-Nbc+c,iy) *= -1.0; // right boundary 
                Y.EMF().Bx()(c,iy) *= -1.0;  // left  boundary
            }
        }


        // if (Input::List().hydromotion)
        // {
        //     // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(c) *= -1.0;
        //         // Y.HYDRO().temperature(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //         // Y.HYDRO().kpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //         // Y.HYDRO().mpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //     }    
    
        //     // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(Nx_local-Nbc+c) *= -1.0;// Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().temperature(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().kpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().mpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);

        //     }
        // }

    }
//--------------------------------------------------------------
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_bound_Y(State2D& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the y direction
//--------------------------------------------------------------
    
        switch (BNDY()) {
            case 0:                   // periodic
                sameNode_periodic_Y(Y);

                break;
            case 1:                   // mirror boundary
                sameNode_mirror_Y(Y);
                break;
            default:
                cout<<"Not a valid boundary condition." << endl;
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_periodic_Y(State2D& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------
        // valarray<double> temp(0.0,Y.SH(0,0,0).numx());
  
       

       // temp = Y.DF(0).getcurrent(0);
       // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
       //  {        
       //      std::cout<<"\n original[" << ix  << "] = " << temp[ix] << "\n";      
            
            
       //  }

        // Harmonics:x0 "Right-Bound ---> Left-Guard" 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i) {
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                    for(size_t c(0); c < Nbc; c++) {
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            
                            (Y.DF(s)(i))(p, ix, c) = (Y.DF(s)(i))(p, ix, Ny_local-2*Nbc+c);
                            // std::cout << "\n 1: DF(" << i << "," << p << "," << c << ") = " << (Y.DF(s)(i))(p, Nx_local-2*Nbc+c) << "\n";
                        }
                    }
                }
            }
        }
        
        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for(size_t i(0); i < Y.EMF().dim(); ++i) {
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t c(0); c < Nbc; c++) {
                    Y.FLD(i)(ix, c) = Y.FLD(i)(ix, Ny_local-2*Nbc+c);
                    // std::cout << "\n 2: FLD(" << i << "," << c << ") = " << Y.FLD(i)(Nx_local-2*Nbc+c) << "\n";
                }
            }
        }
        // std::cout << "13 \n";


        // Harmonics:x0 "Left-Bound ---> Right-Guard" 
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t i(0); i < Y.DF(s).dim(); ++i) {
                for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                    for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                        for(size_t c(0); c < Nbc; c++) {
                            // std::cout << "\n 3: DF(" << i << "," << p << "," << Nx_local-Nbc+c << ") = " << (Y.DF(s)(i))(p, Nx_local-Nbc+c) << "\n";

                            (Y.DF(s)(i))(p, ix, Ny_local-Nbc+c) = (Y.DF(s)(i))(p, ix, Nbc+c);

                            // std::cout << "\n 3n: DF(" << i << "," << p << "," << Nx_local-Nbc+c << ") = " << (Y.DF(s)(i))(p, Nx_local-Nbc+c-100.0) << "\n";
                            // if (c==2)

                            
                        }
                    }
                }
            }
        }
        // std::cout << "14 \n";
        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for(size_t i(0); i < Y.EMF().dim(); ++i) {
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t c(0); c < Nbc; c++) {
                    // if (i == 0)                 std::cout << "\n 4: FLD(" << i << "," << Nx_local-Nbc+c << ") = " << Y.FLD(i)(Nx_local-Nbc+c) << "\n";
                    Y.FLD(i)(ix, Ny_local-Nbc+c) = Y.FLD(i)(ix, Nbc+c);
                    // if (i == 0)                 std::cout << "\n 4n: FLD(" << i << "," << Nx_local-Nbc+c << ") = " << Y.FLD(i)(Nx_local-Nbc+c) << "\n";
                }
            }
            // Y.FLD(i)(Nx_local-1) = 0.0;
        }
        // std::cout << "15 \n";
        // temp = Y.DF(0).getcurrent(0);
        // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
        // {        
        //     std::cout<<"\n      new[" << ix  << "] = " << temp[ix] << "\n";      
            
        // }


   //      if (Input::List().hydromotion)
   //      {
   //          // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
            // for(size_t c(0); c < Nbc; c++) {
            //  Y.HYDRO().density(c) = Y.HYDRO().density(Nx_local-2*Nbc+c);
   //              Y.HYDRO().velocity(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            //  Y.HYDRO().temperature(c) = Y.HYDRO().temperature(Nx_local-2*Nbc+c);
            //  // Y.HYDRO().kpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            //  // Y.HYDRO().mpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
            // }    
    
   //          // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
            // for(size_t c(0); c < Nbc; c++) {
            //  Y.HYDRO().density(Nx_local-Nbc+c) =  Y.HYDRO().density(Nbc+c);
   //              Y.HYDRO().velocity(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
            //  Y.HYDRO().temperature(Nx_local-Nbc+c) =  Y.HYDRO().temperature(Nbc+c);
            //  // Y.HYDRO().kpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
            //  // Y.HYDRO().mpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);

            // }
   //      }

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications_2D::sameNode_mirror_Y(State2D& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
        int sign(1);
        

        // Mirror the harmonics
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
                for(size_t l(0); l < Y.DF(s).l0(); ++l){
                    for(size_t m(0); m < ((Y.DF(s).m0() < l)? Y.DF(s).m0():l)+1; ++m){
                        sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                        // right boundary
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, ix, Ny_local-c-1) = Y.SH(s,l,m)(p, ix, Ny_local-2*Nbc+c);
                                Y.SH(s,l,m)(p, ix, Ny_local-c-1) *= sign;
                            }
                        }

                        // left boundary
                        for(size_t p(0); p < Y.SH(s,0,0).nump(); ++p) {
                            for (size_t c(0); c < Nbc; ++c) {
                                Y.SH(s,l,m)(p, ix, c) = Y.SH(s,l,m)(p, ix, 2*Nbc-c-1);
                                Y.SH(s,l,m)(p, ix, c) *= sign;
                            }
                        }

                    }
                }
            }           
        }

        // Mirror the fields 
        for(size_t i(0); i < Y.EMF().dim(); ++i){
            for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
               // right boundary
               for (size_t c(0); c < Nbc; ++c) 
                   Y.FLD(i)(ix, Ny_local-c-1) = Y.FLD(i)(ix, Ny_local-2*Nbc+c);
               // left boundary
               for (size_t c(0); c < Nbc; ++c) 
                  Y.FLD(i)(ix, c) = Y.FLD(i)(ix, 2*Nbc-c-1);
            }
        }

        for(size_t ix(0); ix < Nx_local; ++ix){  // All the x cells
            for(size_t c(0); c < Nbc; c++) {
                //Ey
                Y.EMF().Ey()(ix, Ny_local-Nbc+c) *= -1.0; // right boundary 
                Y.EMF().Ey()(ix, c) *= -1.0;  // left  boundary

                //Ez
                Y.EMF().Ez()(ix, Ny_local-Nbc+c) *= -1.0; // right boundary 
                Y.EMF().Ez()(ix, c) *= -1.0;  // left  boundary
           
                //Bx
                Y.EMF().Bx()(ix, Ny_local-Nbc+c) *= -1.0; // right boundary 
                Y.EMF().Bx()(ix, c) *= -1.0;  // left  boundary
            }
        }


        // if (Input::List().hydromotion)
        // {
        //     // Hydro Quantities:   x0 "Right-Bound ---> Left-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(c) *= -1.0;
        //         // Y.HYDRO().temperature(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //         // Y.HYDRO().kpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //         // Y.HYDRO().mpressure(c) = Y.HYDRO().velocity(Nx_local-2*Nbc+c);
        //     }    
    
        //     // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        //     for(size_t c(0); c < Nbc; c++) {
        //         Y.HYDRO().velocity(Nx_local-Nbc+c) *= -1.0;// Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().temperature(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().kpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
        //         // Y.HYDRO().mpressure(Nx_local-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);

        //     }
        // }

    }
//--------------------------------------------------------------

//**************************************************************
//**************************************************************
//   Definition of the Parallel Environment
//**************************************************************
//**************************************************************


//--------------------------------------------------------------
    Parallel_Environment_2D:: Parallel_Environment_2D() : 
//--------------------------------------------------------------
//  Constructor, domain decomposition
//--------------------------------------------------------------
        bndX(Input::List().bndX),           // Type of boundary
        bndY(Input::List().bndY),           // Type of boundary
        
        MPI_Processes_X(Input::List().MPI_X[0]),   // Number of processes in X-direction
        MPI_Processes_Y(Input::List().MPI_X[1]),   // Number of processes in Y-direction
        MPI_Procs(MPI_Processes_X*MPI_Processes_Y)
    {
        // Determination of the rank and size of the run
        MPI_Comm_size(MPI_COMM_WORLD, &MPI_Procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        rankx.push_back(rank % MPI_Processes_X);
        rankx.push_back(rank / MPI_Processes_X);

        if (error_check()) {MPI_Finalize(); exit(1);}

        // Determination of the local computational domain (i.e. the x-axis and the y-axis) 
        for(size_t i(0); i < Input::List().xminLocal.size(); ++i) {
            Input::List().xminLocal[i] = Input::List().xminGlobal[i]
                                        + rankx[i] * Input::List().NxLocalnobnd[i] * Input::List().globdx[i] 
                                        - Input::List().BoundaryCells * Input::List().globdx[i];
            Input::List().xmaxLocal[i] = Input::List().xminLocal[i] 
                                        + (Input::List().NxLocal[i]) * Input::List().globdx[i];
        }

//        exit(0);

         // Restart files will be generated when output_step % restart_step == 0 
//          if ( ( (Input::List().n_outsteps + 1) >  Input::List().n_restarts ) &&
//                (Input::List().n_restarts > 0)  ) { 
//               restart_step = Input::List().n_outsteps / Input::List().n_restarts ;
//          } 
//          else restart_step = -1;
//          if ( Input::List().restart_time > Input::List().n_restarts ) {
//              restart_time = 0; 
//          }
//          else restart_time = Input::List().restart_time;

    }
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Parallel_Environment_2D:: ~Parallel_Environment_2D(){ }
//--------------------------------------------------------------

//--------------------------------------------------------------
    bool Parallel_Environment_2D:: error_check() {
//--------------------------------------------------------------
//  Temporary error check
//--------------------------------------------------------------
        if ( MPI_Procs != (MPI_Processes_X * MPI_Processes_Y) ) { 
            if (rank == 0) std:: cout << "the number of nodes in the input deck is "
                                      << (MPI_Processes_X * MPI_Processes_Y) << ", terminating ..." << endl;
            return true; 
        }
        // Test the number of cells
        if (rank == 0){ 
            if (Input::List().NxLocal[0] < 2*Input::List().BoundaryCells+1){
                std::cout<<"Not enough cells per processor in the x direction" << endl;
                return true;
            }
            if (Input::List().NxLocal[1] < 2*Input::List().BoundaryCells+1){ 
                std::cout<<"Not enough cells per processor in the y direction" << endl;
                return true;
            }
         } 
         return false;
    }
//--------------------------------------------------------------
 
//--------------------------------------------------------------
//  Basic information
//--------------------------------------------------------------
    int Parallel_Environment_2D:: RANK()  const {return rank;} 
    int Parallel_Environment_2D:: RANKX() const {return rankx[0];} 
    int Parallel_Environment_2D:: RANKY() const {return rankx[1];} 
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    int Parallel_Environment_2D:: MPI_Processes() const {return MPI_Procs;} 
    int Parallel_Environment_2D:: MPI_X() const {return MPI_Processes_X;} 
    int Parallel_Environment_2D:: MPI_Y() const {return MPI_Processes_Y;} 
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    int Parallel_Environment_2D:: BNDX()  const {return bndX;} 
    int Parallel_Environment_2D:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Environment_2D::Neighbor_ImplicitE_Communications(State2D& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------

        int moduloX(RANKX()%2); 
        int RNx((RANKX()+1)%MPI_X() +RANKY()*MPI_X()),         // This is the right neighbor 
            LNx((RANKX()-1+MPI_X())%MPI_X()+RANKY()*MPI_X()); // This is the left  neighbor 
        
        
        if (MPI_X() > 1) {
            //even nodes 
            if (moduloX==0){
                Bfield_Data.Send_right_X(Y,RNx);                  //   (Send) 0 --> 1               
                if ((RANKX() != 0) || (BNDX()==0)){
                    Bfield_Data.Recv_from_left_X(Y,LNx);          //          1 --> 0 (Receive)   
                    Bfield_Data.Send_left_X(Y,LNx);               //          1 <-- 0 (Send)
                }
                else {
                    if (BNDX()==1) {
                        Bfield_Data.mirror_bound_Xleft(Y);        // Update node "0" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
               } 
               Bfield_Data.Recv_from_right_X(Y,RNx);               // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                Bfield_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)               
                if ((RANKX()!=(MPI_X()-1)) || (BNDX()==0)){
                     Bfield_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0 
                     Bfield_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
                }
                else {
                    if (BNDX()==1) {
                        Bfield_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;;
                    }
                } 
                Bfield_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)              
            } 
        }
        else { Bfield_Data.sameNode_bound_X(Y); }


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

         int moduloY(RANKY()%2);
         int RNy((RANK()+MPI_X())%MPI_Processes()),                  // This is the right neighbor 
             LNy((RANK()-MPI_X()+MPI_Processes())%MPI_Processes());          // This is the left  neighbor 

         if (MPI_Y() > 1) {
            //even nodes 
            if (moduloY==0){
                Bfield_Data.Send_right_Y(Y,RNy);                     //   (Send) 0 --> 1               
                if ((RANKY() != 0) || (BNDY()==0)){
                    Bfield_Data.Recv_from_left_Y(Y,LNy);             //          1 --> 0 (Receive)                
                    Bfield_Data.Send_left_Y(Y,LNy);                  //          1 <-- 0 (Send)
                }
                else {
                    if (BNDY()==1) {
                        Bfield_Data.mirror_bound_Yleft(Y);          //  Update node "0" in the y direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
                }  
                Bfield_Data.Recv_from_right_Y(Y,RNy);                // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                Bfield_Data.Recv_from_left_Y(Y,LNy);                 //           0 --> 1 (Receive)               
                if ((RANKY()!=(MPI_Y()-1)) || (BNDY()==0)){
                     Bfield_Data.Send_right_Y(Y,RNy);                //   (Send)  1 --> 0 
                     Bfield_Data.Recv_from_right_Y(Y,RNy);           // (Receive) 1 <-- 0
                }
                else {
                    if (BNDY()==1) {
                        Bfield_Data.mirror_bound_Yright(Y);          // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
                } 
                Bfield_Data.Send_left_Y(Y,LNy);                      //           0 <-- 1 (Send)              
            }
         }
         else { Bfield_Data.sameNode_bound_Y(Y); }
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Environment_2D::Neighbor_Communications(State2D& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------

        int moduloX(RANKX()%2); 
        int RNx((RANKX()+1)%MPI_X() +RANKY()*MPI_X()),         // This is the right neighbor 
            LNx((RANKX()-1+MPI_X())%MPI_X()+RANKY()*MPI_X()); // This is the left  neighbor 


        // std::cout << "Parallel: \n";
        // std::cout << "MPI_X() = "<< MPI_X() << " \n";
        // std::cout << "MPI_Y() = "<< MPI_Y() << " \n";    

        if (MPI_X() > 1) {
            //even nodes 
            if (moduloX==0){
               X_Data.Send_right_X(Y,RNx);                  //   (Send) 0 --> 1               
                if ((RANKX() != 0) || (BNDX()==0)){
                    X_Data.Recv_from_left_X(Y,LNx);          //          1 --> 0 (Receive)   
                    X_Data.Send_left_X(Y,LNx);               //          1 <-- 0 (Send)
                }
                else {
                    if (BNDX()==1) {
                        X_Data.mirror_bound_Xleft(Y);        // Update node "0" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
               } 
               X_Data.Recv_from_right_X(Y,RNx);               // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                X_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)               
                if ((RANKX()!=(MPI_X()-1)) || (BNDX()==0)){
                     X_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0 
                     X_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
                }
                else {
                    if (BNDX()==1) {
                        X_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
               } 
                X_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)              
            } 
         }
         else { X_Data.sameNode_bound_X(Y); }


//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

         int moduloY(RANKY()%2);
         int RNy((RANK()+MPI_X())%MPI_Processes()),                  // This is the right neighbor 
             LNy((RANK()-MPI_X()+MPI_Processes())%MPI_Processes());          // This is the left  neighbor 

         if (MPI_Y() > 1) {
            //even nodes 
            if (moduloY==0){
                X_Data.Send_right_Y(Y,RNy);                     //   (Send) 0 --> 1               
                if ((RANKY() != 0) || (BNDY()==0)){
                    X_Data.Recv_from_left_Y(Y,LNy);             //          1 --> 0 (Receive)                
                    X_Data.Send_left_Y(Y,LNy);                  //          1 <-- 0 (Send)
                }
                else {
                    if (BNDY()==1) {
                        X_Data.mirror_bound_Yleft(Y);          //  Update node "0" in the y direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
                }  
                X_Data.Recv_from_right_Y(Y,RNy);                // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                X_Data.Recv_from_left_Y(Y,LNy);                 //           0 --> 1 (Receive)               
                if ((RANKY()!=(MPI_Y()-1)) || (BNDY()==0)){
                     X_Data.Send_right_Y(Y,RNy);                //   (Send)  1 --> 0 
                     X_Data.Recv_from_right_Y(Y,RNy);           // (Receive) 1 <-- 0
                }
                else {
                    if (BNDY()==1) {
                        X_Data.mirror_bound_Yright(Y);          // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary." << endl;
                    }
                } 
                X_Data.Send_left_Y(Y,LNy);                      //           0 <-- 1 (Send)              
            }
         }
         else { X_Data.sameNode_bound_Y(Y); }

    }
//--------------------------------------------------------------

//**************************************************************


//**************************************************************
