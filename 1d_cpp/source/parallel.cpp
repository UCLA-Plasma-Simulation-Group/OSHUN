///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras, Benjamin Winjum
//
//	Last Modified:	September 1, 2016
///////////////////////////////////////////////////////////

//   
//   Contains the declerations for the communications
//   between nodes, boundaries and the parallel output
///////////////////////////////////////////////////////////
//
//   This file contains three modules:
//
//   1. class Node_Communications: 
//        Allows the nodes to exchange information in order
//        to update their guard cells. For boundary nodes 
//        it provides the appropriate boundary conditions.
//
//   2. class Parallel_Output: 
//        Collects the output information from all of the
//        nodes and combines them to the final file ready
//        to be exported. For the moments of the distribution
//        function it also performs the integration over 
//        momentum space and for detailed information on 
//        phasespace it calls the appropriate functions from
//        "Output" to convert the spherical harmonics to 
//        cartesian geometry.
// 
//   3. class Parallel_Environment:
//        - It decomposes the computational domain
//        - It controls the node communications
//        - It controls the parallel output
//        - It controls the restart facility
//
///////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

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
#include "H5Cpp.h"
#include "lib-array.h"
#include "lib-algorithms.h"
#include <map>

//  Declerations
#include "input.h"
#include "state.h"
//     #include "decl-output.h"
// #include "setup.h"
//     #include "export.h"
#include "parallel.h"


//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Node_ImplicitE_Communications:: Node_ImplicitE_Communications() :
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
Node_ImplicitE_Communications:: ~Node_ImplicitE_Communications(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] msg_bufX;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
int Node_ImplicitE_Communications:: BNDX()  const {return bndX;}
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_ImplicitE_Communications::Send_right_X(State1D& Y, int dest) {
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
void Node_ImplicitE_Communications::Recv_from_left_X(State1D& Y, int origin) {
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
void Node_ImplicitE_Communications::Send_left_X(State1D& Y, int dest) {
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
void Node_ImplicitE_Communications::Recv_from_right_X(State1D& Y, int origin) {
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
void Node_ImplicitE_Communications::mirror_bound_Xleft(State1D& Y) {
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
void Node_ImplicitE_Communications::mirror_bound_Xright(State1D& Y) {
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
void Node_ImplicitE_Communications::sameNode_bound_X(State1D& Y) {
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
void Node_ImplicitE_Communications::sameNode_periodic_X(State1D& Y) {
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
void Node_ImplicitE_Communications::sameNode_mirror_X(State1D& Y) {
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
Node_Communications:: Node_Communications() :
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
        msg_sizeX += temp*Input::List().ps[s];
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
    msg_bufX = new complex<double>[msg_sizeX];

}
//--------------------------------------------------------------

//--------------------------------------------------------------
Node_Communications:: ~Node_Communications(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    delete[] msg_bufX;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
int Node_Communications:: BNDX()  const {return bndX;}
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_Communications::Send_right_X(State1D& Y, int dest) {
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
    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications::Recv_from_left_X(State1D& Y, int origin) {
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


}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications::Send_left_X(State1D& Y, int dest) {
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

    MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications::Recv_from_right_X(State1D& Y, int origin) {
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
void Node_Communications::mirror_bound_Xleft(State1D& Y) {
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


}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications::mirror_bound_Xright(State1D& Y) {
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

}
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
void Node_Communications::sameNode_bound_X(State1D& Y) {
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
void Node_Communications::sameNode_periodic_X(State1D& Y) {
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

    // temp = Y.DF(0).getcurrent(0);
    // for (size_t ix(0);ix<Y.SH(0,0,0).numx();++ix)
    // {
    //     std::cout<<"\n      new[" << ix  << "] = " << temp[ix] << "\n";

    // }


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
            // Y.HYDRO().kpressure(c) = Y.HYDRO().velocity(Y.EMF().Ex().numx()-2*Nbc+c);
            // Y.HYDRO().mpressure(c) = Y.HYDRO().velocity(Y.EMF().Ex().numx()-2*Nbc+c);
        }

        // Hydro Quantities:   x0 "Left-Bound ---> Right-Guard"
        for(int c(0); c < Nbc; c++) {
            Y.HYDRO().density(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().density(Nbc+c);
            Y.HYDRO().vx(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vx(Nbc+c);
            Y.HYDRO().vy(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vy(Nbc+c);
            Y.HYDRO().vz(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().vz(Nbc+c);
            Y.HYDRO().temperature(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().temperature(Nbc+c);
            Y.HYDRO().Z(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().Z(Nbc+c);
            // Y.HYDRO().kpressure(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);
            // Y.HYDRO().mpressure(Y.EMF().Ex().numx()-Nbc+c) =  Y.HYDRO().velocity(Nbc+c);

        }
    }

}
//--------------------------------------------------------------

//--------------------------------------------------------------
void Node_Communications::sameNode_mirror_X(State1D& Y) {
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
        Nnodes(Input::List().NnodesX)   // Number of nodes in X-direction
{
    // Determination of the rank and size of the run
    MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (error_check()) {MPI_Finalize(); exit(1);}

    // Determination of the local computational domain (i.e. the x-axis and the y-axis)

    for(size_t i(0); i < Input::List().xminLocal.size(); ++i) {

        Input::List().xminLocal[i] = Input::List().xminGlobal[i]
                                     + rank * Input::List().NxLocalnobnd[i] * Input::List().globdx[i]
                                     - Input::List().BoundaryCells * Input::List().globdx[i];
        Input::List().xmaxLocal[i] = Input::List().xminLocal[i]
                                     + (Input::List().NxLocal[i]) * Input::List().globdx[i];
    }

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
        if ( Nnodes != (Input::List().NnodesX) ) {
            std:: cout << "the number of nodes in the input deck is "
                       << (Input::List().NnodesX) << ", terminating ..." << endl;
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
int Parallel_Environment_1D:: NODES() const {return Nnodes;}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int Parallel_Environment_1D:: BNDX()  const {return bndX;}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Parallel_Environment_1D::Neighbor_ImplicitE_Communications(State1D& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------

    int moduloX(RANK()%2);
    int RNx((RANK()+1)%NODES()),         // This is the right neighbor
            LNx((RANK()-1+NODES())%NODES()); // This is the left  neighbor

    if (NODES() > 1) {
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
            if ((RANK()!=(NODES()-1)) || (BNDX()==0)){
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
    int RNx((RANK() + 1) % NODES()),         // This is the right neighbor
            LNx((RANK() - 1 + NODES()) % NODES()); // This is the left  neighbor

    if (NODES() > 1) {
        //even nodes
//        if (moduloX==0){
        if (BNDX() == 0) {
            if (((RANK() != 0) && (RANK() != NODES() - 1))) {
                X_Data.Send_right_X(Y, RNx);                  //   (Send) 0 --> 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
                X_Data.Recv_from_right_X(Y, RNx);               // (Receive) 0 <-- 1
            } else if (RANK() == 0) {                         /// Update node "0" in the x direction
                X_Data.Recv_from_left_X(Y, NODES() - 1);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, NODES() - 1);               //          1 <-- 0 (Send)
                X_Data.Send_right_X(Y, RNx);                ///   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, RNx);           /// (Receive) 0 <-- 1
            } else if (RANK() == NODES() - 1) {               ///        // Update node "NODES()" in the x direction
                X_Data.Send_right_X(Y, 0);                  //   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, 0);           /// (Receive) 0 <-- 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
            }
        } else if (BNDX() == 1) {
            if (((RANK() != 0) && (RANK() != NODES() - 1))) {
                X_Data.Send_right_X(Y, RNx);                  //   (Send) 0 --> 1
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
                X_Data.Recv_from_right_X(Y, RNx);               // (Receive) 0 <-- 1
            } else if (RANK() == 0) {                         /// Update node "0" in the x direction
                X_Data.mirror_bound_Xleft(Y);
                X_Data.Send_right_X(Y, RNx);                ///   (Send) 0 --> 1
                X_Data.Recv_from_right_X(Y, RNx);           /// (Receive) 0 <-- 1
            } else if (RANK() == NODES() - 1) {               ///        // Update node "NODES()" in the x direction
                X_Data.mirror_bound_Xright(Y);
                X_Data.Recv_from_left_X(Y, LNx);          //          1 --> 0 (Receive)
                X_Data.Send_left_X(Y, LNx);               //          1 <-- 0 (Send)
            }

        } else {
            cout << "Invalid Boundary." << endl;
        }

//            //odd nodes
//        else {
//            X_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)
//            if ((RANK()!=(NODES()-1)) || (BNDX()==0)){
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
