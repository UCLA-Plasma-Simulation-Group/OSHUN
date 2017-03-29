///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras, Benjamin Winjum
//
//	Modified:	September 1 2016
///////////////////////////////////////////////////////////

//   
//   Contains the declerations for the communications
//   between nodes, boundaries and the parallel output
///////////////////////////////////////////////////////////
//
// 
//   This file contains three modules:
//
//   1. class Node_ImplicitE_Communications: 
//
//   2. class Node_Communications:
//        Allows the nodes to exchange information in order
//        to update their guard cells. For boundary nodes 
//        it provides the appropriate boundary conditions.
// 
//   3. class Parallel_Environment:
//        - It decomposes the computational domain
//        - It controls the node communications
//        - It controls the parallel output
//        - It controls the restart facility
//
///////////////////////////////////////////////////////////
//

    #ifndef PARALLEL_ENVIRONMENT_H
    #define PARALLEL_ENVIRONMENT_H


//**************************************************************
//--------------------------------------------------------------
        class Node_ImplicitE_Communications {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_ImplicitE_Communications(); 
            ~Node_ImplicitE_Communications();
         
//          Boundary conditions
            int BNDX()   const;

//          Data exchange in x direction
            void Send_right_X(State1D& Y, int dest);
            void Recv_from_left_X(State1D& Y,int origin);
            void Send_left_X(State1D& Y, int dest);
            void Recv_from_right_X(State1D& Y, int origin); 

//          Boundaries 
            void mirror_bound_Xleft(State1D& Y);
            void mirror_bound_Xright(State1D& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(State1D& Y);

        private:
//          Domain information
            int Nbc, bndX;

//          Information exchange
            int  msg_sizeX;
			complex<double> *msg_bufX;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(State1D& Y);
            void sameNode_mirror_X(State1D& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Node_Communications {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_Communications(); 
            ~Node_Communications();
         
//          Boundary conditions
            int BNDX()   const;

//          Data exchange in x direction
            void Send_right_X(State1D& Y, int dest);
            void Recv_from_left_X(State1D& Y,int origin);
            void Send_left_X(State1D& Y, int dest);
            void Recv_from_right_X(State1D& Y, int origin); 

//          Boundaries 
            void mirror_bound_Xleft(State1D& Y);
            void mirror_bound_Xright(State1D& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(State1D& Y);

        private:
//          Domain information
            int Nbc, bndX;
            int numspec, numpmax;

//          Information exchange
            int  msg_sizeX;
            complex<double> *msg_bufX;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(State1D& Y);
            void sameNode_mirror_X(State1D& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Parallel_Environment_1D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Parallel_Environment_1D(); 
            ~Parallel_Environment_1D(); 
         
//          Parallel parameters
            int RANK()  const;
            int NODES() const;
            int BNDX()  const;

//          Restart
//             bool READ_RESTART() const;
//             void Read_Restart(State1D& Y);

//             bool WRITE_RESTART(const size_t step) const;
//             void Write_Restart(State1D& Y, size_t step);

//             size_t T_IN() const;
          
//          Information exchange
            void Neighbor_ImplicitE_Communications(State1D& Y);
            void Neighbor_Communications(State1D& Y);

        private:
//          Parallel parameters
            int rank;
            int Nnodes;

//          Boundaries 
            int bndX;

//          Information Exchange
            Node_ImplicitE_Communications Bfield_Data;
            Node_Communications X_Data;

//          Error Checking of the constructor
            bool error_check();
        };
//--------------------------------------------------------------
//**************************************************************



    #endif
