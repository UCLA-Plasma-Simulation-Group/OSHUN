/*!\brief  Parallelization routines - Definitions
* \author  PICKSC
 * \date   March, 2017
 * \file   parallel.cpp
 *
 * In here are the structures that enable parallelization
 * 
 * Periodic and reflecting have been implemented.
 * 
 * 
 */

    #ifndef PARALLEL_ENVIRONMENT_H
    #define PARALLEL_ENVIRONMENT_H


//**************************************************************
//--------------------------------------------------------------
        class Node_ImplicitE_Communications_1D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_ImplicitE_Communications_1D(); 
            ~Node_ImplicitE_Communications_1D();
         
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
            int  msg_parsizeX;
			
            complex<double> *msg_bufX;

            // std::vector<complex<double> > msg_parX;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(State1D& Y);
            void sameNode_mirror_X(State1D& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Node_Communications_1D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_Communications_1D(); 
            ~Node_Communications_1D();
         
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
            int MPI_Processes() const;
            int BNDX()  const;
          
//          Information exchange
            void Neighbor_ImplicitE_Communications(State1D& Y);
            void Neighbor_Communications(State1D& Y);

        private:
//          Parallel parameters
            int rank;
            int MPI_Procs;

//          Boundaries 
            int bndX;

//          Information Exchange
            Node_ImplicitE_Communications_1D Bfield_Data;
            Node_Communications_1D X_Data;

//          Error Checking of the constructor
            bool error_check();
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Node_ImplicitE_Communications_2D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_ImplicitE_Communications_2D(); 
            ~Node_ImplicitE_Communications_2D();
         
//          Boundary conditions
            int BNDX()   const;
            int BNDY()   const;

//          Data exchange in x direction
            void Send_right_X(State2D& Y, int dest);
            void Recv_from_left_X(State2D& Y,int origin);
            void Send_left_X(State2D& Y, int dest);
            void Recv_from_right_X(State2D& Y, int origin); 

//          Data exchange in y direction
            void Send_right_Y(State2D& Y, int dest);
            void Recv_from_left_Y(State2D& Y,int origin);
            void Send_left_Y(State2D& Y, int dest);
            void Recv_from_right_Y(State2D& Y, int origin);             

//          Boundaries 
            void mirror_bound_Xleft(State2D& Y);
            void mirror_bound_Xright(State2D& Y);

            void mirror_bound_Yleft(State2D& Y);
            void mirror_bound_Yright(State2D& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(State2D& Y);
            void sameNode_bound_Y(State2D& Y);

        private:
//          Domain information
            size_t Nbc, bndX, bndY;
            size_t Nx_local, Ny_local;


//          Information exchange
            int  msg_sizeX, msg_sizeY;
            complex<double> *msg_bufX, *msg_bufY;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(State2D& Y);
            void sameNode_periodic_Y(State2D& Y);
            void sameNode_mirror_X(State2D& Y);
            void sameNode_mirror_Y(State2D& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Node_Communications_2D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_Communications_2D(); 
            ~Node_Communications_2D();
         
//          Boundary conditions
            int BNDX()   const;
            int BNDY()   const;

//          Data exchange in x direction
            void Send_right_X(State2D& Y, int dest);
            void Recv_from_left_X(State2D& Y,int origin);
            void Send_left_X(State2D& Y, int dest);
            void Recv_from_right_X(State2D& Y, int origin); 

            void Send_right_Y(State2D& Y, int dest);
            void Recv_from_left_Y(State2D& Y,int origin);
            void Send_left_Y(State2D& Y, int dest);
            void Recv_from_right_Y(State2D& Y, int origin); 

//          Boundaries 
            void mirror_bound_Xleft(State2D& Y);
            void mirror_bound_Xright(State2D& Y);

            void mirror_bound_Yleft(State2D& Y);
            void mirror_bound_Yright(State2D& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(State2D& Y);
            void sameNode_bound_Y(State2D& Y);

        private:
//          Domain information
            // int Nbc, bndX, bndY;
            int numspec, numpmax;

            size_t Nbc, bndX, bndY;
            size_t Nx_local, Ny_local;

//          Information exchange
            int  msg_sizeX, msg_sizeY; 
            complex<double> *msg_bufX, *msg_bufY;
            
//          Boundaries for single-node configurations
            void sameNode_periodic_X(State2D& Y);
            void sameNode_mirror_X(State2D& Y);

            void sameNode_periodic_Y(State2D& Y);
            void sameNode_mirror_Y(State2D& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Parallel_Environment_2D {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Parallel_Environment_2D(); 
            ~Parallel_Environment_2D(); 
         
//          Parallel parameters
            int RANK()  const;
            int MPI_Processes() const;
            int BNDX()  const;
            int BNDY()  const;
            int RANKX()  const;
            int RANKY()  const;
            int MPI_X() const;
            int MPI_Y() const;

//          Restart 
//             bool READ_RESTART() const;
//             void Read_Restart(State2D& Y); 
// 
//             bool WRITE_RESTART(const size_t step) const;
//             void Write_Restart(size_t step, State2D& Y); 
// 
//             size_t T_IN() const;
          
//          Information exchange
            void Neighbor_ImplicitE_Communications(State2D& Y);
            void Neighbor_Communications(State2D& Y);

        private:
//          Parallel parameters

            int rank;
            vector<int> rankx;
            int MPI_Procs, MPI_Processes_X, MPI_Processes_Y;

//          Boundaries 
            int bndX, bndY;


//          Information Exchange
            Node_ImplicitE_Communications_2D Bfield_Data;
            Node_Communications_2D X_Data;

//          Error Checking of the constructor
            bool error_check();

//          Restart files 
//             int restart_time;
//             int restart_step;
        };
//--------------------------------------------------------------
//**************************************************************

    #endif
