/*! \brief Implicit Electric Field - Definitions
 * \author PICKSC
 * \date   2017
 * \file   implicitE.cpp
 * 
 * Contains:
 * 1) the explicit electric field methods
 * 2) implicit electric field solver -- raison d'etre
 */

    #ifndef DECL_IMPLICITEFIELD_H
    #define DECL_IMPLICITEFIELD_H

//**************************************************************
//**************************************************************
//   Definition for the Electric_Field_Methods namespace
//**************************************************************
//**************************************************************


    namespace Electric_Field_Methods { 
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



//**************************************************************
//--------------------------------------------------------------
        class Current_xyz {
//--------------------------------------------------------------
//      Decleration of the Current 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Current_xyz(); 

//          Choose component
            Field1D& J_1D(int component);  
            void calculate_J_1D(State1D& Yin);  

            Field2D& J_2D(int component);  
            void calculate_J_2D(State2D& Yin);  
            
//          Components
            Field1D& Jx_1D();              
            Field1D& Jy_1D();              
            Field1D& Jz_1D();  
            Field2D& Jx_2D();              
            Field2D& Jy_2D();              
            Field2D& Jz_2D();  
            
        private:

            size_t Nbc, szx, szy;

            Field1D jayx_1D, jayy_1D, jayz_1D;


            Array2D<double>              current1D;


            Field2D jayx_2D, jayy_2D, jayz_2D;
            Array3D<double>              current2D;
        };
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
        class Efield_xyz {
//--------------------------------------------------------------
//      Decleration of the Current 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Efield_xyz(); 

//          Choose component
            Field1D& E_1D(int component);  
            Field2D& E_2D(int component);  

//          Components
            Field1D& Ex_1D();  
            Field1D& Ey_1D();  
            Field1D& Ez_1D();

            Field2D& Ex_2D();  
            Field2D& Ey_2D();  
            Field2D& Ez_2D();  

        private:
            Field1D efieldx_1D, efieldy_1D, efieldz_1D;
            Field2D efieldx_2D, efieldy_2D, efieldz_2D;

//            complex<double> Delta_p, small;
//            Axis< complex<double> >  p3og;

        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Efield_Method {
//--------------------------------------------------------------
//      Abstract class for explicit methods
//--------------------------------------------------------------
        public:
            virtual void advance(Algorithms::RK2<State1D>* rk, State1D& Y, collisions_1D& coll, VlasovFunctor1D_implicitE_p2* rkF)=0;//, double time, double dt) = 0; // "covariant" return
            virtual void advance(Algorithms::RK2<State2D>* rk, State2D& Y, collisions_2D& coll, VlasovFunctor2D_implicitE_p2* rkF)=0;//, double time, double dt) = 0; // "covariant" return
            virtual ~Efield_Method() = 0;
        };
//--------------------------------------------------------------
//**************************************************************


// //**************************************************************
// //--------------------------------------------------------------
//         class Explicit_E_Field: public Efield_Method {
// //--------------------------------------------------------------
// //      class for explicit electric field              
// //--------------------------------------------------------------
//         public:
// //          Constructor
//             Explicit_E_Field(State1D& Yin, int tout_start);

// //          Main function
//             void advance(Algorithms::RK2<State1D>* rk, Euler_Backward* eb, RKFunctor1D* rkF);    
         
// //          Is the electric field implicit?
//             bool implicitE() const;

//         private:

// //          Reference to Y
//             State1D& Y;

// //          Time
//             double t;

// //          This class does not do anything at this stage            

//         };
// //--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Implicit_E_Field: public Efield_Method {
//--------------------------------------------------------------
//      class for ==> implicit electric field              
//--------------------------------------------------------------
        public:
//          Constructor
            Implicit_E_Field(const double& deltat, const Algorithms::AxisBundle<double> axes);

//          Main function
            void advance(Algorithms::RK2<State1D>* rk, State1D& Y, collisions_1D& coll, VlasovFunctor1D_implicitE_p2* rkF);
            void advance(Algorithms::RK2<State2D>* rk, State2D& Y, collisions_2D& coll, VlasovFunctor2D_implicitE_p2* rkF);

        private:
//          Boundary Cells
            size_t      Nbc, szx, szy;

            //          Current and electric field 
            Current_xyz JN, J0, J_Ex, J_Ey, J_Ez;
            Efield_xyz  EN, E0, DE;
            double dt;
            complex<double>   idx, idy;            

//          Ampere's law JN = rot(B)
            void Ampere(EMF1D& emf);
            void FindDE(EMF1D& emf);

            void Ampere(EMF2D& emf);
            void FindDE(EMF2D& emf);

        };
//--------------------------------------------------------------

//**************************************************************

    }
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #endif
