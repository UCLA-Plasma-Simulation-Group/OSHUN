/*! \brief Implicit Electric Field - Definitions
 * \author PICKSC
 * \date   September 23, 2016
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
            Current_xyz(EMF1D& emf); 

//          Choose component
            Field1D& J(int component);  
            void calculate_J(State1D& Yin);  
            

//          Components
            Field1D& Jx();  
            // void calculate_Jx(State1D& Yin);  
            Field1D& Jy();  
            // void calculate_Jy(State1D& Yin);  
            Field1D& Jz();  
            // void calculate_Jz(State1D& Yin);  

        private:
            Field1D jayx, jayy, jayz;

            // complex<double> Delta_p, small;
            complex<double> fourpioverthree;
            // valarray< complex<double> >  p3og;
            // vector <valarray< complex<double> >>  p3ogv;

            valarray< complex<double> >  tempx; 
            valarray< complex<double> >  tempy; 
            valarray< complex<double> >  tempz; 

            size_t Nbc, szx;

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
            Efield_xyz(EMF1D& emf); 

//          Choose component
            Field1D& E(int component);  

//          Components
            Field1D& Ex();  
            Field1D& Ey();  
            Field1D& Ez();  

        private:
            Field1D efieldx, efieldy, efieldz;

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
            virtual void advance(Algorithms::RK3<State1D>* rk, State1D& Y, collisions& coll, VlasovFunctor1D_implicitE_p2* rkF)=0;//, double time, double dt) = 0; // "covariant" return
            // virtual bool implicitE()  const = 0;
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
//             void advance(Algorithms::RK3<State1D>* rk, Euler_Backward* eb, RKFunctor1D* rkF);    
         
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
            Implicit_E_Field(EMF1D& emf, const double& deltat, const double& deltax);

//          Main function
            void advance(Algorithms::RK3<State1D>* rk, State1D& Y, collisions& coll, VlasovFunctor1D_implicitE_p2* rkF);//, double time, double dt);
         
//          Is the electric field implicit?
            // bool implicitE() const;

        private:
//          Reference to Y
            // State1D& Y;

//          Current and electric field 
            Current_xyz JN, J0, J_Ex, J_Ey, J_Ez;
            Efield_xyz  EN, E0, DE;
            double dt;
            complex<double>   idx;//, idy;            

//          Storage of harmonics
            // State1D Ysafe;
            // SHarmonic1D f00; 
            // SHarmonic1D f10, f11;
            // SHarmonic1D f20, f21, f22;

//          Time
            // double t;


//          Ampere's law JN = rot(B)
            void Ampere(EMF1D& emf);
            void FindDE(EMF1D& emf);

//          Boundary Cells
            int               Nbc, szx;//, szy;
//          Derivative constants

            // size_t            l0, m0;

//          Implicit methods
            // bool if_implicitES;
            // bool if_implicit1D;
        };
//--------------------------------------------------------------

//**************************************************************

    }
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #endif
