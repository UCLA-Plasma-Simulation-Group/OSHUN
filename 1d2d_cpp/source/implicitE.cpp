/*! \brief Implicit Electric Field - Definitions
* \author PICKSC
 * \date   April 2017
 * \file   implicitE.cpp
 * 
 * Contains:
 * 1) the explicit electric field methods
 * 2) implicit electric field solver -- raison d'etre
 */


//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>
    #include <float.h>
    #include <math.h>
    #include <map>
    #include <limits>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declarations
    #include "input.h"
    #include "state.h"
    #include "formulary.h"
    #include "nmethods.h"
    #include "collisions.h"
    #include "fluid.h"
    #include "vlasov.h"
    #include "functors.h"
    #include "implicitE.h"
//**************************************************************
//**************************************************************
//**************************************************************
//*******THE CLASS FOR THE CURRENT DEFINITION ******************
//**************************************************************
//**************************************************************
//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Current_xyz::Current_xyz() 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :    Nbc(Input::List().BoundaryCells),
            szx(Input::List().NxLocal[0]),
            szy(Input::List().NxLocal[1]),
            jayx_1D(szx), 
            jayy_1D(szx), 
            jayz_1D(szx), 
            current1D(3,szx),
            jayx_2D(szx,szy), 
            jayy_2D(szx,szy), 
            jayz_2D(szx,szy), 
            current2D(3,szx,szy)
         {}
//--------------------------------------------------------------
//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::J_1D(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Jx_1D(); 
            break;
        }
        case 2: { 
            return Jy_1D(); 
            break;
        }
        case 3: { 
            return Jz_1D(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Jx_1D();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jx_1D() {
       return jayx_1D; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jy_1D() {
        return jayy_1D;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jz_1D() {
        return jayz_1D;
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::J_2D(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Jx_2D(); 
            break;
        }
        case 2: { 
            return Jy_2D(); 
            break;
        }
        case 3: { 
            return Jz_2D(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Jx_2D();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jx_2D() {
       return jayx_2D; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jy_2D() {
        return jayy_2D;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jz_2D() {
        return jayz_2D;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_J_1D(State1D& Yin) {
//--------------------------------------------------------------
//  Update the total current 
//--------------------------------------------------------------
        jayx_1D = static_cast<complex<double> >(0.0);
        jayy_1D = static_cast<complex<double> >(0.0);
        jayz_1D = static_cast<complex<double> >(0.0);

        for (size_t s(0); s < Yin.Species(); ++s)
        {   
            current1D = Yin.DF(s).getcurrent();      
            for (size_t ix(0); ix < szx; ++ix)
            {   
                jayx_1D(ix) +=   static_cast<complex<double> >(current1D(0,ix));
                jayy_1D(ix) +=   static_cast<complex<double> >(current1D(1,ix));
                jayz_1D(ix) +=   static_cast<complex<double> >(current1D(2,ix));
            }
        }
}    
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_J_2D(State2D& Yin) {
//--------------------------------------------------------------
//  Update the total current 
//--------------------------------------------------------------
        jayx_2D = static_cast<complex<double> >(0.0);
        jayy_2D = static_cast<complex<double> >(0.0);
        jayz_2D = static_cast<complex<double> >(0.0);

        for (size_t s(0); s < Yin.Species(); ++s)
        {   
            current2D = Yin.DF(s).getcurrent();      
            for (size_t ix(0); ix < szx; ++ix)
            {   
                for (size_t iy(0); iy < szy; ++iy)
                {   
                    jayx_2D(ix,iy) +=   static_cast<complex<double> >(current2D(0,ix,iy));
                    jayy_2D(ix,iy) +=   static_cast<complex<double> >(current2D(1,ix,iy));
                    jayz_2D(ix,iy) +=   static_cast<complex<double> >(current2D(2,ix,iy));
                }
            }
        }
}    
//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Efield_xyz::Efield_xyz() 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :
            efieldx_1D(Input::List().NxLocal[0]), 
            efieldy_1D(Input::List().NxLocal[0]), 
            efieldz_1D(Input::List().NxLocal[0]),
            efieldx_2D(Input::List().NxLocal[0],Input::List().NxLocal[1]), 
            efieldy_2D(Input::List().NxLocal[0],Input::List().NxLocal[1]), 
            efieldz_2D(Input::List().NxLocal[0],Input::List().NxLocal[1])
          {}
//--------------------------------------------------------------
//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::E_1D(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Ex_1D(); 
            break;
        }
        case 2: { 
            return Ey_1D(); 
            break;
        }
        case 3: { 
            return Ez_1D(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Ex_1D();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ex_1D() {
       return efieldx_1D; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ey_1D() {
        return efieldy_1D;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ez_1D() {
        return efieldz_1D;
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::E_2D(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Ex_2D(); 
            break;
        }
        case 2: { 
            return Ey_2D(); 
            break;
        }
        case 3: { 
            return Ez_2D(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Ex_2D();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ex_2D() {
       return efieldx_2D; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ey_2D() {
        return efieldy_2D;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ez_2D() {
        return efieldz_2D;
    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
//  Definition of the pure virtual destructor for the abstract
//  class Explicit_Method
//--------------------------------------------------------------
    Electric_Field_Methods::Efield_Method::
    ~Efield_Method(){}
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the implicit electric field 
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Implicit_E_Field::
    Implicit_E_Field(const double& deltat, const Algorithms::AxisBundle<double> axes)://, int tout_start)
//--------------------------------------------------------------
//  Constructor for the implicit electric field
//--------------------------------------------------------------
        Nbc(Input::List().BoundaryCells),
        szx(Input::List().NxLocal[0]),
        szy(Input::List().NxLocal[1]),
        JN(),   J0(),                     // New current, Current from E = 0
        J_Ex(), J_Ey(), J_Ez(),        // Current due to the effect of Ex, Ey, Ez
        EN(),   E0(),   DE(),          // New E, old E, perturbation E
   
         dt(deltat),
         idx(static_cast< complex<double> >(0.5/axes.dx(0))),
         idy(static_cast< complex<double> >(0.5/axes.dx(1)))
         {}
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::
    advance(Algorithms::RK2<State1D>* rk, State1D& Yin, collisions_1D& coll, VlasovFunctor1D_implicitE_p2* rkF, const double step_size){//, double time, double dt){
//--------------------------------------------------------------
//  Calculate the implicit electric field
//--------------------------------------------------------------

        int zeros_in_det(1);      // This counts the number of zeros in the determinant 
        int execution_attempt(0); // This counts the number of attempts to find invert the E-field
        State1D Yh(Yin);
        Yh = 0.0;
        
        FindDE(Yin.EMF());                           //  Reset DE

        SHarmonic1D f00(Yin.SH(0,0,0));
        SHarmonic1D f10(Yin.SH(0,1,0)), f11(Yin.SH(0,1,1));
        
// - - - - - - - - - - - - - - - - - - - - - -
        while ( (zeros_in_det > 0) && ( execution_attempt < 4) ) {  // Execute this loop at most twice
            zeros_in_det = 0;                                       // Count the zeros of the determinant
            ++execution_attempt;                                // Count the execusion attempts
// - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - -  
            // Effect of E = 0 on f00, f10, f11
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            
            J0.calculate_J_1D(Yh);
            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;

            // Ex
            Yin.EMF().Ex() = DE.Ex_1D();
            Yin = (*rk)(Yin, dt, rkF, 1);           // Effect of DEx on Y00, Y10, Y11, Y20, Y21, Y22            
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ex.calculate_J_1D(Yh);                   // Evaluate J(DEx)
            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;

            // Ey
            Yin.EMF().Ey() = DE.Ey_1D();               // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 2);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ey.calculate_J_1D(Yh);                   // Evaluate J(DEy)           

            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;

            // Ez
            Yin.EMF().Ez() = DE.Ez_1D();               // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 3);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22            
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ez.calculate_J_1D(Yh);                   // Evaluate J(DEy)
            
            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;


            Ampere(Yin.EMF());                           // Calculate JN
                           
// - - - - - - - - - - - - - - - - - - - - - -
            
            //  calculate new EN 

            Array2D< complex<double> > sgm(3,3);     
            valarray< complex<double> > clm(3);

            for (size_t ix(0); ix < szx; ++ix)
            {
                sgm(0,0) =  ( J_Ex.Jx_1D()(ix) - J0.Jx_1D()(ix) ) / DE.Ex_1D()(ix);  // sxx = dJ(E_x)_x/DE_x
                sgm(0,1) =  ( J_Ey.Jx_1D()(ix) - J0.Jx_1D()(ix) ) / DE.Ey_1D()(ix);  // sxy = dJ(E_y)_x/DE_y
                sgm(0,2) =  ( J_Ez.Jx_1D()(ix) - J0.Jx_1D()(ix) ) / DE.Ez_1D()(ix);  // sxz = dJ(E_z)_x/DE_z

                sgm(1,0) =  ( J_Ex.Jy_1D()(ix) - J0.Jy_1D()(ix) ) / DE.Ex_1D()(ix);  // syx = dJ(E_x)_y/DE_x
                sgm(1,1) =  ( J_Ey.Jy_1D()(ix) - J0.Jy_1D()(ix) ) / DE.Ey_1D()(ix);  // syy = dJ(E_y)_y/DE_y
                sgm(1,2) =  ( J_Ez.Jy_1D()(ix) - J0.Jy_1D()(ix) ) / DE.Ez_1D()(ix);  // syz = dJ(E_z)_y/DE_z

                sgm(2,0) =  ( J_Ex.Jz_1D()(ix) - J0.Jz_1D()(ix) ) / DE.Ex_1D()(ix);  // szx = dJ(E_x)_z/DE_x
                sgm(2,1) =  ( J_Ey.Jz_1D()(ix) - J0.Jz_1D()(ix) ) / DE.Ey_1D()(ix);  // szy = dJ(E_y)_z/DE_y
                sgm(2,2) =  ( J_Ez.Jz_1D()(ix) - J0.Jz_1D()(ix) ) / DE.Ez_1D()(ix);  // szz = dJ(E_z)_z/DE_z

                clm[0] =  JN.Jx_1D()(ix) - J0.Jx_1D()(ix);
                clm[1] =  JN.Jy_1D()(ix) - J0.Jy_1D()(ix);
                clm[2] =  JN.Jz_1D()(ix) - J0.Jz_1D()(ix);

                // std::cout << "\n\n" << sgm(0,0) << "  ,  " << sgm(0,1) << "  ,  " << sgm(0,2) << "\n";
                // std::cout           << sgm(1,0) << "  ,  " << sgm(1,1) << "  ,  " << sgm(1,2) << "\n";
                // std::cout           << sgm(2,0) << "  ,  " << sgm(2,1) << "  ,  " << sgm(2,2) << "\n";

                
            // Solve the 3 by 3 system of equations
                complex<double> D_sgm( Det33(sgm) );            // The Determinant of the conductivity tensor
                if ( abs( D_sgm.real() ) > 6.0*DBL_MIN ) 
                {
                    EN.Ex_1D()(ix) = Detx33(clm, sgm) / D_sgm;
                    EN.Ey_1D()(ix) = Dety33(clm, sgm) / D_sgm;
                    EN.Ez_1D()(ix) = Detz33(clm, sgm) / D_sgm;
                }
                else 
                {
                    ++zeros_in_det;                                 // Try again with a new perturbation as below
                
                    DE.Ex_1D()(ix) *= 2;                        // DEx is 2 times larger than before
                    DE.Ey_1D()(ix) *= 1.3;                      // DEy is 1.3 times larger  than before

                    EN.Ex_1D()(ix) = 0.0; 
                    EN.Ey_1D()(ix) = 0.0;
                    EN.Ez_1D()(ix) = 0.0;
                }
            }
                // }
// - - - - - - - - - - - - - - - - - - - - - -

        } // End the while statement
// - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - -
        // for (size_t ix(0); ix < szx; ++ix)
        // {
        //     Yin.EMF().Ex()(ix+Nbc) = EN.Ex_1D()(ix);
        //     Yin.EMF().Ey()(ix+Nbc) = EN.Ey_1D()(ix);
        //     Yin.EMF().Ez()(ix+Nbc) = EN.Ez_1D()(ix);
        // }

        Yin.EMF().Ex() = EN.Ex_1D();
        Yin.EMF().Ey() = EN.Ey_1D();
        Yin.EMF().Ez() = EN.Ez_1D();

        if ( zeros_in_det > 0) {
                        cout << "WARNING, Det = 0 in "<<zeros_in_det <<"locations" << endl;
                        if ( zeros_in_det > 8) { exit(1);}
        }
    }
// }
//--------------------------------------------------------------


// //--------------------------------------------------------------
//     bool Electric_Field_Methods::Implicit_E_Field::
//     implicitE() const { 
//         return true;
//     } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::FindDE(EMF1D& emf){
//--------------------------------------------------------------
//  Reset DE
//--------------------------------------------------------------
        double Eps(16.0*numeric_limits<double>::epsilon());  
        double LargeEps(sqrt(Eps));  

        // for (size_t ix(0); ix < szx; ++ix)
        // {
        //     DE.Ex_1D()(ix) = emf.Ex()(ix+Nbc);
        //     DE.Ey_1D()(ix) = emf.Ey()(ix+Nbc);
        //     DE.Ez_1D()(ix) = emf.Ez()(ix+Nbc);
        // }
        DE.Ex_1D() = emf.Ex();
        DE.Ey_1D() = emf.Ey();
        DE.Ez_1D() = emf.Ez();
        // Calculate DEx = LargeEps*(|E|)+Eps
        // for (size_t iy(0); iy < szy; ++iy){
        for (size_t ix(0); ix < szx; ++ix){
            DE.Ex_1D()(ix) *=  DE.Ex_1D()(ix);
            DE.Ey_1D()(ix) *=  DE.Ey_1D()(ix);
            DE.Ez_1D()(ix) *=  DE.Ez_1D()(ix);

            DE.Ex_1D()(ix) +=  DE.Ey_1D()(ix);
            DE.Ex_1D()(ix) +=  DE.Ez_1D()(ix);
            DE.Ex_1D()(ix)  =  sqrt(DE.Ex_1D()(ix)); 
            DE.Ex_1D()(ix) *=  LargeEps;
            DE.Ex_1D()(ix) +=  Eps;
        }
        // }
        DE.Ey_1D() = DE.Ex_1D(); 
        DE.Ez_1D() = DE.Ex_1D(); 

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::Ampere(EMF1D& emf){
//--------------------------------------------------------------
//  Use Ampere's law to calculate JN
//--------------------------------------------------------------
        Field1D tmpJi(emf.Bz());

//      Jy = -dBz/dx       
        tmpJi    = emf.Bz(); 
        tmpJi   *= (-1.0) * idx;
        JN.Jy_1D()  = tmpJi.Dx(Input::List().dbydx_order);  

//      Jz += dBy/dx       
        tmpJi    = emf.By(); 
        tmpJi   *=  idx;
        JN.Jz_1D() += tmpJi.Dx(Input::List().dbydx_order);   
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::
    advance(Algorithms::RK2<State2D>* rk, State2D& Yin, collisions_2D& coll, VlasovFunctor2D_implicitE_p2* rkF, const double step_size){//, double time, double dt){
//--------------------------------------------------------------
//  Calculate the implicit electric field
//--------------------------------------------------------------

        int zeros_in_det(1);      // This counts the number of zeros in the determinant 
        int execution_attempt(0); // This counts the number of attempts to find invert the E-field
        State2D Yh(Yin);
        Yh = 0.0;
        
        FindDE(Yin.EMF());                           //  Reset DE

        SHarmonic2D f00(Yin.SH(0,0,0));
        SHarmonic2D f10(Yin.SH(0,1,0)), f11(Yin.SH(0,1,1));
        
// - - - - - - - - - - - - - - - - - - - - - -
        while ( (zeros_in_det > 0) && ( execution_attempt < 10) ) {  // Execute this loop at most twice
            zeros_in_det = 0;                                       // Count the zeros of the determinant
            ++execution_attempt;                                // Count the execusion attempts
// - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - -  
            // Effect of E = 0 on f00, f10, f11
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J0.calculate_J_2D(Yh);
            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;


            // Ex
            Yin.EMF().Ex() = DE.Ex_2D();   // Ex = DEx
            Yin = (*rk)(Yin, dt, rkF, 1);           // Effect of DEx on Y00, Y10, Y11, Y20, Y21, Y22            
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ex.calculate_J_2D(Yh);                   // Evaluate J(DEx)

            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;
            

            // Ey
            Yin.EMF().Ey() = DE.Ey_2D();               // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 2);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ey.calculate_J_2D(Yh);                   // Evaluate J(DEy)           

            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;
            
            // Ez
            Yin.EMF().Ez() = DE.Ez_2D();               // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 3);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22            
            coll.advancef1(Yin,Yh,step_size);                 // Collisions for f10, f11
            J_Ez.calculate_J_2D(Yh);                   // Evaluate J(DEy)

            Yin.SH(0,0,0) = f00;
            Yin.SH(0,1,0) = f10;
            Yin.SH(0,1,1) = f11;
            
            Ampere(Yin.EMF());                           // Calculate JN
                           
// - - - - - - - - - - - - - - - - - - - - - -
            
            //  calculate new EN 

            Array2D< complex<double> > sgm(3,3);     
            valarray< complex<double> > clm(3);

            for (size_t ix(0); ix < szx; ++ix)
            {
                for (size_t iy(0); iy < szy; ++iy)
                {
                    sgm(0,0) =  ( J_Ex.Jx_2D()(ix,iy) - J0.Jx_2D()(ix,iy) ) / DE.Ex_2D()(ix,iy);  // sxx = dJ(E_x)_x/DE_x
                    sgm(0,1) =  ( J_Ey.Jx_2D()(ix,iy) - J0.Jx_2D()(ix,iy) ) / DE.Ey_2D()(ix,iy);  // sxy = dJ(E_y)_x/DE_y
                    sgm(0,2) =  ( J_Ez.Jx_2D()(ix,iy) - J0.Jx_2D()(ix,iy) ) / DE.Ez_2D()(ix,iy);  // sxz = dJ(E_z)_x/DE_z

                    sgm(1,0) =  ( J_Ex.Jy_2D()(ix,iy) - J0.Jy_2D()(ix,iy) ) / DE.Ex_2D()(ix,iy);  // syx = dJ(E_x)_y/DE_x
                    sgm(1,1) =  ( J_Ey.Jy_2D()(ix,iy) - J0.Jy_2D()(ix,iy) ) / DE.Ey_2D()(ix,iy);  // syy = dJ(E_y)_y/DE_y
                    sgm(1,2) =  ( J_Ez.Jy_2D()(ix,iy) - J0.Jy_2D()(ix,iy) ) / DE.Ez_2D()(ix,iy);  // syz = dJ(E_z)_y/DE_z

                    sgm(2,0) =  ( J_Ex.Jz_2D()(ix,iy) - J0.Jz_2D()(ix,iy) ) / DE.Ex_2D()(ix,iy);  // szx = dJ(E_x)_z/DE_x
                    sgm(2,1) =  ( J_Ey.Jz_2D()(ix,iy) - J0.Jz_2D()(ix,iy) ) / DE.Ey_2D()(ix,iy);  // szy = dJ(E_y)_z/DE_y
                    sgm(2,2) =  ( J_Ez.Jz_2D()(ix,iy) - J0.Jz_2D()(ix,iy) ) / DE.Ez_2D()(ix,iy);  // szz = dJ(E_z)_z/DE_z

                    clm[0] =  JN.Jx_2D()(ix,iy) - J0.Jx_2D()(ix,iy);
                    clm[1] =  JN.Jy_2D()(ix,iy) - J0.Jy_2D()(ix,iy);
                    clm[2] =  JN.Jz_2D()(ix,iy) - J0.Jz_2D()(ix,iy);

                // std::cout << "\n\n" << sgm(0,0) << "  ,  " << sgm(0,1) << "  ,  " << sgm(0,2) << "\n";
                // std::cout           << sgm(1,0) << "  ,  " << sgm(1,1) << "  ,  " << sgm(1,2) << "\n";
                // std::cout           << sgm(2,0) << "  ,  " << sgm(2,1) << "  ,  " << sgm(2,2) << "\n";

                
            // Solve the 3 by 3 system of equations
                    complex<double> D_sgm( Det33(sgm) );            // The Determinant of the conductivity tensor
                    if ( abs( D_sgm.real() ) > 6.0*DBL_MIN ) 
                    {
                        EN.Ex_2D()(ix,iy) = Detx33(clm, sgm) / D_sgm;
                        EN.Ey_2D()(ix,iy) = Dety33(clm, sgm) / D_sgm;
                        EN.Ez_2D()(ix,iy) = Detz33(clm, sgm) / D_sgm;
                    }
                    else 
                    {
                        ++zeros_in_det;                                 // Try again with a new perturbation as below
                    
                        DE.Ex_2D()(ix,iy) *= 2;                        // DEx is 2 times larger than before
                        DE.Ey_2D()(ix,iy) *= 1.3;                      // DEy is 1.3 times larger  than before

                        EN.Ex_2D()(ix,iy) = 0.0; 
                        EN.Ey_2D()(ix,iy) = 0.0;
                        EN.Ez_2D()(ix,iy) = 0.0;
                    }
                }
            }
                // }
// - - - - - - - - - - - - - - - - - - - - - -

        } // End the while statement
// - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - -
        // for (size_t ix(0); ix < szx; ++ix)
        // {
        //     for (size_t iy(0); iy < szy; ++iy)
        //     {
        //         Yin.EMF().Ex()(ix+Nbc,iy+Nbc) = EN.Ex_2D()(ix,iy);
        //         Yin.EMF().Ey()(ix+Nbc,iy+Nbc) = EN.Ey_2D()(ix,iy);
        //         Yin.EMF().Ez()(ix+Nbc,iy+Nbc) = EN.Ez_2D()(ix,iy);
        //     }
        // }
        Yin.EMF().Ex() = EN.Ex_2D();
        Yin.EMF().Ey() = EN.Ey_2D();
        Yin.EMF().Ez() = EN.Ez_2D();

        if ( zeros_in_det > 0) {
                        cout << "WARNING, Det = 0 in "<<zeros_in_det <<"locations" << endl;
                        if ( zeros_in_det > 8) { exit(1);}
        }
    }
// }
//--------------------------------------------------------------


// //--------------------------------------------------------------
//     bool Electric_Field_Methods::Implicit_E_Field::
//     implicitE() const { 
//         return true;
//     } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::FindDE(EMF2D& emf){
//--------------------------------------------------------------
//  Reset DE
//--------------------------------------------------------------
        double Eps(16.0*numeric_limits<double>::epsilon());  
        double LargeEps(sqrt(Eps));  

        // for (size_t ix(0); ix < szx; ++ix)
        // {
        //     for (size_t iy(0); iy < szy; ++iy)
        //     {
        //         DE.Ex_2D()(ix,iy) = emf.Ex()(ix+Nbc, iy+Nbc);
        //         DE.Ey_2D()(ix,iy) = emf.Ey()(ix+Nbc, iy+Nbc);
        //         DE.Ez_2D()(ix,iy) = emf.Ez()(ix+Nbc, iy+Nbc);
        //     }
        // }

        DE.Ex_2D() = emf.Ex();
        DE.Ey_2D() = emf.Ey();
        DE.Ez_2D() = emf.Ez();

        // Calculate DEx = LargeEps*(|E|)+Eps
        for (size_t iy(0); iy < szy; ++iy)
        {
            for (size_t ix(0); ix < szx; ++ix)
            {
                DE.Ex_2D()(ix,iy) *=  DE.Ex_2D()(ix,iy);
                DE.Ey_2D()(ix,iy) *=  DE.Ey_2D()(ix,iy);
                DE.Ez_2D()(ix,iy) *=  DE.Ez_2D()(ix,iy);

                DE.Ex_2D()(ix,iy) +=  DE.Ey_2D()(ix,iy);
                DE.Ex_2D()(ix,iy) +=  DE.Ez_2D()(ix,iy);
                DE.Ex_2D()(ix,iy)  =  sqrt(DE.Ex_2D()(ix,iy)); 
                DE.Ex_2D()(ix,iy) *=  LargeEps;
                DE.Ex_2D()(ix,iy) +=  Eps;
            }
        }
        
        DE.Ey_2D() = DE.Ex_2D(); 
        DE.Ez_2D() = DE.Ex_2D(); 

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::Ampere(EMF2D& emf){
//--------------------------------------------------------------
//  Use Ampere's law to calculate JN
//--------------------------------------------------------------
        Field2D tmpJi(emf.Bz());

//      Jx =  dBz/dy       
        // tmpJi    = emf.Bz(); //(same as assignment above)
        tmpJi   *= idy;
        JN.Jx_2D()  = tmpJi.Dy(Input::List().dbydy_order);

//      Jy = -dBz/dx       
        tmpJi    = emf.Bz(); 
        tmpJi   *= (-1.0) * idx;
        JN.Jy_2D()  = tmpJi.Dx(Input::List().dbydx_order);  

//      Jz = -dBx/dy       
        tmpJi    = emf.Bx(); 
        tmpJi   *= (-1.0) * idy;
        JN.Jz_2D()  = tmpJi.Dy(Input::List().dbydy_order);    

//      Jz += dBy/dx       
        tmpJi    = emf.By(); 
        tmpJi   *=  idx;
        JN.Jz_2D() += tmpJi.Dx(Input::List().dbydx_order);   
    }
//--------------------------------------------------------------

//--------------------------------------------------------------

//**************************************************************
