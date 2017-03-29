/*! \brief Implicit Electric Field - Definitions
* \author PICKSC
 * \date   September 23, 2016
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
    #include "fluid.h"
    #include "vlasov.h"
    #include "functors.h"
    #include "formulary.h"
    #include "nmethods.h"
    #include "collisions.h"
    
    #include "implicitE.h"
//**************************************************************


// vector<float>   vfloat(const vector<complex<double> > vDouble); 
// vector<float>   vfloat_complex(const vector<complex<double> > vDouble); 
//**************************************************************
//**************************************************************
//*******THE CLASS FOR THE CURRENT DEFINITION ******************
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Current_xyz::Current_xyz(EMF1D& emf) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : jayx(emf.Ex()), 
         jayy(emf.Ex()), 
         jayz(emf.Ex()), 
         // p3og(Algorithms::MakeAxis(
         //    static_cast< complex<double> >(Input::List().pmin), 
         //    static_cast< complex<double> >(Input::List().pmaxs[0]),
         //    Input::List().ps[0])),
         // Delta_p(static_cast< complex<double> >((Input::List().pmaxs[0]-Input::List().pmin)/Input::List().ps[0])),
         // small(static_cast< complex<double> >(Input::List().pmin)),
         tempx(0.0,emf.Ex().numx()),
         tempy(0.0,emf.Ex().numx()),
         tempz(0.0,emf.Ex().numx()),
         fourpioverthree(4.0*M_PI/3.0,0.0)
         {
            Nbc = Input::List().BoundaryCells;
            szx = Input::List().NxLocal[0];
            // for (size_t s(0); s<Yin.Species();++s)
            // {
            //  for (size_t i(0); i<p3og.size(); ++i) { // calculate p^3/g
            //     p3og[i] = (p3og[i]*p3og[i])*(p3og[i]/sqrt(1.0+p3og[i]*p3og[i])); 
            //  }
            //     p3ogv.push_back(p3og);
            // }
            //  small *= small; small *= small; small *= Input::List().pmin; 
            //  small *= 0.2; small *= 1.0/(Input::List().pmin+Delta_p); 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::J(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Jx(); 
            break;
        }
        case 2: { 
            return Jy(); 
            break;
        }
        case 3: { 
            return Jz(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Jx();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jx() {
       return jayx; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jy() {
        return jayy;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Current_xyz::Jz() {
        return jayz;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_J(State1D& Yin) {
//--------------------------------------------------------------
//  Update the total current 
//--------------------------------------------------------------
        complex<double> c01(0.0,1.0);
        complex<double> pmin, pmax;
        size_t nump;
        Array2D<double> current(3,Yin.EMF().Ex().numx());

        // jayx=0.0;jayy=0.0;jayz=0.0;
        jayx = static_cast<complex<double> >(0.0);
        jayy = static_cast<complex<double> >(0.0);
        jayz = static_cast<complex<double> >(0.0);


        for (size_t s(0); s < Yin.Species(); ++s){  
            current = Yin.DF(s).getcurrent();      
            for (size_t ix(0); ix < Yin.EMF().Ex().numx(); ++ix)
            {
                jayx(ix) +=   static_cast<complex<double>>(current(0,ix));
                jayy(ix) +=   static_cast<complex<double>>(current(1,ix));
                jayz(ix) +=   static_cast<complex<double>>(current(2,ix));

            }
            
            // jayx.array() += static_cast<complex<double>>(Yin.DF(s).getcurrent(0)); 
            // jayy.array() += static_cast<complex<double>>(Yin.DF(s).getcurrent(1)); 
            // jayz.array() += static_cast<complex<double>>(Yin.DF(s).getcurrent(2)); 
        }
}    
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Efield_xyz::Efield_xyz(EMF1D& emf) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : efieldx(emf.Ex()), 
         efieldy(emf.Ey()), 
         efieldz(emf.Ez()) {
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::E(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Ex(); 
            break;
        }
        case 2: { 
            return Ey(); 
            break;
        }
        case 3: { 
            return Ez(); 
            break;
        }
        default: {
                cout << "There is no such component for the current!" << endl;
                exit(1);
                break;
        }
    }
    return Ex();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ex() {
       return efieldx; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ey() {
        return efieldy;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field1D& Electric_Field_Methods::Efield_xyz::Ez() {
        return efieldz;
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
    Implicit_E_Field(EMF1D& emf, const double& deltat, const double& deltax)://, int tout_start)
//--------------------------------------------------------------
//  Constructor for the implicit electric field
//--------------------------------------------------------------
         JN(emf),   J0(emf),                     // New current, Current from E = 0
         J_Ex(emf), J_Ey(emf), J_Ez(emf),        // Current due to the effect of Ex, Ey, Ez
         EN(emf),   E0(emf),   DE(emf),          // New E, old E, perturbation E
   
         dt(deltat),
         idx(static_cast< complex<double> >(0.5/deltax))
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            Nbc = Input::List().BoundaryCells;
            szx = Input::List().NxLocal[0];
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::
    advance(Algorithms::RK3<State1D>* rk, State1D& Yin, collisions& coll, VlasovFunctor1D_implicitE_p2* rkF){//, double time, double dt){
//--------------------------------------------------------------
//  Calculate the implicit electric field
//--------------------------------------------------------------

        int zeros_in_det(1);      // This counts the number of zeros in the determinant 
        int execution_attempt(0); // This counts the number of attempts to find invert the E-field
        State1D Yh(Yin);
        Yh = 0.0;

        // vector<SHarmonic1D> f00;
        // vector<SHarmonic1D> f10;vector<SHarmonic1D> f11;
        // vector<SHarmonic1D> f20;vector<SHarmonic1D> f21;vector<SHarmonic1D> f22;

        // for (size_t s(0); s < Yin.Species(); ++s){
        //     f00.push_back(Yin.SH(s,0,0));
        //     f10.push_back(Yin.SH(s,1,0));f11.push_back(Yin.SH(s,1,1));
        //     f20.push_back(Yin.SH(s,2,0));f21.push_back(Yin.SH(s,2,1));f22.push_back(Yin.SH(s,2,2));
        // }

        FindDE(Yin.EMF());                           //  Reset DE

// - - - - - - - - - - - - - - - - - - - - - -
        while ( (zeros_in_det > 0) && ( execution_attempt < 4) ) {  // Execute this loop at most twice
            zeros_in_det = 0;                                       // Count the zeros of the determinant
            ++execution_attempt;                                // Count the execusion attempts
// - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - -  
                                                        // Effect of E = 0 on f00, f10, f11
            
            coll.advancef1(Yin,Yh);                 // Collisions for f10, f11
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,1,0) = Yh.SH(s,1,0);
            //     Yin.SH(s,1,1) = Yh.SH(s,1,0);
            // }   
            J0.calculate_J(Yh);
            
            
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,0,0) = f00[s];
            //     Yin.SH(s,1,0) = f10[s];Yin.SH(s,1,1) = f11[s];
            // // Yin.SH(s,2,0) = f20[s];Yin.SH(s,2,1) = f21[s];Yin.SH(s,2,2) = f22[s];
            // }       

            Yin.EMF().Ex() = DE.Ex();             // Ex = DEx
            Yin = (*rk)(Yin, dt, rkF, 1);           // Effect of DEx on Y00, Y10, Y11, Y20, Y21, Y22            
            coll.advancef1(Yin,Yh);                 // Collisions for f10, f11
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,1,0) = Yh.SH(s,1,0);
            //     Yin.SH(s,1,1) = Yh.SH(s,1,0);
            // }                    // Collisions for f10, f11
            J_Ex.calculate_J(Yh);                // Evaluate J(DEx)
            
            // Ytemp = Yin;
            // for (size_t s(0); s < Yin.Species(); ++s){
            // Yin.SH(s,0,0) = f00[s];
            // Yin.SH(s,1,0) = f10[s];Yin.SH(s,1,1) = f11[s];
            // Yin.SH(s,2,0) = f20[s];Yin.SH(s,2,1) = f21[s];Yin.SH(s,2,2) = f22[s];
            // }       
           
                // 
            Yin.EMF().Ey() = DE.Ey();             // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 2);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22
            coll.advancef1(Yin,Yh);                 // Collisions for f10, f11
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,1,0) = Yh.SH(s,1,0);
            //     Yin.SH(s,1,1) = Yh.SH(s,1,0);
            // }                    // Collisions for f10, f11
            J_Ey.calculate_J(Yh);                 // Evaluate J(DEy)
            
            // Ytemp = Yin;
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,0,0) = f00[s];
            //     Yin.SH(s,1,0) = f10[s];Yin.SH(s,1,1) = f11[s];
            //     Yin.SH(s,2,0) = f20[s];Yin.SH(s,2,1) = f21[s];Yin.SH(s,2,2) = f22[s];
            // }       
// - - - - - - - - - - - - - - - - - - - - - -

           
            Yin.EMF().Ez() = DE.Ez();             // Ey = DEy
            Yin = (*rk)(Yin, dt, rkF, 3);           // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22
            
            coll.advancef1(Yin,Yh);                 // Collisions for f10, f11
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,1,0) = Yh.SH(s,1,0);
            //     Yin.SH(s,1,1) = Yh.SH(s,1,0);
            // }   

            J_Ez.calculate_J(Yh); //exit(1);               // Evaluate J(DEy)
            
            // Ytemp = Yin;
            // for (size_t s(0); s < Yin.Species(); ++s){
            //     Yin.SH(s,0,0) = f00[s];
            //     Yin.SH(s,1,0) = f10[s];Yin.SH(s,1,1) = f11[s];
            //     Yin.SH(s,2,0) = f20[s];Yin.SH(s,2,1) = f21[s];Yin.SH(s,2,2) = f22[s];
            // }    
// - - - - - - - - - - - - - - - - - - - - - -

                
            Ampere(Yin.EMF());                           // Calculate JN
                
            
// - - - - - - - - - - - - - - - - - - - - - -
            
            //  calculate new EN 

            Array2D< complex<double> > sgm(3,3);     
            valarray< complex<double> > clm(3);

            
                // for (size_t iy(0); iy < szy; ++iy){
                    for (size_t ix(0); ix < szx; ++ix){
                        sgm(0,0) =  ( J_Ex.Jx()(ix) - J0.Jx()(ix) ) / DE.Ex()(ix);  // sxx = dJ(E_x)_x/DE_x
                        sgm(0,1) =  ( J_Ey.Jx()(ix) - J0.Jx()(ix) ) / DE.Ey()(ix);  // sxy = dJ(E_y)_x/DE_y
                        sgm(0,2) =  ( J_Ez.Jx()(ix) - J0.Jx()(ix) ) / DE.Ez()(ix);  // sxz = dJ(E_z)_x/DE_z

                        sgm(1,0) =  ( J_Ex.Jy()(ix) - J0.Jy()(ix) ) / DE.Ex()(ix);  // syx = dJ(E_x)_y/DE_x
                        sgm(1,1) =  ( J_Ey.Jy()(ix) - J0.Jy()(ix) ) / DE.Ey()(ix);  // syy = dJ(E_y)_y/DE_y
                        sgm(1,2) =  ( J_Ez.Jy()(ix) - J0.Jy()(ix) ) / DE.Ez()(ix);  // syz = dJ(E_z)_y/DE_z

                        sgm(2,0) =  ( J_Ex.Jz()(ix) - J0.Jz()(ix) ) / DE.Ex()(ix);  // szx = dJ(E_x)_z/DE_x
                        sgm(2,1) =  ( J_Ey.Jz()(ix) - J0.Jz()(ix) ) / DE.Ey()(ix);  // szy = dJ(E_y)_z/DE_y
                        sgm(2,2) =  ( J_Ez.Jz()(ix) - J0.Jz()(ix) ) / DE.Ez()(ix);  // szz = dJ(E_z)_z/DE_z
  
                        clm[0] =  JN.Jx()(ix) - J0.Jx()(ix);
                        clm[1] =  JN.Jy()(ix) - J0.Jy()(ix);
                        clm[2] =  JN.Jz()(ix) - J0.Jz()(ix);

                        

                        // std::cout << "\n\n" << sgm(0,0) << "  ,  " << sgm(0,1) << "  ,  " << sgm(0,2) << "\n";
                        // std::cout           << sgm(1,0) << "  ,  " << sgm(1,1) << "  ,  " << sgm(1,2) << "\n";
                        // std::cout           << sgm(2,0) << "  ,  " << sgm(2,1) << "  ,  " << sgm(2,2) << "\n";

                        
                    // Solve the 3 by 3 system of equations
                        complex<double> D_sgm( Det33(sgm) );            // The Determinant of the conductivity tensor
                        if ( abs( D_sgm.real() ) > 6.0*DBL_MIN ) {
                            EN.Ex()(ix) = Detx33(clm, sgm) / D_sgm;
                            EN.Ey()(ix) = Dety33(clm, sgm) / D_sgm;
                            EN.Ez()(ix) = Detz33(clm, sgm) / D_sgm;
                        }
                        else {
                            ++zeros_in_det;                                 // Try again with a new perturbation as below
                        // exit(1);
                            DE.Ex()(ix) *= 2;                        // DEx is 2 times larger than before
                            DE.Ey()(ix) *= 1.3;                      // DEy is 1.3 times larger  than before

                            EN.Ex()(ix) = 0.0; 
                            EN.Ey()(ix) = 0.0;
                            EN.Ez()(ix) = 0.0;
                        }
                    }
                // }
// - - - - - - - - - - - - - - - - - - - - - -

        } // End the while statement
// - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - -
        Yin.EMF().Ex() = EN.Ex();
        Yin.EMF().Ey() = EN.Ey();
        Yin.EMF().Ez() = EN.Ez();

        if ( zeros_in_det > 0) {
                        cout << "WARNING, Det = 0 in "<<zeros_in_det <<"locations" << endl;
                        if ( zeros_in_det > 8) { exit(1);}
        }

        // return Ytemp;
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

        DE.Ex() = emf.Ex();
        DE.Ey() = emf.Ey();
        DE.Ez() = emf.Ez();

        // Calculate DEx = LargeEps*(|E|)+Eps
        // for (size_t iy(0); iy < szy; ++iy){
            for (size_t ix(0); ix < szx; ++ix){
                DE.Ex()(ix) *=  DE.Ex()(ix);
                DE.Ey()(ix) *=  DE.Ey()(ix);
                DE.Ez()(ix) *=  DE.Ez()(ix);

                DE.Ex()(ix) +=  DE.Ey()(ix);
                DE.Ex()(ix) +=  DE.Ez()(ix);
                DE.Ex()(ix)  =  sqrt(DE.Ex()(ix)); 
                DE.Ex()(ix) *=  LargeEps;
                DE.Ex()(ix) +=  Eps;
            }
        // }
        DE.Ey() = DE.Ex(); 
        DE.Ez() = DE.Ex(); 

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::Ampere(EMF1D& emf){
//--------------------------------------------------------------
//  Use Ampere's law to calculate JN
//--------------------------------------------------------------
        Field1D tmpJi(emf.Bz());

//      Jx =  dBz/dy       
//      tmpJi    = Yin.EMF().Bz(); (same as assignment above)
        // tmpJi   *= idy;
        // JN.Jx()  = tmpJi.Dy();

//      Jy = -dBz/dx       
        tmpJi    = emf.Bz(); 
        tmpJi   *= (-1.0) * idx;
        tmpJi.Dx();
        // tmpJi(tmpJi.numx()-1) = 0.0;
        JN.Jy()  = tmpJi;        

//      Jz = -dBx/dy       
        // tmpJi    = Y.EMF().Bx(); 
        // tmpJi   *= (-1.0) * idy;
        // JN.Jz()  = tmpJi.Dy();    

//      Jz += dBy/dx       
        tmpJi    = emf.By(); 
        tmpJi   *=  idx;
        tmpJi.Dx();
        // tmpJi(tmpJi.numx()-1) = 0.0;
        JN.Jz() += tmpJi;   
    }
//--------------------------------------------------------------

//--------------------------------------------------------------

//**************************************************************
