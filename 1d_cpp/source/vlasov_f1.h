/*!\brief  Vlasov Equation - Declarations
* \author PICKSC
 * \date   October 10, 2016
 * \file   vlasov_f1.h
 *
 * Includes declarations for spatial advection, electric field advection, current, and 
 * the RK1D Functor responsible for stepping forward
 * 
 */

#ifndef DECL_VLASOVMAXWELLF1_H
#define DECL_VLASOVMAXWELLF1_H

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
class Spatial_Advection_1D_f1 {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Spatial_Advection_1D_f1(size_t Nl, size_t Nm,
                            double pmin, double pmax, size_t Np,
                            double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din, DistFunc1D& Dh);

private:
    SHarmonic1D                   	fd1;
//            Array2D< complex<double> >  	A1, A2;
    valarray< complex<double> >  	vr;
    complex<double>                 A00, A10, A20;

};
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Electric field
class Electric_Field_1D_f1 {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Electric_Field_1D_f1(size_t Nl, size_t Nm,
                         double pmin, double pmax, size_t Np,
                         double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din,
                    const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                    DistFunc1D& Dh);
    void Implicit_Ex(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh);
    void Implicit_Ey(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh);
    void Implicit_Ez(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh);

private:
    void MakeG00(SHarmonic1D& f);
    void MakeGH( SHarmonic1D& f, size_t l);

    SHarmonic1D H, G, TMP;

    complex<double>                 A100, C100, A210, B211, C311, A310;

    //    Array2D< complex<double> >   A1, A2;
//    valarray< complex<double> >  B1, B2;
//    valarray< complex<double> >  C1, C3;
//    Array2D< complex<double> >   C2, C4;
    valarray< complex<double> >  pr, invpr, Hp0;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Magnetic field
class Magnetic_Field_1D_f1 {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Magnetic_Field_1D_f1(size_t Nl, size_t Nm,
                         double pmin, double pmax, size_t Np,
                         double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din,
                    const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                    DistFunc1D& Dh);
    void implicit(DistFunc1D& Din,
                  const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                  double dt);

private:

    SHarmonic1D FLM;

    valarray< complex<double> >		A1, B1;
    Array2D< complex<double> > 		A2;
    complex<double> 				A3;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
class VlasovFunctor1D_f1_explicitE : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------    	
public:
//          Constructor
    VlasovFunctor1D_f1_explicitE(vector<size_t> Nl, vector<size_t> Nm,
                                 vector<double> pmax, vector<size_t> Np,
                                 double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_f1_explicitE(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);
    // void implicit_rest(const State1D& Yin, State1D& Yslope);
    // void implicit_E(const State1D& Yin, State1D& Yslope, size_t dir);



private:
    vector<Spatial_Advection_1D_f1> SA;
    vector<Electric_Field_1D_f1>    EF;
    vector<Current_1D>           JX;
    vector<Ampere_1D>    		 AM;
    vector<Magnetic_Field_1D_f1>    BF;
    vector<Faraday_1D>    		 FA;
//            vector<Hydro_Advection_1D>   HA;

};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_f1_explicitEB : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_f1_explicitEB(vector<size_t> Nl, vector<size_t> Nm,
                                  vector<double> pmax, vector<size_t> Np,
                                  double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_f1_explicitEB(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);
    // void implicit_rest(const State1D& Yin, State1D& Yslope);
    // void implicit_E(const State1D& Yin, State1D& Yslope, size_t dir);



private:
    vector<Spatial_Advection_1D_f1> SA;
    vector<Electric_Field_1D_f1>    EF;
    vector<Current_1D>           JX;
    vector<Ampere_1D>    		 AM;
//    vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>    		 FA;
//            vector<Hydro_Advection_1D>   HA;

};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
class VlasovFunctor1D_f1_implicitE_p1 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------        
public:
//          Constructor
    VlasovFunctor1D_f1_implicitE_p1(vector<size_t> Nl, vector<size_t> Nm,
                                    vector<double> pmax, vector<size_t> Np,
                                    double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_f1_implicitE_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection_1D_f1> SA;
    // vector<Electric_Field_1D>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
    vector<Magnetic_Field_1D_f1>    BF;
//            vector<Faraday_1D>           FA;
//            vector<Hydro_Advection_1D>   HA;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_f1_implicitEB_p1 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_f1_implicitEB_p1(vector<size_t> Nl, vector<size_t> Nm,
                                     vector<double> pmax, vector<size_t> Np,
                                     double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_f1_implicitEB_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection_1D_f1> SA;
    // vector<Electric_Field_1D>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
//    vector<Magnetic_Field_1D>    BF;
//            vector<Faraday_1D>           FA;
//            vector<Hydro_Advection_1D>   HA;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
class VlasovFunctor1D_f1_implicitE_p2 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------        
public:
//          Constructor
    VlasovFunctor1D_f1_implicitE_p2(vector<size_t> Nl, vector<size_t> Nm,
                                    vector<double> pmax, vector<size_t> Np,
                                    double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_f1_implicitE_p2(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);


private:
    // vector<Spatial_Advection_1D> SA;
    vector<Electric_Field_1D_f1>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
    // vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>           FA;
};
//--------------------------------------------------------------



/** @} */


#endif
