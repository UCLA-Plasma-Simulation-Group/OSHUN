/*!\brief  Functors for various time-integration methodds - Declarations
 * \author PICKSC
 * \file   vlasov.h
 *
 * Includes functors for fully explicit, implicit B, implicit E
 *
 */
#ifndef OSHUN1D_FUNCTORS_H
#define OSHUN1D_FUNCTORS_H

#endif //OSHUN1D_FUNCTORS_H


//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_explicitE : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_explicitE(vector<size_t> Nl, vector<size_t> Nm,
                              vector<double> pmax, vector<size_t> Np,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_explicitE(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);
    // void implicit_rest(const State1D& Yin, State1D& Yslope);
    // void implicit_E(const State1D& Yin, State1D& Yslope, size_t dir);



private:
    vector<Spatial_Advection_1D> SA;
    vector<Electric_Field_1D>    EF;
    vector<Current_1D>           JX;
    vector<Ampere_1D>    		 AM;
    vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>    		 FA;
//            vector<Hydro_Advection_1D>   HA;

};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_spatialpush : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_spatialpush(vector<size_t> Nl, vector<size_t> Nm,
                              vector<double> pmax, vector<size_t> Np,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_spatialpush(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection_1D> SA;

};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_momentumpush : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_momentumpush(vector<size_t> Nl, vector<size_t> Nm,
                              vector<double> pmax, vector<size_t> Np,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_momentumpush(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:

    vector<Electric_Field_1D>    EF;
    vector<Current_1D>           JX;
    vector<Ampere_1D>    		 AM;
    vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>    		 FA;
//            vector<Hydro_Advection_1D>   HA;

};
//--------------------------------------------------------------



//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_explicitE_implicitB : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_explicitE_implicitB(vector<size_t> Nl, vector<size_t> Nm,
                               vector<double> pmax, vector<size_t> Np,
                               double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_explicitE_implicitB(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);
    // void implicit_rest(const State1D& Yin, State1D& Yslope);
    // void implicit_E(const State1D& Yin, State1D& Yslope, size_t dir);



private:
    vector<Spatial_Advection_1D> SA;
    vector<Electric_Field_1D>    EF;
    vector<Current_1D>           JX;
    vector<Ampere_1D>    		 AM;
//    vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>    		 FA;
//            vector<Hydro_Advection_1D>   HA;

};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_implicitE_p1 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_implicitE_p1(vector<size_t> Nl, vector<size_t> Nm,
                                 vector<double> pmax, vector<size_t> Np,
                                 double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_implicitE_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection_1D> SA;
    // vector<Electric_Field_1D>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
    vector<Magnetic_Field_1D>    BF;
//            vector<Faraday_1D>           FA;
    vector<Hydro_Advection_1D>   HA;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_implicitE_implicitB_p1 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_implicitE_implicitB_p1(vector<size_t> Nl, vector<size_t> Nm,
                                  vector<double> pmax, vector<size_t> Np,
                                  double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_implicitE_implicitB_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection_1D> SA;
    // vector<Electric_Field_1D>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
    vector<Magnetic_Field_1D>    BF;
//            vector<Faraday_1D>           FA;
    vector<Hydro_Advection_1D>   HA;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_implicitE_p2 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_implicitE_p2(vector<size_t> Nl, vector<size_t> Nm,
                                 vector<double> pmax, vector<size_t> Np,
                                 double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_implicitE_p2(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, const State1D& Y2in, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);


private:
    // vector<Spatial_Advection_1D> SA;
    vector<Electric_Field_1D>    EF;
    // vector<Current_1D>           JX;
    // vector<Ampere_1D>            AM;
    // vector<Magnetic_Field_1D>    BF;
    vector<Faraday_1D>           FA;
};
//--------------------------------------------------------------



/** @} */