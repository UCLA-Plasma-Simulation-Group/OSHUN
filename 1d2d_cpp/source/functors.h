/*!\brief  Functors for various time-integration methodds - Declarations
 * \author PICKSC
 * \file   functors.h
 *
 * Includes functors for fully explicit, implicit B, implicit E
 *
 */
#ifndef OSHUN1D_FUNCTORS_H
#define OSHUN1D_FUNCTORS_H



//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_explicitE : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_explicitE(vector<size_t> Nl, vector<size_t> Nm,
                              // vector<double> pmax, vector<size_t> Np,
                                vector<valarray<double> > dp,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_explicitE(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, double dt);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection> SA;
    vector<Electric_Field>    EF;
    vector<Current>           JX;
    
    vector<Ampere>    		 AM;
    vector<Magnetic_Field>    BF;
    vector<Faraday>    		 FA;

//            vector<Hydro_Advection>   HA;

};
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_spatialAdvection : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_spatialAdvection(vector<size_t> Nl, vector<size_t> Nm,
                              // vector<double> pmax, vector<size_t> Np,
                                vector<valarray<double> > dp,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_spatialAdvection(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection> SA;
};
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_fieldUpdate : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_fieldUpdate(vector<size_t> Nl, vector<size_t> Nm,
                              // vector<double> pmax, vector<size_t> Np,
                                vector<valarray<double> > dp,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_fieldUpdate(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:

    vector<Current>             JX;
    vector<Ampere>              AM;
    vector<Faraday>             FA;

};
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_momentumAdvection : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_momentumAdvection(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                              double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_momentumAdvection(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Electric_Field>    EF;
    vector<Magnetic_Field>    BF;
};


//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor2D_explicitE : public Algorithms::AbstFunctor<State2D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor2D_explicitE(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                              double xmin, double xmax, size_t Nx,
                              double ymin, double ymax, size_t Ny);
    ~VlasovFunctor2D_explicitE(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State2D& Yin, State2D& Yslope);
    void operator()(const State2D& Yin, State2D& Yslope, size_t dir);

private:
    vector<Spatial_Advection> SA;
    vector<Electric_Field>    EF;
    vector<Current>           JX;
    
    vector<Ampere>           AM;
    vector<Magnetic_Field>    BF;
    vector<Faraday>          FA;

//            vector<Hydro_Advection>   HA;

};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_implicitE_p1 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_implicitE_p1(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                                 double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_implicitE_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);

private:
    vector<Spatial_Advection> SA;
    vector<Magnetic_Field>    BF;
    // vector<Hydro_Advection_1D>   HA;
};
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor2D_implicitE_p1 : public Algorithms::AbstFunctor<State2D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor2D_implicitE_p1(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                                 double xmin, double xmax, size_t Nx,
                                 double ymin, double ymax, size_t Ny);
    ~VlasovFunctor2D_implicitE_p1(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State2D& Yin, State2D& Yslope);
    void operator()(const State2D& Yin, State2D& Yslope, size_t dir);

private:
    vector<Spatial_Advection> SA;
    vector<Magnetic_Field>    BF;
    // vector<Hydro_Advection_1D>   HA;
};

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor1D_implicitE_p2 : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor1D_implicitE_p2(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                                 double xmin, double xmax, size_t Nx);
    ~VlasovFunctor1D_implicitE_p2(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State1D& Yin, State1D& Yslope);
    void operator()(const State1D& Yin, State1D& Yslope, size_t dir);


private:
    vector<Electric_Field>    EF;
    vector<Faraday>           FA;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods
class VlasovFunctor2D_implicitE_p2 : public Algorithms::AbstFunctor<State2D> {
//--------------------------------------------------------------
public:
//          Constructor
    VlasovFunctor2D_implicitE_p2(vector<size_t> Nl, vector<size_t> Nm,
                                vector<valarray<double> > dp,
                                 double xmin, double xmax, size_t Nx,
                                 double ymin, double ymax, size_t Ny);
    ~VlasovFunctor2D_implicitE_p2(){ };

//          Collect all the operators and apply on Yin
    void operator()(const State2D& Yin, State2D& Yslope);
    void operator()(const State2D& Yin, State2D& Yslope, size_t dir);


private:
    vector<Electric_Field>    EF;
    vector<Faraday>           FA;
};
//--------------------------------------------------------------

#endif //OSHUN1D_FUNCTORS_H


/** @} */
