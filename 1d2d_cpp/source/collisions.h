/*! \brief Collisions - Declarations
 * \author PICKSC
 * \date   2017
 * \file   collisions.h
 * 
 * Contains:
 * 1) the explicit energy-conserving algorithm for effect of 
 * electron-electron collisions on f_00
 * 2) implicit algorithm used for e-e + e-i collisions on f_lm
 * And all their containers.
 * 
 */

#ifndef DECLARATION_FOKKER_PLANCK_H
#define DECLARATION_FOKKER_PLANCK_H

//-------------------------------------------------------------------
/** \addtogroup vfp1d
 *  @{
 */
class self_f00_implicit_step {
//-------------------------------------------------------------------
private:

    double mass;
    double dt;
    bool   ib;
    ///     Define the velocity axis
    valarray<double>  vr;
    valarray<double>  dvr;
    valarray<double>  vrh;
    valarray<double>  dtoverv2;

    ///     Various coefficients for the integrals
    valarray<double>  p2dp, p2dpm1, phdp, phdpm1, p4dp, laser_Inv_Uav6;

    ///     Rosenbluth Potentials
    valarray<double>  C_RB, D_RB;
    double            I4_Lnee;

    ///     Chang-Cooper weighting delta
    valarray<double>  delta_CC;

    ///     Constants
    double c_kpre;
    double vw_coeff_cube;

    Formulary formulas;

    void   update_C_Rosenbluth(valarray<double>& fin);
    double update_D_Rosenbluth(const size_t& k, valarray<double>& fin, const double& delta);
    void   update_D_and_delta(valarray<double>& fin);
    void   update_D_inversebremsstrahlung(const double& Z0, const double& heatingcoefficient, const double& vos);
    double calc_delta_ChangCooper(const size_t& k, const double& C, const double& D);

public:

    // self_f00_implicit_step(const size_t &nump, const double &pmax, const double &_mass, const double &_deltat, bool& _ib);
    self_f00_implicit_step(const valarray<double> &dp, const double &_mass, const double &_deltat, bool& _ib);

    void takestep(valarray<double> &fin, valarray<double> &fh, const double& Z0, const double& heating, const double& cooling);
    void getleftside(valarray<double> &fin, const double& Z0, const double& heating, const double& cooling, Array2D<double> &LHStemp);
};

//-------------------------------------------------------------------
/**
* \class self_f00_explicit_step
* \brief The top-level container for collisions on l=0.
*
*   The self_f00_explicit_step class describes the object that
*   contains the RK4 algorithm that sets up the explicit solve step.
*   self_f00_explicit_step has a "loop()" function that pulls the relevant
*   distribution function and passes it to the RK4 solver.
*/
class self_f00_implicit_collisions {
//--------------------------------------------------------------
public:
    /// Constructors/Destructors
    self_f00_implicit_collisions(const valarray<double>& dp, const double& charge, const double& mass, const double& deltat);

    /// This loop calls the RK4_f00 private member that is
    /// responsible for setting up the RK4 algorithm to advance
    /// the collision step.
    void loop(SHarmonic1D& SHin, valarray<double>& Zarray, const double time, SHarmonic1D& SHout);
    void loop(SHarmonic2D& SHin, Array2D<double>& Zarray, const double time, SHarmonic2D& SHout);

private:
    //  Variables
    valarray<double>            fin, fout;
    valarray<double>            xgrid, ygrid;
    bool                        ib;
    self_f00_implicit_step      collide;

    ///     Switches for inverse bremsstrahlung and maxwellian cooling
    bool                        IB_heating;
    bool                        MX_cooling;

    valarray<double>            heatingprofile_1d;
    valarray<double>            coolingprofile_1d;

    Array2D<double>            heatingprofile_2d;
    Array2D<double>            coolingprofile_2d;


    size_t                         Nbc; ///< Number of boundary cells in each direction
    size_t                         szx,szy; ///< Total cells including boundary cells in x-direction

};
//--------------------------------------------------------------

//-------------------------------------------------------------------
/** \class  f00_step
 *  \brief  Collisions for l=0
 *  
 *   This class contains all the algebra necessary for the explicit, 
 *   energy-conserving algorithm that describes the effect of 
 *   collisions on the isotropic component of the distribution. An 
 *   object of this class is only constructed once per RK4_f00 object. 
 *   The RK4_f00 object is constructed once per Collide object. This 
 *   means that this is at the base of the heirarchy for collisions 
 *   affecting the isotropic component.
*/        
class self_f00_explicit_step {
//-------------------------------------------------------------------
    private:
    ///     Define the velocity axis
        valarray<double>  vr;

    ///     Various coefficients for the integrals
        valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;

    ///     The integrals
        valarray<double>  J1, I2, I4, U3, Qn, Pn;

    ///     Constants
        double c_kpre;
        size_t    NB;
        Formulary formulas;

        double G(const int& n, const valarray<double>& fin);

    public:    
        self_f00_explicit_step(const valarray<double>& dp); 

        void takestep(const valarray<double>& fin, valarray<double>& fh);    
};
//--------------------------------------------------------------
/** @} */ 

//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
    class self_f00_RKfunctor : public Algorithms::AbstFunctor<valarray<double> > {
//--------------------------------------------------------------        
        public:
//          Constructor
            self_f00_RKfunctor(const valarray<double>& dp);
            ~self_f00_RKfunctor(){ };

//          Collect all the operators and apply on Yin
            void operator()(const valarray<double>& fin, valarray<double>& fslope);
            void operator()(const valarray<double>& fin, valarray<double>& fslope, size_t dir);
        private:
            //          Variables
            self_f00_explicit_step collide;
             
    };
//--------------------------------------------------------------    



/** \addtogroup vfp1d
 *  @{
 */
//-------------------------------------------------------------------
/** 
* \class self_f00_explicit_step
* \brief The top-level container for collisions on l=0.
* 
*   The self_f00_explicit_step class describes the object that 
*   contains the RK4 algorithm that sets up the explicit solve step. 
*   self_f00_explicit_step has a "loop()" function that pulls the relevant
*   distribution function and passes it to the RK4 solver.
*/        
        class self_f00_explicit_collisions {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            self_f00_explicit_collisions(const valarray<double>& dp, const double& deltat); //DistFunc1D& DFin,

        /// This loop calls the RK4_f00 private member that is
        /// responsible for setting up the RK4 algorithm to advance
        /// the collision step. 
            void loop(SHarmonic1D& SHin, SHarmonic1D& SHout);
            void loop(SHarmonic2D& SHin, SHarmonic2D& SHout);

        private:
        //  Variables
            valarray<double>                        fin; //, fout;


            /// This object contains the RK4 algorithm that advances 
            /// the collision operator. Inside of it is the Collide object
            /// that contains all the relevant collision integral algebra.
            
            Algorithms::RK4<valarray<double> >       RK;

            self_f00_RKfunctor                      rkf00;

            size_t                                  num_h;
            double                                  h;

            size_t Nbc; ///< Number of boundary cells in each direction
            size_t szx, szy; ///< Total cells including boundary cells in x-direction

        };
//--------------------------------------------------------------

//--------------------------------------------------------------

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
/** 
 * \class self_flm_implicit_step
 * \brief Collisions for l >= 1.
 * 
 *   This class forms the numerical backbone for quantifying
 *   the effect of collisions on the anisotropic components of
 *   the distribution function.
*/   
//--------------------------------------------------------------
        class self_flm_implicit_step {
//-------------------------------------------------------------------
        private:

//          Define the powers of the velocity
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;

//          Variables for the matrix A, such that A*x = b
            // valarray<double> Scattering_Term, df0, ddf0;
            // double _LOGee;
            // Array2D<double> Alpha_Tri;

            bool if_tridiagonal;
            double mass;

//          Constant

            double Dt, kpre;

            vector<double>              _LOGee_x;
            vector<valarray<double> >   Scattering_Term_x; 
            vector<Array2D<double>  >   Alpha_Tri_x; 
            vector<valarray<double> >   df0_x, ddf0_x;
            
            Formulary formulas;

        public:
//          Constructors/Destructors
            // self_flm_implicit_step(double pmax, size_t nump, double mass); 
            self_flm_implicit_step(const size_t numxtotal, valarray<double> dp); 
         
//          Calculate the coefficients
            void reset_coeff(const valarray<double>& f00, double Zvalue, const double Delta_t, size_t position);

//          Implicit Advance
            void advance(valarray<complex<double> >& fin, const int el, size_t position);    
        };
//-------------------------------------------------------------------
/** @} */ 

/** \addtogroup vfp1d
 *  @{
 */
//-------------------------------------------------------------------
/** 
 * \class Anisotropic_Collisions
 * \brief Middle container for collisions on  l >= 1.
 * 
 *   The Anisotropic_Collisions class describes the object that 
 *   contains the implicit solve step. Note that the main procedures
 *   return a new state. These procedures extract the data from the relevant
 *   harmonic and send it to the implicit_step object.
*/
        class self_flm_implicit_collisions {
        public:
///          Constructors/Destructors
            self_flm_implicit_collisions(const size_t _l0, const size_t _m0,
                             const valarray<double>& dp, const double& deltat); 

            void advancef1(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh);
            void advanceflm(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh);

            void advancef1(DistFunc2D& DF, Array2D<double>& Zarray, DistFunc2D& DFh);
            void advanceflm(DistFunc2D& DF, Array2D<double>& Zarray, DistFunc2D& DFh);

        private:

            bool if_tridiagonal;

            size_t l0; ///< Number of m harmonics.
            size_t m0; ///< Number of m harmonics.
            size_t Nbc; ///< Number of boundary cells in each direction
            size_t szx; ///< Total cells including boundary cells in x-direction
            size_t szy;
            
///         The object that is responsible for performing the algebra required for the integrals.
            self_flm_implicit_step  implicit_step;
            double Dt;

            Formulary formulas;
/// ---------------------------------------------------------------------------- ///
            size_t f1_m_upperlimit;

        };
//-------------------------------------------------------------------

//-------------------------------------------------------------------
/** 
* \class self_collisions
* \brief The top-level container for self collisions over all species on l=0.
* 
*   
*/        
        class self_collisions {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            self_collisions(const size_t _l0, const size_t _m0,
                             const valarray<double>& dp, const double& charge, const double& mass, const double& deltat); 
            // ~self_collisions();
            void advancef00(SHarmonic1D& f00, valarray<double>& Zarray, const double time, SHarmonic1D& f00h);
            void advancef1(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh);
            void advanceflm(DistFunc1D& DF, valarray<double>& Zarray, DistFunc1D& DFh);

            void advancef00(SHarmonic2D& f00, Array2D<double>& Zarray, const double time, SHarmonic2D& f00h);
            void advancef1(DistFunc2D& DF, Array2D<double>& Zarray, DistFunc2D& DFh);
            void advanceflm(DistFunc2D& DF, Array2D<double>& Zarray, DistFunc2D& DFh);


        private:
        //  Variables
            self_f00_explicit_collisions self_f00_exp_collisions;
            self_f00_implicit_collisions self_f00_imp_collisions;
            self_flm_implicit_collisions self_flm_imp_collisions;
        };

//-------------------------------------------------------------------
/** 
* \class self_collisions
* \brief The top-level container for self collisions over all species on l=0.
* 
*   
*/        
        class collisions_1D {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            collisions_1D(const State1D& Yin, const double& deltat); 
            // ~self_collisions();
            void advance(State1D& Y, const Clock& W);
            void advancef0(State1D& Y, const Clock& W, State1D& Yh);
            void advancef1(State1D& Y, State1D& Yh);
            void advanceflm(State1D& Y, State1D& Yh);

            vector<self_collisions> self();
            // void advancef1(State1D& Y);
            // void advanceflm(State1D& Y);

        private:
        //  Variables
            State1D Yh;
            vector<self_collisions> self_coll;
            // vector<interspecies_collisions> unself_coll;
//            vector<interspecies_f00_explicit_collisions> unself_f00_coll;
        };

        class collisions_2D {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            collisions_2D(const State2D& Yin, const double& deltat); 
            // ~self_collisions();
            void advance(State2D& Y, const Clock& W);
            void advancef0(State2D& Y, const Clock& W, State2D& Yh);
            void advancef1(State2D& Y, State2D& Yh);
            void advanceflm(State2D& Y, State2D& Yh);

            vector<self_collisions> self();
            // void advancef1(State1D& Y);
            // void advanceflm(State1D& Y);

        private:
        //  Variables
            State2D Yh;
            vector<self_collisions> self_coll;
            // vector<interspecies_collisions> unself_coll;
//            vector<interspecies_f00_explicit_collisions> unself_f00_coll;
        };


    #endif
