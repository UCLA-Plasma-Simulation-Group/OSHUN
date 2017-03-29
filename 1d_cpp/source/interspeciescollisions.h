/*! \brief Interspecies Collisions - Definitions
  * \date   October 17, 2016
 * \file   interspeciescollisions.h
 * 
 * Contains:
 * 1) the explicit energy-conserving algorithm for effect of 
 * electron-electron collisions on f_00
 * 2) implicit algorithm used for e-e + e-i collisions on f_lm
 * And all their containers.
 * 
 */

#ifndef DECLARATION_INTERSPECIES_FP_H
#define DECLARATION_INTERSPECIES_FP_H


/** \addtogroup vfp1d
 *  @{
 */
//-------------------------------------------------------------------
/** 
* \class explicit_f00
* \brief The top-level container for collisions on l=0.
* 
*   The explicit_f00 class describes the object that 
*   contains the RK4 algorithm that sets up the explicit solve step. 
*   explicit_f00 has a "loop()" function that pulls the relevant
*   distribution function and passes it to the RK4 solver.
*/        
        class interspecies_f00_explicit_step {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            interspecies_f00_explicit_step(const DistFunc1D& DF1, const DistFunc1D& DF2, const double& deltat);//, int tout_start) //DistFunc1D& DFin,

            valarray<double>  takestep(const valarray<double>& f1in, const valarray<double>& f2in);


        private:
        //  Variables

            // SHarmonic1D&        f1, f2;
            valarray<double>    fslope;

            double              dt;
            Formulary           formulary;
            void                calculateintegrals(const valarray<double>& f1in, const valarray<double>& f2in);            
            // vector<double>      f1;
            void                remapintegrals();

            double              m1, m2, z1, z2, n1, n2, T1, T2, Gamma12, kpre;
            
            valarray<double>            pgrid_s1, pgrid_s2;
            valarray<double>            df0;
            valarray<double>            U4, U4m1, U2, U2m1, U1, U1m1, U3, Qn, Pn;

            ///     The integrals
            valarray<double>            J1_s2, I2_s2, I4_s2;

            ///     The integrals
            valarray<double>            J1_s1, I2_s1, I4_s1;

            int                         Nbc; ///< Number of boundary cells in each direction
            int                         szx; ///< Total cells including boundary cells in x-direction

        };
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
    class interspecies_f00_RKfunctor : public Algorithms::AbstFunctor<valarray<double>> {
//--------------------------------------------------------------        
        public:
//          Constructor
            interspecies_f00_RKfunctor(const DistFunc1D& DF1,const DistFunc1D& DF2,const double& smalldt);
            ~interspecies_f00_RKfunctor(){ };

//          Collect all the operators and apply on Yin
            void operator()(const valarray<double>& fin1, valarray<double>& fslope);
            void operator()(const valarray<double>& fin1, const valarray<double>& fin2, valarray<double>& fslope);
            void operator()(const valarray<double>& fin1, valarray<double>& fslope, size_t dir);

        private:
            //          Variables
            interspecies_f00_explicit_step collide;
             
    };
//--------------------------------------------------------------    

/** \addtogroup vfp1d
 *  @{
 */
//-------------------------------------------------------------------
/** 
* \class interspecies_f00_explicit_step
* \brief The top-level container for collisions on l=0.
* 
*   The interspecies_f00_explicit_step class describes the object that 
*   contains the RK4 algorithm that sets up the explicit solve step. 
*   interspecies_f00_explicit_step has a "loop()" function that pulls the relevant
*   distribution function and passes it to the RK4 solver.
*/        
        class interspecies_f00_explicit_collisions {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            interspecies_f00_explicit_collisions(const DistFunc1D& DF1,const DistFunc1D& DF2, const double& deltat);//, int tout_start) //DistFunc1D& DFin,

        /// This loop calls the RK4_f00 private member that is
        /// responsible for setting up the RK4 algorithm to advance
        /// the collision step. 
            void rkloop(SHarmonic1D& SH1, const SHarmonic1D& SH2);

        private:
        //  Variables
            interspecies_f00_RKfunctor              rkf00;
            valarray<double>                        fin1,fin2;
            Algorithms::RK4<valarray<double>>            RK4;
            

            size_t                                  num_h;
            double                                  h;

            int Nbc; ///< Number of boundary cells in each direction
            int szx; ///< Total cells including boundary cells in x-direction

        };
//-------------------------------------------------------------------

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
/** 
 * \class interspecies_flm_implicit_step
 * \brief Collisions for l >= 1.
 * 
 *   This class forms the numerical backbone for quantifying
 *   the effect of collisions on the anisotropic components of
 *   the distribution function.
*/   
//--------------------------------------------------------------
        class interspecies_flm_implicit_step {
//-------------------------------------------------------------------
        private:

//          Define the powers of the velocity
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;

//          Define the integrals
            valarray<double>  J1m, I0, I2;

//          Define the derivatives of f0
            valarray<double>  df0, ddf0;

//          Variables for the matrix A, such that A*x = b
            valarray<double> Scattering_Term;
            Array2D<double> Alpha_Tri;
            bool if_tridiagonal;

//          Constant
            double I0_density, I2_temperature;
            double _ZLOGei, _LOGee, kpre, Dt;
            Formulary formulas;

        public:
//          Constructors/Destructors
            interspecies_flm_implicit_step(double pmax, size_t nump); 
         
//          Calculate the coefficients
            void reset_coeff(const valarray<double>& f00, const double Delta_t);    

//          Explicit Advance
            void advance(valarray<complex<double> >& fin, const int el);    
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
        class interspecies_flm_implicit_collisions {
        public:
///          Constructors/Destructors
            interspecies_flm_implicit_collisions(const DistFunc1D &DFin, const double& deltat); 

///         Advance f1_loop for 2D code.
            void advancef1(DistFunc1D& DF); 

///         Advance flm_loop for 2D code.
            void advanceflm(DistFunc1D& DF);

        private:
            // State1D& Y;
            double Dt;

            int Nbc; ///< Number of boundary cells in each direction
            int szx; ///< Total cells including boundary cells in x-direction

            size_t l0; ///< Number of m harmonics.
            size_t m0; ///< Number of m harmonics.

            valarray<complex<double> >  fc;   ///< Dummy array
            valarray<double> f00;   ///< Array for isotropic component distribution function. Needed for calculating coefficients.

///         The object that is responsible for performing the algebra required for the integrals.
            interspecies_flm_implicit_step  implicit_step;  

        };
//--------------------------------------------------------------
/** @} */ 

//-------------------------------------------------------------------
/** 
* \class interspecies_collisions
* \brief The top-level container for interspecies collisions over all species on l=0.
* 
*   
*/        
        class interspecies_collisions {
//--------------------------------------------------------------            
        public:
        /// Constructors/Destructors
            interspecies_collisions(const State1D& Yin, const size_t& sind, const double& deltat); 

            void advancef00(State1D& Y, const size_t& sind);
            // void advancef1(State1D& Y, const size_t sin);
            // void advanceflm(State1D& Y, const size_t sin);

        private:
        //  Variables
            vector<interspecies_f00_explicit_collisions> interspecies_f00_collisions;
            // vector<implicit_flm_interspecies> flm_collisions;

            // double              t, tout;
            // int Nbc; ///< Number of boundary cells in each direction
            // int szx; ///< Total cells including boundary cells in x-direction


        };


#endif