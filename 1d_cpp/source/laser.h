/*!\brief  Laser source - Declarations
 * \author PICKSC
 * \date   August 9, 2016
 * \file   laser.h
 *
 * The laser package contains:
 * - Inverse bremsstrahlung heating
 * - External field driver (TBD)
 * - Phenomenological heat source (TBD)
 * - Ray tracing (TBD)
 * - Antenna (TBD)
 * 
 */
//**************************************************************

#ifndef DECLERATION_PLASERSOURCE_H
#define DECLERATION_PLASERSOURCE_H

//--------------------------------------------------------------
        class IB_f00{
//--------------------------------------------------------------
//      Decleration of Energy/Number-Conserving Algorithm
//--------------------------------------------------------------
        private:
//          retain a reference to the slope
            valarray<double>&  fh;

//          Define the velocity axis
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1;//, U1, U1m1;
            valarray<double>  Inv_Uav6, gn, Qn, Pn;

//          Define the integrals
            // valarray<double>  I2, I4; // J1, 

//          Constant
            double c_kpre, Qn_coeff, vw_coeff_cube;
            double Inv_Uav6_nm1, Pnm1;

            Formulary formulas;

        public:
//          Constructors/Destructors
            IB_f00(valarray<double>& fslope, double pmax);
         
//          Explicit Advance
            valarray<double>& Getslope(const valarray<double>& fin, const double vos, const double Zval);
        };
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//**************************************************************
//   Decleration for the standard RK4 class for valarrays
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class RK4_IB {
//--------------------------------------------------------------
//      Decleration of the 4rth order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructor
            RK4_IB(valarray<double>& fin, double pmax, int tout_start);

//          Main function
            RK4_IB& advance(const double vosc, const double Zval);
        
//          Access
            double&  tout();
            double&  time();
            size_t&  numh();
            double&  th();

        private:
//          Helper functions
            valarray<double>& F(const valarray<double>& fin, const double vosc, const double Zval);

//          Variables
    	    valarray<double>  f0, f1, fh;
            valarray<double>&  f;


            IB_f00 Inversebremsstrahlung;

            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class InverseBremsstrahlung{
//--------------------------------------------------------------
//      Decleration of the Phenomenological Laser Source 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            // InverseBremsstrahlung(DistFunc1D& DFin, int tout_start, const valarray<double>& grid); 
            InverseBremsstrahlung(double pmax, size_t nump, size_t numx, int tout_start, const valarray<double>& grid);
         
//          Explicit Advance
            void loop(SHarmonic1D& f0, valarray<double>& Zarray, const double& tnew);

            valarray<double>  U2, U2m1;//, U1, U1m1;

//          Laser Temporal Evolution
            // double t_Profile(const double& t2); 

//          Laser polarization
            int IBSOURCE() const; 

        private:

//          Variables
            // DistFunc1D& DF;

            double              t, tout;

//          harmonic data
            valarray<double>    fc;


//          Algorithm
            RK4_IB      rk4_ib;

	        valarray<double>  vr; 

            //          Status 
            int                 InvBremsstrahlung;


            int                 Nbc, szx, szy;

//          Define the velocity axis
                       



//          Laser profile
            valarray<double>     IL_xprofile;//, IL_yprofile;

//          Laser profile
            // Axis< double >  xaxis;

//          Define laser profile
            // const double rise_time, fall_time, flat_time;
            const double vos, omega_0, omega_p, w0overwp;
            // const bool   linear_Ivst, polynomial_Ivst;
        };
//--------------------------------------------------------------
//**************************************************************

#endif
