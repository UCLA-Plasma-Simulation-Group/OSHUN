/*!\brief  Vlasov Equation - Declarations
 * \author Michail Tzoufras, Benjamin Winjum, Archis Joglekar
 * \date   October 10, 2016
 * \file   vlasov.h
 *
 * Includes declarations for spatial advection, electric field advection, current, and 
 * the RK1D Functor responsible for stepping forward
 * 
 * \todo Needs Spatial Advection 2D
 * \todo Electric Field 2D
 * \todo RK2D Functor
 * 
 */

    #ifndef DECL_VLASOVMAXWELL_H
    #define DECL_VLASOVMAXWELL_H

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
    class Spatial_Advection_1D {
//--------------------------------------------------------------
        public:       	
//      Constructors/Destructors
            Spatial_Advection_1D(size_t Nl, size_t Nm, 
                             double pmin, double pmax, size_t Np, 
                             double xmin, double xmax, size_t Nx); 
//          Advance
            void operator()(const DistFunc1D& Din, DistFunc1D& Dh);

        private:
            SHarmonic1D                   	fd1, fd2;
            Array2D< complex<double> >  	A1, A2;
            valarray< complex<double> >  	vr;
    };
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Electric field
    class Electric_Field_1D {
//--------------------------------------------------------------    	
        public:
//      Constructors/Destructors
            Electric_Field_1D(size_t Nl, size_t Nm,
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

            Array2D< complex<double> >   A1, A2;
            valarray< complex<double> >  B1, B2;
            valarray< complex<double> >  C1, C3;
            Array2D< complex<double> >   C2, C4;
            valarray< complex<double> >  pr, invpr, Hp0;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Magnetic field
    class Magnetic_Field_1D {
//--------------------------------------------------------------    	
        public:
//      Constructors/Destructors
            Magnetic_Field_1D(size_t Nl, size_t Nm,
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
//--------------------------------------------------------------
//  Current
    class Current_1D {
//--------------------------------------------------------------
        public:
//      Constructors/Destructors
            Current_1D( double pmin, double pmax, size_t Np, 
                        size_t Nx ); 
//          Advance
            void operator()(const DistFunc1D& Din, 
            				Field1D& FExh, Field1D& FEyh, Field1D& FEzh);

        private:
            Field1D Jx, Jy, Jz;
            valarray< complex<double> >  pr, invg;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update B with E term from Faraday's Law
    class Faraday_1D {
//--------------------------------------------------------------
        public:       	
//      Constructors/Destructors
            Faraday_1D(double xmin, double xmax, size_t Nx); 
//          Advance
            void operator()(EMF1D& EMFin, EMF1D& EMFh);

        private:
			Field1D             tmpE;
			complex<double>		idx;
            size_t              numx;

    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update E with B from Ampere's Law
    class Ampere_1D {
//--------------------------------------------------------------
        public:       	
//      Constructors/Destructors
            Ampere_1D(double xmin, double xmax, size_t Nx); 
//          Advance
            void operator()(EMF1D& EMFin, EMF1D& EMFh);

        private:
            Field1D             tmpB;
			complex<double>		idx;
            size_t              numx;

    };
//--------------------------------------------------------------

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
            vector<Hydro_Advection_1D>   HA;
            
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
            vector<Faraday_1D>           FA;
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
            // vector<Faraday_1D>           FA;
    };
//--------------------------------------------------------------



/** @} */ 


    #endif
