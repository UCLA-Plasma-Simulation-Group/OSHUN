  
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
    class Fluid_Equation_1D {
//--------------------------------------------------------------
        public:       	
            /// Constructors/Destructors
            Fluid_Equation_1D(double xmin, double xmax, size_t Nx); 
            
            /// Momentum (and Energy and Density?) Update
            void density(const Hydro1D& Hin, Hydro1D& Hslope);//State1D& Yin, double deltat);
            void chargefraction(const Hydro1D& Hin, Hydro1D& Hslope);//State1D& Yin, double deltat);
            void velocity(const State1D& Yin, Hydro1D& Hslope,
                            valarray<double>& electrondensity, valarray<double>& electronpressure,
                            Array2D<double>& electroncurrent); //State1D& Yin,double deltat);
            // void energy()(const valarray<double>& energy, valarray<double>& newenergy);
            

            /// Electric Field Update
            void updateE( Hydro1D& HYDRO,
                                EMF1D& EMF);


            /// Evaluate Quantities from Distribution Functions
//            void calcquantities(const State1D& Yin);
            // void evaluatekineticenergyexchange(const State1D& Yin);            

        private:

            double                          idx;
            size_t                          szx, Nbc;

            valarray<double>                dummy;
            
            
    };
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 
    class Hydro_Functor : public Algorithms::AbstFunctor<State1D> {
//--------------------------------------------------------------        
        public:
//          Constructor
            Hydro_Functor(double xmin, double xmax, size_t numx); 
            ~Hydro_Functor(){ };

//          Collect all the operators and apply on Yin
            void operator()(const State1D& Hin, State1D& Hslope);
            void operator()(const State1D& Hin, State1D& Hslope, size_t dir);
            void operator()(const State1D& H1in, const State1D& H2in, State1D& Hslope);
        private:
            //          Variables
            double idx;
            size_t szx;
            Fluid_Equation_1D FEQ;
    };
//--------------------------------------------------------------  
//--------------------------------------------------------------
//  Spatial advection
    class Hydro_Advection_1D {
//--------------------------------------------------------------
        public:         
        //  Constructors/Destructors
            Hydro_Advection_1D(size_t Nl, size_t Nm, 
                             double pmin, double pmax, size_t Np, 
                             double xmin, double xmax, size_t Nx); 

    //          Advance
            void operator()(const DistFunc1D& Din, const Hydro1D& hydro, DistFunc1D& Dh);

        private:
            void MakeGH( SHarmonic1D& f, size_t l);
            void MakeG00( SHarmonic1D& f);

            SHarmonic1D H, G, TMP;
            SHarmonic1D                     fd1, fd2;
            Array2D< complex<double> >      XX1, XX2, XX3, XX4;
            Array2D< complex<double> >      A1, A2;
            valarray< complex<double> >     pr, invpr, Hp0;
            size_t                          szx, Nbc;
            double                          idp, idx;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------