/*! \brief  Export Files - Declarations
 *  \author PICKSC
 *  \date   2017
 *  \file   export.h
 * 
 */ 
//--------------------------------------------------------------

    #ifndef DECL_EXPORT_H
    #define DECL_EXPORT_H
    
//--------------------------------------------------------------
//--------------------------------------------------------------
namespace Export_Files{

    namespace ofconventions {
        const int    ofile_digits = 5;
        const string ofile_extension = ".txt";
        const int    ofile_precision = 6; 

        const int    rfile_digits = 3;
        const int    rank_digits = 6;
        const string rfile_extension = ".dat";

        const string h5file_extension = ".h5";
    }
//--------------------------------------------------------------

//  Folder and filename functions
    template<typename T> inline std::string stringify(T const& x) {
    std::ostringstream out;
    out << x;
    return out.str();
}

int Makefolder(string _name);

void Folders();


//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
class DefaultTags {
public:
//          Contents
    vector<string> time;
    vector<string> space;
    vector<string> fld;
    // vector<string> fld2d;
    vector<string> mom;
    // vector<string> mom2d;
    vector<string> part;
    vector<string> pvsx;
    // vector<string> pvsx2d;
    vector<string> fvsx;
    // vector<string> fvsx2d;
    vector<string> pvspvsx;
    // vector<string> pvspvsx2d;

    DefaultTags(size_t species); 
    ~DefaultTags(){}
};
//--------------------------------------------------------------

//  redefinitions of "<<" operator for data containers
    template <class T> 
ofstream& operator<<(ofstream& s, const vector<T>& v) {
    s << setprecision(ofconventions::ofile_precision);
    s << 1 <<"\n";
    s << v.size()<<"\n";
    for (size_t i(0); i < v.size(); ++i) {
        s << v[i]<<"\n";
    }
    return s;
}

    template <class T> 
ofstream& operator<<(ofstream& s, const valarray<T>& v) {
    s << setprecision(ofconventions::ofile_precision);
    s << 1 <<"\n";
    s << v.size()<<"\n";
    for (size_t i(0); i < v.size(); ++i) {
        s << v[i]<<"\n";
    }
    return s;
}

    template <class T> 
ofstream& operator<<(ofstream& s, const Array2D<T>& array2D) {
    s << setprecision(ofconventions::ofile_precision);
    s << 2 <<"\n";
    s << array2D.dim1()<<"\n";
    s << array2D.dim2()<<"\n";
    for (size_t i(0); i < array2D.dim(); ++i) {
        s << array2D(i)<<"\n";
    }
    return s;
}

    template <class T> 
ofstream& operator<<(ofstream& s, const Array3D<T>& array3D) {
    s << setprecision(ofconventions::ofile_precision);
    s << 3 <<"\n";
    s << array3D.dim1()<<"\n";
    s << array3D.dim2()<<"\n";
    s << array3D.dim3()<<"\n";
    for (size_t i(0); i < array3D.dim(); ++i) {
        s << array3D(i)<<"\n";
    }
    return s;
}

template <class T>
ofstream& operator<<(ofstream& s, const Array4D<T>& array4D) {
    s << setprecision(ofconventions::ofile_precision);
    s << 4 <<"\n";
    s << array4D.dim1()<<"\n";
    s << array4D.dim2()<<"\n";
    s << array4D.dim3()<<"\n";
    s << array4D.dim4()<<"\n";
    for (size_t i(0); i < array4D.dim(); ++i) {
        s << array4D(i)<<"\n";
    }
    return s;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Facilitates the generation of a header
        class Header {
//--------------------------------------------------------------        
        public:
//          Constructor
            Header() { };
            Header(string _axis_units, string _quantity_units, string _outputDir);
            
            double  AxisUnits()                {return formulary().Uconv(axis_units);}
            double  QuantityUnits()            {return formulary().Uconv(quantity_units);}
            string  Directory()                {return outputDir;}

        private:
            string              axis_units;
            string              quantity_units;
            string              outputDir;
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Main facility for exporting data 
        class Xport {
//--------------------------------------------------------------        
        public:
//          Constructor
            Xport(const Algorithms::AxisBundle<double>& _axis, 
              const vector< string > oTags,
              string homedir=""); 

            void Export_h5(const std::string tag, std::vector<double> &axis1, 
                std::vector<double> &data, 
                const size_t step, const double time, const double dt,
                const int spec = -1);

            void Export_h5(const std::string tag, 
                std::vector<double> &axis1, std::vector<double> &axis2, 
                Array2D<double> &dataA, 
                const size_t  step, const double time, const double dt,
                const int spec = -1);

            void Export_h5(const std::string tag, 
                std::vector<double> &axis1, std::vector<double> &axis2, std::vector<double> &axis3, 
                Array3D<double> &dataA, 
                const size_t  step, const double time, const double dt,
                const int spec = -1);

            void Export_h5(const std::string tag, 
                std::vector<double> &axis1, std::vector<double> &axis2, std::vector<double> &axis3, std::vector<double> &axis4, 
                Array4D<double> &dataA, 
                const size_t  step, const double time, const double dt,
                const int spec = -1);
            
            void add_attributes(HighFive::DataSet &dataset, const std::string tag, 
                const double time, const double dt);

        private:
            map< string, Header > Hdr; // Dictionary of headers
            string oH5Fextension(size_t step, int species = -1);

        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        class Restart_Facility {
//--------------------------------------------------------------
        public:
            Restart_Facility(const int rank, string homedir="");

            void Read(const int rank, const size_t re_step, State1D& Y, double time_start);
            void Write(const int rank, const size_t re_step, State1D& Y, double time_dump);

            void Read(const int rank, const size_t re_step, State2D& Y, double time_start);
            void Write(const int rank, const size_t re_step, State2D& Y, double time_dump);

        private:
            string hdir;
            string rFextension(const int rank, const size_t rstep);
        }; 
//--------------------------------------------------------------
    }
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    namespace Output_Data{
//--------------------------------------------------------------
//  LEGENDRE VALUES
//  A 2D space is generated from px, |p| axis. The cos8 for
//  this 2D space is calculated and then the Legendre polynomials
//  for each cos8 are calculated. We end up with a 1D array for l
//  containing 2D matrices of the polynomials for each cos8
        class PLegendre2D {
        public:
//      Constructors/Destructors
            PLegendre2D(size_t Nl, size_t Nm, double pmax,
             valarray<double> px , valarray<double> py);
            PLegendre2D(const PLegendre2D& other);
            ~PLegendre2D();

//      Access
            size_t dim()   const { return (plegendre).size(); }
            
            Array2D<double>& operator()(size_t i)       { return (plegendre)[i]; }
            Array2D<double>  operator()(size_t i) const { return (plegendre)[i]; }   

        private:
            vector< Array2D<double> > plegendre;
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
//--------------------------------------------------------------
        class fulldistvsposition {
//--------------------------------------------------------------        
        public:
            // Constructor/Destructor
            fulldistvsposition(const Grid_Info& _G);
        

            ~fulldistvsposition();

            // 1P at a single point
            valarray<double> p1(DistFunc1D& df, size_t x0, size_t s) ;
            valarray<double> p2(DistFunc1D& df, size_t x0, size_t s) ;
            valarray<double> p3(DistFunc1D& df, size_t x0, size_t s) ;

            // 2P at a single point
            Array2D<double> p1p2(DistFunc1D& df, size_t x0, size_t s) ;
            Array2D<double> p2p3(DistFunc1D& df, size_t x0, size_t s) ;
            Array2D<double> p1p3(DistFunc1D& df, size_t x0, size_t s) ;
            Array3D<double> p1p2p3(DistFunc1D& df, size_t x0, size_t s) ;

            // 1P, integrate over x or y
            // valarray<double> p1(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;
            // valarray<double> p2(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;
            // valarray<double> p3(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;

            // 1P, no spatial integration
            valarray<double> p1(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;
            valarray<double> p2(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;
            valarray<double> p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;

            // 2P at a single point, integrated over the other x dimension
            // Array2D<double> p1p2(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;
            // Array2D<double> p2p3(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;
            // Array2D<double> p1p3(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;
            // Array3D<double> p1p2p3(size_t integrationdimension, DistFunc2D& df, size_t x0, size_t s) ;

            // 2D at a single point, not integrated over the other dimension
            Array2D<double> p1p2(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;
            Array2D<double> p2p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;
            Array2D<double> p1p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;
            // Array3D<double> p1p2p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) ;

            // Access
            Grid_Info           gridinfo()              const   { return grid;}
        

        private:
            Grid_Info           grid;  

            vector<vector<double> >    pvec;
            vector<valarray<double> >  pxvec,pyvec,pzvec;
            

            vector< PLegendre2D     >  PL2D;
            vector< valarray<double>  >  pout1D_p1, pout1D_p2, pout1D_p3;
            vector< Array2D<double>  >  pout2D_p1p2, pout2D_p1p3, pout2D_p2p3;
            vector< Array3D<double>  >  pout3D;

            // Interpolation quantities
            vector< Array3D<double>  >  pradius;
            vector< Array2D<double>  >  phi;
            vector< valarray<double> >            dpx,dpy,dpz;
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
//--------------------------------------------------------------
        class harmonicvsposition {
//--------------------------------------------------------------        
        public:
//      Constructor/Destructor
            harmonicvsposition(const Grid_Info& _G);
            harmonicvsposition(const harmonicvsposition& other);

            ~harmonicvsposition();

//      Methods
            Array2D<double> operator()(DistFunc1D& df, size_t l, size_t m, size_t x0, size_t s) ;
            
            // Array2D<double> operator()(DistFunc2D& df, size_t l, size_t m, size_t x0, size_t s) ;

            // Array2D<double> operator()(size_t integrationdimension, DistFunc2D& df, size_t l, size_t m, size_t x0, size_t s) ;

            Array2D<double> operator()(DistFunc2D& df, size_t l, size_t m, size_t x0, size_t y0, size_t s) ;

//      Access
            size_t Species()         const { return nump.size(); }
            size_t Np(size_t s)     const  { return nump[s]; }
            double  Pmin(size_t s)    const  { return pmin[s]; }
            double  Pmax(size_t s)    const  { return pmax[s]; }
            // valarray<double> Dp(size_t s) const {return deltap[s];}
            valarray<double> paxis(size_t s) const {return pvec[s];}

        private:
           vector<double> pmin, pmax;
           vector<double> nump;
           // vector<valarray<double> > deltap;
           vector<valarray<double> > pvec;

       };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Output Functor
       class Output_Preprocessor {
//--------------------------------------------------------------        
       public:
//      Constructor
        Output_Preprocessor(const Grid_Info& _grid, 
         const vector< string > _oTags, 
         string homedir="")  
        : expo( _grid.axis, _oTags, homedir), p_x( _grid),f_x( _grid),oTags(_oTags) { }

//      Functor
        void operator()(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void distdump(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void bigdistdump(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);

        void operator()(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void distdump(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void bigdistdump(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);

    private:
        size_t                          Nbc;
        Export_Files::Xport             expo;
        fulldistvsposition              p_x;
        harmonicvsposition              f_x;
        vector< string >                oTags;
        
        // Fields
        void Ex(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Ey(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Ez(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Bx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void By(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Bz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);

        void Ex(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Ey(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Ez(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Bx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void By(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Bz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);

        // P-Integrated distributions
        void px(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void py(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void pz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void pxpy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_1D& PE);
        void pypz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_1D& PE);
        void pxpz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_1D& PE);

        void px(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void py(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void pz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void pxpy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_2D& PE);
        void pypz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_2D& PE);
        void pxpz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_2D& PE);

        // Full distribution
        void pxpypz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
         const Parallel_Environment_1D& PE);

        // Raw distributions
        void f0(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void f10(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);        
        void f11(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void f20(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);        
        void fl0(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);

        // For 2D , integrated over space
        // void f0i(size_t integrationdimension, const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
        //     const Parallel_Environment_2D& PE);
        // void f10i(size_t integrationdimension, const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
        //     const Parallel_Environment_2D& PE);        
        // void f11i(size_t integrationdimension, const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
        //     const Parallel_Environment_2D& PE);
        // void f20i(size_t integrationdimension, const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
        //     const Parallel_Environment_2D& PE);        
        // void fl0i(size_t integrationdimension, const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
        //     const Parallel_Environment_2D& PE);

        void f0(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void f10(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);        
        void f11(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void f20(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);        
        void fl0(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);

        // f0 moments
        void n(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void T(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);

        void n(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void T(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);


        // Particle tracker
        void particles_x(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void particles_px(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void particles_py(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void particles_pz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);

        // f1 moments
        void Jx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Jy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Jz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Qx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Qy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Qz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void vNx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void vNy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
           const Parallel_Environment_1D& PE);
        void vNz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
           const Parallel_Environment_1D& PE);

        void Jx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Jy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Jz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Qx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Qy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void Qz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void vNx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_2D& PE);
        void vNy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
           const Parallel_Environment_2D& PE);
        void vNz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
           const Parallel_Environment_2D& PE);

        // Hydro
        void Ux(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Uy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Uz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Z(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void ni(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
        void Ti(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
            const Parallel_Environment_1D& PE);
    };

}

    #endif
