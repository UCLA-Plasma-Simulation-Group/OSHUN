/*! \brief  Export Files - Declarations
 * \author PICKSC
 *  \date   September 1, 2016
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
        vector<string> mom;
        vector<string> pvsx;
        vector<string> fvsx;
        vector<string> pvspvsx;
        
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
//--------------------------------------------------------------

//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class oAxis {
        public:
//          Contents
            string label;           // e.g. label = cm 
            string units;
            float min, max;  
            size_t sz;

//          Constructors
            oAxis();
            oAxis(const float _m, const float _M, const size_t _sz); 
            oAxis(const string _l, const string _u, const float _m, const float _M, 
                  const size_t _sz);

//          Copy constructor 
            oAxis(const oAxis& other);
            ~oAxis(){}
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Facilitates the generation of a header
    class Header {
//--------------------------------------------------------------        
        public:
//          Constructor
            Header() { };
            Header(oAxis _x,                                        // 1D
                   string _Ql, float _Qc, string _tl, string _tu, float _tc, string _oD);
            Header(oAxis _x, oAxis _y,                             // 2D 
                   string _Ql, float _Qc,  string _tl, string _tu, float _tc, string _oD);
            Header(oAxis _x, oAxis _y, oAxis _z,                  // 3D
                   string _Ql, float _Qc, string _tl, string _tu, float _tc, string _oD);
            Header(vector< oAxis > _xyz,                            // xD
                   string _Ql, float _Qc, string _tl, string _tu, float _tc, string _oD);
            size_t dim();    

            valarray<float> axis(const size_t i); // this is axis 0, 1, 2 
            string          label(const size_t i);
            string          units(const size_t i);
            float           conv(const size_t i);

            float           min(const size_t i) {return(xyz_axis[i].min);}
            float           max(const size_t i) {return(xyz_axis[i].max);}
            
            string  Title_label(); 
            float   Title_conv(); 
            string  Time_label();
            float   Time_conv();
            string  Directory();
        
        private:
            vector< oAxis >    xyz_axis;         // axis 0, 1, 2  : size 0-2
            string             title,  time, timeU;
            float              titleC, timeC;
            string         oDir;
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
        
            void Export_h5(const std::string tag, std::valarray<float> ex, 
                            const size_t& step, const int spec = -1);
            void Export_h5(const std::string tag, Array2D<float> ex, 
                            const size_t& step, const int spec = -1);
            void Export_h5(const std::string tag, Array3D<float> ex, 
                            const size_t& step, const int spec = -1);

            H5::H5File hmake_file(string ofilename);
            void hclose_file(H5::H5File &file);
            void hinit_attr(H5::H5File &hfilehandle, const std::string tag, 
                                size_t step, float xmax, float xmin);
            void hinit_attr2(H5::H5File &hfilehandle, const std::string tag, 
                                size_t step, float xmax[], float xmin[]);
            void haxis(H5::Group &hgrouphandle, string axismainname, float axisrange[2], 
                                string axislongname, string axisname, 
                                string axistype, string axisunits);
            void hfile_add_attr(H5::H5File &hfilehandle, string attrname, int attrdata);
            void hfile_add_attr2(H5::H5File &hfilehandle, string attrname, int attrdata[2]);
            void hfile_add_attr(H5::H5File &hfilehandle, string attrname, float attrdata);
            void hfile_add_attr2(H5::H5File &hfilehandle, string attrname, float attrdata[2]);
            void hfile_add_attr(H5::H5File &hfilehandle, string attrname, string attrdata);
            void hfile_add_attr_todataset(H5::DataSet &hdatasethandle, string attrname, string attrdata);


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
        
        void Read(const int rank, const size_t re_step, State1D& Y);
        void Write(const int rank, const size_t re_step, State1D& Y);

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
//--------------------------------------------------------------
//  LEGENDRE VALUES
//  A 2D space is generated from px, |p| axis. The cos8 for
//  this 2D space is calculated and then the Legendre polynomials
//  for each cos8 are calculated. We end up with a 1D array for l
//  containing 2D matrices of the polynomials for each cos8
    class PLegendre1D {
    public:
//      Constructors/Destructors
        PLegendre1D(size_t Nl, size_t Nm, size_t Np,  float pmin, float pmax,
                               size_t Npx );
        PLegendre1D(const PLegendre1D& other);
        ~PLegendre1D();

//      Access
        size_t dim()   const { return (*plegendre).size(); }
        Array2D<float>& operator()(size_t i)       { return (*plegendre)[i]; }
        Array2D<float>  operator()(size_t i) const { return (*plegendre)[i]; }   
    
private:
        vector< Array2D<float> > *plegendre;
        size_t lmax;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  LEGENDRE VALUES
//  A 2D space is generated from px, |p| axis. The cos8 for
//  this 2D space is calculated and then the Legendre polynomials
//  for each cos8 are calculated. We end up with a 1D array for l
//  containing 2D matrices of the polynomials for each cos8
    class PLegendre2D {
    public:
//      Constructors/Destructors
        PLegendre2D(size_t Nl, size_t Nm, size_t Np,  float pmin, float pmax,
                               size_t Npx , size_t Npy);
        PLegendre2D(const PLegendre2D& other);
        ~PLegendre2D();

//      Access
        size_t dim()   const { return (*plegendre).size(); }
        Array2D<float>& operator()(size_t i)       { return (*plegendre)[i]; }
        Array2D<float>  operator()(size_t i) const { return (*plegendre)[i]; }   
    
private:
        vector< Array2D<float> > *plegendre;
        size_t lmax;
    };
//--------------------------------------------------------------


//--------------------------------------------------------------
//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
//--------------------------------------------------------------
    class p1x1_1D {
//--------------------------------------------------------------        
    public:
//      Constructor/Destructor
        p1x1_1D(const Grid_Info& _G);
        p1x1_1D(const p1x1_1D& other);
                        
        ~p1x1_1D();

//      Methods
        valarray<float> operator()(DistFunc1D& df, size_t x0, size_t s) ;

//      Access
        size_t Species()         const { return p1x1.size(); }
//        size_t Nl(size_t s)      const  { return Pl[s].dim(); }
//        size_t Npx(size_t s)     const  { return polarR[s].dim1(); }
//        size_t Np(size_t s)      const  { return polarR[s].dim2(); }
        float  Pmin(size_t s)    const  { return pmin[s]; }
        float  Pmax(size_t s)    const  { return pmax[s]; }

//        PLegendre1D      PL(size_t s)     const { return Pl[s]; }
        valarray<float>  p1_x1(size_t s)  const { return p1x1[s]; }
//        Array2D<float>   PolarR(size_t s) const { return polarR[s]; }

    private:
//        vector< PLegendre1D >    Pl;
//        vector< Array2D<float> > polarR;

        vector< valarray<float> > p1x1;
        vector<float> pmin, pmax;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
//--------------------------------------------------------------
    class p2p1x1_1D {
//--------------------------------------------------------------        
    public:
//      Constructor/Destructor
        p2p1x1_1D(const Grid_Info& _G);
        p2p1x1_1D(const p2p1x1_1D& other);
                        
        ~p2p1x1_1D();

//      Methods
        Array2D<float> operator()(DistFunc1D& df, size_t x0, size_t s) ;

//      Access
        size_t Species()         const { return p2p1x1.size(); }
        size_t Nl(size_t s)      const  { return numl[s]; }
        size_t Nm(size_t s)      const  { return numm[s]; }
        size_t Npx(size_t s)     const  { return Pl[s](0).dim1(); }
        size_t Npy(size_t s)     const  { return Pl[s](0).dim2(); }
        // size_t Np(size_t s)      const  { return polarR[s].dim2(); }
        float  Pmin(size_t s)    const  { return pmin[s]; }
        float  Pmax(size_t s)    const  { return pmax[s]; }


        PLegendre2D      PL(size_t s)     const { return Pl[s]; }
        Array2D<float>  p2p1_x1(size_t s)  const { return p2p1x1[s]; }
        Array2D<float>   prad(size_t s) const { return pr[s]; }
        Array2D<size_t>  npcell(size_t s) const {return nextpcell[s];}
        Array2D<float>  dcell(size_t s) const {return distancetothatcell[s];}

    private:
        vector<float> pmin, pmax;
        vector<size_t> numl, numm;
        vector< PLegendre2D >    Pl;
        vector< Array2D<float> > pr;
        vector< Array2D<size_t> > nextpcell;
        vector< Array2D<float> > distancetothatcell;

        vector< Array2D<float> > p2p1x1;
        
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
//--------------------------------------------------------------
    class fx1_1D {
//--------------------------------------------------------------        
    public:
//      Constructor/Destructor
        fx1_1D(const Grid_Info& _G);
        fx1_1D(const fx1_1D& other);
                        
        ~fx1_1D();

//      Methods
        Array2D<float> operator()(DistFunc1D& df, size_t l, size_t m, size_t x0, size_t s) ;

//      Access
        size_t Species()         const { return nump.size(); }
        // size_t Npx(size_t s)     const  { return f0x1[s].size(); }
        size_t Np(size_t s)     const  { return nump[s]; }
         float  Pmin(size_t s)    const  { return pmin[s]; }
         float  Pmax(size_t s)    const  { return pmax[s]; }

        
        // valarray<float>  f_x1(size_t s)  const { return f0x1[s]; }
        

    private:
        // vector< valarray<float> > f0x1;
         vector<float> pmin, pmax;
        vector<float> nump;
        
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Output Functor
    class Output_Preprocessor_1D {
//--------------------------------------------------------------        
    public:
//      Constructor
        Output_Preprocessor_1D(const Grid_Info& _grid, 
                               const vector< string > _oTags, 
                               string homedir="")  
        : expo( _grid.axis, _oTags, homedir),
          px_x( _grid ), //pxpy_x( _grid),
          f_x( _grid),
          oTags(_oTags) { }

//      Functor
//         void operator()(const State1D& Y, const size_t tout);
//         void operator()(const State1D& Y, const Grid_Info& grid, const size_t tout,
//                                      const Parallel_Environment_1D& PE);
        void operator()(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void distdump(const State1D& Y, const Grid_Info& grid, const size_t tout,
                        const Parallel_Environment_1D& PE);

    private:
        size_t              Nbc;
        Export_Files::Xport expo;
        p1x1_1D             px_x;
        fx1_1D             f_x;
//        p2p1x1_1D           pxpy_x;
        vector< string >    oTags;

        void Ex(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Ey(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Ez(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Bx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void By(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Bz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void px(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
//        void pxpy(const State1D& Y, const Grid_Info& grid, const size_t tout,
//                                        const Parallel_Environment_1D& PE);
        void f0(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void f10(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);        
        void f11(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void f20(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);        
        void fl0(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void n(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void T(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Jx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Jy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Jz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Qx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Qy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                const Parallel_Environment_1D& PE);
        void Qz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                const Parallel_Environment_1D& PE);
        void vNx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void vNy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                 const Parallel_Environment_1D& PE);
        void vNz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                 const Parallel_Environment_1D& PE);

        void Ux(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Uy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                const Parallel_Environment_1D& PE);
        void Uz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                const Parallel_Environment_1D& PE);
        void Z(const State1D& Y, const Grid_Info& grid, const size_t tout,
                const Parallel_Environment_1D& PE);
        void ni(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);
        void Ti(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                        const Parallel_Environment_1D& PE);


    };

}

    #endif
