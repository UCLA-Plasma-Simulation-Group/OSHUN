/*!\brief  Grid Setup - Declaration
 * \author Michail Tzoufras, Benjamin Winjum
 * \date   September 1, 2016
 * \file   setup.h
 *
 * Declares a Grid_Info class
 * 
 * \todo Create a function parser
 * \todo Change plasma input parameters such that they accept a toggle between uniform, piecewise, and functional profiles
 */

    #ifndef DECL_SETUP_H
    #define DECL_SETUP_H

//--------------------------------------------------------------        
//--------------------------------------------------------------
class Grid_Info {
//--------------------------------------------------------------
public:
//  Constructor
    Grid_Info(const vector<size_t> _Nl,   const vector<size_t> _Nm, 
              const vector<double> _mass, const vector<double> _charge,
//            local  spatial axes
              const vector<double> _xmin, const vector<double> _xmax,  const vector<size_t> _Nx,     
//            global spatial axes
              const vector<double> _xgmin,const vector<double> _xgmax, const vector<size_t> _Nxg,   
//                            momentum axes
                                          const vector<double> _pmax,  const vector<size_t> _Np,     
//                                             output momentum axes
              const vector<size_t> _Npx,   const vector<size_t> _Npy, const vector<size_t> _Npz) 
         : Nl(_Nl), Nm(_Nm), 
           Np(_Np), mass(_mass), charge(_charge),
                axis( _xmin, _xmax, _Nx, _xgmin, _xgmax, _Nxg, _pmax, _Np, _Npx, _Npy, _Npz ) {}

//--------------------------------------------------------------
//  Copy constructor
    Grid_Info(const Grid_Info& other): Nl(other.Nl), Nm(other.Nm), 
                                        Np(other.Np),
                                        mass(other.mass), charge(other.charge),
                                        axis(other.axis){}
    ~Grid_Info(){};

    const vector<size_t> Nl;
    const vector<size_t> Nm;
    const vector<size_t> Np;
    const vector<double> mass;
    const vector<double> charge;

    const Algorithms::AxisBundle<double> axis;
};
//--------------------------------------------------------------
//--------------------------------------------------------------


namespace Setup_Y {  
    
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    typedef exprtk::parser_error::type error_t;

    void checkparse(parser_t& parser, std::string& expression_str, expression_t& expression);
    void parseprofile(const valarray<double>& grid, std::string& str_profile, valarray<double>& profile);
    void parseprofile(const double& input, std::string& str_profile, double& ouput);
    // void startmessages(State1D& Y);
    
    void init_f0(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x, 
    valarray<double>& density, valarray<double>& temperature,const double mass);    

    void init_f1(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x, 
    valarray<double>& density, valarray<double>& temperature,const double mass);    

    void init_f2(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x, 
    valarray<double>& density, valarray<double>& temperature,const double mass);    

//      Initialize the appropriate density and temperature profiles (from the list below)
    void initialize(State1D &Y, Grid_Info &grid);

    void applyexternalfields(Grid_Info &grid, State1D &Y, double time);

}



//**************************************************************

    #endif
