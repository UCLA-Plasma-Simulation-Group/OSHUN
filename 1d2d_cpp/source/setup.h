/*!\brief  Grid Setup - Declaration
 * \author PICKSC
 * \date   September 1, 2016
 * \file   setup.h
 *
 * Declares a Grid_Info class
 *
 */

    #ifndef DECL_SETUP_H
    #define DECL_SETUP_H

//--------------------------------------------------------------
//--------------------------------------------------------------
class Grid_Info {
//--------------------------------------------------------------
public:
//  Constructor
    Grid_Info(const vector<size_t> _l0,   const vector<size_t> _m0,              
//            local  spatial axes
              const vector<double> _xmin, const vector<double> _xmax,  const vector<size_t> _Nx,
//            global spatial axes
              const vector<double> _xgmin,const vector<double> _xgmax, const vector<size_t> _Nxg,
//                            momentum axes
                                          // const vector<double> _pmax,  const vector<size_t> _Np,
                               const vector<valarray<double> > _dp,
                               const vector<valarray<double> > _dpx,const vector<valarray<double> > _dpy,const vector<valarray<double> > _dpz)
//                                             output momentum axes
              // const vector<size_t> _Npx,   const vector<size_t> _Npy, const vector<size_t> _Npz)
         : l0(_l0), m0(_m0),
           axis( _xmin, _xmax, _Nx, _xgmin, _xgmax, _Nxg, _dp, _dpx, _dpy, _dpz ) {}

//--------------------------------------------------------------
//  Copy constructor
    Grid_Info(const Grid_Info& other): l0(other.l0), m0(other.m0),
                                        axis(other.axis){}
    ~Grid_Info(){};

    const vector<size_t> l0;
    const vector<size_t> m0;

    const Algorithms::AxisBundle<double> axis;
};
//--------------------------------------------------------------
//--------------------------------------------------------------


namespace Setup_Y {
    // void startmessages(State1D& Y);

    void init_f0(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
    valarray<double>& density, valarray<double>& temperature,const double mass, const valarray<double>& pedestal);

    void init_f1(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
    valarray<double>& density, valarray<double>& temperature, valarray<double>& f10x, const SHarmonic1D& f0, const double mass);

    void init_f2(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
    valarray<double>& density, valarray<double>& temperature, valarray<double>& f20x, const double mass);

//      Initialize the appropriate density and temperature profiles (from the list below)
    void initialize(State1D &Y, Grid_Info &grid);

    void applyexternalfields(Grid_Info &grid, State1D &Y, double time);
    void applytravelingwave(Grid_Info &grid, State1D &Y, double time, double stepsize);


    // 2D

    void init_f0(size_t s, SHarmonic2D& h, const valarray<double>& p, const valarray<double>& x, const valarray<double>& y,
    Array2D<double>& density, Array2D<double>& temperature,const double mass, const Array2D<double>& pedestal);

    // void init_f1(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
    // valarray<double>& density, valarray<double>& temperature, valarray<double>& f10x, const SHarmonic1D& f0, const double mass);

    // void init_f2(size_t s, SHarmonic1D& h, const valarray<double>& p, const valarray<double>& x,
    // valarray<double>& density, valarray<double>& temperature, valarray<double>& f20x, const double mass);

//      Initialize the appropriate density and temperature profiles (from the list below)
    void initialize(State2D &Y, Grid_Info &grid);

    void applyexternalfields(Grid_Info &grid, State2D &Y, double time);
    void applytravelingwave(Grid_Info &grid, State2D &Y, double time, double stepsize);

}

//**************************************************************

    #endif
