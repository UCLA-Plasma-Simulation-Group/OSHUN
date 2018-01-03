/*! \brief %Algorithms
 * \author PICKSC
 * \date   March 2017
 * \file   lib-algorithms.h
 * 
 * This file includes
 * - Algorithms namespace that contains MakeAxis and Axis
 * - AxisBundle which contains global axes information
 * - %algorithms to calculate Legendre polynomials
 * - moments of the distribution function
 * - RK2, RK3, and RK4 definitions
 */
//--------------------------------------------------------------
#ifndef ALGORITHM_LIBRARY_H
#define ALGORITHM_LIBRARY_H

//**************************************************************
namespace Algorithms {



//--------------------------------------------------------------
    template<typename T>
    valarray<T> MakeAxis(const T min, const T max, const size_t N){
        valarray<T> v(N);
        for (size_t i(0); i < N; ++i) {
            v[i] = static_cast<T>(i);
        }
        v *= (max-min)/(static_cast<T>(N-1));
        v += min;
        return v;
    }

    template<typename T>
    valarray<T> MakeCAxis(const T min, const T max, const size_t N){
        valarray<T> v(N);
        for (size_t i(0); i < N; ++i) {
            v[i] = static_cast<T>(i);
        }
        v *= (max-min)/(static_cast<T>(N));
        v += min+0.5*(max-min)/(static_cast<T>(N));
        return v;
    }

    template<typename T>
    valarray<T> MakeCAxis(const T min, const T max, const T delta){
        size_t N = (max-min)/delta;
        valarray<T> v(N);
        for (size_t i(0); i < N; ++i) {
            v[i] = static_cast<T>(i);
        }
        v *= delta;
        v += min+0.5*delta;
        return v;
    }

    template<typename T>
    valarray<T> MakeCAxis(const T min, const valarray<T> delta_grid)
    {
        size_t N = delta_grid.size();
        valarray<T> v(N);
        v[0] = min + 0.5*delta_grid[0];
        for (size_t i(1); i < N; ++i) {
            v[i]  = delta_grid[i];
            v[i] += v[i-1];
            // std::cout << "\n vR[" << i << "] = " << v[i] << std::endl;
        }
        return v;
    }

    template<typename T> class Axis {
    public:
        valarray<T>* v;

//  Constructor
        Axis(const T min=0, const T max=0, const size_t N=1){
            if (N > 1) {
                v = new valarray<T>(MakeAxis(min,max,N));
            }
            else v = NULL;
        }
//  Copy constructor

        Axis(const Axis& other){
            if (other.v != NULL) {
                v = new valarray<T>(*(other.v));
            }
            else v = NULL;
        }

//  Destructor 
        ~Axis() { delete v; }
    };

    template<typename T> class CAxis {
    public:
        valarray<T>* v;

//  Constructor
        CAxis(const T min=0, const T max=0, const size_t N=1){
            if (N > 1) {
                v = new valarray<T>(MakeCAxis(min,max,N));
            }
            else v = NULL;
        }
        CAxis(const T min=0, const T max=0, const T delta=1){
            size_t N = (max-min)/delta;
            if (N > 1) {
                v = new valarray<T>(MakeCAxis(min,max,delta));
            }
            else v = NULL;
        }
        CAxis(const T min=0, const valarray<T> delta_grid=1){
            if (delta_grid.size() > 1) {
                v = new valarray<T>(MakeCAxis(min,delta_grid));
            }
            else v = NULL;
        }

//  Copy constructor

        CAxis(const CAxis& other){
            if (other.v != NULL) {
                v = new valarray<T>(*(other.v));
            }
            else v = NULL;
        }

//  Destructor 
        ~CAxis() { delete v; }
    };

    template<typename T> class AxisBundle {
    public:
//  Local  axes x(0) --> x, x(1) --> y, x(2) --> z            
        valarray<T> x(size_t i)     const { return *(_x[i].v); }
        size_t      Nx(size_t i)    const  { return (*(_x[i].v)).size(); }
        T           xmin(size_t i)  const  { return (*(_x[i].v))[0]-0.5*_dx[i]; }
        T           xmax(size_t i)  const  { return (*(_x[i].v))[Nx(i)-1]+0.5*_dx[i]; }
        size_t      xdim()          const  { return _x.size();  }
        T           dx(size_t i)    const  { return _dx[i]; }

//  Global axes xg(0) --> x, xg(1) --> y, xg(2) --> z            
        valarray<T> xg(size_t i)    const  { return *(_xg[i].v); }
        size_t      Nxg(size_t i)   const  { return (*(_xg[i].v)).size(); }
        T           xgmin(size_t i) const  { return (*(_xg[i].v))[0]-0.5*_dx[i]; }
        T           xgmax(size_t i) const  { return (*(_xg[i].v))[Nxg(i)-1]+0.5*_dx[i]; }
        size_t      xgdim()         const  { return _xg.size();  }

//  Global axes p(0) --> species1, p(1) --> species2 ...
        valarray<T> p(size_t i)     const  { return *(_p[i].v); }
        size_t      Np(size_t i)    const  { return (*(_p[i].v)).size(); }
        T           pmin(size_t i)  const  { return (*(_p[i].v))[0]-0.5*_dp[i][0]; }
        T           pmax(size_t i)  const  { return (*(_p[i].v))[Np(i)-1]+0.5*_dp[i][Np(i)-1]; }
        size_t      pdim()          const  { return _p.size();  }
        // T           dp(size_t i)    const  { return _dp[i]; }
        // Non-uniform grid
        valarray<T> dp(size_t i)    const  { return (_dp[i]); }

//  Global axes px(0) --> species1, px(1) --> species2 ...
        valarray<T> px(size_t i)    const  { return *(_px[i].v); }
        size_t      Npx(size_t i)   const  { return (*(_px[i].v)).size(); }
        T           pxmin(size_t i) const  { return (*(_px[i].v))[0]; }
        T           pxmax(size_t i) const  { return (*(_px[i].v))[Npx(i)-1]+0.5*_dpx[i][Npx(i)-1]; }
        size_t      pxdim()         const  { return _px.size();  }
        valarray<T> dpx(size_t i)    const  { return (_dpx[i]); }

        //  Global axes py(0) --> species1, py(1) --> species2 ...
        valarray<T> py(size_t i)    const  { return *(_py[i].v); }
        size_t      Npy(size_t i)   const  { return (*(_py[i].v)).size(); }
        T           pymin(size_t i) const  { return (*(_py[i].v))[0]; }
        T           pymax(size_t i) const  { return (*(_py[i].v))[Npy(i)-1]+0.5*_dpy[i][Npy(i)-1]; }
        size_t      pydim()         const  { return _py.size();  }
        valarray<T> dpy(size_t i)    const  { return (_dpy[i]); }

        //  Global axes py(0) --> species1, py(1) --> species2 ...
        valarray<T> pz(size_t i)    const  { return *(_pz[i].v); }
        size_t      Npz(size_t i)   const  { return (*(_pz[i].v)).size(); }
        T           pzmin(size_t i) const  { return (*(_pz[i].v))[0]; }
        T           pzmax(size_t i) const  { return (*(_pz[i].v))[Npz(i)-1]+0.5*_dpz[i][Npz(i)-1]; }
        size_t      pzdim()         const  { return _pz.size();  }
        valarray<T> dpz(size_t i)    const  { return (_dpz[i]); }

//  Constructor
        AxisBundle( const vector<T> _xmin,  const vector<T> _xmax,  const vector<size_t> _Nx,     // local  spatial axes
                    const vector<T> _xgmin, const vector<T> _xgmax, const vector<size_t> _Nxg,    // global spatial axes
                    const vector<valarray<T> > _deltap, 
                    const vector<valarray<T> > _deltapx, const vector<valarray<T> > _deltapy, const vector<valarray<T> > _deltapz){     // momentum axis for each species
                    // dimensions in configuration space
            
            
            for (size_t i(0); i < _Nx.size(); ++i) {
                // std::cout << "\n xmin = " << _xmin[i] << " \n ";
                // std::cout << "\n xmax = " << _xmax[i] << " \n ";
                // std::cout << "\n Nx = " << _Nx[i] << " \n ";
                _x.push_back( CAxis<T>( _xmin[i], _xmax[i], _Nx[i]));
            }
            // std::cout << "\n _Nxg.size() = " << _Nxg.size() << " \n ";
            for (size_t i(0); i < _Nxg.size(); ++i) {
                // std::cout << "\n xgmin = " << _xgmin[i] << " \n ";
                // std::cout << "\n xgmax = " << _xgmax[i] << " \n ";
                // std::cout << "\n Nxg = " << _Nxg[i] << " \n ";
                _xg.push_back( CAxis<T>( _xgmin[i], _xgmax[i], _Nxg[i]));
            }
            for (size_t i(0); i < _Nxg.size(); ++i) {

                _dx.push_back((_xmax[i]-_xmin[i])/(static_cast<T>(_Nx[i])));

            }
            // for (size_t i(0); i < _Np.size(); ++i) {
            //     _p.push_back( Axis<T>( _pmax[i]/(static_cast<T>(2 * _Np[i])), _pmax[i], _Np[i]) );
            // }
            for (size_t i(0); i < _deltap.size(); ++i) {
                // _p.push_back( CAxis<T>(static_cast<T>(0.0), _pmax[i], _Np[i]) );
                // Non-uniform grid
                _p.push_back( CAxis<T>(static_cast<T>(0.0), _deltap[i]) );
            }

            for (size_t i(0); i < _deltap.size(); ++i) {
                // _dp.push_back((_pmax[i])/(static_cast<T>(_Np[i])));
                // Non-uniform grid
                _dp.push_back(_deltap[i]);
            }

            for (size_t i(0); i < _deltapx.size(); ++i) {
                // _dp.push_back((_pmax[i])/(static_cast<T>(_Np[i])));
                // Non-uniform grid
                _dpx.push_back(_deltapx[i]);
            }

            for (size_t i(0); i < _deltapy.size(); ++i) {
                // _dp.push_back((_pmax[i])/(static_cast<T>(_Np[i])));
                // Non-uniform grid
                _dpy.push_back(_deltapy[i]);
            }

            for (size_t i(0); i < _deltapz.size(); ++i) {
                // _dp.push_back((_pmax[i])/(static_cast<T>(_Np[i])));
                // Non-uniform grid
                _dpz.push_back(_deltapz[i]);
            }

            for (size_t i(0); i < _deltapx.size(); ++i) {
                // _px.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npx[i]));
                // _px.push_back( Axis<T>( static_cast<T>(-1.0) * pmax(i), pmax(i), _Npx[i]));
                _px.push_back( CAxis<T>( static_cast<T>(-1.0) * pmax(i), _deltapx[i]));
            }
            for (size_t i(0); i < _deltapy.size(); ++i) {
                // _py.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npy[i]));
                // _py.push_back( Axis<T>( static_cast<T>(-1.0) * pmax(i), pmax(i), _Npy[i]));
                _py.push_back( CAxis<T>( static_cast<T>(-1.0) * pmax(i), _deltapy[i]));
            }
            for (size_t i(0); i < _deltapz.size(); ++i) {
                // _pz.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npz[i]));
                // _pz.push_back( Axis<T>( static_cast<T>(-1.0) * pmax(i), pmax(i), _Npz[i]));
                _pz.push_back( CAxis<T>( static_cast<T>(-1.0) * pmax(i), _deltapz[i]));
            }
        }

        AxisBundle(const AxisBundle& a){
            // std::cout << "\n 12 \n ";
            for (size_t i(0); i < a.xdim(); ++i) {
                // std::cout << "\n a.xmin = " << a.xmin(i) << " \n ";
                // std::cout << "\n a.xmax = " << a.xmax(i) << " \n ";
                // std::cout << "\n a.Nx = " << a.Nx(i) << " \n ";
                _x.push_back( CAxis<T>(a.xmin(i), a.xmax(i), a.Nx(i)) );
            }

            for (size_t i(0); i < a.xgdim(); ++i) {
                // std::cout << "\n a.xgmin = " << a.xgmin(i) << " \n ";
                // std::cout << "\n a.xgmax = " << a.xgmax(i) << " \n ";
                // std::cout << "\n a.Nxg = " << a.Nxg(i) << " \n ";
                _xg.push_back( CAxis<T>(a.xgmin(i), a.xgmax(i), a.Nxg(i)) );
            }
            // std::cout << "\n 14 \n ";
            for (size_t i(0); i < a.xdim(); ++i) {
                _dx.push_back( (a.xgmax(i)-a.xgmin(i))/a.Nxg(i) );
            }
            // for (size_t i(0); i < a.pdim(); ++i) {
            //     _p.push_back( Axis<T>(a.pmax(i)/static_cast<T>(2 * a.Np(i)),a.pmax(i),a.Np(i)) );
            // }
            // std::cout << "\n 15 \n ";
            for (size_t i(0); i < a.pdim(); ++i) {
                // _p.push_back( CAxis<T>(a.pmin(i),a.pmax(i),a.Np(i)) );
                _p.push_back( CAxis<T>(a.pmin(i),a.dp(i)) );
            }
            for (size_t i(0); i < a.pdim(); ++i) {
                // _dp.push_back( a.pmax(i)/a.Np(i) );
                _dp.push_back( a.dp(i) );
            }

            for (size_t i(0); i < a.pxdim(); ++i) {
                // _px.push_back( Axis<T>( static_cast<T>(-1.0) * a.pmax(i), a.pmax(i), a.Npx(i)) );
                _px.push_back( CAxis<T>( static_cast<T>(-1.0) * a.pmax(i), a.dpx(i)) );
            }
            for (size_t i(0); i < a.pydim(); ++i) {
                _py.push_back( CAxis<T>( static_cast<T>(-1.0) * a.pmax(i), a.dpy(i)) );
            }
            for (size_t i(0); i < a.pzdim(); ++i) {
                _pz.push_back( CAxis<T>( static_cast<T>(-1.0) * a.pmax(i), a.dpz(i)) );
            }
        }

        ~AxisBundle(){ }

    private:
        // vector< Axis<T> > _x, _xg, _p, _px;
        // 
        vector< CAxis<T> >  _px, _py, _pz;
        vector< CAxis<T> > _x, _xg, _p;
        vector<T>           _dx;//, _dp;
        // Non-uniform grid
        vector<valarray<T> > _dp;
        vector<valarray<T> > _dpx,_dpy,_dpz;
        // vector<valarray<T>> _dx, _dp;
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
// LEGENDRE POLYNOMIALS
// Calculate the legendre polynomials using the recurrance relations
// For m0 = 0
    template<class T>
    valarray<T> Legendre(const T x, const size_t Nl){ // where Nl is to denote l0+1

        valarray<T> P_Legendre(0.0, Nl);
        if ( abs(x) > 1 ) return P_Legendre;

        T r1, sqrtx = sqrt(1.0-x*x);

//      Initialization 
        P_Legendre[0] = 1.0;
        P_Legendre[1] = x;

        for (size_t l(1); l < Nl - 1; ++l){
            r1 = 1.0 / double(l + 1);
            P_Legendre[l+1] = P_Legendre[l] * (x*(2.0*l+1.0) * r1) -
                              P_Legendre[l-1]*(double(l) * r1);
        }

        return P_Legendre;
    }

// For m0 > 0
    template<class T>
    Array2D<T> Legendre(const T x, const size_t Nl, const size_t Nm){

//      Local variables
        if (Nl < Nm) {
            cout << "ERROR: " << "l0+1 = " << Nl << " < " << Nm << " = m0+1\n";
            exit(1);
        }

        Array2D<T> P_Legendre(Nl+1,Nm+1);

        if ( abs(x) > 1 ) {
//             cout << "ERROR: " << "|" << x << "| > 1 is not a valid cosine\n";
//      cout << "setting values to 0.0";
            P_Legendre = 0.0;

        }
        else
        {

            T r1, sqrtx = sqrt(1.0-x*x), fact = 1.0;

            //      Initialization
            P_Legendre = 0.0;
            P_Legendre(0,0) = 1.0;

            for (size_t l = 1; l < Nm+1; ++l){
                // std::cout << "\n l=" << l;
                P_Legendre(l,l) = - P_Legendre(l-1,l-1)*(fact*sqrtx);
                fact += 2.0;
            }
            // std::cout << "\n ";

            for (size_t l = 0; l < ((Nm < Nl) ? Nm+1 : Nl); ++l){
                // std::cout << "\n l=" << l;
                P_Legendre(l+1,l) = P_Legendre(l,l)*(x*(2.0*l+1.0));
            }
            // std::cout << "\n ";

            for (size_t m = 0; m < Nm+1; ++m){
                for (size_t l = m+1; l < Nl ; ++l){
                    // std::cout << "\n l=" << l << ", m=" << m;
                    r1 = 1.0 / double(l - m + 1);
                    P_Legendre(l+1,m) = P_Legendre(l,m) * (x*(2.0*l+1.0) * r1) -
                                        P_Legendre(l-1,m)*(double(l+m) * r1);
                }
            }
        }

        // std::cout << "\n leaving \n ";

        // exit(1);

        return P_Legendre;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
// MOMENTS
// p-th moment of a quantity q(x)
    template<class T>
    T moment(const vector<T> q, const vector<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * (pow(x[0], p))                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2]);
      return integral*0.5;

        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (0.5*((x[1]-x[0])+x[0]));
        // for (size_t i(1); i < q.size(); ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += 0.5*(q[i] * pow(x[i], p) + q[i-1] * pow(x[i-1], p))
        //                 *(x[i] - x[i-1]);
        // }
        // return integral;
        
// 
        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (x[1]-x[0]);
        // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += q[i] * pow(x[i], p)
        //                 * 0.5*(x[i+1] - x[i-1]);
        // }
        // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
        //             * (x[q.size()-1]-x[q.size()-2]);
        // return integral;
    }

    template<class T>
    T moment(const vector<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * (pow(x[0], p))                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2]);
      return integral*0.5;


      // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
      //               * (0.5*((x[1]-x[0])+x[0]));
      //   for (size_t i(1); i < q.size(); ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
      //       integral += 0.5*(q[i] * pow(x[i], p) + q[i-1] * pow(x[i-1], p))
      //                   *(x[i] - x[i-1]);
      //   }
      //   return integral;

        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (x[1]-x[0]);
        // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += q[i] * pow(x[i], p)
        //                 * 0.5*(x[i+1] - x[i-1]);
        // }
        // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
        //             * (x[q.size()-1]-x[q.size()-2]);
        // return integral;
    }
    template<class T>
    T moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * (pow(x[0], p))                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2]);
      return integral*0.5;


      // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
      //               * (0.5*((x[1]-x[0])+x[0]));
      //   for (size_t i(1); i < q.size(); ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
      //       integral += 0.5*(q[i] * pow(x[i], p) + q[i-1] * pow(x[i-1], p))
      //                   *(x[i] - x[i-1]);
      //   }
      //   return integral;

    //     integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
    //                 * (x[1]-x[0]);
    //     for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
    //         integral += q[i] * pow(x[i], p)
    //                     * 0.5*(x[i+1] - x[i-1]);
    //     }
    //     integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
    //                 * (x[q.size()-1]-x[q.size()-2]);
    //     return integral;
    }
//--------------------------------------------------------------
// relativistic moments for inverse gamma and gamma
//--------------------------------------------------------------
    template<class T>
    T relativistic_invg_moment(const vector<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0])
                  / sqrt(1.0+x[0]*x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1])
                      / sqrt(1.0+x[i]*x[i]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2])
                      / sqrt(1.0+x[q.size()-1]*x[q.size()-1]);
      return integral*0.5;

        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (x[1]-x[0])
        //             / sqrt(1.0+x[0]*x[0]);
        // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += q[i] * pow(x[i], p)
        //                 * 0.5*(x[i+1] - x[i-1])
        //                 / sqrt(1.0+x[i]*x[i]);
        // }
        // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
        //             * (x[q.size()-1]-x[q.size()-2])
        //             / sqrt(1.0+x[q.size()-1]*x[q.size()-1]);
        // return integral;
    }

    template<class T>
    T relativistic_invg_moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0])
                  / sqrt(1.0+x[0]*x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1])
                      / sqrt(1.0+x[i]*x[i]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2])
                      / sqrt(1.0+x[q.size-1]*x[q.size-1]);
      return integral*0.5;
        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (x[1]-x[0])
        //             / sqrt(1.0+x[0]*x[0]);
        // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += q[i] * pow(x[i], p)
        //                 * 0.5*(x[i+1] - x[i-1])
        //                 / sqrt(1.0+x[i]*x[i]);
        // }
        // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
        //             * (x[q.size()-1]-x[q.size()-2])
        //             / sqrt(1.0+x[q.size-1]*x[q.size-1]);
        // return integral;
    }
//--------------------------------------------------------------
    template<class T>
    T relativistic_gamma_moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                  * (x[1]-x[0])
                  * sqrt(1.0+x[0]*x[0]);
      for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
          integral += q[i] * pow(x[i], p)
                      * (x[i+1] - x[i-1])
                      * sqrt(1.0+x[i]*x[i]);
      }
      integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                      * (x[q.size()-1]-x[q.size()-2])
                      * sqrt(1.0+x[q.size-1]*x[q.size-1]);
      return integral*0.5;
        // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
        //             * (x[1]-x[0])
        //             * sqrt(1.0+x[0]*x[0]);
        // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
        //     integral += q[i] * pow(x[i], p)
        //                 * 0.5*(x[i+1] - x[i-1])
        //                 * sqrt(1.0+x[i]*x[i]);
        // }
        // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
        //             * (x[q.size()-1]-x[q.size()-2])
        //             * sqrt(1.0+x[q.size-1]*x[q.size-1]);
        // return integral;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
//  RUNGE-KUTTA METHODS
//--------------------------------------------------------------
    template<typename T> class AbstFunctor {
//  abstract functor interface 
    public:
        virtual void operator()(const T& Yin, T& Yslope)=0;  // call using operator
        virtual void operator()(const T& Yin, T& Yslope, size_t dir)=0;  // call using operator
    };


//  RK2 
    template<class T> class RK2 {
    public:
//      Constructor
        RK2(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F);
        T& operator()(T& Y, double h, AbstFunctor<T>* F, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& RK2<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F) {
//      Take a step using RK2

//      Initialization
        Y0 = Y;

//      Step 1
        (*F)(Y0,Yh); Yh *= h;     // Yh = h*F(Y0)
        Y0 += Yh;                 // Y0 = Y0 + h*Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2)*F(Y0)      
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y0,Yh); Yh *= 0.5*h; // Yh = (h/2)*F(Y0)
        Y += Yh;                  // Y = Y + (h/2)*F(Y0)
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
    template<class T> T& RK2<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F, size_t dir) {
//      Take a step using RK2

//      Initialization
        Y0 = Y;

//      Step 1
        (*F)(Y0,Yh,dir); Yh *= h;     // Yh = h*F(Y0)
        Y0 += Yh;                 // Y0 = Y0 + h*Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2)*F(Y0)      
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y0,Yh,dir); Yh *= 0.5*h; // Yh = (h/2)*F(Y0)
        Y += Yh;                  // Y = Y + (h/2)*F(Y0)
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  RK3 (Osher&Shu Method) 
    template<class T> class RK3 {
    public:
//      Constructor
        RK3(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F);
        T& operator()(T& Y, double h, AbstFunctor<T>* F, size_t dir);


    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& RK3<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F) {
//      Take a step using RK3

//      Initialization
        Y0 = Y;

//      Step 1
        (*F)(Y0,Yh); Yh *= h;                         // Yh = h*F(Y0)
        Y0 += Yh;
//      Y0 = Y0 + h*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y0,Yh); Yh *= h;                      // Yh = h*F(Y0)
        Y0 += Yh; Y *= 3.0; Y0 += Y;  Y0 *= 0.25;  // Changed Y to 3*Y!
//      Y0 = 1/4 * ( 3*Y + (Y0 + h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        (*F)(Y0,Yh); Yh *= h;                        // Yh = h*F(Y0)
        Y0 += Yh; Y0 *= 6.0; Y += Y0; Y *= (1.0/9.0);
//      Y  = 1/3 * ( Y + 2 * (Y0+h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        return Y;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    template<class T> T& RK3<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F, size_t dir) {
//      Take a step using RK3

//      Initialization
        Y0 = Y;

//      Step 1
        (*F)(Y0,Yh,dir); Yh *= h;                         // Yh = h*F(Y0)
        Y0 += Yh;
//      Y0 = Y0 + h*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y0,Yh,dir); Yh *= h;                      // Yh = h*F(Y0)
        Y0 += Yh; Y *= 3.0; Y0 += Y;  Y0 *= 0.25;  // Changed Y to 3*Y!
//      Y0 = 1/4 * ( 3*Y + (Y0 + h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        (*F)(Y0,Yh,dir); Yh *= h;                        // Yh = h*F(Y0)
        Y0 += Yh; Y0 *= 6.0; Y += Y0; Y *= (1.0/9.0);
//      Y  = 1/3 * ( Y + 2 * (Y0+h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        return Y;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


//  RK4 
    template<class T> class RK4 {
    public:
//      Constructor
        RK4(T& Yin): Y0(Yin), Y1(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F);
        T& operator()(T& Y, double h, AbstFunctor<T>* F, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Y1, Yh;
    };

    template<class T> T& RK4<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F) {
//      Take a step using RK4

//      Initialization
        Y0 = Y; Y1 = Y;

//      Step 1
        (*F)(Y1,Yh);                    // slope in the beginning
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y1 + (h/2)*Yh
        Yh *= (1.0/3.0); Y  += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y1,Yh);     Y1  = Y0;      // slope in the middle
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y0 + (h/2)*Yh
        Yh *= (2.0/3.0); Y  += Yh;      // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        (*F)(Y1,Yh);                    // slope in the middle again
        Yh *= h;          Y0 += Yh;     // Y0 = Y0 + h*Yh
        Yh *= (1.0/3.0);  Y  += Yh;     // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 4
        (*F)(Y0,Yh);                    // slope at the end
        Yh *= (h/6.0);    Y += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        return Y;
    }
    template<class T> T& RK4<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F, size_t dir) {
//      Take a step using RK4

//      Initialization
        Y0 = Y; Y1 = Y;

//      Step 1
        (*F)(Y1,Yh,dir);                    // slope in the beginning
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y1 + (h/2)*Yh
        Yh *= (1.0/3.0); Y  += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        (*F)(Y1,Yh,dir);     Y1  = Y0;      // slope in the middle
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y0 + (h/2)*Yh
        Yh *= (2.0/3.0); Y  += Yh;      // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        (*F)(Y1,Yh,dir);                    // slope in the middle again
        Yh *= h;          Y0 += Yh;     // Y0 = Y0 + h*Yh
        Yh *= (1.0/3.0);  Y  += Yh;     // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 4
        (*F)(Y0,Yh,dir);                    // slope at the end
        Yh *= (h/6.0);    Y += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        return Y;
    }


//  RK4 
    template<class T> class RKCK54 {
    public:
//      Constructor
        RKCK54(T& Yin): Yh1(Yin), Yh3(Yin), Yh4(Yin), Yh5(Yin), Yh6(Yin), Yt(Yin) { }

//      Main function
        T& operator()(T& Y5, T&Y4, double h, AbstFunctor<T>* F);
        // T& operator()(T& Y5, T&Y4, double h, AbstFunctor<T>* F, size_t dir);

    private:
//      R-K copies for the data
        T  Yh1, Yh3, Yh4, Yh5, Yh6, Yt;
        double a21 = 0.2;
        double a31 = 3./40., a32 = 9./40.;
        double a41 = .3, a42 = -.9, a43 = 1.2;
        double a51 = -11./54., a52 = 2.5 ,a53 = -70./27. ,a54 = 35./27.;
        double a61 = 1631./55296., a62 = 175./512., a63 = 575./13824., a64 = 44275./110592., a65 = 253./4096.;

        double b1_5 = 37./378., b3_5 = 250./621., b4_5 = 125./594., b6_5 = 512./1771.;
        double b1_4 = 2825./27648., b3_4 = 18575./48384., b4_4 = 13525./55296., b5_4 = 277./14336., b6_4 = 0.25;
    };

    template<class T> T& RKCK54<T>::operator()
            (T& Y5, T& Y4, double h, AbstFunctor<T>* F) {
//      Take a step using RKCK54

//      Initialization
        // Yh1 = Y5;   Yh3 = Y5;   Yh4 = Y5;   Yh5 = Y5;   Yh6 = Y5;
        // Yt = Y5;

//      Step 1
        (*F)(Y4,Yh1); Yh1 *= h;
        Yh1 *= a21;
        Yt = Y4;    Yt  += Yh1;                              // Y1 = Y1 + (h/5)*Yh

        //      Step 2
        // (*F)(Yt,Yh2);                                   // f(Y1)
        (*F)(Yt,Y5); Y5 *= h;                                   // f(Y1)
        Yh1 *= a31/a21; Y5 *= a32;  
        Yt = Y4;    Yt += Yh1;  Yt += Y5;

        //      Step 3
        (*F)(Yt,Yh3); Yh3 *= h;
        Yh1 *= a41/a31; Y5 *= a42/a32;   Yh3 *= a43;
        Yt = Y4;    Yt += Yh1; Yt += Y5; Yt += Yh3;
        
        //      Step 4
        (*F)(Yt,Yh4); Yh4 *= h;
        Yh1 *= a51/a41; Y5 *= a52/a42;   Yh3 *= a53/a43;    Yh4 *= a54;
        Yt = Y4;    Yt += Yh1; Yt += Y5; Yt += Yh3; Yt += Yh4;
        
        //      Step 5
        (*F)(Yt,Yh5); Yh5 *= h;
        Yh1 *= a61/a51; Y5 *= a62/a52;   Yh3 *= a63/a53;    Yh4 *= a64/a54; Yh5 *= a65;
        Yt = Y4;    Yt += Yh1; Yt += Y5; Yt += Yh3; Yt += Yh4; Yt += Yh5;
            

        //      Step 6
        (*F)(Yt,Yh6); Yh6 *= h;


        //      Assemble 5th order solution
        Y5 = Y4;
        Yh1 *= b1_5/a61;    Y5 += Yh1;
        Yh3 *= b3_5/a63;    Y5 += Yh3;
        Yh4 *= b4_5/a64;    Y5 += Yh4;
        Yh6 *= b6_5;        Y5 += Yh6;

        //      Assemble 4th order solution
        Yh1 *= b1_4/b1_5;   Y4 += Yh1;
        Yh3 *= b3_4/b3_5;   Y4 += Yh3;
        Yh4 *= b4_4/b4_5;   Y4 += Yh4;
        Yh5 *= b5_4/a65;    Y4 += Yh5;
        Yh6 *= b6_4/b6_5;   Y4 += Yh6;

//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // return Y;
    }
    
//     template<class T> T& RKCK54<T>::operator()
//             (T& Y5, T&Y4, double h, AbstFunctor<T>* F, size_t dir) {
// //      Take a step using RKCK54

// //      Initialization
//         Y0 = Y; Y1 = Y;

// //      Step 1
//         (*F)(Y1,Yh,dir);                    // slope in the beginning
//         Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y1 + (h/2)*Yh
//         Yh *= (1.0/3.0); Y  += Yh;      // Y  = Y  + (h/6)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 2
//         (*F)(Y1,Yh,dir);     Y1  = Y0;      // slope in the middle
//         Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y0 + (h/2)*Yh
//         Yh *= (2.0/3.0); Y  += Yh;      // Y  = Y  + (h/3)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 3
//         (*F)(Y1,Yh,dir);                    // slope in the middle again
//         Yh *= h;          Y0 += Yh;     // Y0 = Y0 + h*Yh
//         Yh *= (1.0/3.0);  Y  += Yh;     // Y  = Y  + (h/3)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 4
//         (*F)(Y0,Yh,dir);                    // slope at the end
//         Yh *= (h/6.0);    Y += Yh;      // Y  = Y  + (h/6)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

//         return Y;
//     }


//--------------------------------------------------------------

    //  Leapfrog space (Position verlet)
    template<class T> class LEAPs {
    public:
//      Constructor
        LEAPs(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, AbstFunctor<T>* F_field);
        // T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& LEAPs<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum,
                AbstFunctor<T>* F_field) {
//      Take a step using LEAPspace

//      Initialization
        Y0 = Y;

        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x*  = h/2 * F_space(Y0)
        Y0  += Yh;
        (*F_field)(Y0,Yh);      Yh *= 0.5*h;            //  E*  = h/2 * F_field(Y0(x*))
        Y0  += Yh;
        (*F_momentum)(Y0,Yh);   Yh *= h;                //  p  = h * F_momentum(Y0(x*,E*))
        Y0  += Yh;                                   
        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h/2 * F_space(Y0(p))
        Y0  += Yh;
        (*F_field)(Y0,Yh);      Yh *= 0.5*h;            //  E*  = h/2 * F_field(Y0(x*))
        Y0  += Yh;
        
        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
//     template<class T> T& LEAPs<T>::operator()
//             (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir) {
// //      Take a step using LEAPspace

// //      Initialization
//         Y0 = Y;

// //      First step momentum
//         (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
//         Y0  += Yh;
//         (*F_momentum)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
//         Y0 += Yh;                                   //  Y0 = Y0 + h*Yh
//         (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
//         Y0  += Yh;

//         Y = Y0;
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//         return Y;
//     }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//--------------------------------------------------------------

    //  Leapfrog space (velocity verlet)
    template<class T> class LEAPv {
    public:
//      Constructor
        LEAPv(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, AbstFunctor<T>* F_field);
        // T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& LEAPv<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum,
                AbstFunctor<T>* F_field) {
//      Take a step using LEAPmomentum

//      Initialization
        Y0 = Y;

        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        (*F_space)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh; 
        (*F_field)(Y0,Yh);      Yh *= h;            //  E*  = h/2 * F_field(Y0(x*))
        Y0  += Yh;
        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;

        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
//     template<class T> T& LEAPv<T>::operator()
//             (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir) {
// //      Take a step using LEAPmomentum

// //      Initialization
//         Y0 = Y;

        
//         (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  p*  = h/2 * F_momentum(Y0)
//         Y0  += Yh;
        
//         (*F_space)(Y0,Yh);   Yh *= h;                   //  x  = h * F_space(Y0)
//         Y0 += Yh;

//         (*F_field)(Y0,Yh);      Yh *= 0.5*h;            //  E*  = h/2 * F_field(Y0(x*))
//         Y0  += Yh;                                   
//         (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;         //  p  = h * F_space(Y0)
//         Y0  += Yh;
//         (*F_field)(Y0,Yh);      Yh *= 0.5*h;            //  E*  = h/2 * F_field(Y0(x*))
//         Y0  += Yh;

//         Y = Y0;
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//         return Y;
//     }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //  Position Extended Forest-Ruth Like
    template<class T> class PEFRL {
    public:
//      Constructor
        PEFRL(T& Yin): Y0(Yin), Yh(Yin),
                       xsi(0.1786178958448091),
                       lambda(-0.2123418310626054),
                       chi(-0.06626458266981849)
        {}

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, AbstFunctor<T>* F_field);
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
        double xsi, lambda, chi;
    };

    template<class T> T& PEFRL<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum,  AbstFunctor<T>* F_field) {
//      Take a step using PEFRL

//      Initialization
        Y0 = Y;

//      First step space
        (*F_space)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        // (*F_field)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        // Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= (1.0-2.0*lambda)*0.5*h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        // (*F_field)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        // Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= lambda*h;             //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                           //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= (1.0-2.0*(chi+xsi))*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        // (*F_field)(Y0,Yh);      Yh *= (1.0-2.0*(chi+xsi))*h;                              //  x  = h * F_space(Y0)
        // Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= lambda*h;             //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                           //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        // (*F_field)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        // Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= (1.0-2.0*lambda)*0.5*h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;

        (*F_space)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        // (*F_field)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        // Y0 += Yh;        

        (*F_field)(Y0,Yh);      Yh *= h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;        
        
        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
    template<class T> T& PEFRL<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir) {
//      Take a step using PEFRL

//      Initialization
        Y0 = Y;

//      First step space
        (*F_space)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= (1.0-2.0*lambda)*0.5*h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= lambda*h;             //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                           //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= (1.0-2.0*(chi+xsi))*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= lambda*h;             //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                           //  Y0 = Y0 + h*Yh

        (*F_space)(Y0,Yh);      Yh *= chi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;

        (*F_momentum)(Y0,Yh);   Yh *= (1.0-2.0*lambda)*0.5*h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;

        (*F_space)(Y0,Yh);      Yh *= xsi*h;                              //  x  = h * F_space(Y0)
        Y0 += Yh;        
        
        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return Y;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
//**************************************************************



#endif
