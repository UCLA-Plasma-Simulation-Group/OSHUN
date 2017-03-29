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
        T           pmin(size_t i)  const  { return (*(_p[i].v))[0]-0.5*_dp[i]; }
        T           pmax(size_t i)  const  { return (*(_p[i].v))[Np(i)-1]+0.5*_dp[i]; }
        size_t      pdim()          const  { return _p.size();  }
        T           dp(size_t i)    const  { return _dp[i]; }

//  Global axes px(0) --> species1, px(1) --> species2 ...
        valarray<T> px(size_t i)    const  { return *(_px[i].v); }
        size_t      Npx(size_t i)   const  { return (*(_px[i].v)).size(); }
        T           pxmin(size_t i) const  { return (*(_px[i].v))[0]; }
        T           pxmax(size_t i) const  { return (*(_px[i].v))[Npx(i)-1]; }
        size_t      pxdim()         const  { return _px.size();  }

        //  Global axes py(0) --> species1, py(1) --> species2 ...
        valarray<T> py(size_t i)    const  { return *(_py[i].v); }
        size_t      Npy(size_t i)   const  { return (*(_py[i].v)).size(); }
        T           pymin(size_t i) const  { return (*(_py[i].v))[0]; }
        T           pymax(size_t i) const  { return (*(_py[i].v))[Npy(i)-1]; }
        size_t      pydim()         const  { return _py.size();  }

        //  Global axes py(0) --> species1, py(1) --> species2 ...
        valarray<T> pz(size_t i)    const  { return *(_pz[i].v); }
        size_t      Npz(size_t i)   const  { return (*(_pz[i].v)).size(); }
        T           pzmin(size_t i) const  { return (*(_pz[i].v))[0]; }
        T           pzmax(size_t i) const  { return (*(_pz[i].v))[Npz(i)-1]; }
        size_t      pzdim()         const  { return _pz.size();  }

//  Constructor
        AxisBundle(const vector<T> _xmin,  const vector<T> _xmax,  const vector<size_t> _Nx,     // local  spatial axes
                   const vector<T> _xgmin, const vector<T> _xgmax, const vector<size_t> _Nxg,    // global spatial axes
                   const vector<T> _pmax,  const vector<size_t> _Np,     // momentum axis for each species
                   const vector<size_t> _Npx, const vector<size_t> _Npy, const vector<size_t> _Npz  ){  // dimensions in configuration space
            for (size_t i(0); i < _Nx.size(); ++i) {
                _x.push_back( CAxis<T>( _xmin[i], _xmax[i], _Nx[i]));
            }
            for (size_t i(0); i < _Nxg.size(); ++i) {
                _xg.push_back( CAxis<T>( _xgmin[i], _xgmax[i], _Nxg[i]));
            }
            for (size_t i(0); i < _Nxg.size(); ++i) {
                _dx.push_back((_xmax[i]-_xmin[i])/(static_cast<T>(_Nx[i])));
            }
           // for (size_t i(0); i < _Np.size(); ++i) {
           //     _p.push_back( Axis<T>( _pmax[i]/(static_cast<T>(2 * _Np[i])), _pmax[i], _Np[i]) );
           // }
            for (size_t i(0); i < _Np.size(); ++i) {
                _p.push_back( CAxis<T>(static_cast<T>(0.0), _pmax[i], _Np[i]) );
            }

            for (size_t i(0); i < _Np.size(); ++i) {
                _dp.push_back((_pmax[i])/(static_cast<T>(_Np[i])));
            }

            for (size_t i(0); i < _Npx.size(); ++i) {
                _px.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npx[i]));
            }
            for (size_t i(0); i < _Npy.size(); ++i) {
                _py.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npy[i]));
            }
            for (size_t i(0); i < _Npz.size(); ++i) {
                _pz.push_back( Axis<T>( static_cast<T>(-1.0) * _pmax[i], _pmax[i], _Npz[i]));
            }
        }

        AxisBundle(const AxisBundle& a){
            for (size_t i(0); i < a.xdim(); ++i) {
                _x.push_back( CAxis<T>(a.xmin(i), a.xmax(i), a.Nx(i)) );
            }
            for (size_t i(0); i < a.xgdim(); ++i) {
                _xg.push_back( CAxis<T>(a.xgmin(i), a.xgmax(i), a.Nxg(i)) );
            }
            for (size_t i(0); i < a.xdim(); ++i) {
                _dx.push_back( (a.xgmax(i)-a.xgmin(i))/a.Nxg(i) );
            }
            // for (size_t i(0); i < a.pdim(); ++i) {
            //     _p.push_back( Axis<T>(a.pmax(i)/static_cast<T>(2 * a.Np(i)),a.pmax(i),a.Np(i)) );
            // }
            for (size_t i(0); i < a.pdim(); ++i) {
                _p.push_back( CAxis<T>(a.pmin(i),a.pmax(i),a.Np(i)) );
            }
            for (size_t i(0); i < a.pdim(); ++i) {
                _dp.push_back( a.pmax(i)/a.Np(i) );
            }

            for (size_t i(0); i < a.pxdim(); ++i) {
                _px.push_back( Axis<T>( a.pxmin(i), a.pxmin(i), a.Npx(i)) );
            }
            for (size_t i(0); i < a.pydim(); ++i) {
                _py.push_back( Axis<T>( a.pymin(i), a.pymin(i), a.Npy(i)) );
            }
            for (size_t i(0); i < a.pzdim(); ++i) {
                _pz.push_back( Axis<T>( a.pzmin(i), a.pzmin(i), a.Npz(i)) );
            }
        }

        ~AxisBundle(){ }

    private:
        // vector< Axis<T> > _x, _xg, _p, _px;
        // 
        vector< Axis<T> >  _px, _py, _pz;
        vector< CAxis<T> > _x, _xg, _p;
        vector<T>           _dx, _dp;




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

            for (size_t l = 0; l < ((Nm < Nl) ? Nm+1 : Nl+1); ++l){
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
//       integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
//                   * (x[1]-x[0]);
//       for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
//           integral += q[i] * pow(x[i], p)
//                       * (x[i+1] - x[i-1]);
//       }
//       integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
//                       * (x[q.size()-1]-x[q.size()-2]);
//       return integral*0.5;

        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2]);
        return integral;
    }

    template<class T>
    T moment(const vector<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
//       integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
//                   * (x[1]-x[0]);
//       for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
//           integral += q[i] * pow(x[i], p)
//                       * (x[i+1] - x[i-1]);
//       }
//       integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
//                       * (x[q.size()-1]-x[q.size()-2]);
//       return integral*0.5;
        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2]);
        return integral;
    }
    template<class T>
    T moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
//       integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
//                   * (x[1]-x[0]);
//       for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
//           integral += q[i] * pow(x[i], p)
//                       * (x[i+1] - x[i-1]);
//       }
//       integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
//                       * (x[q.size()-1]-x[q.size()-2]);
//       return integral*0.5;
        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2]);
        return integral;
    }
//--------------------------------------------------------------
// relativistic moments for inverse gamma and gamma
//--------------------------------------------------------------
    template<class T>
    T relativistic_invg_moment(const vector<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
      //             * (x[1]-x[0])
      //             / sqrt(1.0+x[0]*x[0]);
      // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
      //     integral += q[i] * pow(x[i], p)
      //                 * (x[i+1] - x[i-1])
      //                 / sqrt(1.0+x[i]*x[i]);
      // }
      // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
      //                 * (x[q.size()-1]-x[q.size()-2])
      //                 / sqrt(1.0+x[q.size()-1]*x[q.size()-1]);
      // return integral*0.5;

        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0])
                    / sqrt(1.0+x[0]*x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1])
                        / sqrt(1.0+x[i]*x[i]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2])
                    / sqrt(1.0+x[q.size()-1]*x[q.size()-1]);
        return integral;
    }

    template<class T>
    T relativistic_invg_moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
      // integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
      //             * (x[1]-x[0])
      //             / sqrt(1.0+x[0]*x[0]);
      // for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
      //     integral += q[i] * pow(x[i], p)
      //                 * (x[i+1] - x[i-1])
      //                 / sqrt(1.0+x[i]*x[i]);
      // }
      // integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
      //                 * (x[q.size()-1]-x[q.size()-2])
      //                 / sqrt(1.0+x[q.size-1]*x[q.size-1]);
      // return integral*0.5;
        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0])
                    / sqrt(1.0+x[0]*x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1])
                        / sqrt(1.0+x[i]*x[i]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2])
                    / sqrt(1.0+x[q.size-1]*x[q.size-1]);
        return integral;
    }
//--------------------------------------------------------------
    template<class T>
    T relativistic_gamma_moment(const valarray<T> q, const valarray<T> x, const int p){

        T integral(0.0);
// TODO:   Integral values up to the zeroth cell and above the last cell
//       integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
//                   * (x[1]-x[0])
//                   * sqrt(1.0+x[0]*x[0]);
//       for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
//           integral += q[i] * pow(x[i], p)
//                       * (x[i+1] - x[i-1])
//                       * sqrt(1.0+x[i]*x[i]);
//       }
//       integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
//                       * (x[q.size()-1]-x[q.size()-2])
//                       * sqrt(1.0+x[q.size-1]*x[q.size-1]);
//       return integral*0.5;
        integral += q[0] * pow(x[0], p)                   // += Q_0*x_0^p * (x_1-x_0)
                    * (x[1]-x[0])
                    * sqrt(1.0+x[0]*x[0]);
        for (size_t i(1); i < q.size()-1; ++i){           // += Q_i*x_i^p * (x_{i+1}-x_{i-1})
            integral += q[i] * pow(x[i], p)
                        * 0.5*(x[i+1] - x[i-1])
                        * sqrt(1.0+x[i]*x[i]);
        }
        integral += q[q.size()-1] * pow(x[q.size()-1], p) // += Q_n*x_n^p * (x_{n}-x_{n-1})
                    * (x[q.size()-1]-x[q.size()-2])
                    * sqrt(1.0+x[q.size-1]*x[q.size-1]);
        return integral;
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


//--------------------------------------------------------------

    //  Leapfrog space (Position verlet)
    template<class T> class LEAPs {
    public:
//      Constructor
        LEAPs(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum);
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& LEAPs<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum) {
//      Take a step using LEAPspace

//      Initialization
        Y0 = Y;

        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        (*F_momentum)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh
        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        
        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
    template<class T> T& LEAPs<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir) {
//      Take a step using LEAPspace

//      Initialization
        Y0 = Y;

//      First step momentum
        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        (*F_momentum)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh
        (*F_space)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;

        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//--------------------------------------------------------------

    //  Leapfrog space (velocity verlet)
    template<class T> class LEAPv {
    public:
//      Constructor
        LEAPv(T& Yin): Y0(Yin), Yh(Yin) { }

//      Main function
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum);
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
    };

    template<class T> T& LEAPv<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum) {
//      Take a step using LEAPmomentum

//      Initialization
        Y0 = Y;

        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        (*F_space)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh
        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;

        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return Y;
    }
    template<class T> T& LEAPv<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir) {
//      Take a step using LEAPmomentum

//      Initialization
        Y0 = Y;

        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;
        (*F_space)(Y0,Yh);   Yh *= h;            //  p  = h * F_momentum(Y0)
        Y0 += Yh;                                   //  Y0 = Y0 + h*Yh
        (*F_momentum)(Y0,Yh);      Yh *= 0.5*h;            //  x  = h * F_space(Y0)
        Y0  += Yh;

        Y = Y0;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        return Y;
    }
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
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum);
        T& operator()(T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum, size_t dir);

    private:
//      R-K copies for the data
        T  Y0, Yh;
        double xsi, lambda, chi;
    };

    template<class T> T& PEFRL<T>::operator()
            (T& Y, double h, AbstFunctor<T>* F_space, AbstFunctor<T>* F_momentum) {
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



//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can
//  be incremented;

class Clock {
public:
//      Constructor
    Clock(int tout_start, double dt_out, double CFL) {
        _hstep = 0;
        _t_start = double(tout_start)*dt_out;
        _numh  = size_t(static_cast<int>(dt_out/CFL))+1;
        _h     = dt_out/static_cast<double>(_numh);
    }

//      Clock readings
    double h()     const {return _h;}                    // Step size
    size_t numh()  const {return _numh;}                 // # steps
    size_t tick()  const {return _hstep;}                // Current step
    double time()  const {return tick()*h() + _t_start;} // Current time

//      Increment time
    Clock& operator++() { ++_hstep; return *this;}

private:
    double _t_start;   // Initial output timestep
    size_t _numh;      // # of steps derived from this CFL
    double _h;         // resulting time-step from CFL
    size_t _hstep;
};
//--------------------------------------------------------------
//**************************************************************


#endif
