/*! \brief Numerical Methods - Declarations
* \author PICKSC
 * \date   September 1, 2016
 * \file   nmethods.h
 * 
 * This cpp file contains the definitions for the functions
 * required for the numerical methods. 
 */

#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H


//------------------------------------------------------------------------------
/// @brief      Performs Gauss-Seidel method on Ax = b
///
/// @param      A     Input matrix that contains collisional integrals and coefficients
/// @param      b     Vector representing distribution function at current time-step 
/// @param      xk    The solution vector 
///
/// @return     The solution vector representing the distribution function at the next time-step
/// 
///  Fills solution into xk. The other matrices are not modified
///  The function returns "false" if the matrix A is not diagonally
///  dominant
bool Gauss_Seidel(Array2D<double>& A,
                  valarray<complex<double> >& b,
                  valarray<complex<double> >& xk);
bool Gauss_Seidel(Array2D<double>& A,
                  valarray<double >& b,
                  valarray<double >& xk);
//-------------------------------------------------------------------

//-------------------------------------------------------------------

/**
 * @brief      Set up tridiagonal solver
 *
 * @param[in]  a     ???
 * @param[in]  b     ???
 * @param      c     ???
 * @param[in]  d     right side
 * @param      x     solution
 */
void TridiagonalSolve ( size_t calculations_per_loop,
                        const valarray<double>& a,
                        const valarray<double>& b,
                        valarray<double>& c,
                        valarray<complex<double> >  d,
                        valarray<complex<double> >& x);

void TridiagonalSolve ( size_t calculations_per_loop,
                        const valarray<double>& a,
                        const valarray<double>& b,
                        valarray<double>& c,
                        valarray<double>  d,
                        valarray<double>& x);

void TridiagonalSolve (const valarray<double>& a,
                       const valarray<double>& b,
                       valarray<double>& c,
                       valarray<complex<double> >  d,
                       valarray<complex<double> >& x);

void TridiagonalSolve (const valarray<double>& a,
                       const valarray<double>& b,
                       valarray<double>& c,
                       valarray<double>  d,
                       valarray<double>& x);

void TridiagonalSolve (const valarray<complex<double> >& a,
                       const valarray<complex<double> >& b,
                       valarray<complex<double> >& c,
                       valarray<complex<double> >  d,
                       valarray<complex<double> >& x);

//-------------------------------------------------------------------


//------------------------------------------------------------------------------
/// @brief      The tridiagonal solver for implicit collisions
///
/// @param      A     Input matrix 
/// @param      d     Right Side
/// @param      xk    Solution vector
///
/// @return     xk is returned as the solution vector
///
///
bool Thomas_Tridiagonal(size_t calculations_per_loop,
                        Array2D<double>& A,
                        valarray<double>& d,
                        valarray<double>& xk);
bool Thomas_Tridiagonal(size_t calculations_per_loop,
                        Array2D<double>& A,
                        valarray<complex<double> >& d,
                        valarray<complex<double> >& xk);
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<double>& d,
                        valarray<double>& xk);
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<complex<double> >& d,
                        valarray<complex<double> >& xk);
bool Thomas_Tridiagonal(Array2D<complex<double> >& A,
                        valarray<complex<double> >& d,
                        valarray<complex<double> >& xk);

//-------------------------------------------------------------------




//-------------------------------------------------------------------
complex<double> Det33(/*const valarray<double>& D, */
        Array2D<complex <double> >& A);
//-------------------------------------------------------------------
complex<double> Detx33(valarray<complex <double> >& D,
                       Array2D<complex <double> >& A);
//-------------------------------------------------------------------

//-------------------------------------------------------------------
complex<double> Dety33(valarray<complex <double> >& D,
                       Array2D<complex <double> >& A);
//-------------------------------------------------------------------

//-------------------------------------------------------------------
complex<double> Detz33(valarray<complex <double> >& D,
                       Array2D<complex <double> >& A);
//-------------------------------------------------------------------

//  Convert double structure to float structure
valarray<float>             vfloat(const valarray<double>& vDouble);
valarray<complex<float> >   vfloat(const valarray<complex<double> >& vDouble);

valarray<float>             vfloat_real(const valarray<complex<double> >& vDouble);
valarray<float>             vfloat_complex(const valarray<complex<double> >& vDouble);

// Vectors
vector<float>               vfloat(const vector<double>& vDouble);
vector<complex<float> >     vfloat(const vector<complex<double> >& vDouble);

vector<float>               vfloat_real(const vector<complex<double> >& vDouble);
vector<float>               vfloat_complex(const vector<complex<double> >& vDouble);



// Conver float to double
// 
valarray<double>            vdouble(const valarray<float>& vFloat);
valarray<complex<double> >  vdouble(const valarray<complex<float> >& vFloat);

vector<double>              vdouble(const vector<float>& vFloat);
vector<complex<double> >    vdouble(const vector<complex<float> >& vFloat);


// Convert double to double
vector<double>              valtovec(const valarray<double>& vDouble);
vector<double>               vdouble_real(const vector<complex<double> >& vDouble);
vector<double>               vdouble_imag(const vector<complex<double> >& vDouble);


valarray<complex<double> >  vdoubletocomplex(const valarray<double>& vDouble);
//// Gradients


valarray<double> df_4thorder(const valarray<double>& f);
valarray<double> df_4thorder(valarray<double>& f);

Array2D<complex<double> > df1_4thorder(Array2D<complex<double> >& f);
Array2D<complex<double> > df2_4thorder(Array2D<complex<double> >& f);

#endif
