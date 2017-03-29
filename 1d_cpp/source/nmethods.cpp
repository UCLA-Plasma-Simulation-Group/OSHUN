/*! \brief Numerical Methods - Declarations
* \author PICKSC
 * \date   September 1, 2016
 * \file   nmethods.cpp
 * 
 * This cpp file contains the definitions for the functions
 * required for the numerical methods. 
 */  


//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>
    #include <float.h>

    #include <math.h>
    #include <map>

//  My libraries
    #include "lib-array.h"



//-------------------------------------------------------------------
     bool Gauss_Seidel(Array2D<double>& A, 
                       valarray<complex<double> >& b,
                       valarray<complex<double> >& xk) {
//-------------------------------------------------------------------
        double     tol(1.0e-1);   //< Tolerance for absolute error
        int        MAXiter(5);    //< Maximum iteration allowed


//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
        if ( ( A.dim1() != A.dim2()  ) || 
             ( A.dim1() != b.size()  ) ||
             ( A.dim1() != xk.size() )    )  {
            cout << "Error: The Matrices don't have the right dimensions!" << endl;
            exit(1);
        }
//      Check if the matrix A is diagonally dominant
//      -------------------------------------------------------------
        for (int i(0); i < A.dim1(); ++i){
            double rowi(0.0);
            for (int j(0); j < A.dim2(); ++j){
                rowi += A(i,j);
            }
            if (!(rowi < 2.0*A(i,i))) return false;
        }
//      -------------------------------------------------------------


//      Calculate and invert the diagonal elements once and for all
//      -------------------------------------------------------------      
        valarray<double> invDIAG(A.dim1());
        for (int i(0); i < invDIAG.size(); ++i) invDIAG[i] = 1.0 / A(i,i);

        valarray<complex<double> > xold(xk);
        int iteration(0);         // used to count iterations
        int conv(0);              // used to test convergence

//      Start the iteration loop
//      -------------------------------------------------------------      
        while ( (iteration++ < MAXiter) && (conv < b.size()) ) {

            xold = xk;
            for (int i(0); i < A.dim1(); ++i){
                complex<double> sigma(0.0);    // Temporary sum
                for (int j(0); j < i; ++j){
                    sigma += A(i,j)*xk[j];   
                }
                for (int j(i+1); j < A.dim2(); ++j){
                    sigma += A(i,j)*xk[j];   
                }
                xk[i] = invDIAG[i] * (b[i] - sigma);
            }

            // Calculate Dx = x_old - x_new
            xold -= xk;

//          If the relative error < prescribed tolerance everywhere the method has converged
//          |Dx(i)| < t*|x(i)| + eps 
            conv = 0;
            while ( ( conv < b.size() ) && 
                    ( abs(xold[conv]) < (tol*abs(xk[conv] + 20.0*DBL_MIN)) ) ){ 
                ++conv;
            } 

            //----> Output for testing
            //--------------------------------
            //cout << "iteration = " << iteration << "    ";
            //for (int i(0); i < b.size(); ++i){
            //    cout << xk[i] << "        ";
            //}
            //cout << "\n";
            //--------------------------------

        }

        //----> Output for testing
        //--------------------------------
        // cout << "Iterations = " << iteration-1  <<"\n";
        //for (int i(0); i < b.size(); ++i) {
        //    cout << "Error |Dx| = " << abs(xold[i]) 
        //         << ",    " 
        //         << "Tolerance * |x| = " << tol*abs(xk[i]) <<"\n";
        //}
        //--------------------------------

        return true;
    }
//-------------------------------------------------------------------
//*******************************************************************

//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve (const valarray<double>& a, 
                            const valarray<double>& b, 
                                  valarray<double>& c,      
                                  valarray<complex<double> >  d,
                                  valarray<complex<double> >& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d! 
//-------------------------------------------------------------------
        size_t n(x.size());
	// Modify the coefficients. 
	c[0] /= b[0];                            // Division by zero risk. 
	d[0] /= b[0];                            // Division by zero would imply a singular matrix. 
	for (int i(1); i < n; ++i){
            double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk. 
	    c[i] *= id;	                         // Last value calculated is redundant.
	    d[i] -= d[i-1] * a[i];
	    d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id 
	}
 
	// Now back substitute. 
	x[n-1] = d[n-1];
	for (int i(n-2); i > -1; --i){
	    x[i]  = d[i];
            x[i] -= c[i] * x[i+1];               // x[i] = d[i] - c[i] * x[i + 1];
        }
    }
//-------------------------------------------------------------------

//*******************************************************************
//-------------------------------------------------------------------
void TridiagonalSolve (const valarray<double>& a,
                       const valarray<double>& b,
                       valarray<double>& c,
                       valarray<double>  d,
                       valarray<double>& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d!
//-------------------------------------------------------------------
    size_t n(x.size());
    // Modify the coefficients.
    c[0] /= b[0];                            // Division by zero risk.
    d[0] /= b[0];                            // Division by zero would imply a singular matrix.
    for (int i(1); i < n; ++i){
        double id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk.
        c[i] *= id;	                         // Last value calculated is redundant.
        d[i] -= d[i-1] * a[i];
        d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id
    }

    // Now back substitute.
    x[n-1] = d[n-1];
    for (int i(n-2); i > -1; --i){
        x[i]  = d[i];
        x[i] -= c[i] * x[i+1];               // x[i] = d[i] - c[i] * x[i + 1];
    }
}
//-------------------------------------------------------------------



//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve (const valarray<complex<double>>& a, 
                            const valarray<complex<double>>& b, 
                                  valarray<complex<double>>& c,      
                                  valarray<complex<double>>  d,
                                  valarray<complex<double>>& x) {
//-------------------------------------------------------------------
//   Fills solution into x. Warning: will modify c and d! 
//-------------------------------------------------------------------
        size_t n(x.size());
    // Modify the coefficients. 
    c[0] /= b[0];                            // Division by zero risk. 
    d[0] /= b[0];                            // Division by zero would imply a singular matrix. 
    for (int i(1); i < n; ++i){
           complex<double> id(1.0/(b[i]-c[i-1]*a[i]));   // Division by zero risk. 
        c[i] *= id;                          // Last value calculated is redundant.
        d[i] -= d[i-1] * a[i];
        d[i] *= id;                          // d[i] = (d[i] - d[i-1] * a[i]) * id 
    }
 
    // Now back substitute. 
    x[n-1] = d[n-1];
    for (int i(n-2); i > -1; --i){
        x[i]  = d[i];
            x[i] -= c[i] * x[i+1];               // x[i] = d[i] - c[i] * x[i + 1];
        }
    }
//-------------------------------------------------------------------


//-------------------------------------------------------------------
//*******************************************************************

//-------------------------------------------------------------------
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<double> & d,
                        valarray<double> & xk) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
    if ( ( A.dim1() != A.dim2()  ) ||
         ( A.dim1() != d.size()  ) ||
         ( A.dim1() != xk.size() )    )  {
        cout << "Error: The Matrices don't have the right dimensions!" << endl;
        exit(1);
    }
//      -------------------------------------------------------------

    valarray<double> a(d.size()), b(d.size()), c(d.size());

    for (int i(0); i < A.dim1()-1; ++i){
        a[i+1] = A(i+1,i);
    }
    for (int i(0); i < A.dim1(); ++i){
        b[i] = A(i,i);
    }
    for (int i(0); i < A.dim1()-1; ++i){
        c[i] = A(i,i+1);
    }

//        valarray< double > dcopy(d);
    TridiagonalSolve(a,b,c,d,xk);

    return true;
}

//-------------------------------------------------------------------
bool Thomas_Tridiagonal(Array2D<double>& A,
                        valarray<complex<double> >& d,
                        valarray<complex<double> >& xk) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
    if ( ( A.dim1() != A.dim2()  ) ||
         ( A.dim1() != d.size()  ) ||
         ( A.dim1() != xk.size() )    )  {
        cout << "Error: The Matrices don't have the right dimensions!" << endl;
        exit(1);
    }
//      -------------------------------------------------------------

    valarray<double> a(d.size()), b(d.size()), c(d.size());

    for (int i(0); i < A.dim1()-1; ++i){
        a[i+1] = A(i+1,i);
    }
    for (int i(0); i < A.dim1(); ++i){
        b[i] = A(i,i);
    }
    for (int i(0); i < A.dim1()-1; ++i){
        c[i] = A(i,i+1);
    }

//        valarray< double > dcopy(d);
    TridiagonalSolve(a,b,c,d,xk);

    return true;
}

//-------------------------------------------------------------------
//*******************************************************************
//-------------------------------------------------------------------
     bool Thomas_Tridiagonal(Array2D<complex<double>>& A, 
                       valarray<complex<double> >& d,
                       valarray<complex<double> >& xk) {
//-------------------------------------------------------------------
//   Fills solution into xk. The other matrices are not modified
//   The function returns "false" if the matrix A is not diagonally
//   dominant
//-------------------------------------------------------------------

//      The Matrices all have the right dimensions
//      -------------------------------------------------------------
        if ( ( A.dim1() != A.dim2()  ) || 
             ( A.dim1() != d.size()  ) ||
             ( A.dim1() != xk.size() )    )  {
            cout << "Error: The Matrices don't have the right dimensions!" << endl;
            exit(1);
        }
//      -------------------------------------------------------------

        valarray<complex<double>> a(d.size()), b(d.size()), c(d.size());

        for (int i(0); i < A.dim1()-1; ++i){
           a[i+1] = A(i+1,i);
        }
        for (int i(0); i < A.dim1(); ++i){
           b[i] = A(i,i);
        }
        for (int i(0); i < A.dim1()-1; ++i){
           c[i] = A(i,i+1);
        }
     
//        valarray< double > dcopy(d);
        TridiagonalSolve(a,b,c,d,xk);

        return true;
    }

//-------------------------------------------------------------------
//*******************************************************************
//
//*******************************************************************
//-------------------------------------------------------------------
    complex <double> Det33(/*const valarray<double>& D, */
                          Array2D<complex <double>>& A) {           // Determinant for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( A(1,1)*A(2,2) - A(2,1)*A(1,2) ) -
               A(1,0) * ( A(0,1)*A(2,2) - A(2,1)*A(0,2) ) +
               A(2,0) * ( A(0,1)*A(1,2) - A(1,1)*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Detx33(valarray<complex <double>>& D, 
                           Array2D<complex <double>>& A) {         // Determinant x for a 3*3 system
//-------------------------------------------------------------------
        return D[0] * ( A(1,1)*A(2,2) - A(2,1)*A(1,2) ) -
               D[1] * ( A(0,1)*A(2,2) - A(2,1)*A(0,2) ) +
               D[2] * ( A(0,1)*A(1,2) - A(1,1)*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Dety33(valarray<complex <double>>& D, 
                           Array2D<complex <double>>& A) {         // Determinant y for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( D[1]*A(2,2) - D[2]*A(1,2) ) -
               A(1,0) * ( D[0]*A(2,2) - D[2]*A(0,2) ) +
               A(2,0) * ( D[0]*A(1,2) - D[1]*A(0,2) );
    }
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex <double> Detz33(valarray<complex <double>>& D,
                           Array2D<complex <double>>& A) {         // Determinant z for a 3*3 system
//-------------------------------------------------------------------
        return A(0,0) * ( A(1,1)*D[2] - A(2,1)*D[1] ) -
               A(1,0) * ( A(0,1)*D[2] - A(2,1)*D[0] ) +
               A(2,0) * ( A(0,1)*D[1] - A(1,1)*D[0] );
    }
//-------------------------------------------------------------------

// Convert data structure to float structure
//
// @param[in]  vDouble  The v double
//
// @return     float structure
//
valarray<float> vfloat(const valarray<double>& vDouble) {
    valarray<float> vf(vDouble.size());
    for (size_t i(0); i < vf.size(); ++i) {
        vf[i] = static_cast<float>(vDouble[i]);
    }
    return vf;
}


/**
 * @brief      Convert data structure to float structure
 *
 * @param[in]  vDouble  The v double
 *
 * @return     float structure
 */
vector<float> vfloat(const vector<double> vDouble) {
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i]));
    }
    return vf;
}
vector<float> vfloat(const vector<complex<double> > vDouble) {
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i].real()));
    }
    return vf;
}

vector<float> vfloat_complex(const vector<complex<double> > vDouble){
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i].imag()));
    }
    return vf;
}
//--------------------------------------------------------------

valarray<double> df_4thorder(const valarray<double>& f) {
    valarray<double> df(f.size());

    df[0] = f[1]-f[0];
    df[1] = 1.0/12.0*(f[4]-6.0*f[3]+18.0*f[2]-10.0*f[1]-3.0*f[0]);

    for (size_t i(2); i < df.size()-2; ++i) {
        df[i] = 1.0/12.0*(-f[i+2]+8.0*f[i+1]-8.0*f[i-1]+f[i-2]);
    }

    df[df.size()-2] = 1.0/12.0*(3.0*f[df.size()-1]+10.0*f[df.size()-2]-18.0*f[df.size()-3]+6.0*f[df.size()-4]-f[df.size()-5]);
    df[df.size()-1] = f[df.size()-1]-f[df.size()-2];

    return df;
}

valarray<double> df_4thorder(valarray<double>& f) {
    valarray<double> df(f.size());

    df[0] = f[1]-f[0];
    df[1] = 1.0/12.0*(f[4]-6.0*f[3]+18.0*f[2]-10.0*f[1]-3.0*f[0]);

    for (size_t i(2); i < df.size()-2; ++i) {
        df[i] = 1.0/12.0*(-f[i+2]+8.0*f[i+1]-8.0*f[i-1]+f[i-2]);
    }

    df[df.size()-2] = 1.0/12.0*(3.0*f[df.size()-1]+10.0*f[df.size()-2]-18.0*f[df.size()-3]+6.0*f[df.size()-4]-f[df.size()-5]);
    df[df.size()-1] = f[df.size()-1]-f[df.size()-2];

    return df;
}

valarray<complex<double> > df_4thorder(const valarray<complex<double> >& f) {
    valarray<complex<double> > df(f.size());

    df[0] = f[1]-f[0];
    df[1] = 1.0/12.0*(f[4]-6.0*f[3]+18.0*f[2]-10.0*f[1]-3.0*f[0]);

    for (size_t i(2); i < df.size()-2; ++i) {
        df[i] = 1.0/12.0*(-f[i+2]+8.0*f[i+1]-8.0*f[i-1]+f[i-2]);
    }

    df[df.size()-2] = 1.0/12.0*(3.0*f[df.size()-1]+10.0*f[df.size()-2]-18.0*f[df.size()-3]+6.0*f[df.size()-4]-f[df.size()-5]);
    df[df.size()-1] = f[df.size()-1]-f[df.size()-2];

    return df;
}

    Array2D<complex<double> > df1_4thorder(Array2D<complex<double> >& f) {
    Array2D<complex<double> > df(f.dim1(),f.dim2());

    for (long i2(0); i2<f.dim1();++i2){

        df(0,i2) = f(1,i2)-f(0,i2);
        df(1,i2) = 1.0/12.0*(f(4,i2)-6.0*f(3,i2)+18.0*f(2,i2)-10.0*f(1,i2)-3.0*f(0,i2));

        for (long i1(2); i1<f.dim1()-2;++i1){
            df(i1,i2) = 1.0/12.0*(-f(i1+2,i2)+8.0*f(i1+1,i2)-8.0*f(i1-1,i2)+f(i1-2,i2));
        }

        df(f.dim1()-2,i2) = 1.0/12.0*(3.0*f(f.dim2()-1,i2)+10.0*f(f.dim2()-2,i2)-18.0*f(f.dim2()-3,i2)+6.0*f(f.dim2()-4,i2)-f(f.dim2()-5,i2));
        df(f.dim1()-1,i2) = f(f.dim2()-1,i2)-f(f.dim2()-2,i2);

    }

    return df;
}

    Array2D<complex<double> > df2_4thorder(Array2D<complex<double> >& f) {
    Array2D<complex<double> > df(f.dim1(),f.dim2());

    for (long i1(0); i1<f.dim1();++i1){

        df(i1,0) = f(i1,1)-f(i1,0);
        df(i1,1) = 1.0/12.0*(f(i1,4)-6.0*f(i1,3)+18.0*f(i1,2)-10.0*f(i1,1)-3.0*f(i1,0));

        for (long i2(2); i2<f.dim2()-2;++i2){
            df(i1,i2) = 1.0/12.0*(-f(i1,i2+2)+8.0*f(i1,i2+1)-8.0*f(i1,i2-1)+f(i1,i2-2));
        }

        df(i1,f.dim2()-2) = 1.0/12.0*(3.0*f(i1,f.dim2()-1)+10.0*f(i1,f.dim2()-2)-18.0*f(i1,f.dim2()-3)+6.0*f(i1,f.dim2()-4)-f(i1,f.dim2()-5));
        df(i1,f.dim2()-1) = f(i1,f.dim2()-1)-f(i1,f.dim2()-2);

    }

    return df;
}

