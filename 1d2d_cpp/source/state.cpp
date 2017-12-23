/*! \brief Fields, Distributions, Harmonics, States - Definitions
* \author PICKSC
 * \date   2017
 * \file   state.cpp
 * 
 *   
 *   This cpp file contains the definitions for the data
 *   structures that characterize the state of the 
 *   system:
 *
 *   1.A. SHarmonic1D :
 *       A wrapper class for Array2D<...>. It is used to
 *       describe a single spherical harmonic in 1D.
 *
 *   1.B. Field1D :
 *       This is a class declaration for Ex since it's the only field component in the 1D code
 *
 *   1.C. Distribution1D:
 *       The collection of 1D spherical harmonics. 
 *
 *   1.D. State1D:
 *       The collection of 1D distribution and 1D Fields. 
 */
//--------------------------------------------------------------

// Standard Libraries
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>


// My Libraries
#include "lib-array.h"
#include "lib-algorithms.h"

// Declerations
#include "state.h"

//--------------------------------------------------------------
//  Definition of the 1D spherical harmonic
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
SHarmonic1D::SHarmonic1D(size_t nump, size_t numx) {
    sh = new Array2D<complex<double> >(nump,numx);
}
//  Copy constructor
SHarmonic1D::SHarmonic1D(const SHarmonic1D& other){
    sh = new Array2D<complex<double> >(other.nump(),other.numx());
    *sh = other.array();
}
//  Destructor
SHarmonic1D:: ~SHarmonic1D(){
    delete sh;
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
SHarmonic1D& SHarmonic1D::operator=(const complex<double> & d){
    *sh = d;
    return *this;
}
SHarmonic1D& SHarmonic1D::operator=(const SHarmonic1D& other){
    if (this != &other) {   //self-assignment
        *sh = other.array();
    }
    return *this;
}
//  *= 
SHarmonic1D& SHarmonic1D::operator*=(const complex<double> & d){
    (*sh) *=d;
    return *this;
}
SHarmonic1D& SHarmonic1D::operator*=(const SHarmonic1D& shmulti){
    (*sh) *= shmulti.array();
    return *this;
}
//  +=
SHarmonic1D& SHarmonic1D::operator+=(const complex<double> & d){
    (*sh) +=d;
    return *this;
}
SHarmonic1D& SHarmonic1D::operator+=(const SHarmonic1D& shadd){
    (*sh) += shadd.array();
    return *this;
}
//  -= 
SHarmonic1D& SHarmonic1D::operator-=(const complex<double> & d){
    (*sh) -=d;
    return *this;
}
SHarmonic1D& SHarmonic1D::operator-=(const SHarmonic1D& shmin){
    (*sh) -= shmin.array();
    return *this;
}

//--------------------------------------------------------------
//   Other Algebra
//--------------------------------------------------------------
SHarmonic1D& SHarmonic1D::mpaxis(const valarray <complex<double> > & shmulti){
    (*sh).multid1(shmulti);
    return *this;
}
SHarmonic1D& SHarmonic1D::mxaxis(const valarray <complex<double> > & shmulti){
    (*sh).multid2(shmulti);
    return *this;
}
SHarmonic1D& SHarmonic1D::Re(){
    for (size_t i(0); i < dim(); ++i) {
        (*sh)(i) = (*sh)(i).real();
    }
    return *this;
}
//--------------------------------------------------------------

//  P-difference
SHarmonic1D& SHarmonic1D::Dp(){
    //--------------------------------------------------------//
    //--------------------------------------------------------//
    /// 2nd order
        valarray  <complex<double> >  plast(this->numx());

        for (size_t i(0); i < plast.size(); ++i) {
            plast[i] = (*sh)(nump()-2,i) - (*sh)(nump()-1,i);
        }
        *sh = (*sh).Dd1();
        for (size_t i(0); i < plast.size(); ++i) {
          // TODO                The Dp at the zeroth cell is taken care off
          (*sh)(0,i) = 0.0;     //separately, both for the E-field and the collisions.
            (*sh)(nump()-1,i) = 2.0*plast[i];
        }

    //--------------------------------------------------------//
    //--------------------------------------------------------//
    /// 4th order
    // Array2D<complex<double> > df((*this).nump(),(*this).numx());

    // for (long i2(0); i2<numx();++i2){

    //     df(0,i2) = (*this)(1,i2)-(*this)(0,i2);
    //     df(1,i2) = 1.0/12.0*((*this)(4,i2)-6.0*(*this)(3,i2)+18.0*(*this)(2,i2)-10.0*(*this)(1,i2)-3.0*(*this)(0,i2));

    //     for (long i1(2); i1<nump()-2;++i1){
    //         df(i1,i2) = 1.0/12.0*(-(*this)(i1+2,i2)+8.0*(*this)(i1+1,i2)-8.0*(*this)(i1-1,i2)+(*this)(i1-2,i2));
    //     }

    //     df(nump()-2,i2) = 0.0; //1.0/12.0*(3.0*(*this)(nump()-1,i2)+10.0*(*this)(nump()-2,i2)-18.0*(*this)(nump()-3,i2)+6.0*(*this)(nump()-4,i2)-(*this)(nump()-5,i2));
    //     df(nump()-1,i2) = 0.0; //(*this)(nump()-1,i2)-(*this)(nump()-2,i2);

    // }

    // for (long i2(0); i2<numx();++i2) {
    //     for (long i1(0); i1 < nump(); ++i1) {
    //         (*this)(i1, i2) = -2.0*df(i1, i2);
    //     }
    // }

    return *this;
}
//--------------------------------------------------------------

//  X-difference
SHarmonic1D& SHarmonic1D::Dx(size_t order){

    //--------------------------------------------------------//
    //--------------------------------------------------------//
    /// 2nd order
    // for (size_t ix(0); ix < this->numx(); ++ix) {
    //     for (size_t ip(0); ip < this->nump(); ++ip) {
    //         (*sh)(ip,ix) = ix;
    //     }
    // }
    if (order == 2) *sh = (*sh).Dd2_2nd_order();                          // Worry about boundaries elsewhere
    if (order == 4) *sh = (*sh).Dd2_4th_order();                          // Worry about boundaries elsewhere


    
    
    // for (size_t ix(0); ix < this->numx(); ++ix) {
    //     for (size_t ip(0); ip < this->nump(); ++ip) {
    //         std::cout << "*sh(" << ip << "," << ix << ")" <<  (*sh)(ip,ix) << "\n";
    //     }
    // }

                              
    //--------------------------------------------------------//
    //--------------------------------------------------------//

    return *this;
}
//--------------------------------------------------------------

//  Filter Pcells
SHarmonic1D& SHarmonic1D::Filterp(size_t N){
    *sh = (*sh).Filterd1(N);
    return *this;
}

//  Debug
void SHarmonic1D::checknan(){

    for (size_t i(0); i<numx();++i){
        for (size_t p(0); p<nump();++p){
            if ((isnan((*sh)(p,i).real())) || (isnan((*sh)(p,i).imag())))
            {
                std::cout << "NaN @ (" << p << "," << i << ")\n";
                // MPI_Bcast(&error, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                // if (error != 0) 
                {
                    int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
                    if (rank == 0)
                    {
                            fprintf(stderr, "Error: Program terminated with error code %d\n", 1);
                    }
                    MPI_Finalize();
                    exit(1);
                } 
            }
        }
    }
//    std::cout << "SH OK! \n";
    return;


}
//**************************************************************
//  Definition of the 2D spherical harmonic
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    SHarmonic2D::SHarmonic2D(size_t nump, size_t numx, size_t numy) {
        sh = new Array3D <complex <double> >(nump,numx,numy);
    }
//  Copy constructor
    SHarmonic2D::SHarmonic2D(const SHarmonic2D& other){
        sh = new Array3D < complex <double> >(other.nump(),other.numx(),other.numy());
        *sh = other.array();
    }
//  Destructor
    SHarmonic2D:: ~SHarmonic2D(){
        delete sh; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    SHarmonic2D& SHarmonic2D::operator=(const complex<double>& d){
        *sh = d;
        return *this;
    }
    SHarmonic2D& SHarmonic2D::operator=(const SHarmonic2D& other){
        if (this != &other) {   //self-assignment
            
            *sh = other.array(); 
            
        }
        return *this;
    }
//  *= 
    SHarmonic2D& SHarmonic2D::operator*=(const complex<double>& d){
        (*sh) *=d;
        return *this;
    }
    SHarmonic2D& SHarmonic2D::operator*=(const SHarmonic2D& shmulti){
        (*sh) *= shmulti.array();
        return *this;
    }
//  +=
    SHarmonic2D& SHarmonic2D::operator+=(const complex<double>& d){
        (*sh) +=d;
        return *this;
    }
    SHarmonic2D& SHarmonic2D::operator+=(const SHarmonic2D& shadd){
        (*sh) += shadd.array(); 
        return *this;
    }
//  -= 
    SHarmonic2D& SHarmonic2D::operator-=(const complex<double>& d){
        (*sh) -=d;
        return *this;
    }
    SHarmonic2D& SHarmonic2D::operator-=(const SHarmonic2D& shmin){
        (*sh) -= shmin.array();
        return *this;
    }

//--------------------------------------------------------------
//   Other Algebra
//--------------------------------------------------------------
    SHarmonic2D& SHarmonic2D::mpaxis(const valarray < complex <double> > & shmulti){
        (*sh).multid1(shmulti);
        return *this;
    }
    SHarmonic2D& SHarmonic2D::mxaxis(const valarray < complex <double> > & shmulti){
        (*sh).multid2(shmulti);
        return *this;
    }
    SHarmonic2D& SHarmonic2D::myaxis(const valarray < complex <double> > & shmulti){
        (*sh).multid3(shmulti);
        return *this;
    }

    SHarmonic2D& SHarmonic2D::Re(){
        for (int i(0); i < dim(); ++i) {
            (*sh)(i) = (*sh)(i).real();
        }
        return *this;
    }
//--------------------------------------------------------------
//--------------------------------------------------------------
    SHarmonic2D& SHarmonic2D::Dp(){
//--------------------------------------------------------------
//  P-difference with derivative at #0 set to 0, and f at #np equal to #np-1  
//--------------------------------------------------------------
        Array2D< complex<double> > plast(numx(), numy());
        plast = 0.0;

        for (size_t ix(0); ix < numx(); ++ix) {
            for (size_t iy(0); iy < numy(); ++iy) {
                plast(ix,iy) = (*sh)(nump()-2,ix,iy) - (*sh)(nump()-1,ix,iy); 
            }
        }
        *sh = (*sh).Dd1();
        for (size_t ix(0); ix < numx(); ++ix) {
            for (size_t iy(0); iy < numy(); ++iy) {
           // TODO                The Dp at the zeroth cell is taken care off
                (*sh)(0,ix,iy) = 0.0;   //separately, both for the E-field and the collisions.                     
                (*sh)(nump()-1,ix,iy) = 2.0*plast(ix,iy); 
            }
        }
        


        // GSlice_iter< complex<double> > it1(p0(nump()-2)), it2(p0(nump()-1));  
        // for(int i=0; i< numx()*numy(); ++i){ 
        //     ma(i)  = *it1; ma(i) -= *it2;
        //     ++it1; ++it2;
        // }

        // *sh = (*sh).Dd1();

        // it1 = p0(0), it2 = p0(nump()-1);  
        // for(int i=0; i< numx()*numy(); ++i){
        //     *it1 = 0.0; *it2 = ma(i);
        //      ++it1; ++it2;
        // }
        return *this;
    }
//--------------------------------------------------------------
//  X-difference 
    SHarmonic2D& SHarmonic2D::Dx(size_t order){

    if (order == 2) *sh = (*sh).Dd2_2nd_order();                          // Worry about boundaries elsewhere
    if (order == 4) *sh = (*sh).Dd2_4th_order();                          // Worry about boundaries elsewhere
        // *sh = (*sh).Dd2();                          // Worry about boundaries elsewhere
    return *this;
    }
//  y-difference 
    SHarmonic2D& SHarmonic2D::Dy(size_t order){

        // *sh = (*sh).Dd3();                          // Worry about boundaries elsewhere

        if (order == 2) *sh = (*sh).Dd3_2nd_order();                          // Worry about boundaries elsewhere
        if (order == 4) *sh = (*sh).Dd3_4th_order();                          // Worry about boundaries elsewhere
        return *this;
    }
    SHarmonic2D& SHarmonic2D::mxy_matrix(Array2D< complex<double> >& shmultiM){
        int st(0), nxt(nump());
        for (int im(0); im < shmultiM.dim(); ++im) {
            for (int ip(st); ip < nxt; ++ip)
                (*sh)(ip) *= shmultiM(im);
            st += nump(); nxt += nump();
        }
        return *this;
    }    
//--------------------------------------------------------------

//  Filter Pcells
    SHarmonic2D& SHarmonic2D::Filterp(size_t N){
        *sh = (*sh).Filterd1(N);
        return *this;
    }


    void SHarmonic2D::checknan(){

    for (size_t ix(0); ix<numx();++ix)
    {
        for (size_t iy(0); iy<numx();++iy)
        {
            for (size_t p(0); p<nump();++p)
            {
                if ((isnan((*sh)(p,ix,iy).real())) || (isnan((*sh)(p,ix,iy).imag())))
                {
                    std::cout << "NaN @ (" << p << "," << ix << "," << iy << ")\n";
                    // MPI_Bcast(&error, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                    // if (error != 0) 
                    {
                        int rank;
                        MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
                        if (rank == 0)
                        {
                                fprintf(stderr, "Error: Program terminated with error code %d\n", 1);
                        }
                        MPI_Finalize();
                        exit(1);
                    } 
                }
            }
        }
        return;
    }
}
//**************************************************************

//**************************************************************    

//**************************************************************
//  Definition of the "Field1D" Class
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
Field1D:: Field1D(size_t numx) {
    fi = new  valarray<complex<double> >(numx);
}
//  Copy constructor
Field1D:: Field1D(const Field1D& other){
    fi = new valarray<complex<double> >(other.numx());
    *fi = other.array();
}
//  Destructor
Field1D:: ~Field1D(){
    delete fi;
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
Field1D& Field1D::operator=(const complex<double> & d){
    *fi = d;
    return *this;
}
Field1D& Field1D::operator=(const valarray<complex<double> >& other){
    *fi = other;
    return *this;
}
Field1D& Field1D::operator=(const Field1D& other){
    if (this != &other) {   //self-assignment
        *fi = other.array();
    }
    return *this;
}
//  *= 
Field1D& Field1D::operator*=(const complex<double> & d){
    (*fi) *=d;
    return *this;
}
Field1D& Field1D::operator*=(const valarray<complex<double> >& fimulti){
    (*fi) *= fimulti;
    return *this;
}
Field1D& Field1D::operator*=(const Field1D& fimulti){
    (*fi) *= fimulti.array();
    return *this;
}
//  +=
Field1D& Field1D::operator+=(const complex<double> & d){
    (*fi) +=d;
    return *this;
}
Field1D& Field1D::operator+=(const Field1D& fiadd){
    (*fi) += fiadd.array();
    return *this;
}
//  -= 
Field1D& Field1D::operator-=(const complex<double> & d){
    (*fi) -=d;
    return *this;
}
Field1D& Field1D::operator-=(const Field1D& fimin){
    (*fi) -= fimin.array();
    return *this;
}

//--------------------------------------------------------------
Field1D& Field1D::Re(){
//--------------------------------------------------------------
    for (size_t i(0); i < numx(); ++i) {
        (*fi)[i] = (*fi)[i].real();
    }
    return *this;
}

//--------------------------------------------------------------
Field1D& Field1D::Dx(size_t order){
//--------------------------------------------------------------
    //--------------------------------------------------------//
    //--------------------------------------------------------//
    /// 4th order
    if (order == 4)
    {
        valarray<complex<double> > df(numx());

        df[0] = -2.0*((*fi)[1]-(*fi)[0]);
        df[1] = -1.0*((*fi)[2]-(*fi)[0]);

        for (long i(2); i < numx()-2; ++i) {
            df[i] = -2.0/12.0*(-(*fi)[i+2]+8.0*(*fi)[i+1]-8.0*(*fi)[i-1]+(*fi)[i-2]);
        }

        df[numx()-2] = -1.0*((*fi)[numx()-1]-(*fi)[numx()-3]);
        df[numx()-1] = -2.0*((*fi)[numx()-1]-(*fi)[numx()-2]);

    }
    //--------------------------------------------------------//
    //--------------------------------------------------------//
    // 2nd order
        // df[0] = 2.0*((*fi)[0]-(*fi)[1]);
        // for(long i(0); i < long(numx())-1; ++i) {
        //     df[i] = (*fi)[i-1]-(*fi)[i+1];
        // }
        // df[numx()-1] = 2.0*((*fi)[numx()-2]-(*fi)[numx()-1]);
        // *fi = df;
    else
    {   
        for(long i(0); i< long(numx())-2; ++i) {
            (*fi)[i] -= (*fi)[i+2];
        }
       
        for(long i(numx()-3); i>-1; --i) {
            (*fi)[i+1] = (*fi)[i];
        }
    }
    //--------------------------------------------------------//
    //--------------------------------------------------------//
    return *this;
}

//**************************************************************
//  Definition of the "Field2D" Class
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    Field2D:: Field2D(size_t numx, size_t numy) {
        fi = new  Array2D < complex <double> >(numx,numy);
    }
//  Copy constructor
    Field2D:: Field2D(const Field2D& other){
        fi = new Array2D < complex < double > >(other.numx(),other.numy());
        *fi = other.array();
    }
//  Destructor
    Field2D:: ~Field2D(){
        delete fi; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    Field2D& Field2D::operator=(const complex<double>& d){
        *fi = d;
        return *this;
    }
    Field2D& Field2D::operator=(const Field2D& other){
        if (this != &other) {   //self-assignment
           *fi = other.array(); 
        }
        return *this;
    }
//  *= 
    Field2D& Field2D::operator*=(const complex<double>& d){
        (*fi) *=d;
        return *this;
    }

    Field2D& Field2D::operator*=(const Field2D& fimulti){
        (*fi) *= fimulti.array();
        return *this;
    }
//  +=
    Field2D& Field2D::operator+=(const complex<double>& d){
        (*fi) +=d;
        return *this;
    }
    Field2D& Field2D::operator+=(const Field2D& fiadd){
        (*fi) += fiadd.array(); 
        return *this;
    }
//  -= 
    Field2D& Field2D::operator-=(const complex<double>& d){
        (*fi) -=d;
        return *this;
    }
    Field2D& Field2D::operator-=(const Field2D& fimin){
        (*fi) -= fimin.array();
        return *this;
    }

//--------------------------------------------------------------
    Field2D& Field2D::Dx(size_t order){
//--------------------------------------------------------------
        // *fi = (*fi).Dd1();
        if (order == 2) *fi = (*fi).Dd1_2nd_order();                          // Worry about boundaries elsewhere
        if (order == 4) *fi = (*fi).Dd1_4th_order();                          // Worry about boundaries elsewhere

        return *this;
    }
//--------------------------------------------------------------
    Field2D& Field2D::Dy(size_t order){
//--------------------------------------------------------------
        // *fi = (*fi).Dd2();

        if (order == 2) *fi = (*fi).Dd2_2nd_order();                          // Worry about boundaries elsewhere
        if (order == 4) *fi = (*fi).Dd2_4th_order();                          // Worry about boundaries elsewhere

        return *this;
    } 
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
//**************************************************************
//**************************************************************
//**************************************************************
//  Definition of the Electromagnetic Fields Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor for 1D
//--------------------------------------------------------------
//  Constructor
EMF1D:: EMF1D(size_t nx) {
    fie = new vector<Field1D> (6,Field1D(nx));

}
//  Copy constructor
EMF1D:: EMF1D(const EMF1D& other){
    fie = new vector<Field1D>(6,Field1D(other(1).numx()));
    for (size_t i=0; i < other.dim() ; ++i) (*fie)[i] = other(i);
}
//  Destructor
EMF1D:: ~EMF1D(){
    delete fie;
}
//--------------------------------------------------------------
//  Operators for 1D
//--------------------------------------------------------------

//  Copy assignment operator
EMF1D& EMF1D::operator=(const complex<double>& d){
    for (size_t i=0; i < dim() ; ++i)
        (*fie)[i] = d;
    return *this;
}
EMF1D& EMF1D::operator=(const Field1D& h){
    for (size_t i=0; i < dim() ; ++i){
        if (&((*fie)[i]) != &h) {   //self-assignment
            (*fie)[i] = h;
        }
    }
    return *this;
}
EMF1D& EMF1D::operator=(const EMF1D& other){
    if (this != &other) {   //self-assignment
        for (size_t i=0; i < dim() ; ++i)
            (*fie)[i] = other(i);
    }
    return *this;
}
//  *=
EMF1D& EMF1D::operator*=(const complex<double>& d){
    for (size_t i=0; i < dim() ; ++i)
        (*fie)[i] *= d;
    return *this;
}
EMF1D& EMF1D::operator*=(const EMF1D& other){
    if (this != &other) {   //self-assignment
        for (size_t i=0; i < dim() ; ++i)
            (*fie)[i] *= other(i);
    }
    return *this;
}
//  +=
EMF1D& EMF1D::operator+=(const complex<double>& d){
    for (size_t i=0; i < dim() ; ++i)
        (*fie)[i] += d;
    return *this;
}
EMF1D& EMF1D::operator+=(const EMF1D& other){
    if (this != &other) {   //self-assignment
        for (size_t i=0; i < dim() ; ++i)
            (*fie)[i] += other(i);
    }
    return *this;
}
//  -=
EMF1D& EMF1D::operator-=(const complex<double>& d){
    for (size_t i=0; i < dim() ; ++i)
        (*fie)[i] -= d;
    return *this;
}
EMF1D& EMF1D::operator-=(const EMF1D& other){
    if (this != &other) {   //self-assignment
        for (size_t i=0; i < dim() ; ++i)
            (*fie)[i] -= other(i);
    }
    return *this;
}
//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor for 2D
//--------------------------------------------------------------
//  Constructor
    EMF2D:: EMF2D(size_t nx, size_t ny) {
        fie = new vector<Field2D> (6,Field2D(nx,ny)); 

    }
//  Copy constructor
    EMF2D:: EMF2D(const EMF2D& other){
        fie = new vector<Field2D>(6,Field2D(other(1).numx(),other(1).numy())); 
        for(int i=0; i < other.dim() ; ++i) (*fie)[i] = other(i); 
    }
//  Destructor
    EMF2D:: ~EMF2D(){
        delete fie;
    }
//--------------------------------------------------------------
//  Operators for 2D
//--------------------------------------------------------------

//  Copy assignment operator
    EMF2D& EMF2D::operator=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] = d;
        return *this;
    }
    EMF2D& EMF2D::operator=(const Field2D& h){
        for(int i=0; i < dim() ; ++i){  
            if (&((*fie)[i]) != &h) {   //self-assignment
                (*fie)[i] = h;
            }
        }
        return *this;
    }
    EMF2D& EMF2D::operator=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] = other(i);
        }
        return *this;
    }
//  *=
    EMF2D& EMF2D::operator*=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] *= d;
        return *this;
    }
    EMF2D& EMF2D::operator*=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] *= other(i);
        }
        return *this;
    }
//  +=
    EMF2D& EMF2D::operator+=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] += d;
        return *this;
    }
    EMF2D& EMF2D::operator+=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] += other(i);
        }
        return *this;
    }
//  -=
    EMF2D& EMF2D::operator-=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] -= d;
        return *this;
    }
    EMF2D& EMF2D::operator-=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] -= other(i);
        }
        return *this;
    }    


//**************************************************************
//  Definition of the 1D distribution function
//--------------------------------------------------------------
//  Constructors and Destructor
//--------------------------------------------------------------
//  Constructor
DistFunc1D:: DistFunc1D(size_t l, size_t m, 
                        valarray<double> _dp, 
                        size_t nx, 
                        double q, double _ma)
        : lmax(l), mmax(m), sz(((m+1)*(2*l-m+2))/2),
        dp(_dp), 
        charge(q), ma(_ma), ind(l+1,m+1) {

//      Initialize the array of the harmonics
    if (lmax < 1) {
        cout << "l0 < 1 is not acceptable.\n";
        exit(1);
    }

//      Generate container for the harmonics
    // sz = ((mmax+1)*(2*lmax-mmax+2))/2;
    df = new vector<SHarmonic1D>(sz,SHarmonic1D(_dp.size(),nx));
    
//      Define the index for the triangular array 
    ind = -1;

    if (mmax == 0)
    {
        for (size_t il(0); il < lmax+1 ; ++il)
        {
            ind(il,0) = il;
        }

    }
    else
    {
        for (size_t il=0; il < lmax+1 ; ++il)
        {
            for (size_t im=0; im < ((mmax < il)? mmax:il)+1; ++im)
            {
                ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
                              (il*(mmax+1)-(mmax*(mmax+1))/2 + im));
            }
        }
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Copy constructor
DistFunc1D:: DistFunc1D(const DistFunc1D& other)
        : lmax(other.l0()), mmax(other.m0()),
        sz(((other.m0()+1)*(2*other.l0()-other.m0()+2))/2),
            dp(other.getdp()),
          charge(other.q()), ma(other.mass()), 
          ind(other.l0()+1,other.m0()+1)
          {

//      Generate container for the harmonics
    // sz = ((mmax+1)*(2*lmax-mmax+2))/2;
    df = new vector<SHarmonic1D>(sz,SHarmonic1D(other(0).nump(),other(0).numx()));
    

    for(size_t i(0); i < sz ; ++i){
        (*df).push_back(other(i));
    }

     // Define the index for the triangular array 
    ind = -1;
    if (mmax == 0)
    {
        for (size_t il(0); il < lmax+1 ; ++il)
        {
            ind(il,0) = il;
        }

    }
    else
    {
        for (size_t il(0); il < lmax+1 ; ++il)
        {
            for (size_t im(0); im < ((mmax < il)? mmax:il)+1; ++im)
            {
                ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
                              (il*(mmax+1)-(mmax*(mmax+1))/2 + im));
        
            }
        }
    }

}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Destructor
DistFunc1D:: ~DistFunc1D(){
    delete df;
}
//--------------------------------------------------------------
//  Access 
//--------------------------------------------------------------
//  Pointer to the l-th harmonic
// SHarmonic1D& DistFunc1D::operator()(int i) {
// //         if ((l < 0) || (l> lmax)) return NULL;
//     return (*df)[size_t(i)];
// }

// SHarmonic1D& DistFunc1D::operator()(int i) const {
// //         if ((l < 0) || (l> lmax)) return NULL;
//     return (*df)[size_t(i)];
// }

// SHarmonic1D& DistFunc1D::operator()(size_t l, size_t m) {
// //         if ((l < 0) || (l> lmax)) return NULL;
//     return (*df)[ind(l,m)];
// }

// SHarmonic1D& DistFunc1D::operator()(size_t l, size_t m) const {
// //         if ((l < 0) || (l> lmax)) return NULL;
//     return (*df)[ind(l,m)];
// }

//  Pointer to the "n" neighbor of the l harmonic
//     SHarmonic1D* DistFunc1D::Compus(size_t l, _Compus1D n) const {
//         return _Neighbors[2*l+n];
//     }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
DistFunc1D& DistFunc1D::operator=(const complex<double> & d){
    for(size_t i(0); i < dim() ; ++i){
        (*df)[i] = d;
    }
    return *this;
}
DistFunc1D& DistFunc1D::operator=(const SHarmonic1D& h){
    for(size_t i(0); i < dim() ; ++i){
        if (&((*df)[i]) != &h) {   //self-assignment
            (*df)[i] = h;
        }
    }
    return *this;
}
DistFunc1D& DistFunc1D::operator=(const DistFunc1D& other){
    if (this != &other) {   //self-assignment
        for(size_t i(0); i < dim() ; ++i) {
            (*df)[i] = other(i);
        }
    }
    return *this;
}
//  *=
DistFunc1D& DistFunc1D::operator*=(const complex<double> & d){
    for(size_t i(0); i < dim() ; ++i) {
        (*df)[i] *= d;
    }
    return *this;
}
DistFunc1D& DistFunc1D::operator*=(const DistFunc1D& other){
    if (this != &other) {   //self-assignment
        for(size_t i(0); i < dim() ; ++i) {
            (*df)[i] *= other(i);
        }
    }
    return *this;
}
//  +=
DistFunc1D& DistFunc1D::operator+=(const complex<double> & d){
    for(size_t i(0); i < dim() ; ++i) {
        (*df)[i] += d;
    }
    return *this;
}
DistFunc1D& DistFunc1D::operator+=(const DistFunc1D& other){
    if (this != &other) {   //self-assignment
        for(size_t i(0); i < dim() ; ++i) {
            (*df)[i] += other(i);
        }
    }
    return *this;
}
//  -=
DistFunc1D& DistFunc1D::operator-=(const complex<double> & d){
    for(size_t i(0); i < dim() ; ++i) {
        (*df)[i] -= d;
    }
    return *this;
}
DistFunc1D& DistFunc1D::operator-=(const DistFunc1D& other){
    if (this != &other) {   //self-assignment
        for(size_t i(0); i < dim() ; ++i) {
            (*df)[i] -= other(i);
        }
    }
    return *this;
}

DistFunc1D& DistFunc1D::Filterp(){
    for(size_t i(1); i < dim() ; ++i) {
        // (*df)[i].Filterp(i);
        // (*df)[i].Filterp(filter_ceiling[i]);
    }
    return *this;
}
//--------------------------------------------------------------------------------------------------------------------------
//  Moments for Hydro
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getdensity(){

    valarray<double> out((*df)[0].numx());
    valarray<complex<double> > vr(Algorithms::MakeCAxis(
    static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    for (size_t i(0); i<(*df)[0].numx();++i){
        out[i] = (4.0*M_PI*Algorithms::moment((*df)[0].xVec(i),vr,2.0)).real();
    }

    return out;
}
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getdensity() const {

    valarray<double> out((*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    for (size_t i(0); i<(*df)[0].numx();++i){
        out[i] = (4.0*M_PI*Algorithms::moment((*df)[0].xVec(i),vr,2.0)).real();
    }

    return out;
}
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getcurrent(size_t dir){
    valarray<double> out((*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
        static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
        
    }

    if (dir == 0)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (4.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[1].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 1)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (8.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 2)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (-8.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).imag();
        }
    }

    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;


}
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getcurrent(size_t dir) const{
    valarray<double> out((*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }
    if (dir == 0)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (4.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[1].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 1)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (8.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 2)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (-8.0/3.0*M_PI*charge/ma)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).imag();
        }
    }

    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;


}
//--------------------------------------------------------------------------------------------------------------------------    
Array2D<double> DistFunc1D::getcurrent() const{
    
    Array2D<double> out(3,(*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
        static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    double current_c1(4.0/3.0*M_PI*charge/ma);
    double current_c2(2.0*current_c1);
    double current_c3(-1.0*current_c2);

    for (size_t i(0); i<(*df)[0].numx();++i)
    {
        out(0,i) = (current_c1)*(Algorithms::moment((*df)[1].xVec(i),vr,3.0)).real();

        out(1,i) = (current_c2)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).real();

        out(2,i) = (current_c3)*(Algorithms::moment((*df)[2].xVec(i),vr,3.0)).imag();
    }

    return out;
}
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getrelativisticcurrent(size_t dir){
    
    valarray<double> out((*df)[0].numx());
    valarray<complex<double> > vr(Algorithms::MakeCAxis(
        static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    if (dir == 0)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (4.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 1)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (8.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 2)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (-8.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).imag();
        }
    }
    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;
}
//--------------------------------------------------------------------------------------------------------------------------
valarray<double> DistFunc1D::getrelativisticcurrent(size_t dir) const{
   
    valarray<double> out((*df)[0].numx());
    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    if (dir == 0)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (4.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 1)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (8.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).real();
        }
    }
    else if (dir == 2)
    {
        for (size_t i(0); i<(*df)[0].numx();++i){
            out[i] = (-8.0/3.0*M_PI*charge/ma)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).imag();
        }
    }

    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;


}
//--------------------------------------------------------------------------------------------------------------------------    
Array2D<double> DistFunc1D::getrelativisticcurrent() const{
    Array2D<double> out(3,(*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
        static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    double current_c1(4.0/3.0*M_PI*charge/ma);
    double current_c2(2.0*current_c1);
    double current_c3(-1.0*current_c2);

    for (size_t i(0); i<(*df)[0].numx();++i)
    {
        out(0,i) = (current_c1)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(i),vr,3.0)).real();

        out(1,i) = (current_c2)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).real();

        out(2,i) = (current_c3)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(i),vr,3.0)).imag();
    }


    return out;


}
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//  Moments for Hydro
valarray<double> DistFunc1D::getpressure(){

    valarray<double> out((*df)[0].numx());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    for (size_t i(0); i<(*df)[0].numx();++i){
        out[i] = (4.0/3.0/ma*M_PI*Algorithms::moment((*df)[0].xVec(i),vr,4.0)).real();
    }

    return out;
}

//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//  Debug
void DistFunc1D::checknan(){

    for (size_t indx(0); indx<dim();++indx){
        for (size_t i(0); i<(*df)[indx].numx();++i){
            for (size_t p(0); p<(*df)[indx].nump();++p){
                if (  isnan((*df)[indx](p,i).real()) || isnan((*df)[indx](p,i).imag())   )
                {
                    std::cout << "NaN @ (" << indx << "," << p << "," << i << ")\n";
                    int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
                    if (rank == 0)
                    {
                            fprintf(stderr, "Error: Program terminated with error code %d\n", 1);
                    }
                    MPI_Finalize();
                    exit(1);
                } 
            }
        }
    }
//    std::cout << "OK! \n";
    return;


}
//*********************************************************************************************************************
//  Definition of the 2D distribution function
//--------------------------------------------------------------
//  Constructors and Destructor
//--------------------------------------------------------------
//  Constructor
    DistFunc2D:: DistFunc2D(size_t l, size_t m, 
                            valarray<double> _dp, 
                            size_t nx, size_t ny,
                            double q=1, double _ma=1) 
             : lmax(l), mmax(m), sz(((m+1)*(2*l-m+2))/2),
             dp(_dp), charge(q), ma(_ma), ind(l+1,m+1) {
             
//      Initialize the array of the harmonics
        if (lmax < 1 || mmax < 1) {
            cout << "l0 < 1 or m0 < 1 is not acceptable.\n";
            exit(1);
        }

        // sz = ((mmax+1)*(2*lmax-mmax+2))/2;
        //      Generate container for the harmonics
        df = new vector<SHarmonic2D>(sz,SHarmonic2D(_dp.size(),nx,ny)); 
        
//      Define the index for the triangular array 
        ind = -1;

        if (mmax == 0)
        {
            for(size_t il(0); il < lmax+1 ; ++il){
                ind(il,0) = il;

                
            }

        }
        else
        {
            for(size_t il(0); il < lmax+1 ; ++il){
                for(size_t im(0); im < ((mmax < il)? mmax:il)+1; ++im){
                    ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
                                  (il*(mmax+1)-(mmax*(mmax+1))/2 + im));
                 
                }
            }
        }        

     }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Copy constructor
    DistFunc2D:: DistFunc2D(const DistFunc2D& other)
              : lmax(other.l0()), mmax(other.m0()), 
              sz(((other.m0()+1)*(2*other.l0()-other.m0()+2))/2),
              dp(other.getdp()),
              charge(other.q()), ma(other.mass()), ind(other.l0()+1,other.m0()+1)
    {
        // sz = ((mmax+1)*(2*lmax-mmax+2))/2;

//      Generate container for the harmonics
        df = new vector<SHarmonic2D>(sz,SHarmonic2D(other(0).nump(),other(0).numx(),other(0).numy()));


         for(size_t i(0); i < sz ; ++i){  
            (*df).push_back(other(i));
        }

        //      Define the index for the triangular array 
        ind = -1;
        if (mmax == 0)
        {
            for(size_t il(0); il < lmax+1 ; ++il){
                ind(il,0) = il;
                
            }
        }
        else
        {
            for(size_t il(0); il < lmax+1 ; ++il){
                for(size_t im(0); im < ((mmax < il)? mmax:il)+1; ++im){
                    ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
                                  (il*(mmax+1)-(mmax*(mmax+1))/2 + im));
                    
             
                }
            }
        }
        
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Destructor
    DistFunc2D:: ~DistFunc2D(){
        delete df;
    }
//--------------------------------------------------------------
//  Access THESE WERE CHANGED IS BE (L,M) RATHER THAN (I)... WHY?
//--------------------------------------------------------------
////  Pointer to the (l-th, m-th) harmonic

   // SHarmonic2D& DistFunc2D::operator()(size_t l, size_t m) {   
   //     // if ((l < 0) || (l> lmax) || (m < 0) || (m> mmax)) return NULL;
   //     return (*df)[ind(l,m)];
   //  }
//  Pointer to the l-th harmonic
    // SHarmonic2D& DistFunc2D::operator()(size_t i) {   
    //     return ((*df)[size_t(i)]);
    // }      

   // SHarmonic2D& DistFunc2D::operator()(size_t l, size_t m) const {
   //      // if ((l < 0) || (l> lmax) || (m < 0) || (m> mmax)) return NULL;
   //      return (*df)[ind(l,m)];
   //  }
   //  SHarmonic2D& DistFunc2D::operator()(size_t i) const {
   //      // if ((l < 0) || (l> lmax) || (m < 0) || (m > mmax)) return NULL;
   //      return ((*df)[size_t(i)]);
   //  }

// //  Pointer to the "n" neighbor of the l harmonic
//     SHarmonic1D* DistFunc2D::Compus(size_t l, _Compus1D n) const {
//         return _Neighbors[2*l+n];
//     }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    DistFunc2D& DistFunc2D::operator=(const complex <double>& d){
        for(size_t i(0); i < dim() ; ++i){  
            (*df)[i] = d;
        }
        return *this;
    }
    DistFunc2D& DistFunc2D::operator=(const SHarmonic2D& h){
        for(size_t i(0); i < dim() ; ++i){  
            if (&((*df)[i]) != &h) {   //self-assignment
                (*df)[i] = h;
            }
        }
        return *this;
    }
    DistFunc2D& DistFunc2D::operator=(const DistFunc2D& other){
        if (this != &other) 
        {   //self-assignment
            for(size_t i(0); i < other.dim() ; ++i) {  
        
                (*df)[i] = other(i);

            }
        }

        return *this;
    }
//  *=
    DistFunc2D& DistFunc2D::operator*=(const complex <double>& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] *= d;
        }
        return *this;
    }
    DistFunc2D& DistFunc2D::operator*=(const DistFunc2D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i)  
            {
                    (*df)[i] *= (other(i));
            }
        }
        return *this;
    }
//  +=
    DistFunc2D& DistFunc2D::operator+=(const complex <double>& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] += d;
        }
        return *this;
    }
    DistFunc2D& DistFunc2D::operator+=(const DistFunc2D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) {  
                (*df)[i] += (other(i));
            }
        }
        return *this;
    }
//  -=
    DistFunc2D& DistFunc2D::operator-=(const complex <double>& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] -= d;
        }
        return *this;
    }
    DistFunc2D& DistFunc2D::operator-=(const DistFunc2D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) {  
                (*df)[i] -= (other(i));
            }
        }
        return *this;
    }

    DistFunc2D& DistFunc2D::Filterp(){
        for(size_t i(0); i < dim() ; ++i) { 
            // (*df)[i].Filterp(i);
            // (*df)[i].Filterp(filter_ceiling[i]);
        }
        return *this;
    }

    //--------------------------------------------------------------------------------------------------------------------------
    //  Moments for Hydro
    //--------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------
    Array2D<double> DistFunc2D::getdensity()
    {
        
        Array2D<double> out((*df)[0].numx(),(*df)[0].numy());

        valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

        vr[0] = 0.5*dp[0];
        for (size_t ip(1); ip < dp.size(); ++ip)
        {
            vr[ip] = vr[ip-1] + dp[ip];
        }

          
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (4.0*M_PI*Algorithms::moment((*df)[0].xVec(ix,iy),vr,2.0)).real();
                }
            }
        return out;
    }    

//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
   Array2D<double> DistFunc2D::getcurrent(size_t dir){
        
        Array2D<double> out((*df)[0].numx(),(*df)[0].numy());
        

        valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

        vr[0] = 0.5*dp[0];
        for (size_t ip(1); ip < dp.size(); ++ip)
        {
            vr[ip] = vr[ip-1] + dp[ip];
        }

        
        double current_c1(4.0/3.0*M_PI*charge/ma);
        double current_c2(2.0*current_c1);
        double current_c3(-1.0*current_c2);

        if (dir == 0)
        {    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c1*Algorithms::moment((*df)[1].xVec(ix,iy),vr,3.0)).real();
                }
            }
        }   
        else if (dir == 1)
        {    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c2*Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).real();
                }
            }
        }
        else if (dir == 2)
        {    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c3*Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
                }
            }
        }

        else
        {
            std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
            exit(1);
        }

        return out;


    }
//--------------------------------------------------------------------------------------------------------------------------
    Array2D<double> DistFunc2D::getcurrent(size_t dir) const{
        

        Array2D<double> out((*df)[0].numx(),(*df)[0].numy());
        valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

        vr[0] = 0.5*dp[0];
        for (size_t ip(1); ip < dp.size(); ++ip)
        {
            vr[ip] = vr[ip-1] + dp[ip];
        }

        double current_c1(4.0/3.0*M_PI*charge/ma);
        double current_c2(2.0*current_c1);
        double current_c3(-1.0*current_c2);

        if (dir == 0)
        {    
            
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c1*Algorithms::moment((*df)[1].xVec(ix,iy),vr,3.0)).real();
                }
            }
        }   
        else if (dir == 1)
        {    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c2*Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).real();
                }
            }
        }
        else if (dir == 2)
        {    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (current_c3*Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
                }
            }
        }

        else
        {
            std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
            exit(1);
        }

        return out;


    }
//--------------------------------------------------------------------------------------------------------------------------    
    Array3D<double> DistFunc2D::getcurrent() const{
        
        Array3D<double> out(3,(*df)[0].numx(),(*df)[0].numy());
        valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

        vr[0] = 0.5*dp[0];
        for (size_t ip(1); ip < dp.size(); ++ip)
        {
            vr[ip] = vr[ip-1] + dp[ip];
        }

        
        double current_c1(4.0/3.0*M_PI*charge/ma);
        double current_c2(2.0*current_c1);
        double current_c3(-1.0*current_c2);

            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(0,ix,iy) = (current_c1)*(Algorithms::moment((*df)[1].xVec(ix,iy),vr,3.0)).real();
            
                    out(1,ix,iy) = (current_c2)*(Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).real();
            
                    out(2,ix,iy) = (current_c3)*(Algorithms::moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
                }
            }

        return out;


    }
    //--------------------------------------------------------------------------------------------------------------------------

Array2D<double> DistFunc2D::getrelativisticcurrent(size_t dir){
    
    Array2D<double> out((*df)[0].numx(),(*df)[0].numy());
    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    double current_c1(4.0/3.0*M_PI*charge/ma);
    double current_c2(2.0*current_c1);
    double current_c3(-1.0*current_c2);

    if (dir == 0)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c1)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(ix,iy),vr,3.0)).real();
            }
        }
    }
    else if (dir == 1)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c2)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).real();
            }
        }
    }
    else if (dir == 2)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c3)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
            }
        }
    }

    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;
}
//--------------------------------------------------------------------------------------------------------------------------
Array2D<double> DistFunc2D::getrelativisticcurrent(size_t dir) const{
   
    Array2D<double> out((*df)[0].numx(),(*df)[0].numy());
    valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    double current_c1(4.0/3.0*M_PI*charge/ma);
    double current_c2(2.0*current_c1);
    double current_c3(-1.0*current_c2);

    if (dir == 0)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c1)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(ix,iy),vr,3.0)).real();
            }
        }
    }
    else if (dir == 1)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c2)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).real();
            }
        }
    }
    else if (dir == 2)
    {
        for (size_t ix(0); ix<(*df)[0].numx();++ix)
        {
            for (size_t iy(0); iy<(*df)[0].numy();++iy)
            {
                out(ix,iy) = (current_c3)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
            }
        }
    }

    else
    {
        std::cout << "\n\n ERROR: Wrong current Direction \n\n ";
        exit(1);
    }

    return out;


}
//--------------------------------------------------------------------------------------------------------------------------    
Array3D<double> DistFunc2D::getrelativisticcurrent() const{
    Array3D<double> out(3,(*df)[0].numx(),(*df)[0].numy());

    valarray<complex<double> > vr(Algorithms::MakeCAxis(
        static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

    vr[0] = 0.5*dp[0];
    for (size_t ip(1); ip < dp.size(); ++ip)
    {
        vr[ip] = vr[ip-1] + dp[ip];
    }

    double current_c1(4.0/3.0*M_PI*charge/ma);
    double current_c2(2.0*current_c1);
    double current_c3(-1.0*current_c2);

    for (size_t ix(0); ix<(*df)[0].numx();++ix)
    {
        for (size_t iy(0); iy<(*df)[0].numy();++iy)
        {
            out(0,ix,iy) = (current_c1)*(Algorithms::relativistic_invg_moment((*df)[1].xVec(ix,iy),vr,3.0)).real();

            out(1,ix,iy) = (current_c2)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).real();

            out(2,ix,iy) = (current_c3)*(Algorithms::relativistic_invg_moment((*df)[2].xVec(ix,iy),vr,3.0)).imag();
        }
    }


    return out;


}
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
    //  Moments for Hydro
    Array2D<double> DistFunc2D::getpressure(){
        
        Array2D<double> out((*df)[0].numx(),(*df)[0].numy());
        
        valarray<complex<double> > vr(Algorithms::MakeCAxis(
            static_cast<complex<double> > (0.0),static_cast<complex<double> >(1.0),(*df)[0].nump()));

        vr[0] = 0.5*dp[0];
        for (size_t ip(1); ip < dp.size(); ++ip)
        {
            vr[ip] = vr[ip-1] + dp[ip];
        }
          
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    out(ix,iy) = (4.0/3.0/ma*M_PI*Algorithms::moment((*df)[0].xVec(ix,iy),vr,4.0)).real();
                }
            }
        
        return out;
    }         

//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
    //  Debug

    void DistFunc2D::checknan(){
        for (int indx(0); indx<dim();++indx){    
            for (size_t iy(0); iy<(*df)[0].numy();++iy){
                for (size_t ix(0); ix<(*df)[0].numx();++ix){
                    for (size_t p(0); p<(*df)[indx].nump();++p){
                        if (  isnan((*df)[indx](p,ix,iy).real()) || isnan((*df)[indx](p,ix,iy).imag())   )
                        {   
                            std::cout << "NaN @ (" << indx << "," << p << "," << ix << "," << iy << ")\n";
                            exit(1);
                        }
                    }
                }
            }
        }
        // std::cout << "OK! \n";
        return;


    }    

//--------------------------------------------------------------------------------------------------------------------------
//*********************************************************************************************************************

//--------------------------------------------------------------
//**************************************************************    
//**************************************************************
//  Definition of the "Hydro1D" Class
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
Hydro1D:: Hydro1D(size_t nx, double _mass, double _charge): hydromass(_mass), hydrocharge(_charge) {
    hn  = new  valarray<double >(nx);
    hvx = new  valarray<double >(nx);
    hvy = new  valarray<double >(nx);
    hvz = new  valarray<double >(nx);
    ht  = new  valarray<double >(nx);
    hz  = new  valarray<double >(nx);
}
//  Copy constructor
Hydro1D:: Hydro1D(const Hydro1D& other){
    hn = new valarray<double >(other.numx());
    *hn = other.densityarray();

    hvx = new valarray<double >(other.numx());
    *hvx = other.vxarray();

    hvy = new valarray<double >(other.numx());
    *hvy = other.vyarray();

    hvz = new valarray<double >(other.numx());
    *hvz = other.vzarray();

    ht = new valarray<double >(other.numx());
    *ht = other.temperaturearray();

    hz = new valarray<double >(other.numx());
    *hz = other.Zarray();

}
//  Destructor
Hydro1D:: ~Hydro1D(){
    delete hn;
    delete hvx;
    delete hvy;
    delete hvz;
    delete ht;
    delete hz;
}

//  Copy assignment operator
Hydro1D& Hydro1D::operator=(const double & d){
    *hn = d;
    *hvx = d;
    *hvy = d;
    *hvz = d;
    *ht = d;
    *hz = d;
    return *this;
}
Hydro1D& Hydro1D::operator=(const valarray<double >& other){
    *hn  = other;
    *hvx = other;
    *hvy = other;
    *hvz = other;
    *ht  = other;
    *hz  = other;
    return *this;
}
Hydro1D& Hydro1D::operator=(const Hydro1D& other){

    if (this != &other) {   //self-assignment
        *hn = other.densityarray();
        *hvx = other.vxarray();
        *hvy = other.vyarray();
        *hvz = other.vzarray();
        *ht = other.temperaturearray();
        *hz = other.Zarray();
    }
    return *this;
}

//  Copy assignment operator
Hydro1D& Hydro1D::operator*=(const double & d){
    *hn *= d;
    *hvx *= d;
    *hvy *= d;
    *hvz *= d;
    *ht *= d;
    *hz *= d;
    return *this;
}
Hydro1D& Hydro1D::operator*=(const valarray<double >& other){
    *hn *= other;
    *hvx *= other;
    *hvy *= other;
    *hvz *= other;
    *ht *= other;
    *hz *= other;
    return *this;
}
Hydro1D& Hydro1D::operator*=(const Hydro1D& other){
    if (this != &other) {   //self-assignment
        *hn *= other.densityarray();
        *hvx *= other.vxarray();
        *hvy *= other.vyarray();
        *hvz *= other.vzarray();
        *ht *= other.temperaturearray();
        *hz *= other.Zarray();

    }
    return *this;
}

//  Copy assignment operator
Hydro1D& Hydro1D::operator+=(const double & d){
    *hn += d;
    *hvx += d;
    *hvy += d;
    *hvz += d;
    *ht += d;
    *hz += d;

    return *this;
}
Hydro1D& Hydro1D::operator+=(const valarray<double >& other){
    *hn += other;
    *hvx += other;
    *hvy += other;
    *hvz += other;
    *ht += other;
    *hz += other;
    return *this;
}
Hydro1D& Hydro1D::operator+=(const Hydro1D& other){
    if (this != &other) {   //self-assignment
        *hn += other.densityarray();
        *hvx += other.vxarray();
        *hvy += other.vyarray();
        *hvz += other.vzarray();
        *ht += other.temperaturearray();
        *hz += other.Zarray();

    }
    return *this;
}

//  Copy assignment operator
Hydro1D& Hydro1D::operator-=(const double & d){
    *hn -= d;
    *hvx -= d;
    *hvy -= d;
    *hvz -= d;
    *ht -= d;
    *hz -= d;
    return *this;
}
Hydro1D& Hydro1D::operator-=(const valarray<double >& other){
    *hn -= other;
    *hvx -= other;
    *hvy -= other;
    *hvz -= other;
    *ht -= other;
    *hz -= other;

    return *this;
}
Hydro1D& Hydro1D::operator-=(const Hydro1D& other){
    if (this != &other) {   //self-assignment
        *hn -= other.densityarray();
        *hvx -= other.vxarray();
        *hvy -= other.vyarray();
        *hvz -= other.vzarray();
        *ht -= other.temperaturearray();
        *hz -= other.Zarray();

    }
    return *this;
}
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
Hydro2D:: Hydro2D(size_t nx, size_t ny, double _mass, double _charge): hydromass(_mass), hydrocharge(_charge) {
    hn  = new  Array2D<double >(nx,ny);
    hvx = new  Array2D<double >(nx,ny);
    hvy = new  Array2D<double >(nx,ny);
    hvz = new  Array2D<double >(nx,ny);
    ht  = new  Array2D<double >(nx,ny);
    hz  = new  Array2D<double >(nx,ny);
}
//  Copy constructor
Hydro2D:: Hydro2D(const Hydro2D& other){
    hn = new Array2D<double >(other.numx(), other.numy());
    *hn = other.densityarray();

    hvx = new Array2D<double >(other.numx(), other.numy());
    *hvx = other.vxarray();

    hvy = new Array2D<double >(other.numx(), other.numy());
    *hvy = other.vyarray();

    hvz = new Array2D<double >(other.numx(), other.numy());
    *hvz = other.vzarray();

    ht = new Array2D<double >(other.numx(), other.numy());
    *ht = other.temperaturearray();

    hz = new Array2D<double >(other.numx(), other.numy());
    *hz = other.Zarray();

}
//  Destructor
Hydro2D:: ~Hydro2D(){
    delete hn;
    delete hvx;
    delete hvy;
    delete hvz;
    delete ht;
    delete hz;
}

//  Copy assignment operator
Hydro2D& Hydro2D::operator=(const double & d){
    *hn = d;
    *hvx = d;
    *hvy = d;
    *hvz = d;
    *ht = d;
    *hz = d;
    return *this;
}
Hydro2D& Hydro2D::operator=(const Array2D<double >& other){
    *hn  = other;
    *hvx = other;
    *hvy = other;
    *hvz = other;
    *ht  = other;
    *hz  = other;
    return *this;
}
Hydro2D& Hydro2D::operator=(const Hydro2D& other){

    if (this != &other) {   //self-assignment
        *hn = other.densityarray();
        *hvx = other.vxarray();
        *hvy = other.vyarray();
        *hvz = other.vzarray();
        *ht = other.temperaturearray();
        *hz = other.Zarray();
    }
    return *this;
}

//  Copy assignment operator
Hydro2D& Hydro2D::operator*=(const double & d){
    *hn *= d;
    *hvx *= d;
    *hvy *= d;
    *hvz *= d;
    *ht *= d;
    *hz *= d;
    return *this;
}
Hydro2D& Hydro2D::operator*=(const Array2D<double >& other){
    *hn *= other;
    *hvx *= other;
    *hvy *= other;
    *hvz *= other;
    *ht *= other;
    *hz *= other;
    return *this;
}
Hydro2D& Hydro2D::operator*=(const Hydro2D& other){
    if (this != &other) {   //self-assignment
        *hn *= other.densityarray();
        *hvx *= other.vxarray();
        *hvy *= other.vyarray();
        *hvz *= other.vzarray();
        *ht *= other.temperaturearray();
        *hz *= other.Zarray();

    }
    return *this;
}

//  Copy assignment operator
Hydro2D& Hydro2D::operator+=(const double & d){
    *hn += d;
    *hvx += d;
    *hvy += d;
    *hvz += d;
    *ht += d;
    *hz += d;

    return *this;
}
Hydro2D& Hydro2D::operator+=(const Array2D<double >& other){
    *hn += other;
    *hvx += other;
    *hvy += other;
    *hvz += other;
    *ht += other;
    *hz += other;
    return *this;
}
Hydro2D& Hydro2D::operator+=(const Hydro2D& other){
    if (this != &other) {   //self-assignment
        *hn += other.densityarray();
        *hvx += other.vxarray();
        *hvy += other.vyarray();
        *hvz += other.vzarray();
        *ht += other.temperaturearray();
        *hz += other.Zarray();

    }
    return *this;
}

//  Copy assignment operator
Hydro2D& Hydro2D::operator-=(const double & d){
    *hn -= d;
    *hvx -= d;
    *hvy -= d;
    *hvz -= d;
    *ht -= d;
    *hz -= d;
    return *this;
}
Hydro2D& Hydro2D::operator-=(const Array2D<double >& other){
    *hn -= other;
    *hvx -= other;
    *hvy -= other;
    *hvz -= other;
    *ht -= other;
    *hz -= other;

    return *this;
}
Hydro2D& Hydro2D::operator-=(const Hydro2D& other){
    if (this != &other) {   //self-assignment
        *hn -= other.densityarray();
        *hvx -= other.vxarray();
        *hvy -= other.vyarray();
        *hvz -= other.vzarray();
        *ht -= other.temperaturearray();
        *hz -= other.Zarray();

    }
    return *this;
}


//**************************************************************
//  Particle tracker
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
Particle1D::Particle1D(size_t numparticles, double _mass, double _charge): particlemass(_mass), particlecharge(_charge) {
    
    par_posX   = new valarray<double>(numparticles);
    par_momX   = new valarray<double>(numparticles);
    par_momY   = new valarray<double>(numparticles);
    par_momZ   = new valarray<double>(numparticles);
    par_ishere = new valarray<bool>(numparticles);
    par_goingright = new valarray<int>(numparticles);
}

Particle1D::Particle1D(const Particle1D& other){
    par_posX   = new valarray<double>(other.numpar());
    *par_posX   = other.par_posX_array();

    par_momX   = new valarray<double>(other.numpar());
    *par_momX   = other.par_momX_array();

    par_momY   = new valarray<double>(other.numpar());
    *par_momY   = other.par_momY_array();

    par_momZ   = new valarray<double>(other.numpar());
    *par_momZ   = other.par_momZ_array();

    par_ishere = new valarray<bool>(other.numpar());
    *par_ishere = other.par_ishere_array();

    par_goingright = new valarray<int>(other.numpar());
    *par_goingright = other.par_goingright_array();
}

Particle1D:: ~Particle1D(){
    delete par_posX;
    delete par_momX;
    delete par_momY;
    delete par_momZ;
    delete par_ishere;
    delete par_goingright;
}

//  Copy assignment operator
Particle1D& Particle1D::operator=(const double & d){
    // *par_posX = d;
    // *par_momX = d;
    // *par_momY = d;
    // *par_momZ = d;
    // *par_ishere = 0;
    return *this;
}
Particle1D& Particle1D::operator=(const valarray<double >& other){

    // *par_posX   = other;
    // *par_momX   = other;
    // *par_momY   = other;
    // *par_momZ   = other;
    // *par_ishere = 0;
    return *this;
}
Particle1D& Particle1D::operator=(const Particle1D& other){

    if (this != &other) {   //self-assignment
        // *par_posX = other.par_posX_array();
        // *par_momX = other.par_momX_array();
        // *par_momY = other.par_momY_array();
        // *par_momZ = other.par_momZ_array();
        // *par_ishere = other.par_ishere_array();
    }
    return *this;
}


//--------------------------------------------------------------
//  Copy assignment operator
//--------------------------------------------------------------
Particle1D& Particle1D::operator*=(const double & d){
    // *par_posX *= d;
    // *par_momX *= d;
    // *par_momY *= d;
    // *par_momZ *= d;
    // *par_ishere *= d;
    return *this;
}
Particle1D& Particle1D::operator*=(const valarray<double >& other){
    // *par_posX   *= other;
    // *par_momX   *= other;
    // *par_momY   *= other;
    // *par_momZ   *= other;
    // *par_ishere *= other;
    return *this;
}
Particle1D& Particle1D::operator*=(const Particle1D& other){

    if (this != &other) {   //self-assignment
        // *par_posX *= other.par_posX_array();
        // *par_momX *= other.par_momX_array();
        // *par_momY *= other.par_momY_array();
        // *par_momZ *= other.par_momZ_array();
        // *par_ishere *= other.par_ishere_array();
    }
    return *this;
}

//--------------------------------------------------------------
//  Copy assignment operator
//--------------------------------------------------------------
Particle1D& Particle1D::operator+=(const double & d){
    // *par_posX += d;
    // *par_momX += d;
    // *par_momY += d;
    // *par_momZ += d;
    // *par_ishere += d;
    return *this;
}
Particle1D& Particle1D::operator+=(const valarray<double >& other){
    // *par_posX   += other;
    // *par_momX   += other;
    // *par_momY   += other;
    // *par_momZ   += other;
    // *par_ishere += other;
    return *this;
}
Particle1D& Particle1D::operator+=(const Particle1D& other){

    if (this != &other) {   //self-assignment
        // *par_posX += other.par_posX_array();
        // *par_momX += other.par_momX_array();
        // *par_momY += other.par_momY_array();
        // *par_momZ += other.par_momZ_array();
        // *par_ishere += other.par_ishere_array();
    }
    return *this;
}

//--------------------------------------------------------------
//  Copy assignment operator
//--------------------------------------------------------------
Particle1D& Particle1D::operator-=(const double & d){
    // *par_posX -= d;
    // *par_momX -= d;
    // *par_momY -= d;
    // *par_momZ -= d;
    // *par_ishere -= d;
    return *this;
}
Particle1D& Particle1D::operator-=(const valarray<double >& other){
    // *par_posX   -= other;
    // *par_momX   -= other;
    // *par_momY   -= other;
    // *par_momZ   -= other;
    // *par_ishere -= other;
    return *this;
}
Particle1D& Particle1D::operator-=(const Particle1D& other){

    if (this != &other) {   //self-assignment
        // *par_posX -= other.par_posX_array();
        // *par_momX -= other.par_momX_array();
        // *par_momY -= other.par_momY_array();
        // *par_momZ -= other.par_momZ_array();
        // *par_ishere -= other.par_ishere_array();
    }
    return *this;
}

//**************************************************************
//  State for 1D electrostatic code
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
State1D:: State1D( size_t nx, vector<size_t> l0, vector<size_t> m0,
                   // vector<size_t> np, vector<double> pmax, 
                    vector<valarray<double> > dp,
                   vector<double> q, vector<double> ma,
                   double _hydromass, double _hydrocharge, //double filter_dp, double filter_pmax,
                   size_t numparticles, double particlemass, double particlecharge)
        : ns(l0.size()) {
    // if (ns != np.size()) {   // Non uniform
    if (ns != dp.size()) {  
        cout << "ERROR:Overdetermined number of species\n";
        exit(1);
    }
    sp = new vector<DistFunc1D>;
    for(size_t s(0); s < ns; ++s){
        // (*sp).push_back(DistFunc1D(l0[s],m0[s],np[s],pmax[s],nx,q[s],ma[s]));//,filter_dp,filter_pmax));
        (*sp).push_back(DistFunc1D(l0[s],m0[s],dp[s],nx,q[s],ma[s]));//,filter_dp,filter_pmax));
    }
    flds = new EMF1D(nx);

    hydro = new Hydro1D(nx,_hydromass,_hydrocharge);

    prtcls = new Particle1D(numparticles,particlemass,particlecharge);
}

//  Copy constructor
State1D:: State1D(const State1D& other)
        : ns(other.Species()) {
    sp = new vector<DistFunc1D>;
    for(size_t s(0); s < ns; ++s){
        (*sp).push_back(DistFunc1D(other.DF(s)));
    }
    flds = new EMF1D(other.FLD(0).numx());
    *flds = other.EMF();

    hydro = new Hydro1D(other.FLD(0).numx(),other.HYDRO().mass(),other.HYDRO().charge());
    *hydro = other.HYDRO();

    prtcls = new Particle1D(other.particles().numpar(), other.particles().mass(), other.particles().charge());
    *prtcls = other.particles();
}

//  Destructor
State1D:: ~State1D(){
    delete sp;
    delete flds;
    delete hydro;
    delete prtcls;
}
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
State1D& State1D::operator=(const State1D& other){
    if (this != &other) {   //self-assignment
        ns = other.Species();
        for(size_t s(0); s < ns; ++s){
            (*sp)[s] = other.DF(s);
        }
        *flds = other.EMF();
        *hydro = other.HYDRO();
        *prtcls = other.particles();

    }
    return *this;
}
//  =
State1D& State1D::operator=(const complex<double> & d){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] = d;
    }
    *flds = d;
    *hydro = d.real();
    // *prtcls = d.real();
    return *this;
}
//  *=
State1D& State1D::operator*=(const State1D& other){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] *= other.DF(s);
    }
    *flds *= other.EMF();
    *hydro *= other.HYDRO();
    // *prtcls *= other.particles();

    return *this;
}
//  *=
State1D& State1D::operator*=(const complex<double> & d){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] *= d;
    }
    (*flds) *= d;
    *hydro *= d.real();
    // *prtcls *= d.real();
    return *this;
}
//  +=
State1D& State1D::operator+=(const State1D& other){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] += other.DF(s);
    }
    *flds += other.EMF();
    *hydro += other.HYDRO();
    // *prtcls += other.particles();
    return *this;
}
//  +=
State1D& State1D::operator+=(const complex<double> & d){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] += d;
    }
    (*flds) += d;
    *hydro  += d.real();
    // *prtcls += d.real();
    return *this;
}
//  +=
State1D& State1D::operator-=(const State1D& other){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] -= other.DF(s);
    }
    *flds -= other.EMF();
    *hydro -= other.HYDRO();
    // *prtcls -= other.particles();

    return *this;
}
//  +=
State1D& State1D::operator-=(const complex<double> & d){
    for(size_t s(0); s < ns; ++s){
        (*sp)[s] -= d;
    }
    (*flds) -= d;
    *hydro  -= d.real();
    // *prtcls -= d.real();

    return *this;
}
//   //  Debug
void State1D::checknan(){

    for(size_t s(0); s < ns; ++s){
        (*sp)[s].checknan();
    }
//    std::cout << "Y OK! \n";
    return;
}
//**************************************************************
//**************************************************************
//  Definition of the State2D Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    State2D:: State2D(size_t nx, size_t ny, 
        vector<size_t> l0, vector<size_t> m0, 
        vector<valarray<double> > dp, 
        vector<double> q, vector<double> ma,
                        double _hydromass, double _hydrocharge)
         : ns(l0.size()) {
        if (ns != dp.size()) {
            cout << "ERROR:Overdetermined number of species\n";
            exit(1);
        }
        sp = new vector<DistFunc2D>;

        for(size_t s(0); s < ns; ++s){  

            (*sp).push_back(DistFunc2D(l0[s],m0[s],dp[s],nx,ny,q[s],ma[s]));

        }
        
        flds = new EMF2D(nx,ny);
        hydro = new Hydro2D(nx,ny,_hydromass,_hydrocharge);

        
    }

//  Copy constructor
    State2D:: State2D(const State2D& other)
         : ns(other.Species()) {

        sp = new vector<DistFunc2D>;
        for(size_t s(0); s < ns; ++s){  
            (*sp).push_back(DistFunc2D(other.DF(s)));
        }
        flds = new EMF2D(other.FLD(0).numx(),other.FLD(0).numy());        
        *flds = other.EMF();

        hydro = new Hydro2D(other.FLD(0).numx(),other.FLD(0).numy(),other.HYDRO().mass(),other.HYDRO().charge());
        *hydro = other.HYDRO();

    }
    
//  Destructor
    State2D:: ~State2D(){
        delete sp;
        delete flds;
        delete hydro;

    }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    State2D& State2D::operator=(const State2D& other){
        if (this != &other) {   //self-assignment
            ns = other.Species();


            for(size_t s(0); s < ns; ++s){  
                // std::cout << "\n\n ns = " <<(*this).Species() << " \n";
                (*sp)[s] = other.DF(s);

            }
            *flds = other.EMF();
            *hydro = other.HYDRO();


        }
        return *this;
    }
//  =
    State2D& State2D::operator=(const complex<double>& d){
        for(size_t s(0); s < ns; ++s) (*sp)[s] = d;
        *flds = d;
        *hydro = d.real();
        return *this;
    }
//  *=
    State2D& State2D::operator*=(const State2D& other){
        for(size_t s(0); s < ns; ++s) (*sp)[s] *= other.DF(s);
        *flds *= other.EMF();
        *hydro *= other.HYDRO();

        return *this;
    }
//  *=
    State2D& State2D::operator*=(const complex<double>& d){
        for(size_t s(0); s < ns; ++s) (*sp)[s] *= d;
        *flds *= d;
        *hydro *= d.real();
        return *this;
    }
//  +=
    State2D& State2D::operator+=(const State2D& other){
        for(size_t s(0); s < ns; ++s) (*sp)[s] += other.DF(s);
        *flds += other.EMF();
        *hydro += other.HYDRO();
        return *this;
    }
//  +=
    State2D& State2D::operator+=(const complex<double>& d){
        for(size_t s(0); s < ns; ++s) (*sp)[s] += d;
        *flds  += d;
        *hydro += d.real();
        return *this;
    }
//  -=
    State2D& State2D::operator-=(const State2D& other){
        for(size_t s(0); s < ns; ++s) (*sp)[s] -= other.DF(s);
        *flds -= other.EMF();
        *hydro -= other.HYDRO();
        return *this;
    }
//  -=
    State2D& State2D::operator-=(const complex<double>& d){
        for(size_t s(0); s < ns; ++s) (*sp)[s] -= d;
        *flds -= d;
        *hydro  -= d.real();    
        return *this;
    }

    void State2D::checknan(){
        
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s].checknan();
        }
        // std::cout << "Y OK! \n";
        return;


    }       
//--------------------------------------------------------------
//--------------------------------------------------------------

//--------------------------------------------------------------
//**************************************************************    
