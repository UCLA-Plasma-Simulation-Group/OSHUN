/*! \brief Underlying data structures
 * \author PICKSC
 * \date   September 1, 2016
 * \file   lib-array.h
 *
 *   This header file contains the definitions for multidimensional (2D to 4D) array types.
 *   These types define the following: 
 *   --> Fortan-style access to the underlying valarray (_,_,_)
 *   
 *   --> "=, *=, +=, -=" 
 *   
 *   --> Slices of the array given a range of indices in a specific dimension
 *   
 *   --> Multiplication with a valarray.
 *   
 *   --> Central difference
 *   
 *   --> Subsets of the given array 
 *   
 *   Notice that these array-containers are closer to Fortran-style arrays than to standard 
 *   Matrix types, particularly because they are not intended to facilitate standard Matrix
 *   algebra. No error-checking is included.  
 *
 *   The following classes are defined:
 *      
 *   1.  template<class T> class GSlice_iter:
 *        an iterator for the elements of a regular array type.
 * 
 *   2.  template<class T> class CGSlice_iter :
 *        same as above for constant values
 *
 *   3.a.template<class T> class Array2D :
 *        a 2D container with basic access and algebra
 *
 *   3.b.template<class T> class Array2D_cmplx : 
 *        a 2D container of complex with basic access and algebra
 *
 *   4.a.template<class T> class Array3D :
 *        a 3D container 
 *
 *   4.b.template<class T> class Array3D_cmplx :
 *        a 3D container of complex with basic access and algebra
 *
 *   5.a.template<class T> class Array4D :
 *        a 4D container 
 *
 *   5.b.template<class T> class Array4D_cmplx :
 *        a 4D container of complex with basic access and algebra
 *
 */
//--------------------------------------------------------------
//--------------------------------------------------------------

#ifndef ARRAY_LIBRARY_N_AXIS_H
#define ARRAY_LIBRARY_N_AXIS_H

using namespace std;

//**************************************************************
//  Non-Constant Generalized Slice Iterator
//**************************************************************
//--------------------------------------------------------------
// forward Declarations to allow friend declarations
//--------------------------------------------------------------
template<class T> class GSlice_iter;
template<class T> bool operator==(const GSlice_iter<T>&, const GSlice_iter<T>&);
template<class T> bool operator!=(const GSlice_iter<T>&, const GSlice_iter<T>&);
template<class T> bool operator< (const GSlice_iter<T>&, const GSlice_iter<T>&);
//--------------------------------------------------------------

//--------------------------------------------------------------
template<class T> class GSlice_iter : public iterator<forward_iterator_tag, T> {
//--------------------------------------------------------------
//  Iterator definition Stroustrup p552-p553
//  Generalized slice Stroustrup p677-p678
//--------------------------------------------------------------
private:
    valarray<T>* v;
    gslice gs;
    size_t curr;     // index of current element
    valarray<size_t> gsizes;
    T& ref(size_t i) const;

public:
    // Constructor
    GSlice_iter(valarray<T>* vv,gslice gss);

    // Pointer to the end
    GSlice_iter end() const {
        GSlice_iter t = *this;
        t.curr = gsizes[0];
        return t;
    }

    GSlice_iter& operator++() {++curr; return *this;} //prefix
    GSlice_iter  operator++(int) {GSlice_iter t = *this; ++curr; return t;} //postfix

    T& operator[](size_t i) {return ref(i);}  // C style subscript
    T& operator()(size_t i) {return ref(i);}  // Fortran style subscript
    T& operator*() {return ref(curr);}        // Current element

    friend bool operator==<>(const GSlice_iter& p, const GSlice_iter& q);
    friend bool operator!=<>(const GSlice_iter& p, const GSlice_iter& q);
    friend bool operator< <>(const GSlice_iter& p, const GSlice_iter& q);
};

//--------------------------------------------------------------
//  Referencing this iterator
//--------------------------------------------------------------
template<class T>
T& GSlice_iter<T>::ref(size_t i) const {
    size_t loc(0);
    for (size_t ic = 0; ic < gsizes.size()-1; ++ic){
        loc += (i/gsizes[ic+1]) * gs.stride()[ic];
        i %= gsizes[ic+1];
    }
    loc += i * gs.stride()[gsizes.size()-1];
    return (*v)[gs.start()+loc];
}

//--------------------------------------------------------------
//  Generalized slice iterator constructor
//--------------------------------------------------------------
template<class T>
GSlice_iter<T>::GSlice_iter(valarray<T>* vv,gslice gss) :
        v(vv), gs(gss), curr(0), gsizes(gs.size()){
    for (size_t ic(1); ic < gsizes.size(); ++ic) {
        gsizes[gss.size().size()-ic-1] *= gsizes[gss.size().size()-ic];
    }
}

//--------------------------------------------------------------
//  Comparison operators for this iterator:
//--------------------------------------------------------------
template<class T>
bool operator==(const GSlice_iter<T>& p, const GSlice_iter<T>& q){
    size_t count(0);
    if (p.curr==q.curr && p.gs.start() == q.gs.start())
        while (  (p.gs.stride()[count] == p.gs.stride()[count])
                 && (p.gs.size()[count]   == p.gs.size()[count]  ))
            if (count++ == q.gs.size().size()-1) return true;
    return false;
}

template<class T>
bool operator!=(const GSlice_iter<T>& p, const GSlice_iter<T>& q){
    return !(p==q);
}

template<class T>
bool operator< (const GSlice_iter<T>& p, const GSlice_iter<T>& q){
    size_t count(0);
    if (p.curr<q.curr && p.gs.start() == q.gs.start())
        while (  (p.gs.stride()[count] == p.gs.stride()[count])
                 && (p.gs.size()[count]   == p.gs.size()[count]  ))
            if (count++ == q.gs.size().size()-1) return true;
    return false;
}
//--------------------------------------------------------------

//**************************************************************
// Constant Generalized Slice Iterator
//**************************************************************
//--------------------------------------------------------------
// forward Declarations to allow friend declarations
//--------------------------------------------------------------
template<class T> class CGSlice_iter;
template<class T> bool operator==(const CGSlice_iter<T>&, const CGSlice_iter<T>&);
template<class T> bool operator!=(const CGSlice_iter<T>&, const CGSlice_iter<T>&);
template<class T> bool operator< (const CGSlice_iter<T>&, const CGSlice_iter<T>&);
//--------------------------------------------------------------

//--------------------------------------------------------------
template<class T> class CGSlice_iter : public iterator<forward_iterator_tag, T,const T*, const T&> {
//--------------------------------------------------------------
//  Iterator definition: Stroustrup p552-p553 
//  Generalized Slice: Stroustrup p677-p678yy 
//--------------------------------------------------------------
private:
    valarray<T>* v;
    gslice gs;
    size_t curr;     // index of current element
    valarray<size_t> gsizes;
    const T& ref(size_t i) const;

public:
    // Constructor
    CGSlice_iter(valarray<T>* vv,gslice gss);

    // Pointer to the end
    CGSlice_iter end() const {
        CGSlice_iter t = *this;
        t.curr = gsizes[0];
        return t;
    }

    CGSlice_iter& operator++() {++curr; return *this;} //prefix
    CGSlice_iter  operator++(int) {CGSlice_iter t = *this; ++curr; return t;} //postfix

    const T& operator[](size_t i) const {return ref(i);}  // C style subscript
    const T& operator()(size_t i) const {return ref(i);}  // Fortran style subscript
    const T& operator*() const {return ref(curr);}        // Current element

    friend bool operator==<>(const CGSlice_iter& p, const CGSlice_iter& q);
    friend bool operator!=<>(const CGSlice_iter& p, const CGSlice_iter& q);
    friend bool operator< <>(const CGSlice_iter& p, const CGSlice_iter& q);
};

//--------------------------------------------------------------
//  Referencing this iterator
//--------------------------------------------------------------
template<class T>
const T& CGSlice_iter<T>::ref(size_t i) const {
    size_t loc(0);
    for (size_t ic = 0; ic < gsizes.size()-1; ++ic){
        loc += (i/gsizes[ic+1]) * gs.stride()[ic];
        i %= gsizes[ic+1];
    }
    loc += i * gs.stride()[gsizes.size()-1];
    return (*v)[gs.start()+loc];
}

//--------------------------------------------------------------
//  Generalized slice iterator constructor
//--------------------------------------------------------------
template<class T>
CGSlice_iter<T>::CGSlice_iter(valarray<T>* vv,gslice gss) :
        v(vv), gs(gss), curr(0), gsizes(gs.size()){
    for (int ic=1; ic < gsizes.size(); ++ic) {
        gsizes[gss.size().size()-ic-1] *= gsizes[gss.size().size()-ic];
    }
}

//--------------------------------------------------------------
//  Comparison operators for this iterator:
//--------------------------------------------------------------
template<class T>
bool operator==(const CGSlice_iter<T>& p, const CGSlice_iter<T>& q){
    size_t count(0);
    if (p.curr==q.curr && p.gs.start() == q.gs.start())
        while (  (p.gs.stride()[count] == p.gs.stride()[count])
                 && (p.gs.size()[count]   == p.gs.size()[count]  ))
            if (count++ == q.gs.size().size()-1) return true;
    return false;
}

template<class T>
bool operator!=(const CGSlice_iter<T>& p, const CGSlice_iter<T>& q){
    return !(p==q);
}

template<class T>
bool operator< (const CGSlice_iter<T>& p, const CGSlice_iter<T>& q){
    size_t count(0);
    if (p.curr<q.curr && p.gs.start() == q.gs.start())
        while (  (p.gs.stride()[count] == p.gs.stride()[count])
                 && (p.gs.size()[count]   == p.gs.size()[count]  ))
            if (count++ == q.gs.size().size()-1) return true;
    return false;
}
//--------------------------------------------------------------
//**************************************************************

/**************************************************************
 *   2D Array Class
 *
 *   Using Stroustrup's matrices p672-p673
 *   This container is similar to Fortran-style arrays.  
 *   Custom operations, no standard Matrix algebra.
 *   No error-checking.
 *   Compile with -ftree-vectorize    
 */
template<class T> class Array2D {
//--------------------------------------------------------------
//  2D Array decleration
//--------------------------------------------------------------
private:
    valarray<T> *v;
    size_t  d1, d2;            // for VFP: p, x

public:
//      Constructors/Destructors
    Array2D(size_t x, size_t y);
    Array2D(const Array2D& other);
    ~Array2D();

//      Basic Info
    size_t dim()  const {return d1*d2;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    valarray<T>& array() const {return *v;}

//      Access
    T& operator()(size_t i, size_t j); // Fortran-style
    T  operator()(size_t i, size_t j) const;
    T& operator()(size_t i);           // 1D-style
    T  operator()(size_t i) const;
    vector<T>  d2c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);    // b: beginning, e: end (non-contiguous)
    GSlice_iter<T>  d2c(size_t b, size_t e);    // b: beginning, e: end (contiguous)
    GSlice_iter<T> SubArray2D(size_t st, size_t nx, size_t ny); // st: starting cell

//      Constant slice iterators
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray2D(size_t st, size_t nx, size_t ny) const;

//      Operators
    Array2D& operator=(const T& d);
    Array2D& operator=(const Array2D& other);
    Array2D& operator*=(const T& d);
    Array2D& operator*=(const Array2D& vmulti);
    Array2D& operator+=(const T& d);
    Array2D& operator+=(const Array2D& vadd);
    Array2D& operator-=(const T& d);
    Array2D& operator-=(const Array2D& vmin);

//      Array * Vector
    Array2D& multid1(const valarray<T>& vmulti); // M*valarray(d1)
    Array2D& multid2(const valarray<T>& vmulti); // M*valarray(d2)

//      Central difference
    Array2D& Dd1(); // in the direction d1 (requires dim1() > 2)
    Array2D& Dd2(); // in the direction d2 (requires dim2() > 2)

//      Filter first N-cells in d1 direction 
    Array2D& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
template<class T> Array2D<T>:: Array2D(size_t x,size_t y) : d1(x), d2(y) {
    v = new valarray<T>(d1*d2);
}
//  Copy constructor
template<class T> Array2D<T>:: Array2D(const Array2D& other){
    d1  = other.dim1();
    d2  = other.dim2();
    v = new valarray<T>(d1*d2);
    (*v) = other.array();
}
//  Destructor
template<class T> Array2D<T>:: ~Array2D(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  Access Fortan-style
template<class T> inline T& Array2D<T>:: operator()(size_t i, size_t j){
    return (*v)[i+j*d1];
}
//  Constant access Fortan-style
template<class T> inline T Array2D<T>:: operator()(size_t i, size_t j) const {
    return (*v)[i+j*d1];
}
//  1D style access
template<class T> inline T& Array2D<T>:: operator() (size_t i){
    return (*v)[i];
}
//  1D style const access
template<class T> inline T  Array2D<T>:: operator() (size_t i) const {
    return (*v)[i];
}
//  Consider Array2D a list of d2-vectors and multiply with vmulti
template<class T> vector<T> Array2D<T>::d2c(size_t j){
    vector<T> d1vec(d1);
    for (size_t i(0); i < d1; ++i ){
        d1vec[i] = (*v)[i+j*d1];
    }
    return d1vec;
}
//--------------------------------------------------------------
//  Slicers
//--------------------------------------------------------------
//  slices for given d1 value range 
template<class T>
inline GSlice_iter<T> Array2D<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(2), str(2);
    str[1] = d1; str[0] = 1;
    sz[1]  = d2;  sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b,sz,str));
}
//  const slices for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array2D<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(2), str(2);
    str[1] = d1; str[0] = 1;
    sz[1]  = d2;  sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b,sz,str));
}
//  slices for given d2 value range 
template<class T>
inline GSlice_iter<T> Array2D<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(1), str(1);     // sz --> size, str --> stride
    str[0] = 1; sz[0]  = (e-b+1)*d1;
    return GSlice_iter<T>(v,gslice(b*d1,sz,str));
}
//  const slices for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array2D<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(1), str(1);     // sz --> size, str --> stride
    str[0] = 1; sz[0]  = (e-b+1)*d1;
    return CGSlice_iter<T>(v,gslice(b*d1,sz,str));
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  scan Subarray
template<class T>
inline GSlice_iter<T> Array2D<T>::SubArray2D(size_t st, size_t nx, size_t ny) {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;  str[0] = dim1();
    sz[1]  = nx;  sz[0] = ny;
    return GSlice_iter<T>(v,gslice(st,sz,str));
}
//  scan const Subarray
template<class T>
inline CGSlice_iter<T> Array2D<T>::SubArray2D(size_t st, size_t nx, size_t ny) const{
    valarray<size_t> sz(2), str(2);
    str[1] = 1;  str[0] = dim1();
    sz[1]  = nx;  sz[0] = ny;
    return CGSlice_iter<T>(v,gslice(st,sz,str));
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
template<class T> Array2D<T>& Array2D<T>::operator=(const T& d){
    (*v) = d;
    return *this;
}
template<class T> Array2D<T>& Array2D<T>::operator=(const Array2D& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}

//  *= 
template<class T> Array2D<T>& Array2D<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
template<class T> Array2D<T>& Array2D<T>::operator*=(const Array2D& vmulti){
    (*v) *= vmulti.array();
    return *this;
}

//  +=
template<class T> Array2D<T>& Array2D<T>::operator+=(const T& d){
    (*v) +=d;
    return *this;
}
template<class T> Array2D<T>& Array2D<T>::operator+=(const Array2D& vadd){
    (*v) += vadd.array();
    return *this;
}

//  -= 
template<class T> Array2D<T>& Array2D<T>::operator-=(const T& d){
    (*v) -=d;
    return *this;
}
template<class T> Array2D<T>& Array2D<T>::operator-=(const Array2D& vmin){
    (*v) -= vmin.array();
    return *this;
}

//--------------------------------------------------------------
//  Array2D * Vector
//--------------------------------------------------------------
//  Consider Array2D a list of d1-vectors and multiply with vmulti
template<class T> Array2D<T>& Array2D<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< d1*d2; j+=d1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[i+j] *= vmulti[i];
        }
    }
    return *this;
}

//  Consider Array2D a list of d2-vectors and multiply with vmulti
template<class T> Array2D<T>& Array2D<T>::multid2(const valarray<T>& vmulti){
    for (size_t j(0); j< d2; ++j ){
        for (size_t i(j*d1); i< (j+1)*d1; ++i ){
            (*v)[i] *= vmulti[j];
        }
    }
    return *this;
}
//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
// Example: 0 4  8 12 16       -2 -2 -2 -2 -2 
//          1 5  9 13 17  -->  -2 -2 -2 -2 -2 
//          2 6 10 14 18       -2 -2 -2 -2 -2  
//          3 7 11 15 19       -2 -2 -2 -2 19  
// Requires at least 3 elements in d1 
template<class T> Array2D<T>& Array2D<T>::Dd1(){
    for(long i(0); i< long(d1*d2)-2; ++i) {
        (*v)[i] -= (*v)[i+2];
    }
    for(long i(d1*d2-3); i>-1; --i) {
        (*v)[i+1] = (*v)[i];
    }


   //  Array2D<T> temp(*this);

   //  for (long i2(0); i2<long(d2);++i2){


   //      /// Second Order
   //      temp(0,i2) = 2.0*((*this)(0,i2)-(*this)(1,i2));


        
   //      // temp(i1,1) = 1.0/12.0*((*this)(i1,4)-6.0*(*this)(i1,3)+18.0*(*this)(i1,2)-10.0*(*this)(i1,1)-3.0*(*this)(i1,0));
        
   //      for (long i1(1); i1<long(d1)-1;++i1){
   //          // temp(i1,i2) = 1.0/12.0*(-(*this)(i1,i2+2)+8.0*(*this)(i1,i2+1)-8.0*(*this)(i1,i2-1)+(*this)(i1,i2-2));
            
   //          /// Second Order
   //          temp(i1,i2) = (*this)(i1-1,i2) - (*this)(i1+1,i2);
   //      }




   //      // temp(i1,long(d2)-2) = 1.0/12.0*(3.0*(*this)(i1,long(d2)-1)+10.0*(*this)(i1,long(d2)-2)-18.0*(*this)(i1,long(d2)-3)+6.0*(*this)(i1,long(d2)-4)-(*this)(i1,long(d2)-5));
   //      // temp(i1,long(d2)-1) = (*this)(i1,long(d2)-1)-(*this)(i1,long(d2)-2);

   //      /// Second Order
   //      temp(long(d1)-1,i2) = 2.0*((*this)(long(d1)-2,i2)-(*this)(long(d1)-1,i2));
   // }

   // *this = temp;
    return *this;

}

// (minus) Central difference for elements with distance d1 (row-wise)
// Example: 0 4  8 12 16       -8 -8 -8 -8 16 
//          1 5  9 13 17  -->  -8 -8 -8 -8 17 
//          2 6 10 14 18       -8 -8 -8 -8 18
//          3 7 11 15 19       -8 -8 -8 -8 19
// Requires at least 3 elements in d2 
template<class T> Array2D<T>& Array2D<T>::Dd2(){
   //  Array2D<T> temp(*this);

   //  for (long i1(0); i1<long(d1);++i1){
   //      /// Second Order
   //      temp(i1,0) = 2.0*((*this)(i1,0)-(*this)(i1,1));
   //      // std::cout << "this(" << i1 << ")" << (*this)(i1,0) << "\n";
   //      // std::cout << "temp(" << i1 << ")" << temp(i1,0) << "\n";

   //      // temp(i1,1) = 1.0/12.0*((*this)(i1,4)-6.0*(*this)(i1,3)+18.0*(*this)(i1,2)-10.0*(*this)(i1,1)-3.0*(*this)(i1,0));
   //      for (long i2(1); i2<long(d2)-1;++i2){
   //          // temp(i1,i2) = 1.0/12.0*(-(*this)(i1,i2+2)+8.0*(*this)(i1,i2+1)-8.0*(*this)(i1,i2-1)+(*this)(i1,i2-2));
   //          /// Second Order
   //          temp(i1,i2) = (*this)(i1,i2-1) - (*this)(i1,i2+1);
   //          // std::cout << "this(" << i1 << "," << i2 << ")" << (*this)(i1,i2) << "\n";
   //      }
   //      // temp(i1,long(d2)-2) = 1.0/12.0*(3.0*(*this)(i1,long(d2)-1)+10.0*(*this)(i1,long(d2)-2)-18.0*(*this)(i1,long(d2)-3)+6.0*(*this)(i1,long(d2)-4)-(*this)(i1,long(d2)-5));
   //      // temp(i1,long(d2)-1) = (*this)(i1,long(d2)-1)-(*this)(i1,long(d2)-2);
   //      /// Second Order
   //      temp(i1,long(d2)-1) = 2.0*((*this)(i1,long(d2)-2)-(*this)(i1,long(d2)-1));
   //      // std::cout << "this(" << i1 << "," << long(d2)-1 << ")" << (*this)(i1,long(d2)-1) << "\n";
   // }

   // *this = temp;
   //  return *this;
   

    ////////////////// Has boundary errors unless boundary cells increased for each RK level
    /// 
   // std::cout << "\n d1*d2-twod1 = " << long(d1*d2)-2*d1 << "\n";

   long twod1(2*d1);
    for(long i(0); i< long(d1*d2)-twod1; ++i) {
        // std::cout << "v[" << i << "] = " << (*v)[i] 
        // <<  ", v[" << i+twod1 << "] = " << (*v)[i+twod1] << "\n";

        (*v)[i] -= (*v)[i+twod1];
    }
    
    for(long i(d1*d2-twod1-1); i>-1; --i) {
        (*v)[i+d1] = (*v)[i];
    }
    return *this;
   ////////////////// ////////////////// //////////////////
}
//--------------------------------------------------------------

//  Remove data for N cells from dimension 1  
template<class T> Array2D<T>& Array2D<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2; ++j ){
        for (size_t i(j*d1); i< j*d1+N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}
//--------------------------------------------------------------
//**************************************************************

/**************************************************************
 *   2D Array Complex Class
 *
 *   Using Stroustrup's matrices p672-p673
 *   This container is similar to Fortran-style arrays.  
 *   Custom operations, no standard Matrix algebra.
 *   No error-checking.
 *   Compile with -ftree-vectorize    
 */
template<class T> class Array2D_cmplx {
//--------------------------------------------------------------
//  2D Array decleration
//--------------------------------------------------------------
private:
    valarray<T> *v;
    size_t  d1, d2;            // for VFP: p, x
    size_t td1;

public:
//      Constructors/Destructors
    Array2D_cmplx(size_t x, size_t y);
    Array2D_cmplx(const Array2D_cmplx& other);
    ~Array2D_cmplx();

//      Basic Info
    size_t dim()  const {return d1*d2;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    valarray<T>& array() const {return *v;}

//      Access
    T& operator[](size_t i);           // 1D-style
    T  operator[](size_t i) const;

    T& real(size_t i,size_t j);           // 1D-style
    T  real(size_t i, size_t j) const;
    T& imag(size_t i,size_t j);           // 1D-style
    T  imag(size_t i, size_t j) const;

//TODO   vector<T>  d2c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);    // b: beginning, e: end (non-contiguous)
    GSlice_iter<T>  d2c(size_t b, size_t e);    // b: beginning, e: end (contiguous)
    GSlice_iter<T> SubArray2D(size_t st, size_t nx, size_t ny); // st: starting cell

//      Constant slice iterators
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray2D(size_t st, size_t nx, size_t ny) const;

//      Operators
    Array2D_cmplx& operator=(const T& d);
    Array2D_cmplx& operator=(const complex<T>& c);
    Array2D_cmplx& operator=(const Array2D<T>& other);
    Array2D_cmplx& operator=(const Array2D_cmplx& other);
    Array2D_cmplx& operator*=(const T& d);
    Array2D_cmplx& operator*=(const complex<T>& c);
    Array2D_cmplx& operator*=(const Array2D<T>& vmulti);
    Array2D_cmplx& operator*=(const Array2D_cmplx& vmulti);
    Array2D_cmplx& operator+=(const T& d);
    Array2D_cmplx& operator+=(const complex<T>& c);
    Array2D_cmplx& operator+=(const Array2D<T>& vadd);
    Array2D_cmplx& operator+=(const Array2D_cmplx& vadd);
    Array2D_cmplx& operator-=(const T& d);
    Array2D_cmplx& operator-=(const complex<T>& c);
    Array2D_cmplx& operator-=(const Array2D<T>& vmin);
    Array2D_cmplx& operator-=(const Array2D_cmplx& vmin);

//      Array * Vector
    Array2D_cmplx& multid1(const valarray<T>& vmulti); 		// M*valarray(d1)
    Array2D_cmplx& multid1(const valarray< complex<T> >& vmulti);   // M*valarray(d1)
    Array2D_cmplx& multid2(const valarray<T>& vmulti); 		// M*valarray(d2)
    Array2D_cmplx& multid2(const valarray< complex<T> >& vmulti);   // M*valarray(d2)

//      Central difference
    Array2D_cmplx& Dd1(); // in the direction d1 (requires dim1() > 2)
    Array2D_cmplx& Dd2(); // in the direction d2 (requires dim2() > 2)

//      Filter first N-cells in d1 direction 
    Array2D_cmplx& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
template<class T> Array2D_cmplx<T>:: Array2D_cmplx(size_t x,size_t y) : d1(x), d2(y) {
    td1 = 2*d1;
    v = new valarray<T>(2*d1*d2);
}
//  Copy constructor
template<class T> Array2D_cmplx<T>:: Array2D_cmplx(const Array2D_cmplx& other){
    d1  = other.dim1();
    d2  = other.dim2();
    td1 = 2*d1;
    v = new valarray<T>(2*d1*d2);
    (*v) = other.array();
}
//  Destructor
template<class T> Array2D_cmplx<T>:: ~Array2D_cmplx(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  1D style access
template<class T> inline T& Array2D_cmplx<T>::operator[] (size_t i){
    return (*v)[i];
}
//  1D style const access
template<class T> inline T  Array2D_cmplx<T>::operator[] (size_t i) const {
    return (*v)[i];
}
//  Access real part by reference
template<class T> inline T& Array2D_cmplx<T>::real(size_t i,size_t j){
    return (*v)[2*i+j*td1];
}
//  Const access of real part 
template<class T> inline T Array2D_cmplx<T>::real(size_t i,size_t j) const{
    return (*v)[2*i+j*td1];
}
//  access of imag part 
template<class T> inline T& Array2D_cmplx<T>::imag(size_t i,size_t j){
    return (*v)[2*i+1+j*td1];
}
//  Const access of imag part 
template<class T> inline T Array2D_cmplx<T>::imag(size_t i,size_t j) const{
    return (*v)[2*i+1+j*td1];
}
//--------------------------------------------------------------
//  Slicers
//--------------------------------------------------------------
//  slices for given d1 value range 
template<class T>
inline GSlice_iter<T> Array2D_cmplx<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = 1, str[1] = td1; str[0] = 2;
    sz[2] = 2,  sz[1]  = d2;  sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  const slices for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array2D_cmplx<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1, str[1] = td1; str[0] = 2;
    sz[2] = 2,  sz[1]  = d2;  sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  slices for given d2 value range 
template<class T>
inline GSlice_iter<T> Array2D_cmplx<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(1), str(1);     // sz --> size, str --> stride
    str[0] = 1; sz[0]  = (e-b+1)*td1;
    return GSlice_iter<T>(v,gslice(b*td1,sz,str));
}
//  const slices for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array2D_cmplx<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(1), str(1);     // sz --> size, str --> stride
    str[0] = 1; sz[0]  = (e-b+1)*td1;
    return CGSlice_iter<T>(v,gslice(b*td1,sz,str));
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  scan Subarray
template<class T>
inline GSlice_iter<T> Array2D_cmplx<T>::SubArray2D(size_t st, size_t nx, size_t ny) {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;    str[0] = td1;
    sz[1]  = 2*nx;  sz[0] = ny;
    return GSlice_iter<T>(v,gslice(2*st,sz,str));
}
//  scan const Subarray
template<class T>
inline CGSlice_iter<T> Array2D_cmplx<T>::SubArray2D(size_t st, size_t nx, size_t ny) const{
    valarray<size_t> sz(2), str(2);
    str[1] = 1;    str[0] = td1;
    sz[1]  = 2*nx;  sz[0] = ny;
    return CGSlice_iter<T>(v,gslice(2*st,sz,str));
}


//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator=(const T& d){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i] = d;
        (*v)[i+1] = 0.0;
    }
    return *this;
}
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator=(const complex<T>& c){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i]   =c.real();
        (*v)[i+1] =c.imag();
    }
    return *this;
}
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator=(const Array2D<T>& other){
    for (size_t i(0); i< d1*d2; ++i) {
        (*v)[2*i]   = other.array()[i];
        (*v)[2*i+1] = 0.0;
    }
    return *this;
}
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator=(const Array2D_cmplx& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  *= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
//  *= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator*=(const complex<T>& c){
    for (size_t i(0); i< td1*d2; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   c=a+ib
        T vi((*v)[i]*c.imag());
        (*v)[i]    *= c.real();                     // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*c.imag();           // v[i]   = Aa - Bb
        (*v)[i+1]  *= c.real();                     // v[i+1] = Ba
        (*v)[i+1]  += vi;                           // v[i+1] = Ba + Ab
    }
    return *this;
}
//  *= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator*=(const Array2D<T>& vmulti){
    for (size_t i(0); i< d1*d2; ++i) {
        (*v)[2*i]   *= vmulti.array()[i];
        (*v)[2*i+1] *= vmulti.array()[i];
    }
    return *this;
}
//  *= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator*=(const Array2D_cmplx& vmulti){
    for (size_t i(0); i< td1*d2; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   vmulti[i] + ivmulti[i+1] = a+ib
        T vi((*v)[i] * vmulti.array()[i+1]);
        (*v)[i]    *=  vmulti.array()[i];              // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*vmulti.array()[i+1];   // v[i]   = Aa - Bb
        (*v)[i+1]  *= vmulti.array()[i];               // v[i+1] = Ba
        (*v)[i+1]  += vi;                              // v[i+1] = Ba + Ab
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  += 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator+=(const T& d){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i] +=d;
    }
    return *this;
}
//  += 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator+=(const complex<T>& c){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i]   +=c.real();
        (*v)[i+1] +=c.imag();
    }
    return *this;
}
//  += 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator+=(const Array2D<T>& vadd){
    for (size_t i(0); i< d1*d2; ++i) {
        (*v)[2*i]   += vadd.array()[i];
    }
    return *this;
}
//  += 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator+=(const Array2D_cmplx& vadd){
    (*v) += vadd.array();
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  -= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator-=(const T& d){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i] -=d;
    }
    return *this;
}
//  -= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator-=(const complex<T>& c){
    for (size_t i(0); i< td1*d2; i+=2) {
        (*v)[i]   -=c.real();
        (*v)[i+1] -=c.imag();
    }
    return *this;
}
//  -= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator-=(const Array2D<T>& vmin){
    for (size_t i(0); i< d1*d2; ++i) {
        (*v)[2*i]   -= vmin.array()[i];
    }
    return *this;
}
//  -= 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::operator-=(const Array2D_cmplx& vmin){
    (*v) -= vmin.array();
    return *this;
}

//--------------------------------------------------------------
//  Array2D * Vector
//--------------------------------------------------------------
//  Consider Array2D a list of d1-vectors and multiply with vmulti
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< td1*d2; j+= td1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[2*i+j]   *= vmulti[i];
            (*v)[2*i+1+j] *= vmulti[i];
        }
    }
    return *this;
}
//  Consider Array2D a list of d1-vectors and multiply with vmulti
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::multid1(const valarray< complex<T> >& vmulti){
    for (size_t j(0); j< td1*d2; j+= td1 ){
        for (size_t i(0); i< d1; ++i ){
            T vi((*v)[2*i+j]*vmulti[i].imag());
            (*v)[2*i+j]    *= vmulti[i].real();                     // v[i]   = Aa
            (*v)[2*i+j]    -= (*v)[2*i+1+j] * vmulti[i].imag();     // v[i]   = Aa - Bb
            (*v)[2*i+1+j]  *= vmulti[i].real();                     // v[i+1] = Ba
            (*v)[2*i+1+j]  += vi;           			// v[i+1] = Ba + Ab
        }
    }
    return *this;
}

//  Consider Array2D a list of d2-vectors and multiply with vmulti
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::multid2(const valarray<T>& vmulti){
    for (size_t j(0); j< d2; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; ++i ){
            (*v)[i] *= vmulti[j];
        }
    }
    return *this;
}

//  Consider Array2D a list of d2-vectors and multiply with vmulti
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::multid2(const valarray< complex<T> >& vmulti){
    for (size_t j(0); j< d2; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; i+=2 ){
            T vi((*v)[i]* vmulti[j].imag());
            (*v)[i]    *= vmulti[j].real();                 // v[i]   = Aa
            (*v)[i]    -= (*v)[i+1] * vmulti[j].imag();     // v[i]   = Aa - Bb
            (*v)[i+1]  *= vmulti[j].real();                 // v[i+1] = Ba
            (*v)[i+1]  += vi;           			// v[i+1] = Ba + Ab
        }
    }
    return *this;
}

//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
// Example: 0 4  8 12 16       -2 -2 -2 -2 -2 
//          1 5  9 13 17  -->  -2 -2 -2 -2 -2 
//          2 6 10 14 18       -2 -2 -2 -2 -2  
//          3 7 11 15 19       -2 -2 -2 -2 19  
// Requires at least 3 elements in d1 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::Dd1(){
    for(long i(0); i< td1*d2-4; ++i) {
        (*v)[i] -= (*v)[i+4];
    }
    for(long i(td1*d2-5); i>-1; --i) {
        (*v)[i+2] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference for elements with distance d1 (row-wise)
// Example: 0 4  8 12 16       -8 -8 -8 -8 16 
//          1 5  9 13 17  -->  -8 -8 -8 -8 17 
//          2 6 10 14 18       -8 -8 -8 -8 18
//          3 7 11 15 19       -8 -8 -8 -8 19
// Requires at least 3 elements in d2 
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::Dd2(){
    int twotd1 = 2*td1;
    for(long i(0); i< td1*d2-twotd1; ++i) {
        (*v)[i] -= (*v)[i+twotd1];
    }
    for(long i(td1*d2-twotd1-1); i>-1; --i) {
        (*v)[i+td1] = (*v)[i];
    }
    return *this;
}
//--------------------------------------------------------------

//  Remove data for N cells from dimension 1 (if N=0 don't remove anything)  
template<class T> Array2D_cmplx<T>& Array2D_cmplx<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2; ++j ){
        for (size_t i(j*td1); i< j*td1+2*N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}
//--------------------------------------------------------------
//**************************************************************


/**************************************************************
 *   3D Array Class
 *   Using Stroustrup's matrices p672-p673
 *   This container is similar to Fortran-style arrays.  
 *   Custom operations, no standard Matrix algebra.
 *   No error-checking
 *   Compile with -ftree-vectorize
 */
template<class T> class Array3D {
//--------------------------------------------------------------
//  3D Array decleration 
//--------------------------------------------------------------
private:
    valarray<T> *v;
    size_t  d1, d2, d3;            // for 2D VFP: p, x, y
    size_t d1d2;

public:
//      Constructors/Destructors
    Array3D(size_t x, size_t y, size_t z);
    Array3D(const Array3D& other);
    ~Array3D();

//      Basic info
    size_t dim()  const {return d1*d2*d3;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    size_t dim3() const {return d3;}
    valarray<T>&  array() const {return *v;};

//      Access
    T& operator()(size_t i, size_t j, size_t k);       // Fortran-style
    T  operator()(size_t i, size_t j, size_t k) const;
    T& operator()(size_t i);                           // 1D-style
    T  operator()(size_t i) const;

//TODO   vector<T>  d2d3c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);
    GSlice_iter<T>  d2c(size_t b, size_t e);
    GSlice_iter<T>  d3c(size_t b, size_t e);
    GSlice_iter<T>  SubArray3D(size_t st, size_t nx, size_t ny, size_t nz);

//      Constant slice iterators
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d3c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray3D(size_t st, size_t nx, size_t ny, size_t nz) const;

//      Operators
    Array3D& operator=(const T& d);
    Array3D& operator=(const Array3D& other);
    Array3D& operator*=(const T& d);
    Array3D& operator*=(const Array3D& vmulti);
    Array3D& operator+=(const T& d);
    Array3D& operator+=(const Array3D& vadd);
    Array3D& operator-=(const T& d);
    Array3D& operator-=(const Array3D& vmin);

//      Array * Vector
    Array3D& multid1(const valarray<T>& vmulti);// M*valarray(d1)
    Array3D& multid2(const valarray<T>& vmulti);// M*valarray(d2)
    Array3D& multid3(const valarray<T>& vmulti);// M*valarray(d3)

//      Array3D * Array2D
    Array3D& multid2d3(const Array2D<T>& vd2d3);// M*valarray(d3)

//      Central difference
    Array3D& Dd1(); // in the direction d1 (requires dim1() > 2)
    Array3D& Dd2(); // in the direction d2 (requires dim2() > 2)
    Array3D& Dd3(); // in the direction d3 (requires dim3() > 2)

//      Filter first N-cells in d1 direction 
    Array3D& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
template<class T> Array3D<T>::
Array3D(size_t x, size_t y, size_t z) : d1(x), d2(y), d3(z) {
    v = new valarray<T>(d1*d2*d3);
    d1d2 = d1*d2;
}
//  Copy constructor
template<class T> Array3D<T>:: Array3D(const Array3D& other){
    d1  = other.dim1();
    d2  = other.dim2();
    d3  = other.dim3();
    d1d2 = d1*d2;
    v = new valarray<T>(d1*d2*d3);
    (*v) = other.array();
}
//  Destructor
template<class T> Array3D<T>:: ~Array3D(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------

//  Access Fortan-style
template<class T>
inline T& Array3D<T>:: operator()(size_t i, size_t j, size_t k){
    return (*v)[i+j*d1+k*d1d2];
}

//  Access Fortan-style
template<class T>
inline T Array3D<T>:: operator()(size_t i, size_t j, size_t k) const {
    return (*v)[i+j*d1+k*d1d2];
}

//  1D style access
template<class T>
inline T& Array3D<T>:: operator()(size_t i){
    return (*v)[i];
}

//  1D style access
template<class T>
inline T Array3D<T>:: operator()(size_t i) const {
    return (*v)[i];
}
//--------------------------------------------------------------
//  Slicers
//--------------------------------------------------------------
//  slices (surfaces) for given d1 value range 
template<class T>
inline GSlice_iter<T> Array3D<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = d1; str[1] = d1d2; str[0] = 1;
    sz[2]  = d2; sz[1]  = d3;    sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b,sz,str));
}
//  const slices (surfaces) for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array3D<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = d1; str[1] = d1d2; str[0] = 1;
    sz[2]  = d2; sz[1]  = d3;    sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b,sz,str));
}
//  slices (surfaces) for given d2 value range 
template<class T>
inline GSlice_iter<T> Array3D<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = d1d2;  str[0] = d1;
    sz[2]  = d1; sz[1]  = d3;    sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(d1*b,sz,str));
}
//  const slices (surfaces) for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array3D<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = d1d2;  str[0] = d1;
    sz[2]  = d1; sz[1]  = d3;    sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(d1*b,sz,str));
}
//  slices (surfaces) for given d3 value range 
template<class T>
inline GSlice_iter<T> Array3D<T>::d3c(size_t b, size_t e){
    valarray<size_t> sz(2), str(2);
    str[1] = 1; str[0] = d1d2;
    sz[1]  = d1d2; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(d1d2*b,sz,str));
}
//  const slices (surfaces) for given d3 value range 
template<class T>
inline CGSlice_iter<T> Array3D<T>::d3c(size_t b, size_t e) const{
    valarray<size_t> sz(2), str(2);
    str[1] = 1; str[0] = d1d2;
    sz[1]  = d1d2; sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(d1d2*b,sz,str));
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  scan Subarray
template<class T>
inline GSlice_iter<T> Array3D<T>::SubArray3D(size_t st, size_t nx, size_t ny, size_t nz){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = dim1(); str[0] = d1d2;
    sz[2]  = nx;  sz[1] = ny;      sz[0] = nz;
    return GSlice_iter<T>(v,gslice(st,sz,str));
}

//  scan const Subarray
template<class T>
inline CGSlice_iter<T> Array3D<T>::SubArray3D(size_t st, size_t nx, size_t ny, size_t nz) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = dim1(); str[0] = d1d2;
    sz[2]  = nx;  sz[1] = ny;      sz[0] = nz;
    return CGSlice_iter<T>(v,gslice(st,sz,str));
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
template<class T> Array3D<T>& Array3D<T>::operator=(const T& d){
    (*v) = d;
    return *this;
}
//  Copy assignment operator
template<class T> Array3D<T>& Array3D<T>::operator=(const Array3D& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}

//  *= 
template<class T> Array3D<T>& Array3D<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
template<class T> Array3D<T>& Array3D<T>::operator*=(const Array3D& vmulti){
    (*v) *= vmulti.array();
    return *this;
}

//  +=
template<class T> Array3D<T>& Array3D<T>::operator+=(const T& d){
    (*v) +=d;
    return *this;
}
template<class T> Array3D<T>& Array3D<T>::operator+=(const Array3D& vadd){
    (*v) += vadd.array();
    return *this;
}

//  -= 
template<class T> Array3D<T>& Array3D<T>::operator-=(const T& d){
    (*v) -=d;
    return *this;
}
template<class T> Array3D<T>& Array3D<T>::operator-=(const Array3D& vmin){
    (*v) -= vmin.array();
    return *this;
}
//--------------------------------------------------------------
//  Array3D * Vector
//--------------------------------------------------------------
//  Vector multiplies every d1-length of Array3D
template<class T> Array3D<T>& Array3D<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< d1*d2*d3; j+=d1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[i+j] *= vmulti[i];
        }
    }
    return *this;
}

//  Vector multiplies every d2-length of Array3D
template<class T> Array3D<T>& Array3D<T>::multid2(const valarray<T>& vmulti){
    for (size_t k(0); k< d3*d2*d1; k+=d1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*d1); i< (j+1)*d1; ++i ){
                (*v)[i+k] *= vmulti[j];
            }
        }
    }
    return *this;
}

//  Vector multiplies every d3-length of Array3D
template<class T> Array3D<T>& Array3D<T>::multid3(const valarray<T>& vmulti){
    for (size_t k(0); k< d3; ++k) {
        for (size_t i(k*d1*d2); i< (k+1)*d1*d2; ++i ){
            (*v)[i] *= vmulti[k];
        }
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  B(j,k) multiplies every A(i,j,k) 
template<class T> Array3D<T>& Array3D<T>::multid2d3(const Array2D<T>& vd2d3){
    for (size_t j(0); j < d2*d3; ++j ){
        for (size_t i(j*d1); i< (j+1)*d1; ++i ){
            (*v)[i] *= vd2d3.array()[j];
        }
    }
    return *this;
}


//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
template<class T> Array3D<T>& Array3D<T>::Dd1(){
    for(long i(0); i< d1*d2*d3-2; ++i) {
        (*v)[i] -= (*v)[i+2];
    }
    for(long i(d1*d2*d3-3); i>-1; --i) {
        (*v)[i+1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d2 direction
template<class T> Array3D<T>& Array3D<T>::Dd2(){
    long twod1 = 2*d1;
    for(long i(0); i< d1*d2*d3-twod1; ++i) {
        (*v)[i] -= (*v)[i+twod1];
    }
    for(long i(d1*d2*d3-twod1-1); i>-1; --i) {
        (*v)[i+d1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d3 direction
template<class T> Array3D<T>& Array3D<T>::Dd3() {
    long twod1d2 = 2*d1d2;
    for(long i(0); i< d1*d2*d3-twod1d2; ++i) {
        (*v)[i] -= (*v)[i+twod1d2];
    }
    for(long i(d1*d2*d3-twod1d2-1); i>-1; --i) {
        (*v)[i+d1d2] = (*v)[i];
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Remove data for N cells from dimension 1  
template<class T> Array3D<T>& Array3D<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2*d3; ++j ){
        for (size_t i(j*d1); i< j*d1+N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}
//**************************************************************

/**************************************************************
*   3D Complex Array Class
*   Using Stroustrup's matrices p672-p673
*   This container is similar to Fortran-style arrays.
*   Custom operations, no standard Matrix algebra.
*   No error-checking
*   Compile with -ftree-vectorize
*/
template<class T> class Array3D_cmplx {
//--------------------------------------------------------------
//  3D Array decleration 
//--------------------------------------------------------------
private:
    valarray<T> *v;
    size_t  d1, d2, d3;            // for 2D VFP: p, x, y
    size_t td1;
    size_t td1d2;

public:
//      Constructors/Destructors
    Array3D_cmplx(size_t x, size_t y, size_t z);
    Array3D_cmplx(const Array3D_cmplx& other);
    ~Array3D_cmplx();

//      Basic info
    size_t dim()  const {return d1*d2*d3;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    size_t dim3() const {return d3;}
    valarray<T>&  array() const {return *v;};

//      Access
    T& operator[](size_t i);           // 1D-style
    T  operator[](size_t i) const;

    T& real(size_t i,size_t j,size_t k);           // 1D-style
    T  real(size_t i,size_t j,size_t k) const;
    T& imag(size_t i,size_t j,size_t k);           // 1D-style
    T  imag(size_t i,size_t j,size_t k) const;

//TODO   vector<T>  d2d3c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);
    GSlice_iter<T>  d2c(size_t b, size_t e);
    GSlice_iter<T>  d3c(size_t b, size_t e);
    GSlice_iter<T>  SubArray3D(size_t st, size_t nx, size_t ny, size_t nz);

//      Constant slice iterators
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d3c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray3D(size_t st, size_t nx, size_t ny, size_t nz) const;

//      Operators
    Array3D_cmplx& operator=(const T& d);
    Array3D_cmplx& operator=(const complex<T>& c);
    Array3D_cmplx& operator=(const Array3D<T>& other);
    Array3D_cmplx& operator=(const Array3D_cmplx& other);
    Array3D_cmplx& operator*=(const T& d);
    Array3D_cmplx& operator*=(const complex<T>& c);
    Array3D_cmplx& operator*=(const Array3D<T>& vmulti);
    Array3D_cmplx& operator*=(const Array3D_cmplx& vmulti);
    Array3D_cmplx& operator+=(const T& d);
    Array3D_cmplx& operator+=(const complex<T>& c);
    Array3D_cmplx& operator+=(const Array3D<T>& vadd);
    Array3D_cmplx& operator+=(const Array3D_cmplx& vadd);
    Array3D_cmplx& operator-=(const T& d);
    Array3D_cmplx& operator-=(const complex<T>& c);
    Array3D_cmplx& operator-=(const Array3D<T>& vmin);
    Array3D_cmplx& operator-=(const Array3D_cmplx& vmin);

//      Array * Vector
    Array3D_cmplx& multid1(const valarray<T>& vmulti);// M*valarray(d1)
    Array3D_cmplx& multid1(const valarray< complex<T> >& vmulti);// M*valarray(d1)
    Array3D_cmplx& multid2(const valarray<T>& vmulti);// M*valarray(d2)
    Array3D_cmplx& multid2(const valarray< complex<T> >& vmulti);// M*valarray(d2)
    Array3D_cmplx& multid3(const valarray<T>& vmulti);// M*valarray(d3)
    Array3D_cmplx& multid3(const valarray< complex<T>  >& vmulti);// M*valarray(d3)

//      Array3D * Array2D
    Array3D_cmplx& multid2d3(const Array2D<T>& vd2d3);
    Array3D_cmplx& multid2d3(const Array2D_cmplx<T>& vd2d3);

//      Central difference
    Array3D_cmplx& Dd1(); // in the direction d1 (requires dim1() > 2)
    Array3D_cmplx& Dd2(); // in the direction d2 (requires dim2() > 2)
    Array3D_cmplx& Dd3(); // in the direction d3 (requires dim3() > 2)

//      Filter first N-cells in d1 direction 
    Array3D_cmplx& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
template<class T> Array3D_cmplx<T>:: Array3D_cmplx(size_t x,size_t y,size_t z) : d1(x), d2(y), d3(z) {
    v = new valarray<T>(2*d1*d2*d3);
    td1 = 2*d1;
    td1d2 = 2*d1*d2;
}
//  Copy constructor
template<class T> Array3D_cmplx<T>:: Array3D_cmplx(const Array3D_cmplx& other){
    d1  = other.dim1();
    d2  = other.dim2();
    d3  = other.dim3();
    td1 = 2*d1;
    td1d2 = 2*d1*d2;
    v = new valarray<T>(2*d1*d2*d3);
    (*v) = other.array();
}
//  Destructor
template<class T> Array3D_cmplx<T>:: ~Array3D_cmplx(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  1D style access
template<class T> inline T& Array3D_cmplx<T>::operator[] (size_t i){
    return (*v)[i];
}
//  1D style const access
template<class T> inline T  Array3D_cmplx<T>::operator[] (size_t i) const {
    return (*v)[i];
}
//  Access real part by reference
template<class T> inline T& Array3D_cmplx<T>::real(size_t i,size_t j,size_t k){
    return (*v)[2*i+j*td1+k*td1d2];
}
//  Const access of real part 
template<class T> inline T Array3D_cmplx<T>::real(size_t i,size_t j,size_t k) const{
    return (*v)[2*i+j*td1+k*td1d2];
}
//  access of imag part 
template<class T> inline T& Array3D_cmplx<T>::imag(size_t i,size_t j,size_t k){
    return (*v)[2*i+1+j*td1+k*td1d2];
}
//  Const access of imag part 
template<class T> inline T Array3D_cmplx<T>::imag(size_t i,size_t j,size_t k) const{
    return (*v)[2*i+1+j*td1+k*td1d2];
}
//--------------------------------------------------------------
//  Slicers
//--------------------------------------------------------------
//  slices (surfaces) for given d1 value range 
template<class T>
inline GSlice_iter<T> Array3D_cmplx<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(4), str(4);
    str[3] = 1, str[2] = td1; str[1] = td1d2; str[0] = 2;
    sz[3]  = 2, sz[2]  = d2;  sz[1]  = d3;    sz[0]  = e-b+1;
    return GSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  const slices (surfaces) for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array3D_cmplx<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(4), str(4);
    str[3] = 1, str[2] = td1; str[1] = td1d2; str[0] = 2;
    sz[3]  = 2, sz[2]  = d2;  sz[1]  = d3;    sz[0]  = e-b+1;
    return CGSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  slices (surfaces) for given d2 value range 
template<class T>
inline GSlice_iter<T> Array3D_cmplx<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = td1d2;  str[0] = td1;
    sz[2]  = td1; sz[1]  = d3;    sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(td1*b,sz,str));
}
//  const slices (surfaces) for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array3D_cmplx<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;  str[1] = td1d2;  str[0] = td1;
    sz[2]  = td1; sz[1]  = d3;    sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(td1*b,sz,str));
}
//  slices (surfaces) for given d3 value range 
template<class T>
inline GSlice_iter<T> Array3D_cmplx<T>::d3c(size_t b, size_t e){
    valarray<size_t> sz(2), str(2);
    str[1] = 1; str[0] = td1d2;
    sz[1]  = td1d2; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(td1d2*b,sz,str));
}
//  const slices (surfaces) for given d3 value range 
template<class T>
inline CGSlice_iter<T> Array3D_cmplx<T>::d3c(size_t b, size_t e) const{
    valarray<size_t> sz(2), str(2);
    str[1] = 1; str[0] = td1d2;
    sz[1]  = td1d2; sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(td1d2*b,sz,str));
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  scan Subarray
template<class T>
inline GSlice_iter<T> Array3D_cmplx<T>::SubArray3D(size_t st, size_t nx, size_t ny, size_t nz){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;    str[1] = td1;  str[0] = td1d2;
    sz[2]  = 2*nx;  sz[1] = ny;    sz[0] = nz;
    return GSlice_iter<T>(v,gslice(2*st,sz,str));
}

//  scan const Subarray
template<class T>
inline CGSlice_iter<T> Array3D_cmplx<T>::SubArray3D(size_t st, size_t nx, size_t ny, size_t nz) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;    str[1] = td1;  str[0] = td1d2;
    sz[2]  = 2*nx;  sz[1] = ny;    sz[0] = nz;
    return CGSlice_iter<T>(v,gslice(2*st,sz,str));
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator=(const T& d){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i] = d;
        (*v)[i+1] = 0.0;
    }
    return *this;
}
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i]   =c.real();
        (*v)[i+1] =c.imag();
    }
    return *this;
}
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator=(const Array3D<T>& other){
    for (size_t i(0); i< d1*d2*d3; ++i) {
        (*v)[2*i]   = other.array()[i];
        (*v)[2*i+1] = 0.0;
    }
    return *this;
}
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator=(const Array3D_cmplx& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  *= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
//  *= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator*=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   c=a+ib
        T vi((*v)[i]*c.imag());
        (*v)[i]    *= c.real();                     // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*c.imag();           // v[i]   = Aa - Bb
        (*v)[i+1]  *= c.real();                     // v[i+1] = Ba
        (*v)[i+1]  += vi;                           // v[i+1] = Ba + Ab
    }
    return *this;
}
//  *= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator*=(const Array3D<T>& vmulti){
    for (size_t i(0); i< d1*d2*d3; ++i) {
        (*v)[2*i]   *= vmulti.array()[i];
        (*v)[2*i+1] *= vmulti.array()[i];
    }
    return *this;
}
//  *= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator*=(const Array3D_cmplx& vmulti){
    for (size_t i(0); i< td1*d2*d3; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   vmulti[i] + ivmulti[i+1] = a+ib
        T vi((*v)[i] * vmulti.array()[i+1]);
        (*v)[i]    *=  vmulti.array()[i];              // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*vmulti.array()[i+1];   // v[i]   = Aa - Bb
        (*v)[i+1]  *= vmulti.array()[i];               // v[i+1] = Ba
        (*v)[i+1]  += vi;                              // v[i+1] = Ba + Ab
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  += 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator+=(const T& d){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i] +=d;
    }
    return *this;
}
//  += 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator+=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i]   +=c.real();
        (*v)[i+1] +=c.imag();
    }
    return *this;
}
//  += 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator+=(const Array3D<T>& vadd){
    for (size_t i(0); i< d1*d2*d3; ++i) {
        (*v)[2*i]   += vadd.array()[i];
    }
    return *this;
}
//  += 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator+=(const Array3D_cmplx& vadd){
    (*v) += vadd.array();
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  -= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator-=(const T& d){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i] -=d;
    }
    return *this;
}
//  -= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator-=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3; i+=2) {
        (*v)[i]   -=c.real();
        (*v)[i+1] -=c.imag();
    }
    return *this;
}
//  -= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator-=(const Array3D<T>& vmin){
    for (size_t i(0); i< d1*d2*d3; ++i) {
        (*v)[2*i]   -= vmin.array()[i];
    }
    return *this;
}
//  -= 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::operator-=(const Array3D_cmplx& vmin){
    (*v) -= vmin.array();
    return *this;
}

//--------------------------------------------------------------
//  Array3D * Vector
//--------------------------------------------------------------
//  Vector multiplies every d1-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< td1*d2*d3; j+=td1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[2*i+j]   *= vmulti[i];
            (*v)[2*i+1+j] *= vmulti[i];
        }
    }
    return *this;
}

//  Vector multiplies every d1-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid1(const valarray< complex<T> >& vmulti){
    for (size_t j(0); j< td1*d2*d3; j+=td1 ){
        for (size_t i(0); i< d1; ++i ){
            T vi((*v)[2*i+j]*vmulti[i].imag());
            (*v)[2*i+j]    *= vmulti[i].real();                     // v[i]   = Aa
            (*v)[2*i+j]    -= (*v)[2*i+1+j] * vmulti[i].imag();     // v[i]   = Aa - Bb
            (*v)[2*i+1+j]  *= vmulti[i].real();                     // v[i+1] = Ba
            (*v)[2*i+1+j]  += vi;           			// v[i+1] = Ba + Ab
        }
    }
    return *this;
}

//  Vector multiplies every d2-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid2(const valarray<T>& vmulti){
    for (size_t k(0); k< d3*d2*td1; k+=td1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*td1); i< (j+1)*td1; ++i ){
                (*v)[i+k] *= vmulti[j];
            }
        }
    }
    return *this;
}

//  Vector multiplies every d2-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid2(const valarray< complex<T> >& vmulti){
    for (size_t k(0); k< d3*d2*td1; k+=td1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*td1); i< (j+1)*td1; i+=2 ){
                T vi((*v)[i+k]* vmulti[j].imag());
                (*v)[i+k]    *= vmulti[j].real();                 // v[i]   = Aa
                (*v)[i+k]    -= (*v)[i+k+1] * vmulti[j].imag();   // v[i]   = Aa - Bb
                (*v)[i+k+1]  *= vmulti[j].real();                 // v[i+1] = Ba
                (*v)[i+k+1]  += vi;           		       // v[i+1] = Ba + Ab
            }
        }
    }
    return *this;
}

//  Vector multiplies every d3-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid3(const valarray<T>& vmulti){
    for (size_t k(0); k< d3; ++k) {
        for (size_t i(k*td1*d2); i< (k+1)*td1*d2; ++i ){
            (*v)[i] *= vmulti[k];
        }
    }
    return *this;
}

//  Vector multiplies every d3-length of Array3D
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid3(const valarray< complex<T> >& vmulti){
    for (size_t k(0); k< d3; ++k) {
        for (size_t i(k*td1*d2); i< (k+1)*td1*d2; i+=2 ){
            T vi((*v)[i]* vmulti[k].imag());
            (*v)[i]    *= vmulti[k].real();               // v[i]   = Aa
            (*v)[i]    -= (*v)[i+1] * vmulti[k].imag();   // v[i]   = Aa - Bb
            (*v)[i+1]  *= vmulti[k].real();               // v[i+1] = Ba
            (*v)[i+1]  += vi;           		       // v[i+1] = Ba + Ab
        }
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  B(j,k) multiplies every A(i,j,k) 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid2d3(const Array2D<T>& vd2d3){
    for (size_t j(0); j < d2*d3; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; ++i ){
            (*v)[i] *= vd2d3.array()[j];
        }
    }
    return *this;
}


//  B(j,k) multiplies every A(i,j,k) 
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::multid2d3(const Array2D_cmplx<T>& vd2d3){
    for (size_t j(0); j < d2*d3; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; i+=2 ){
            T vi((*v)[i]*vd2d3.array()[2*j+1]);
            (*v)[i]    *= vd2d3.array()[2*j];               // v[i]   = Aa
            (*v)[i]    -= (*v)[i+1]*vd2d3.array()[2*j+1];   // v[i]   = Aa - Bb
            (*v)[i+1]  *= vd2d3.array()[2*j];               // v[i+1] = Ba
            (*v)[i+1]  += vi;           		         // v[i+1] = Ba + Ab
        }
    }
    return *this;
}
//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::Dd1(){
    for(long i(0); i< td1*d2*d3-4; ++i) {
        (*v)[i] -= (*v)[i+4];
    }
    for(long i(td1*d2*d3-5); i>-1; --i) {
        (*v)[i+2] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d2 direction
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::Dd2(){
    long twotd1 = 2*td1;
    for(long i(0); i< td1*d2*d3-twotd1; ++i) {
        (*v)[i] -= (*v)[i+twotd1];
    }
    for(long i(td1*d2*d3-twotd1-1); i>-1; --i) {
        (*v)[i+td1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d3 direction
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::Dd3() {
    long twotd1d2 = 2*td1d2;
    for(long i(0); i< long(td1*d2*d3)-twotd1d2; ++i) {
        (*v)[i] -= (*v)[i+twotd1d2];
    }
    for(long i(td1*d2*d3-twotd1d2-1); i>-1; --i) {
        (*v)[i+td1d2] = (*v)[i];
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Remove data for N cells from dimension 1  
template<class T> Array3D_cmplx<T>& Array3D_cmplx<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2*d3; ++j ){
        for (size_t i(j*td1); i< j*td1+2*N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}



/**************************************************************
*   4D Array Class
***************************************************************
*   Using Stroustrup's matrices p672-p673
*   This container is similar to Fortran-style arrays.
*   Custom operations, no standard Matrix algebra.
*   No error-checking
*   Compile with -ftree-vectorize
*/
template<class T> class Array4D {
//--------------------------------------------------------------
//  4D Array decleration 
//--------------------------------------------------------------

private:
    valarray<T> *v;
    size_t  d1, d2, d3, d4;            // For 3d VFP: p, x, y, z
    size_t d1d2;
    size_t d1d2d3;

public:
//      Constructors/Destructors
    Array4D(size_t x, size_t y, size_t z, size_t w);
    Array4D(const Array4D& other);
    ~Array4D();

//      Basic Info
    size_t dim()  const {return d1*d2*d3*d4;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    size_t dim3() const {return d3;}
    size_t dim4() const {return d4;}
    valarray<T>&  array() const {return *v;};

//      Access
    T& operator()(size_t i, size_t j, size_t k, size_t l);       // Fortran-style
    T  operator()(size_t i, size_t j, size_t k, size_t l) const;
    T& operator()(size_t i);                           // 1D-style
    T  operator()(size_t i) const;

//TODO   vector<T>  d2d3d4c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);
    GSlice_iter<T>  d2c(size_t b, size_t e);
    GSlice_iter<T>  d3c(size_t b, size_t e);
    GSlice_iter<T>  d4c(size_t b, size_t e);
    GSlice_iter<T>  SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw);

//      Constant slice iterators
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d3c(size_t b, size_t e) const;
    CGSlice_iter<T> d4c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw) const;

//      Operators
    Array4D& operator=(const T& d);
    Array4D& operator=(const Array4D& other);
    Array4D& operator*=(const T& d);
    Array4D& operator*=(const Array4D& vmulti);
    Array4D& operator+=(const T& d);
    Array4D& operator+=(const Array4D& vadd);
    Array4D& operator-=(const T& d);
    Array4D& operator-=(const Array4D& vmin);

//      Array * Vector
    Array4D& multid1(const valarray<T>& vmulti);// M*valarray(d1)
    Array4D& multid2(const valarray<T>& vmulti);// M*valarray(d2)
    Array4D& multid3(const valarray<T>& vmulti);// M*valarray(d3)
    Array4D& multid4(const valarray<T>& vmulti);// M*valarray(d3)

//      Array4D * Array3D 
    Array4D& multid2d3d4(const Array3D<T>& vd2d3d4);

//      Central difference
    Array4D& Dd1(); // in the direction d1
    Array4D& Dd2(); // in the direction d2
    Array4D& Dd3(); // in the direction d3
    Array4D& Dd4(); // in the direction d3

//      Filter first N-cells in d1 direction 
    Array4D& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
template<class T> Array4D<T>::
Array4D(size_t x, size_t y, size_t z, size_t w) : d1(x), d2(y), d3(z), d4(w) {
    v = new valarray<T>(d1*d2*d3*d4);
    d1d2   = d1*d2;
    d1d2d3 = d1*d2*d3;
}
//  Copy constructor
template<class T> Array4D<T>:: Array4D(const Array4D& other){
    d1  = other.dim1();
    d2  = other.dim2();
    d3  = other.dim3();
    d4  = other.dim4();
    d1d2 = d1*d2;
    d1d2d3 = d1*d2*d3;
    v = new valarray<T>(d1*d2*d3*d4);
    (*v) = other.array();
}
//  Destructor
template<class T> Array4D<T>:: ~Array4D(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------

//  Access Fortan-style
template<class T>
inline T& Array4D<T>:: operator()(size_t i, size_t j, size_t k, size_t l){
    return (*v)[i+j*d1+k*d1d2+l*d1d2d3];
}

//  Access Fortan-style
template<class T>
inline T Array4D<T>:: operator()(size_t i, size_t j, size_t k, size_t l) const {
    return (*v)[i+j*d1+k*d1d2+l*d1d2d3];
}

//  1D style access
template<class T>
inline T& Array4D<T>:: operator()(size_t i){
    return (*v)[i];
}

//  1D style access
template<class T>
inline T Array4D<T>:: operator()(size_t i) const {
    return (*v)[i];
}
//--------------------------------------------------------------

//  slices (cubes) for given d1 value range 
template<class T>
inline GSlice_iter<T> Array4D<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(4), str(4);
    str[3] = d1; str[2] = d1d2; str[1] = d1d2d3; str[0] = 1;
    sz[3]  = d2; sz[2]  = d3;   sz[1] = d4; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b,sz,str));
}
//  const slices (cubes) for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array4D<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(4), str(4);
    str[3] = d1; str[2] = d1d2; str[1] = d1d2d3; str[0] = 1;
    sz[3]  = d2; sz[2]  = d3;   sz[1] = d4; sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b,sz,str));
}
//  slices (cubes) for given d2 value range 
template<class T>
inline GSlice_iter<T> Array4D<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(4), str(4);
    str[3] = 1;  str[2] = d1d2; str[1] = d1d2d3; str[0] = d1;
    sz[3]  = d1; sz[2]  = d3;   sz[1] = d4; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b*d1,sz,str));
}
//  const slices (cubes) for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array4D<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(4), str(4);
    str[3] = 1;  str[2] = d1d2; str[1] = d1d2d3; str[0] = d1;
    sz[3]  = d1; sz[2]  = d3;   sz[1] = d4; sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*d1,sz,str));
}
//  slices (cubes) for given d3 value range 
template<class T>
inline GSlice_iter<T> Array4D<T>::d3c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;    str[1] = d1d2d3; str[0] = d1d2;
    sz[2]  = d1d2; sz[1]  = d4;     sz[0]  = e-b+1;
    return GSlice_iter<T>(v,gslice(b*d1d2,sz,str));
}
//  const slices (cubes) for given d3 value range 
template<class T>
inline CGSlice_iter<T> Array4D<T>::d3c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;    str[1] = d1d2d3; str[0] = d1d2;
    sz[2]  = d1d2; sz[1]  = d4;     sz[0]  = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*d1d2,sz,str));
}
//  slices (cubes) for given d4 value range 
template<class T>
inline GSlice_iter<T> Array4D<T>::d4c(size_t b, size_t e) {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;     str[0] = d1d2d3;
    sz[1]  = d1d2d3; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b*d1d2d3,sz,str));
}
//  const slices (cubes) for given d4 value range 
template<class T>
inline CGSlice_iter<T> Array4D<T>::d4c(size_t b, size_t e) const {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;     str[0] = d1d2d3;
    sz[1]  = d1d2d3; sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*d1d2d3,sz,str));
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ---------------------

//  Generate Subarray
template<class T>
inline GSlice_iter<T> Array4D<T>::SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw){
    valarray<size_t> sz(4), str(4);
    str[3] = 1;   str[2] = d1;  str[1] = d1d2; str[0] = d1d2d3;
    sz[3]  = nx;  sz[2]  = ny;  sz[1]  = nz;   sz[0]  = nw;
    return GSlice_iter<T>(v,gslice(st,sz,str));
}

//  Generate Subarray
template<class T>
inline CGSlice_iter<T> Array4D<T>::SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw) const{
    valarray<size_t> sz(3), str(3);
    str[3] = 1;   str[2] = d1;  str[1] = d1d2; str[0] = d1d2d3;
    sz[3]  = nx;  sz[2]  = ny;  sz[1]  = nz;   sz[0]  = nw;
    return CGSlice_iter<T>(v,gslice(st,sz,str));
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
template<class T> Array4D<T>& Array4D<T>::operator=(const T& d){
    (*v) = d;
    return *this;
}
//  Copy assignment operator
template<class T> Array4D<T>& Array4D<T>::operator=(const Array4D& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}

//  *= 
template<class T> Array4D<T>& Array4D<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
template<class T> Array4D<T>& Array4D<T>::operator*=(const Array4D& vmulti){
    (*v) *= vmulti.array();
    return *this;
}

//  +=
template<class T> Array4D<T>& Array4D<T>::operator+=(const T& d){
    (*v) +=d;
    return *this;
}
template<class T> Array4D<T>& Array4D<T>::operator+=(const Array4D& vadd){
    (*v) += vadd.array();
    return *this;
}

//  -= 
template<class T> Array4D<T>& Array4D<T>::operator-=(const T& d){
    (*v) -=d;
    return *this;
}
template<class T> Array4D<T>& Array4D<T>::operator-=(const Array4D& vmin){
    (*v) -= vmin.array();
    return *this;
}

//--------------------------------------------------------------
//  Array4D * Vector
//--------------------------------------------------------------
//  Vector multiplies every d1-length of Array4D
template<class T> Array4D<T>& Array4D<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< dim(); j+=d1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[i+j] *= vmulti[i];
        }
    }
    //for(size_t i(0); i< d1; ++i) {
    //    for(GSlice_iter<T> it(d1c(i,i)); it!=it.end(); ++it) *it *= vmulti[i];
    // }
    return *this;
}

//  Vector multiplies every d2-length of Array4D
template<class T> Array4D<T>& Array4D<T>::multid2(const valarray<T>& vmulti){
    for (size_t k(0); k< dim(); k+=d1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*d1); i< (j+1)*d1; ++i ){
                (*v)[i+k] *= vmulti[j];
            }
        }
    }
    /*for(size_t j(0); j< d2; ++j) {
        for(GSlice_iter<T> it(d2c(j,j)); it!=it.end(); ++it) *it *= vmulti[j];
    }*/
    return *this;
}

//  Vector multiplies every d3-length of Array4D
template<class T> Array4D<T>& Array4D<T>::multid3(const valarray<T>& vmulti){
    for (size_t l(0); l< dim(); l+=d1*d2*d3) {
        for (size_t k(0); k< d3; ++k) {
            for (size_t i(k*d1*d2); i< (k+1)*d1*d2; ++i ){
                (*v)[i+l] *= vmulti[k];
            }
        }
    }
    //for(size_t k(0); k< d3; ++k) {
    //    for(GSlice_iter<T> it(d3c(k,k)); it!=it.end(); ++it) *it *= vmulti[k];
    //}
    return *this;
}

//  Vector multiplies every d3-length of Array4D
template<class T> Array4D<T>& Array4D<T>::multid4(const valarray<T>& vmulti){
    for (size_t k(0); k< d4; ++k) {
        for (size_t i(k*d1*d2*d3); i< (k+1)*d1*d2*d3; ++i ){
            (*v)[i] *= vmulti[k];
        }
    }
    //for(size_t l(0); l< d4; ++l) {
    //    for(GSlice_iter<T> it(d4c(l,l)); it!=it.end(); ++it) *it *= vmulti[l];
    //}
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  B(j,k,l) multiplies every A(i,j,k,l) 
template<class T> Array4D<T>& Array4D<T>::multid2d3d4(const Array3D<T>& vd2d3d4){
    for (size_t j(0); j < d2*d3*d4; ++j ){
        for (size_t i(j*d1); i< (j+1)*d1; ++i ){
            (*v)[i] *= vd2d3d4.array()[j];
        }
    }
    return *this;
}

//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
template<class T> Array4D<T>& Array4D<T>::Dd1(){
    for(long i(0); i< dim()-2; ++i) {
        (*v)[i] -= (*v)[i+2];
    }
    for(long i(dim()-3); i>-1; --i) {
        (*v)[i+1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d2 direction
template<class T> Array4D<T>& Array4D<T>::Dd2(){
    long twod1 = 2*d1;
    for(long i(0); i< dim()-twod1; ++i) {
        (*v)[i] -= (*v)[i+twod1];
    }
    for(long i(dim()-twod1-1); i>-1; --i) {
        (*v)[i+d1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d3 direction
template<class T> Array4D<T>& Array4D<T>::Dd3() {
    long twod1d2 = 2*d1d2;
    for(long i(0); i< dim()-twod1d2; ++i) {
        (*v)[i] -= (*v)[i+twod1d2];
    }
    for(long i(dim()-twod1d2-1); i>-1; --i) {
        (*v)[i+d1d2] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d4 direction
template<class T> Array4D<T>& Array4D<T>::Dd4() {
    long twod1d2d3 = 2*d1d2d3;
    for(long i(0); i< dim()-twod1d2d3; ++i) {
        (*v)[i] -= (*v)[i+twod1d2d3];
    }
    for(long i(dim()-twod1d2d3-1); i>-1; --i) {
        (*v)[i+d1d2d3] = (*v)[i];
    }
    return *this;
}


//  Remove data for N cells from dimension 1  
template<class T> Array4D<T>& Array4D<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2*d3*d4; ++j ){
        for (size_t i(j*d1); i< j*d1+N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//   4D Array Class
//**************************************************************
//   Using Stroustrup's matrices p672-p673
//   This container is similar to Fortran-style arrays.  
//   Custom operations, no standard Matrix algebra.
//   No error-checking
//   Compile with -ftree-vectorize
//--------------------------------------------------------------
template<class T> class Array4D_cmplx {
//--------------------------------------------------------------
//  4D Array decleration 
//--------------------------------------------------------------

private:
    valarray<T> *v;
    size_t  d1, d2, d3, d4;            // For 3d VFP: p, x, y, z
    size_t td1;
    size_t td1d2;
    size_t td1d2d3;

public:
//      Constructors/Destructors
    Array4D_cmplx(size_t x, size_t y, size_t z, size_t w);
    Array4D_cmplx(const Array4D_cmplx& other);
    ~Array4D_cmplx();

//      Basic Info
    size_t dim()  const {return d1*d2*d3*d4;}
    size_t dim1() const {return d1;}
    size_t dim2() const {return d2;}
    size_t dim3() const {return d3;}
    size_t dim4() const {return d4;}
    valarray<T>&  array() const {return *v;};

//      Access
    T& operator[](size_t i);           // 1D-style
    T  operator[](size_t i) const;

    T& real(size_t i,size_t j,size_t k,size_t l);           // 1D-style
    T  real(size_t i,size_t j,size_t k,size_t l) const;
    T& imag(size_t i,size_t j,size_t k,size_t l);           // 1D-style
    T  imag(size_t i,size_t j,size_t k,size_t l) const;

//TODO   vector<T>  d2d3d4c(size_t j);

//      Slice iterators
    GSlice_iter<T>  d1c(size_t b, size_t e);
    GSlice_iter<T>  d2c(size_t b, size_t e);
    GSlice_iter<T>  d3c(size_t b, size_t e);
    GSlice_iter<T>  d4c(size_t b, size_t e);
    GSlice_iter<T>  SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw);

//      Constant slice iterators
    CGSlice_iter<T> d1c(size_t b, size_t e) const;
    CGSlice_iter<T> d2c(size_t b, size_t e) const;
    CGSlice_iter<T> d3c(size_t b, size_t e) const;
    CGSlice_iter<T> d4c(size_t b, size_t e) const;
    CGSlice_iter<T> SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw) const;

//      Operators
    Array4D_cmplx& operator=(const T& d);
    Array4D_cmplx& operator=(const complex<T>& c);
    Array4D_cmplx& operator=(const Array4D<T>& other);
    Array4D_cmplx& operator=(const Array4D_cmplx& other);
    Array4D_cmplx& operator*=(const T& d);
    Array4D_cmplx& operator*=(const complex<T>& c);
    Array4D_cmplx& operator*=(const Array4D<T>& vmulti);
    Array4D_cmplx& operator*=(const Array4D_cmplx& vmulti);
    Array4D_cmplx& operator+=(const T& d);
    Array4D_cmplx& operator+=(const complex<T>& c);
    Array4D_cmplx& operator+=(const Array4D<T>& vadd);
    Array4D_cmplx& operator+=(const Array4D_cmplx& vadd);
    Array4D_cmplx& operator-=(const T& d);
    Array4D_cmplx& operator-=(const complex<T>& c);
    Array4D_cmplx& operator-=(const Array4D<T>& vmin);
    Array4D_cmplx& operator-=(const Array4D_cmplx& vmin);

//      Array * Vector
    Array4D_cmplx& multid1(const valarray<T>& vmulti);            // M*valarray(d1)
    Array4D_cmplx& multid1(const valarray< complex<T> >& vmulti); // M*valarray(d1)
    Array4D_cmplx& multid2(const valarray<T>& vmulti);            // M*valarray(d2)
    Array4D_cmplx& multid2(const valarray< complex<T> >& vmulti); // M*valarray(d2)
    Array4D_cmplx& multid3(const valarray<T>& vmulti);            // M*valarray(d3)
    Array4D_cmplx& multid3(const valarray< complex<T>  >& vmulti);// M*valarray(d3)
    Array4D_cmplx& multid4(const valarray<T>& vmulti);            // M*valarray(d4)
    Array4D_cmplx& multid4(const valarray< complex<T>  >& vmulti);// M*valarray(d4)

//      Array4D * Array3D
    Array4D_cmplx& multid2d3d4(const Array3D<T>& vd2d3d4);
    Array4D_cmplx& multid2d3d4(const Array3D_cmplx<T>& vd2d3d4);

//      Central difference
    Array4D_cmplx& Dd1(); // in the direction d1
    Array4D_cmplx& Dd2(); // in the direction d2
    Array4D_cmplx& Dd3(); // in the direction d3
    Array4D_cmplx& Dd4(); // in the direction d3

//      Filter first N-cells in d1 direction 
    Array4D_cmplx& Filterd1(size_t N);
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
template<class T> Array4D_cmplx<T>::
Array4D_cmplx(size_t x,size_t y,size_t z,size_t w) : d1(x), d2(y), d3(z), d4(w) {
    v = new valarray<T>(2*d1*d2*d3*d4);
    td1 = 2*d1;
    td1d2 = 2*d1*d2;
    td1d2d3 = 2*d1*d2*d3;
}
//  Copy constructor
template<class T> Array4D_cmplx<T>:: Array4D_cmplx(const Array4D_cmplx& other){
    d1  = other.dim1();
    d2  = other.dim2();
    d3  = other.dim3();
    d4  = other.dim4();
    td1 = 2*d1;
    td1d2 = 2*d1*d2;
    td1d2d3 = 2*d1*d2*d3;
    v = new valarray<T>(2*d1*d2*d3*d4);
    (*v) = other.array();
}
//  Destructor
template<class T> Array4D_cmplx<T>:: ~Array4D_cmplx(){
    delete v;
}

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  1D style access
template<class T> inline T& Array4D_cmplx<T>::operator[] (size_t i){
    return (*v)[i];
}
//  1D style const access
template<class T> inline T  Array4D_cmplx<T>::operator[] (size_t i) const {
    return (*v)[i];
}
//  Access real part by reference
template<class T> inline T& Array4D_cmplx<T>::real(size_t i,size_t j,size_t k,size_t l){
    return (*v)[2*i+j*td1+k*td1d2+l*td1d2d3];
}
//  Const access of real part 
template<class T> inline T Array4D_cmplx<T>::real(size_t i,size_t j,size_t k,size_t l) const{
    return (*v)[2*i+j*td1+k*td1d2+l*td1d2d3];
}
//  access of imag part 
template<class T> inline T& Array4D_cmplx<T>::imag(size_t i,size_t j,size_t k,size_t l){
    return (*v)[2*i+1+j*td1+k*td1d2+l*td1d2d3];
}
//  Const access of imag part 
template<class T> inline T Array4D_cmplx<T>::imag(size_t i,size_t j,size_t k,size_t l) const{
    return (*v)[2*i+1+j*td1+k*td1d2+l*td1d2d3];
}
//--------------------------------------------------------------
//  Slicers
//--------------------------------------------------------------

//  slices (cubes) for given d1 value range 
template<class T>
inline GSlice_iter<T> Array4D_cmplx<T>::d1c(size_t b, size_t e){
    valarray<size_t> sz(5), str(5);
    str[4] = 1, str[3] = td1; str[2] = td1d2; str[1] = td1d2d3; str[0] = 2;
    sz[4] = 2,  sz[3]  = d2;  sz[2]  = d3;    sz[1]  = d4;       sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  const slices (cubes) for given d1 value range 
template<class T>
inline CGSlice_iter<T> Array4D_cmplx<T>::d1c(size_t b, size_t e) const{
    valarray<size_t> sz(5), str(5);
    str[4] = 1, str[3] = td1; str[2] = td1d2; str[1] = td1d2d3; str[0] = 2;
    sz[4] = 2,  sz[3]  = d2;  sz[2]  = d3;    sz[1]  = d4;       sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(2*b,sz,str));
}
//  slices (cubes) for given d2 value range 
template<class T>
inline GSlice_iter<T> Array4D_cmplx<T>::d2c(size_t b, size_t e){
    valarray<size_t> sz(4), str(4);
    str[3] = 1;  str[2] = td1d2; str[1] = td1d2d3; str[0] = td1;
    sz[3]  = td1; sz[2]  = d3;   sz[1] = d4; sz[0] = e-b+1;
    return GSlice_iter<T>(v,gslice(b*td1,sz,str));
}
//  const slices (cubes) for given d2 value range 
template<class T>
inline CGSlice_iter<T> Array4D_cmplx<T>::d2c(size_t b, size_t e) const{
    valarray<size_t> sz(4), str(4);
    str[3] = 1;  str[2] = td1d2; str[1] = td1d2d3; str[0] = td1;
    sz[3]  = td1; sz[2] = d3;    sz[1]  = d4;       sz[0] = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*td1,sz,str));
}
//  slices (cubes) for given d3 value range 
template<class T>
inline GSlice_iter<T> Array4D_cmplx<T>::d3c(size_t b, size_t e){
    valarray<size_t> sz(3), str(3);
    str[2] = 1;     str[1] = td1d2d3; str[0] = td1d2;
    sz[2]  = td1d2; sz[1]  = d4;      sz[0]  = e-b+1;
    return GSlice_iter<T>(v,gslice(b*td1d2,sz,str));
}
//  const slices (cubes) for given d3 value range 
template<class T>
inline CGSlice_iter<T> Array4D_cmplx<T>::d3c(size_t b, size_t e) const{
    valarray<size_t> sz(3), str(3);
    str[2] = 1;     str[1] = td1d2d3; str[0] = td1d2;
    sz[2]  = td1d2; sz[1]  = d4;      sz[0]  = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*td1d2,sz,str));
}
//  slices (cubes) for given d4 value range 
template<class T>
inline GSlice_iter<T> Array4D_cmplx<T>::d4c(size_t b, size_t e) {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;       str[0] = td1d2d3;
    sz[1]  = td1d2d3; sz[0]  = e-b+1;
    return GSlice_iter<T>(v,gslice(b*td1d2d3,sz,str));
}
//  const slices (cubes) for given d4 value range 
template<class T>
inline CGSlice_iter<T> Array4D_cmplx<T>::d4c(size_t b, size_t e) const {
    valarray<size_t> sz(2), str(2);
    str[1] = 1;       str[0] = td1d2d3;
    sz[1]  = td1d2d3; sz[0]  = e-b+1;
    return CGSlice_iter<T>(v,gslice(b*td1d2d3,sz,str));
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ---------------------

//  Generate Subarray
template<class T>
inline GSlice_iter<T> Array4D_cmplx<T>::
SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw){
    valarray<size_t> sz(4), str(4);
    str[3] = 1;   str[2] = td1;  str[1] = td1d2; str[0] = td1d2d3;
    sz[3]  = 2*nx;  sz[2]  = ny;  sz[1]  = nz;   sz[0]  = nw;
    return GSlice_iter<T>(v,gslice(2*st,sz,str));
}

//  Generate Subarray
template<class T>
inline CGSlice_iter<T> Array4D_cmplx<T>::
SubArray4D(size_t st, size_t nx, size_t ny, size_t nz, size_t nw) const{
    valarray<size_t> sz(4), str(4);
    str[3] = 1;   str[2] = td1;  str[1] = td1d2; str[0] = td1d2d3;
    sz[3]  = 2*nx;  sz[2]  = ny;  sz[1]  = nz;   sz[0]  = nw;
    return CGSlice_iter<T>(v,gslice(2*st,sz,str));
}

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator=(const T& d){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i] = d;
        (*v)[i+1] = 0.0;
    }
    return *this;
}
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i]   =c.real();
        (*v)[i+1] =c.imag();
    }
    return *this;
}
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator=(const Array4D<T>& other){
    for (size_t i(0); i< d1*d2*d3*d4; ++i) {
        (*v)[2*i]   = other.array()[i];
        (*v)[2*i+1] = 0.0;
    }
    return *this;
}
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator=(const Array4D_cmplx& other){
    if (this != &other) {   //self-assignment
        (*v) = other.array();
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  *= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator*=(const T& d){
    (*v) *=d;
    return *this;
}
//  *= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator*=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   c=a+ib
        T vi((*v)[i]*c.imag());
        (*v)[i]    *= c.real();                     // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*c.imag();           // v[i]   = Aa - Bb
        (*v)[i+1]  *= c.real();                     // v[i+1] = Ba
        (*v)[i+1]  += vi;                           // v[i+1] = Ba + Ab
    }
    return *this;
}
//  *= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator*=(const Array4D<T>& vmulti){
    for (size_t i(0); i< d1*d2*d3*d4; ++i) {
        (*v)[2*i]   *= vmulti.array()[i];
        (*v)[2*i+1] *= vmulti.array()[i];
    }
    return *this;
}
//  *= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator*=(const Array4D_cmplx& vmulti){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {            // suppose v[i]+iv[i+1]=A+iB,   vmulti[i] + ivmulti[i+1] = a+ib
        T vi((*v)[i] * vmulti.array()[i+1]);
        (*v)[i]    *=  vmulti.array()[i];              // v[i]   = Aa
        (*v)[i]    -= (*v)[i+1]*vmulti.array()[i+1];   // v[i]   = Aa - Bb
        (*v)[i+1]  *= vmulti.array()[i];               // v[i+1] = Ba
        (*v)[i+1]  += vi;                              // v[i+1] = Ba + Ab
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  += 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator+=(const T& d){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i] +=d;
    }
    return *this;
}
//  += 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator+=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i]   +=c.real();
        (*v)[i+1] +=c.imag();
    }
    return *this;
}
//  += 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator+=(const Array4D<T>& vadd){
    for (size_t i(0); i< d1*d2*d3*d4; ++i) {
        (*v)[2*i]   += vadd.array()[i];
    }
    return *this;
}
//  += 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator+=(const Array4D_cmplx& vadd){
    (*v) += vadd.array();
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  -= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator-=(const T& d){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i] -=d;
    }
    return *this;
}
//  -= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator-=(const complex<T>& c){
    for (size_t i(0); i< td1*d2*d3*d4; i+=2) {
        (*v)[i]   -=c.real();
        (*v)[i+1] -=c.imag();
    }
    return *this;
}
//  -= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator-=(const Array4D<T>& vmin){
    for (size_t i(0); i< d1*d2*d3*d4; ++i) {
        (*v)[2*i]   -= vmin.array()[i];
    }
    return *this;
}
//  -= 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::operator-=(const Array4D_cmplx& vmin){
    (*v) -= vmin.array();
    return *this;
}

//--------------------------------------------------------------
//  Array4D * Vector
//--------------------------------------------------------------
//  Vector multiplies every d1-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid1(const valarray<T>& vmulti){
    for (size_t j(0); j< 2*dim(); j+= td1 ){
        for (size_t i(0); i< d1; ++i ){
            (*v)[2*i+j]   *= vmulti[i];
            (*v)[2*i+1+j] *= vmulti[i];
        }
    }
    return *this;
}
//  Vector multiplies every d1-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid1(const valarray< complex<T> >& vmulti){
    for (size_t j(0); j< 2*dim(); j+=td1 ){
        for (size_t i(0); i< d1; ++i ){
            T vi((*v)[2*i+j]*vmulti[i].imag());
            (*v)[2*i+j]    *= vmulti[i].real();                     // v[i]   = Aa
            (*v)[2*i+j]    -= (*v)[2*i+1+j] * vmulti[i].imag();     // v[i]   = Aa - Bb
            (*v)[2*i+1+j]  *= vmulti[i].real();                     // v[i+1] = Ba
            (*v)[2*i+1+j]  += vi;           			// v[i+1] = Ba + Ab
        }
    }
    return *this;
}

//  Vector multiplies every d2-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid2(const valarray<T>& vmulti){
    for (size_t k(0); k< 2*dim(); k+=td1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*td1); i< (j+1)*td1; ++i ){
                (*v)[i+k] *= vmulti[j];
            }
        }
    }
    return *this;
}
//  Vector multiplies every d2-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid2(const valarray< complex<T> >& vmulti){
    for (size_t k(0); k< 2*dim(); k+=td1*d2 ) {
        for (size_t j(0); j< d2; ++j ) {
            for (size_t i(j*td1); i< (j+1)*td1; i+=2 ){
                T vi((*v)[i+k]* vmulti[j].imag());
                (*v)[i+k]    *= vmulti[j].real();                 // v[i]   = Aa
                (*v)[i+k]    -= (*v)[i+k+1] * vmulti[j].imag();   // v[i]   = Aa - Bb
                (*v)[i+k+1]  *= vmulti[j].real();                 // v[i+1] = Ba
                (*v)[i+k+1]  += vi;           		       // v[i+1] = Ba + Ab
            }
        }
    }
    return *this;
}

//  Vector multiplies every d3-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid3(const valarray<T>& vmulti){
    for (size_t l(0); l< 2*dim(); l+=td1*d2*d3) {
        for (size_t k(0); k< d3; ++k) {
            for (size_t i(k*td1*d2); i< (k+1)*td1*d2; ++i ){
                (*v)[i+l] *= vmulti[k];
            }
        }
    }
    return *this;
}
//  Vector multiplies every d3-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid3(const valarray< complex<T> >& vmulti){
    for (size_t l(0); l< 2*dim(); l+=td1*d2*d3) {
        for (size_t k(0); k< d3; ++k) {
            for (size_t i(k*td1*d2); i< (k+1)*td1*d2; i+=2 ){
                T vi((*v)[i+l]* vmulti[k].imag());
                (*v)[i+l]    *= vmulti[k].real();                 // v[i]   = Aa
                (*v)[i+l]    -= (*v)[i+l+1] * vmulti[k].imag();   // v[i]   = Aa - Bb
                (*v)[i+l+1]  *= vmulti[k].real();                 // v[i+1] = Ba
                (*v)[i+l+1]  += vi;           		      // v[i+1] = Ba + Ab
            }
        }
    }
    return *this;
}

//  Vector multiplies every d4-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid4(const valarray<T>& vmulti){
    for (size_t k(0); k< d4; ++k) {
        for (size_t i(k*td1*d2*d3); i< (k+1)*td1*d2*d3; ++i ){
            (*v)[i] *= vmulti[k];
        }
    }
    return *this;
}
//  Vector multiplies every d4-length of Array4D
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid4(const valarray< complex<T> >& vmulti){
    for (size_t k(0); k< d4; ++k) {
        for (size_t i(k*td1*d2*d3); i< (k+1)*td1*d2*d3; i+=2 ){
            T vi((*v)[i]* vmulti[k].imag());
            (*v)[i]    *= vmulti[k].real();               // v[i]   = Aa
            (*v)[i]    -= (*v)[i+1] * vmulti[k].imag();   // v[i]   = Aa - Bb
            (*v)[i+1]  *= vmulti[k].real();               // v[i+1] = Ba
            (*v)[i+1]  += vi;           		       // v[i+1] = Ba + Ab
        }
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  B(j,k) multiplies every A(i,j,k) 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid2d3d4(const Array3D<T>& vd2d3d4){
    for (size_t j(0); j < d2*d3*d4; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; ++i ){
            (*v)[i] *= vd2d3d4.array()[j];
        }
    }
    return *this;
}


//  B(j,k) multiplies every A(i,j,k) 
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::multid2d3d4(const Array3D_cmplx<T>& vd2d3d4){
    for (size_t j(0); j < d2*d3*d4; ++j ){
        for (size_t i(j*td1); i< (j+1)*td1; i+=2 ){
            T vi((*v)[i]*vd2d3d4.array()[2*j+1]);
            (*v)[i]    *= vd2d3d4.array()[2*j];               // v[i]   = Aa
            (*v)[i]    -= (*v)[i+1]*vd2d3d4.array()[2*j+1];   // v[i]   = Aa - Bb
            (*v)[i+1]  *= vd2d3d4.array()[2*j];               // v[i+1] = Ba
            (*v)[i+1]  += vi;           		           // v[i+1] = Ba + Ab
        }
    }
    return *this;
}

//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::Dd1(){
    for(long i(0); i< 2*dim()-4; ++i) {
        (*v)[i] -= (*v)[i+4];
    }
    for(long i(2*dim()-5); i>-1; --i) {
        (*v)[i+2] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d2 direction
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::Dd2(){
    long twotd1 = 2*td1;
    for(long i(0); i< 2*dim()-twotd1; ++i) {
        (*v)[i] -= (*v)[i+twotd1];
    }
    for(long i(2*dim()-twotd1-1); i>-1; --i) {
        (*v)[i+td1] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d3 direction
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::Dd3() {
    long twotd1d2 = 2*td1d2;
    for(long i(0); i< 2*dim()-twotd1d2; ++i) {
        (*v)[i] -= (*v)[i+twotd1d2];
    }
    for(long i(2*dim()-twotd1d2-1); i>-1; --i) {
        (*v)[i+td1d2] = (*v)[i];
    }
    return *this;
}

// (minus) Central difference in the d3 direction
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::Dd4() {
    long twotd1d2d3 = 2*td1d2d3;
    for(long i(0); i< long(2*dim())-twotd1d2d3; ++i) {
        (*v)[i] -= (*v)[i+twotd1d2d3];
    }
    for(long i(2*dim()-twotd1d2d3-1); i>-1; --i) {
        (*v)[i+td1d2d3] = (*v)[i];
    }
    return *this;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Remove data for N cells from dimension 1  
template<class T> Array4D_cmplx<T>& Array4D_cmplx<T>::Filterd1(size_t N) {
    for (size_t j(0); j < d2*d3*d4; ++j ){
        for (size_t i(j*td1); i< j*td1+2*N; ++i ){
            (*v)[i] = 0.0;
        }
    }
    return *this;
}
//--------------------------------------------------------------
//**************************************************************


#endif
