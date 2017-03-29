/*! \brief Fields, Distributions, Harmonics, States - Declarations
* \author PICKSC
 * \date   September 1, 2016
 * \file   state.h
 * 
 *   
 *   This header file contains the Declarations for the data
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

 * 	 2.A. SHarmonic2D : 
 *		 A wrapper class for Array3D<...>. It is used to
 *       describe a single spherical harmonic in 2D.
 *
 *   2.B. Field2D :
 *       Similar to Field1D but made of a Array2D<complex>
 * 
 *   2.C. EMF2D:
 *       The 6 2D Fields in a 2D-3P code. 
 *
 *   2.D. DistFunc2D:
 *       The collection of 2D spherical harmonics. 
 *
 *   2.E. State2D:
 *       The collection of 2D distributions and an EMF2D object describing all the fields. 
 */
//-------------------------------------------------------------------

#ifndef DECL_STATE_H
#define DECL_STATE_H

/** \addtogroup st1d
 *  @{
 */
//-------------------------------------------------------------------
//-------------------------------------------------------------------
/** \class  SHarmonic1D
 *  \brief  A 1D Spherical Harmonic
 *  
 *   A 1D Spherical Harmonic is composed of an Array2D <double> where the two dimensions are 
 *   the 1D momentum and 1D configuration space.
 *   
*/
class SHarmonic1D {
//-------------------------------------------------------------------    	
private:
    Array2D<complex<double> > *sh;

public:
///     The constructor requires nump, and numx as inputs
    SHarmonic1D(size_t nump, size_t numx);
    SHarmonic1D(const SHarmonic1D& other);
    ~SHarmonic1D();

///     To retrieve the the array that stores the information
    Array2D<complex<double> >& array() const {return (*sh);}
    size_t dim()  const {return (*sh).dim();} //< Total number of values in the array
    size_t nump() const {return (*sh).dim1();} //< Number of momentum cells
    size_t numx() const {return (*sh).dim2();} //< Number of spatial cells

//      Access to the underlying array 
    complex<double> & operator()(size_t i, size_t j) {return (*sh)(i,j);}
    complex<double>   operator()(size_t i, size_t j) const {return (*sh)(i,j);}
    complex<double> & operator()(size_t i) {return (*sh)(i);}             //1D-style
    complex<double>   operator()(size_t i) const {return (*sh)(i);}       //1D-style

    vector<complex<double> > xVec(size_t j) const {return (*sh).d2c(j);}       //1D-style

//      Operators
    SHarmonic1D& operator=(const complex<double> & d);
    SHarmonic1D& operator=(const SHarmonic1D& other);
    SHarmonic1D& operator*=(const complex<double> & d);
    SHarmonic1D& operator*=(const SHarmonic1D& shmulti);
    SHarmonic1D& operator+=(const complex<double> & d);
    SHarmonic1D& operator+=(const SHarmonic1D& shadd);
    SHarmonic1D& operator-=(const complex<double> & d);
    SHarmonic1D& operator-=(const SHarmonic1D& shmin);

//      Other Algebra
    SHarmonic1D& mpaxis(const valarray<complex<double> >& shmulti);
    SHarmonic1D& mxaxis(const valarray<complex<double> >& shmulti);
    SHarmonic1D& Re();

//      Derivatives
    SHarmonic1D& Dp();
    SHarmonic1D& Dx();

//      FilterP
    SHarmonic1D& Filterp(size_t N);


//      Debug
    void checknan();
};
//--------------------------------------------------------------
/** @} */

/** \addtogroup st2d
 *  @{
 */
//-------------------------------------------------------------------
/** \class  SHarmonic2D
 *  \brief  A 2D Spherical Harmonic
 *  
 *   A 2D Spherical Harmonic is composed of an Array3D <complex> where the first dimension is 
 *   the momentum space and the other two are the 2D configuration space.
 *   
*/
class SHarmonic2D {
//--------------------------------------------------------------    	
private:
    Array3D <complex <double> > *sh;

public:
//      Constructors/Destructors
    SHarmonic2D(size_t nump, size_t numx, size_t numy);
    SHarmonic2D(const SHarmonic2D& other);
    ~SHarmonic2D();

//      Basic information
    Array3D < complex <double> >& array() const {return (*sh);}
    size_t dim()  const {return (*sh).dim();}
    size_t nump() const {return (*sh).dim1();}
    size_t numx() const {return (*sh).dim2();}
    size_t numy() const {return (*sh).dim3();}

//      Access to the underlying array 
    complex<double>& operator()(size_t i, size_t j, size_t k) {return (*sh)(i,j,k);}
    complex<double>& operator()(size_t i) {return (*sh)(i);}             //1D-style

    // vector<double> xVec(size_t j) const {return (*sh).d2c(j);}       //1D-style

//      Operators
    SHarmonic2D& operator=(const complex<double>& d);
    SHarmonic2D& operator=(const SHarmonic2D& other);
    SHarmonic2D& operator*=(const complex<double>& d);
    SHarmonic2D& operator*=(const SHarmonic2D& shmulti);
    SHarmonic2D& operator+=(const complex<double>& d);
    SHarmonic2D& operator+=(const SHarmonic2D& shadd);
    SHarmonic2D& operator-=(const complex<double>& d);
    SHarmonic2D& operator-=(const SHarmonic2D& shmin);

//      Other Algebra
    SHarmonic2D& mpaxis(const valarray <complex <double> >& shmulti);
    SHarmonic2D& mxaxis(const valarray <complex <double> >& shmulti);
    SHarmonic2D& myaxis(const valarray <complex <double> >& shmulti);
    SHarmonic2D& mxy_matrix(Array2D <complex <double> >& shmultiM);


//      Derivatives
    // SHarmonic2D& Dp();
    SHarmonic2D& Dx();
    SHarmonic2D& Dy();

//      FilterP
    SHarmonic2D& Filterp(size_t N);
};
//--------------------------------------------------------------    
/** @} */

/** \addtogroup st1d
 *  @{
 */
//-------------------------------------------------------------------
/** \class  Field1D
 *  \brief  A 1D Field
 *  
 *   A 1D Field is composed of a valarray describing the field over 1 dimension.
 *   
*/
class Field1D {
//-------------------------------------------------------------------
private:
    valarray<complex<double> > *fi;

public:
//      Constructors/Destructors
    Field1D(size_t numx);
    Field1D(const Field1D& other);
    ~Field1D();

//      Access to the underlying matrix
    valarray<complex<double> >& array() const {return (*fi);}
    size_t numx() const {return (*fi).size();}
    complex<double> & operator()(size_t i){ return (*fi)[i];}             //1D-style
    complex<double>   operator()(size_t i) const {return (*fi)[i];}

//      Operators
    Field1D& operator=(const complex<double> & d);
    Field1D& operator=(const valarray<complex<double> >& other);
    Field1D& operator=(const Field1D& other);
    Field1D& operator*=(const complex<double> & d);
    Field1D& operator*=(const valarray<complex<double> >& fimulti);
    Field1D& operator*=(const Field1D& fimulti);
    Field1D& operator+=(const complex<double> & d);
    Field1D& operator+=(const Field1D& fiadd);
    Field1D& operator-=(const complex<double> & d);
    Field1D& operator-=(const Field1D& fimin);

//      Other Algebra
    Field1D& Re();

//      Derivatives
    Field1D& Dx();
};
//-------------------------------------------------------------------
/** @} */


/** \addtogroup st2d
 *  @{
 */
//-------------------------------------------------------------------
/** \class  Field2D
 *  \brief  A 2D Field
 *  
 *   A 2D Field is composed of an Array2D <complex> that contains the field over 2 dimensions.
 *   
*/
class Field2D {
//-------------------------------------------------------------------
private:
    Array2D < complex < double > > *fi;

public:
//      Constructors/Destructors
    Field2D(size_t numx, size_t numy);
    Field2D(const Field2D& other);
    ~Field2D();

//      Access to the underlying matrix
    Array2D < complex <double> >& array() const {return (*fi);}
    size_t numx() const {return (*fi).dim1();}
    size_t numy() const {return (*fi).dim2();}
    complex<double>& operator()(size_t i, size_t j){return (*fi)(i,j);} //Fortran-style
    complex<double>& operator()(size_t i){return (*fi)(i);}             //1D-style

//      Operators
    Field2D& operator=(const complex<double>& d);
    // Field2D& operator=(const Array2D< complex <double> >& other);
    Field2D& operator=(const Field2D& other);
    Field2D& operator*=(const complex<double>& d);
    // Field2D& operator*=(const Array2D< complex <double> >& fimulti);
    Field2D& operator*=(const Field2D& fimulti);
    Field2D& operator+=(const complex<double>& d);
    Field2D& operator+=(const Field2D& fiadd);
    Field2D& operator-=(const complex<double>& d);
    Field2D& operator-=(const Field2D& fimin);

//      Derivatives
    Field2D& Dx();
    Field2D& Dy();
};
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
class EMF1D {
//--------------------------------------------------------------
//  Electromagnetic fields decleration
//--------------------------------------------------------------

private:
    vector<Field1D> *fie;

public:
//      Constructors/Destructors
    EMF1D(size_t nx);
    EMF1D(const EMF1D& other);
    ~EMF1D();

//      Access
    size_t dim()  const {return (*fie).size();}
    Field1D& operator()(size_t i) {return (*fie)[i];}
    Field1D  operator()(size_t i) const {return (*fie)[i];}

    Field1D& Ex() {return (*fie)[0];}
    Field1D& Ey() {return (*fie)[1];}
    Field1D& Ez() {return (*fie)[2];}
    Field1D& Bx() {return (*fie)[3];}
    Field1D& By() {return (*fie)[4];}
    Field1D& Bz() {return (*fie)[5];}

//      Operators
    EMF1D& operator=(const complex<double>& d);
    EMF1D& operator=(const Field1D& h);
    EMF1D& operator=(const EMF1D& other);
    EMF1D& operator*=(const complex<double>& d);
    EMF1D& operator*=(const EMF1D& other);
    EMF1D& operator+=(const complex<double>& d);
    EMF1D& operator+=(const EMF1D& other);
    EMF1D& operator-=(const complex<double>& d);
    EMF1D& operator-=(const EMF1D& other);

};
//--------------------------------------------------------------


//--------------------------------------------------------------
/** @} */

/** \addtogroup st2d
 *  @{
 */
//-------------------------------------------------------------------
/** \class  EMF2D
 *  \brief  An EMF2D is the container for the 6 EM fields in a 2D-3P code
 *  
 *   The object is composed of a vector of Field2D. Along with its own member functions,
 *   it also inherits those of Field2D.
 *   
*/
class EMF2D {
//-------------------------------------------------------------------
private:
    vector<Field2D> *fie;

public:
//      Constructors/Destructors
    EMF2D(size_t nx, size_t ny);
    EMF2D(const EMF2D& other);
    ~EMF2D();

//      Access
    size_t dim()  const {return (*fie).size();}
    Field2D& operator()(size_t i) {return (*fie)[i];}
    Field2D  operator()(size_t i) const {return (*fie)[i];}

    Field2D& Ex() {return (*fie)[0];}
    Field2D& Ey() {return (*fie)[1];}
    Field2D& Ez() {return (*fie)[2];}
    Field2D& Bx() {return (*fie)[3];}
    Field2D& By() {return (*fie)[4];}
    Field2D& Bz() {return (*fie)[5];}

//      Operators
    EMF2D& operator=(const complex<double>& d);
    EMF2D& operator=(const Field2D& h);
    EMF2D& operator=(const EMF2D& other);
    EMF2D& operator*=(const complex<double>& d);
    EMF2D& operator*=(const EMF2D& other);
    EMF2D& operator+=(const complex<double>& d);
    EMF2D& operator+=(const EMF2D& other);
    EMF2D& operator-=(const complex<double>& d);
    EMF2D& operator-=(const EMF2D& other);

};
//--------------------------------------------------------------
/** @} */

/** \addtogroup st1d
 *  @{
 */
//-------------------------------------------------------------------
/** \class  DistFunc1D
 *  \brief  The 1D distribution function is the container for all SHarmonic1D per species.
 *  
 *   The object is composed of a vector of SHarmonic1D, the length of which is specified by l0 and m0
 *   Along with its own member functions, it also inherits those of SHarmonic1D. Since each species requires
 *   a DistFunc object, this class also contains information about the l_max, num_p, p_max, charge, and mass.
 *   
*/
class DistFunc1D {
//-------------------------------------------------------------------    	
private:

    vector<SHarmonic1D> *df;
    size_t lmax, mmax, sz;
    double pmx, charge, ma;

    Array2D<int> ind;
    valarray<int> filter_ceiling;

public:

//      Constructors/Destructors
    DistFunc1D(size_t l, size_t m, size_t np, double pma, size_t nx, double q, double _ma);
    DistFunc1D(const DistFunc1D& other);
    ~DistFunc1D();

//      Basic info
    size_t dim()  const {return sz;}
    size_t l0()   const {return lmax;  }
    size_t m0()   const {return mmax;  }
    double pmax() const {return pmx;   }
    double q()    const {return charge;}
    double mass() const {return ma;    }

    Array2D<int> indx() const {return ind;}

    valarray<double> getdensity();
    valarray<double> getcurrent(size_t dir);
    valarray<double> getcurrent(size_t dir) const;
    Array2D<double>  getcurrent() const;
    valarray<double> getpressure();

//      Access
    SHarmonic1D& operator()(int i);                   // Returns the address of a harmonic, NULL if out-of-bounds
    SHarmonic1D& operator()(int i) const;
    SHarmonic1D& operator()(size_t l, size_t m);
    SHarmonic1D& operator()(size_t l, size_t m) const;                  // Returns the address of a harmonic, NULL if out-of-bounds

//      Operators
    DistFunc1D& operator=(const complex<double> & d);
    DistFunc1D& operator=(const SHarmonic1D& h);
    DistFunc1D& operator=(const DistFunc1D& other);
    DistFunc1D& operator*=(const complex<double> & d);
    DistFunc1D& operator*=(const DistFunc1D& other);
    DistFunc1D& operator+=(const complex<double> & d);
    DistFunc1D& operator+=(const DistFunc1D& other);
    DistFunc1D& operator-=(const complex<double> & d);
    DistFunc1D& operator-=(const DistFunc1D& other);

//      Filter
    DistFunc1D& Filterp();

//      Debug
    void checknan();
};
//--------------------------------------------------------------------------------
/** @} */


/** \addtogroup st2d
 *  @{
 */
//--------------------------------------------------------------------------------
/** \class  DistFunc2D
 *  \brief  The 2D distribution function is the container for all SHarmonic2D per species.
 *   
 *   The object is composed of a vector of SHarmonic1D, the length of which is specified by l0 and m0
 *   Along with its own member functions, it also inherits those of SHarmonic1D. Since each species requires
 *   a DistFunc object, this class also contains information about the l_max, num_p, p_max, charge, and mass.
 *   
*/
class DistFunc2D {
//--------------------------------------------------------------------------------    	
private:
    vector<SHarmonic2D> *df;
    size_t lmax, mmax, sz;
    double charge, ma;

//  Index for "trapezoid matrix"
    Array2D<int> ind;
    Array2D<int> indx() const {return ind;}

public:

//      Constructors/Destructors
    DistFunc2D(size_t l, size_t m, size_t np, size_t nx, size_t ny, double q, double _ma);
    DistFunc2D(const DistFunc2D& other);
    ~DistFunc2D();

//      Basic info
    size_t dim()  const {return sz;}
    size_t ldim()  const {return lmax+1;}
    size_t mdim()  const {return mmax+1;}
    size_t l0()   const {return lmax;  }
    size_t m0()   const {return mmax;  }
    double q()    const {return charge;}
    double mass() const {return ma;}

//      Access
    SHarmonic2D* operator()(size_t i); //{return &(*df)[i];}
    SHarmonic2D* operator()(size_t i) const; // const {return &(*df)[i];}
    SHarmonic2D* operator()(size_t l, size_t m) {return &(*df)[ind(l,m)];}
    SHarmonic2D* operator()(size_t l, size_t m) const {return &(*df)[ind(l,m)];}                   // Returns the address of a harmonic, NULL if out-of-bounds

//      Operators
    DistFunc2D& operator=(const complex <double>& d);
    DistFunc2D& operator=(const SHarmonic2D& h);
    DistFunc2D& operator=(const DistFunc2D& other);
    DistFunc2D& operator*=(const complex <double>& d);
    DistFunc2D& operator*=(const DistFunc2D& other);
    DistFunc2D& operator+=(const complex <double>& d);
    DistFunc2D& operator+=(const DistFunc2D& other);
    DistFunc2D& operator-=(const complex <double>& d);
    DistFunc2D& operator-=(const DistFunc2D& other);

//      Filter
    // DistFunc2D& Filterp();
};
//--------------------------------------------------------------    
/** @} */


//-------------------------------------------------------------------
/** \class  Hydro1D
 *  \brief  A Collection of relevant 1D Hydrodynamic Quantities
 *  
 *   A 1D Field is composed of a valarray describing the field over 1 dimension.
 *   
*/
class Hydro1D {
//-------------------------------------------------------------------
private:
    valarray<double>  *hn, *hvx, *hvy, *hvz, *ht, *hz;

    double hydromass, hydrocharge;

public:
//      Constructors/Destructors
    Hydro1D(size_t numx, double _mass, double _charge);
    Hydro1D(const Hydro1D& other);
    ~Hydro1D();

//      Access to the underlying matrix
    // valarray<complex<double> >& array() const {return (*fi);}
    size_t numx() const {return (*hn).size();}
    double mass() const {return hydromass;}
    double charge() const {return hydrocharge;}

    // void Dx_vel(valarray<complex<double>>& vin);
    double & density(size_t i){ return (*hn)[i];}             //1D-style
    double   density(size_t i) const {return (*hn)[i];}

    double & vx(size_t i){ return (*hvx)[i];}
    double   vx(size_t i) const {return (*hvx)[i];}

    double & vy(size_t i){ return (*hvy)[i];}
    double   vy(size_t i) const {return (*hvy)[i];}

    double & vz(size_t i){ return (*hvz)[i];}
    double   vz(size_t i) const {return (*hvz)[i];}

    double & temperature(size_t i){ return (*ht)[i];}             //1D-style
    double   temperature(size_t i) const {return (*ht)[i];}

    double & Z(size_t i){ return (*hz)[i];}             //1D-style
    double   Z(size_t i) const {return (*hz)[i];}

    valarray<double >& densityarray()       const {return (*hn);}
    valarray<double >& vxarray()            const {return (*hvx);}
    valarray<double >& vyarray()            const {return (*hvy);}
    valarray<double >& vzarray()            const {return (*hvz);}
    valarray<double >& temperaturearray()   const {return (*ht);}
    valarray<double >& Zarray()   const {return (*hz);}

//      Operators
    Hydro1D& operator=(const double & d);
    Hydro1D& operator=(const valarray<double >& other);
    Hydro1D& operator=(const Hydro1D& other);

    Hydro1D& operator*=(const double & d);
    Hydro1D& operator*=(const valarray<double >& other);
    Hydro1D& operator*=(const Hydro1D& other);

    Hydro1D& operator+=(const double & d);
    Hydro1D& operator+=(const valarray<double >& other);
    Hydro1D& operator+=(const Hydro1D& other);

    Hydro1D& operator-=(const double & d);
    Hydro1D& operator-=(const valarray<double >& other);
    Hydro1D& operator-=(const Hydro1D& other);

};

/** \addtogroup st1d
 *  @{
 */
//--------------------------------------------------------------
//  Collection of fields and distribution functions decleration
//--------------------------------------------------------------    
class State1D {
private:
    vector<DistFunc1D> *sp;
    EMF1D *flds;
    Hydro1D *hydro;

    // valarray<complex<double>>       *_hydrovelocity, *_hydrotemperature;
    // valarray<complex<double>>       *_kineticspeciespressure, *_magneticpressure;

    size_t ns;


public:
//      Constructors/Destructors
    State1D(size_t nx, vector<size_t> l0, vector<size_t> m0, vector<size_t> np, vector<double> pmax, vector<double> q, vector<double> ma, double hydromass, double hydrocharge);
    State1D(const State1D& other);
    ~State1D();

//      Basic information
    size_t Species() const {return ns;}
    size_t Fields() const  {return 6;}
    void   checknan();

//      Access to underlying structures
//          Distributions
    DistFunc1D&  DF(size_t s)         {return (*sp)[s];}
    DistFunc1D&  DF(size_t s) const   {return (*sp)[s];}
//          Individual Harmonics
    SHarmonic1D& SH(size_t s, size_t lh, size_t mh)        {return ((*sp)[s])(lh,mh);} // Reference to spherical harmonic
    SHarmonic1D& SH(size_t s, size_t lh, size_t mh)  const {return ((*sp)[s])(lh,mh);}
//         SHarmonic1D* SHp(size_t s, size_t lh, size_t mh)       {return   (sp[s])(lh,mh); } // Pointer to spherical harmonic
//         SHarmonic1D* SHp(size_t s, size_t lh, size_t mh) const {return   (sp[s])(lh,mh); }
//          Fields
    EMF1D& EMF() const {return (*flds);}
    Field1D& FLD(size_t ip) const {return (*flds)(ip);}
//          Hydro
    Hydro1D&  HYDRO()         {return (*hydro);}
    Hydro1D&  HYDRO() const   {return (*hydro);}


//      Copy assignment Operator
    State1D& operator=(const State1D& other);
    State1D& operator=(const complex<double> & d);
    State1D& operator*=(const State1D& other);
    State1D& operator*=(const complex<double> & d);
    State1D& operator+=(const State1D& other);
    State1D& operator+=(const complex<double> & d);
    State1D& operator-=(const State1D& other);
    State1D& operator-=(const complex<double> & d);

};
//--------------------------------------------------------------
/** @} */

/** \addtogroup st2d
 *  @{
 */
//--------------------------------------------------------------
//  Collection of fields and distribution functions decleration
//--------------------------------------------------------------    
class State2D {
private:
    vector<DistFunc2D> *sp;
    EMF2D *flds;

    size_t ns;

public:
//      Constructors/Destructors
    State2D(size_t nx, size_t ny, vector<size_t> l0, vector<size_t> m0, vector<size_t> np, vector<double> q, vector<double> ma);
    State2D(const State2D& other);
    ~State2D();

//      Basic information
    size_t Species() const {return ns;}
    size_t Fields() const  {return 6;}

//      Access to underlying structures
//          Distributions
    DistFunc2D&  DF(size_t s)         {return (*sp)[s];}
    DistFunc2D&  DF(size_t s) const   {return (*sp)[s];}
//          Individual Harmonics
    SHarmonic2D& SH(size_t s, size_t lh, size_t mh)        {return *(((*sp)[s])(lh,mh));} // Reference to spherical harmonic
    SHarmonic2D& SH(size_t s, size_t lh, size_t mh)  const {return *(((*sp)[s])(lh,mh));}
    SHarmonic2D* SHp(size_t s, size_t lh, size_t mh)       {return   ((*sp)[s])(lh,mh); } // Pointer to spherical harmonic
    SHarmonic2D* SHp(size_t s, size_t lh, size_t mh) const {return   ((*sp)[s])(lh,mh); }
//          Fields
    EMF2D& EMF() const {return (*flds);}
    Field2D& FLD(size_t ip) const {return (*flds)(ip);}

//      Copy assignment Operator
    State2D& operator=(const State2D& other);
    State2D& operator=(const complex<double>& d);
    State2D& operator*=(const State2D& other);
    State2D& operator*=(const complex<double>& d);
    State2D& operator+=(const State2D& other);
    State2D& operator+=(const complex<double>& d);
    State2D& operator-=(const State2D& other);
    State2D& operator-=(const complex<double>& d);
};
// --------------------------------------------------------------


#endif
