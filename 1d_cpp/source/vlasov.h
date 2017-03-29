/*!\brief  Vlasov Equation - Declarations
 * \author PICKSC
 * \file   vlasov.h
 *
 * Includes declarations for spatial advection, electric field advection, bfield and current
 *
 */

#ifndef DECL_VLASOVMAXWELL_H
#define DECL_VLASOVMAXWELL_H

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
class Spatial_Advection_1D {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Spatial_Advection_1D(size_t Nl, size_t Nm,
                         double pmin, double pmax, size_t Np,
                         double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din, DistFunc1D& Dh);
    void es1d(const DistFunc1D& Din, DistFunc1D& Dh);
    void f1only(const DistFunc1D& Din, DistFunc1D& Dh);

private:
    SHarmonic1D                   	fd1, fd2;
    complex<double>                 A00, A10, A20;
    Array2D< complex<double> >  	A1, A2;
    valarray< complex<double> >  	vr;
};
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Electric field
class Electric_Field_1D {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Electric_Field_1D(size_t Nl, size_t Nm,
                      double pmin, double pmax, size_t Np,
                      double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din,
                    const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                    DistFunc1D& Dh);
    void es1d(const DistFunc1D& Din,
                    const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                    DistFunc1D& Dh);
    void f1only(const DistFunc1D& Din,
                const Field1D& FEx, const Field1D& FEy, const Field1D& FEz,
                DistFunc1D& Dh);
    void Implicit_Ex(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh);
    void Implicit_Ey(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh);
    void Implicit_Ez(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh);

    void Implicit_Ex_f1only(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh);
    void Implicit_Ey_f1only(const DistFunc1D& Din, const Field1D& FEy, DistFunc1D& Dh);
    void Implicit_Ez_f1only(const DistFunc1D& Din, const Field1D& FEz, DistFunc1D& Dh);

private:
    void MakeG00(SHarmonic1D& f);
    void MakeGH( SHarmonic1D& f, size_t l);

    SHarmonic1D H, G, TMP;

    complex<double>                 A100, C100, A210, B211, C311, A310;

    Array2D< complex<double> >   A1, A2;
    valarray< complex<double> >  B1, B2;
    valarray< complex<double> >  C1, C3;
    Array2D< complex<double> >   C2, C4;
    valarray< complex<double> >  pr, invpr, Hp0;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Magnetic field
class Magnetic_Field_1D {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Magnetic_Field_1D(size_t Nl, size_t Nm,
                      double pmin, double pmax, size_t Np,
                      double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(const DistFunc1D& Din,
                    const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                    DistFunc1D& Dh);
    void f1only(const DistFunc1D& Din,
                const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                DistFunc1D& Dh);
    void implicit(DistFunc1D& Din,
                  const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                  double dt);

private:

    SHarmonic1D FLM;

    valarray< complex<double> >		A1, B1;
    Array2D< complex<double> > 		A2;
    complex<double> 				A3;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//--------------------------------------------------------------
//  Current
class Current_1D {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Current_1D( double pmin, double pmax, size_t Np,
                size_t Nx );
//          Advance
    void operator()(const DistFunc1D& Din,
                    Field1D& FExh, Field1D& FEyh, Field1D& FEzh);
    void es1d(const DistFunc1D& Din,
                    Field1D& FExh);

private:
    Field1D Jx, Jy, Jz;
    double small;
//            valarray< complex<double> >  pr, invg;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update B with E term from Faraday's Law
class Faraday_1D {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Faraday_1D(double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(EMF1D& EMFin, EMF1D& EMFh);

private:
    Field1D             tmpE;
    complex<double>		idx;
    size_t              numx;

};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update E with B from Ampere's Law
class Ampere_1D {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Ampere_1D(double xmin, double xmax, size_t Nx);
//          Advance
    void operator()(EMF1D& EMFin, EMF1D& EMFh);

private:
    Field1D             tmpB;
    complex<double>		idx;
    size_t              numx;

};
//--------------------------------------------------------------


#endif
