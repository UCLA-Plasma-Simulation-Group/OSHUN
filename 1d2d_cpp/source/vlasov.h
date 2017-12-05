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
class Spatial_Advection {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Spatial_Advection(size_t Nl, size_t Nm,
                        valarray<double> dp,
                        double xmin, double xmax, size_t Nx,
                        double ymin, double ymax, size_t Ny);
//          Advance
    void operator()(const DistFunc1D& Din, DistFunc1D& Dh);
    void operator()(const DistFunc2D& Din, DistFunc2D& Dh);
    void es1d(const DistFunc1D& Din, DistFunc1D& Dh);
    void f1only(const DistFunc1D& Din, DistFunc1D& Dh);
    void f1only(const DistFunc2D& Din, DistFunc2D& Dh);

private:
    Array2D< complex<double> >      A1, A2, C2, C4;

    valarray< complex<double> >     B1, B2, C1, C3;
    valarray< complex<double> >  	vr;
    
    valarray<size_t>                f_start, f_end;
    valarray<size_t>                dist_il, dist_im;
    valarray<size_t>                nwsediag_il, nwsediag_im;
    valarray<size_t>                neswdiag_il, neswdiag_im;


    complex<double>                 A00, A10, A20;


    
};
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Electric field
class Electric_Field {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Electric_Field(size_t Nl, size_t Nm,                      
                    valarray<double> dp);
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



    void operator()(const DistFunc2D& Din,
                    const Field2D& FEx, const Field2D& FEy, const Field2D& FEz,
                    DistFunc2D& Dh);
    void f1only(const DistFunc2D& Din,
                const Field2D& FEx, const Field2D& FEy, const Field2D& FEz,
                DistFunc2D& Dh);
    void Implicit_Ex(const DistFunc2D& Din, const Field2D& FEx, DistFunc2D& Dh);
    void Implicit_Ey(const DistFunc2D& Din, const Field2D& FEy, DistFunc2D& Dh);
    void Implicit_Ez(const DistFunc2D& Din, const Field2D& FEz, DistFunc2D& Dh);
    void Implicit_Ex_f1only(const DistFunc2D& Din, const Field2D& FEx, DistFunc2D& Dh);
    void Implicit_Ey_f1only(const DistFunc2D& Din, const Field2D& FEy, DistFunc2D& Dh);
    void Implicit_Ez_f1only(const DistFunc2D& Din, const Field2D& FEz, DistFunc2D& Dh);

private:
    // void MakeG00(SHarmonic1D& f);
    void MakeG00(const SHarmonic1D& f, SHarmonic1D& G);
    void MakeG00(const SHarmonic2D& f, SHarmonic2D& G);
    // void MakeGH( SHarmonic1D& f, size_t l);
    void MakeGH(const SHarmonic1D& f, SHarmonic1D& G, SHarmonic1D& H, size_t l);   // OMP version
    void MakeGH(const SHarmonic2D& f, SHarmonic2D& G, SHarmonic2D& H, size_t l);   // OMP version

    

    complex<double>                 A100, C100, A210, B211, C311, A310;

    Array2D< complex<double> >      A1, A2;
    valarray< complex<double> >     B1, B2;
    valarray< complex<double> >     C1, C3;
    Array2D< complex<double> >      C2, C4;
    valarray<complex<double> >      Hp0;


    valarray< complex<double> >     pr, invdp, invpr;

    valarray<size_t>                f_start, f_end;
    valarray<size_t>                dist_il, dist_im;
    valarray<size_t>                nwsediag_il, nwsediag_im;
    valarray<size_t>                neswdiag_il, neswdiag_im;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Magnetic field
class Magnetic_Field {
//--------------------------------------------------------------    	
public:
//      Constructors/Destructors
    Magnetic_Field(size_t Nl, size_t Nm,valarray<double> dp);
//          Advance
    void operator()(const DistFunc1D& Din,
                    const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                    DistFunc1D& Dh);
    void operator()(const DistFunc2D& Din,
                    const Field2D& FBx, const Field2D& FBy, const Field2D& FBz,
                    DistFunc2D& Dh);

    void f1only(const DistFunc1D& Din,
                const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
                DistFunc1D& Dh);

    void f1only(const DistFunc2D& Din,
                const Field2D& FBx, const Field2D& FBy, const Field2D& FBz,
                DistFunc2D& Dh);

    // void implicit(DistFunc1D& Din,
    //               const Field1D& FBx, const Field1D& FBy, const Field1D& FBz,
    //               double dt);

private:
    valarray< complex<double> >		A1, B1;
    Array2D< complex<double> > 		A2;
    complex<double> 				A3;

    valarray<size_t>                f_start, f_end;
    valarray<size_t>                dist_il, dist_im;
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//--------------------------------------------------------------
//  Current
class Current {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Current();
//          Advance
    void operator()(const DistFunc1D& Din,
                    Field1D& FExh, Field1D& FEyh, Field1D& FEzh);

    void operator()(const DistFunc2D& Din,
                    Field2D& FExh, Field2D& FEyh, Field2D& FEzh);

    void es1d(const DistFunc1D& Din,
                    Field1D& FExh);

private:
    
};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update B with E term from Faraday's Law
class Faraday {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Faraday(double xmin, double xmax, size_t Nx,
        double ymin, double ymax, size_t Ny);
//          Advance
    void operator()(EMF1D& EMFin, EMF1D& EMFh);
    void operator()(EMF2D& EMFin, EMF2D& EMFh);

private:

    complex<double>		idx,idy;

};
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Update E with B from Ampere's Law
class Ampere {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Ampere(double xmin, double xmax, size_t Nx,
            double ymin, double ymax, size_t Ny);
//          Advance
    void operator()(EMF1D& EMFin, EMF1D& EMFh);
    void operator()(EMF2D& EMFin, EMF2D& EMFh);

private:
    
    complex<double>		idx, idy;

};
//--------------------------------------------------------------

#endif
