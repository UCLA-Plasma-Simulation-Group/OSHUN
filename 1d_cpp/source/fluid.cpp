/*! \brief  Fluid Equation - Definitions
 * \author PICKSC
 * \date   October 10, 2016
 * \file   fluid.cpp
 *
 * Includes momentum equation, advection in vlasov, and and electric field update routines
 * 
 */
//--------------------------------------------------------------
//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>

    #include <math.h>
    #include <map>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declerations
    #include "state.h"
    #include "formulary.h"
    #include "input.h"
    #include "nmethods.h"
    // #include "vlasov.h"
    #include "fluid.h"

//*******************************************************************************************************************************************************
/// Constructor for hydro equations
//-------------------------------------------------------------------------------------------------------------------------------------------------------
Hydro_Functor::Hydro_Functor(double xmin, double xmax, size_t numx):
    idx(numx/(xmax-xmin)), szx(numx), FEQ(xmin,xmax,numx) {}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------
void Hydro_Functor::operator()(const State1D& Yin, State1D& Yslope){
    Yslope = 0.0;

    valarray<double> electronpressure(0.0,szx), electrondensity(0.0, szx);
    Array2D<double> electroncurrent(3,szx);

    electrondensity = Yin.DF(0).getdensity();
    electronpressure = Yin.DF(0).getpressure();

    FEQ.density(         Yin.HYDRO(), Yslope.HYDRO());
    FEQ.chargefraction(  Yin.HYDRO(), Yslope.HYDRO());
    FEQ.velocity(        Yin        , Yslope.HYDRO(), electrondensity, electronpressure, electroncurrent);
    FEQ.updateE(         Yin.HYDRO(), Yslope.EMF());
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------
void Hydro_Functor::operator()(const State1D& Yin, State1D& Yslope, size_t dir){}
void Hydro_Functor::operator()(const State1D& Y1in, const State1D& Y2in, State1D& Yslope){}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------

//******************************************************************************************
//------------------------------------------------------------------------------------------
//  Fluid Equation Class Constructor

Fluid_Equation_1D::Fluid_Equation_1D(double xmin, double xmax, size_t Nx)
  : idx(Nx/(xmax-xmin)), szx(Nx), dummy(0.0,Nx)
    {
        Nbc = Input::List().BoundaryCells;

    }

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
void Fluid_Equation_1D::density(const Hydro1D& Hin, Hydro1D& Hslope){
    valarray<double> d_nvx(Hin.densityarray());
    for (size_t ix(0);ix<szx;++ix)
        d_nvx[ix] = d_nvx[ix] * Hin.vx(ix);

    d_nvx = df_4thorder(d_nvx);

    Hslope.density(0)       -= idx * d_nvx[0];

    for (size_t ix(1);ix<szx-1;++ix)
    {
        Hslope.density(ix)  -= idx * d_nvx[ix];
    }

    Hslope.density(szx - 1) -= idx * d_nvx[szx-1];

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void Fluid_Equation_1D::chargefraction(const Hydro1D& Hin, Hydro1D& Hslope){
    valarray<double> d_nZx(Hin.densityarray());
    for (size_t ix(0);ix<szx;++ix)
        d_nZx[ix] = d_nZx[ix] * Hin.Z(ix);

    d_nZx = df_4thorder(d_nZx);    

    Hslope.Z(0)       -= idx * d_nZx[0];

    for (size_t ix(1);ix<szx-1;++ix)
    {
        Hslope.Z(ix)  -= idx * d_nZx[ix];
    }

    Hslope.Z(szx - 1) -= idx *  d_nZx[szx-1];

}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
void Fluid_Equation_1D::velocity(const State1D& Yin , Hydro1D& Hslope, 
    valarray<double>& electrondensity, valarray<double>& electronpressure, Array2D<double>& electroncurrent){

    double Zeovermi = Yin.HYDRO().charge()/Yin.HYDRO().mass();
    valarray<double> eventuallynetPovern(Yin.HYDRO().densityarray()), d_vx(Yin.HYDRO().vxarray()), d_vy(Yin.HYDRO().vxarray()), d_vz(Yin.HYDRO().vxarray());
    valarray<double> eventuallychargeseparation(electronpressure);
    Array2D<double> netcurrent(3,szx);


    for (size_t ix(0);ix<szx;++ix){
        eventuallynetPovern[ix] = eventuallynetPovern[ix] * Yin.HYDRO().temperature(ix);
    }

    eventuallynetPovern = df_4thorder(eventuallynetPovern);
    eventuallychargeseparation = df_4thorder(eventuallychargeseparation);

    for (size_t ix(0);ix<szx;++ix){

        // Net pressure
        eventuallynetPovern[ix] = eventuallynetPovern[ix] / Yin.HYDRO().density(ix) + eventuallychargeseparation[ix] / electrondensity[ix];

        // Net current
        netcurrent(0,ix) = Yin.HYDRO().Z(ix)*Yin.HYDRO().charge()*Yin.HYDRO().density(ix)*Yin.HYDRO().vx(ix) - electroncurrent(0,ix);
        netcurrent(1,ix) = Yin.HYDRO().Z(ix)*Yin.HYDRO().charge()*Yin.HYDRO().density(ix)*Yin.HYDRO().vy(ix) - electroncurrent(1,ix);
        netcurrent(2,ix) = Yin.HYDRO().Z(ix)*Yin.HYDRO().charge()*Yin.HYDRO().density(ix)*Yin.HYDRO().vz(ix) - electroncurrent(2,ix);

        // Net charge separation
        eventuallychargeseparation[ix] = Yin.HYDRO().Z(ix)*Yin.HYDRO().charge()*Yin.HYDRO().density(ix)-electrondensity[ix];
    }
    
    d_vx = df_4thorder(d_vx);
    d_vy = df_4thorder(d_vy);
    d_vz = df_4thorder(d_vz);


    Hslope.vx(0)  += - idx * eventuallynetPovern[0] / Yin.HYDRO().mass()
                                    +   Zeovermi * eventuallychargeseparation[0] * (Yin.EMF().Ex()(0).real()
                                    +   netcurrent(1,0) * Yin.EMF().Bz()(0).real()
                                    -   netcurrent(2,0) * Yin.EMF().By()(0).real()); //+ Rie/ni;

    Hslope.vx(0) -= netcurrent(0,0) * (idx * d_vx[0]);

    for (size_t ix(1);ix<szx-1;++ix)
    {

        Hslope.vx(ix)  += -0.5 * idx * eventuallynetPovern[ix] / Yin.HYDRO().mass()
                                    +   Zeovermi * eventuallychargeseparation[ix] * (Yin.EMF().Ex()(ix).real()
                                    +   netcurrent(1,ix) * Yin.EMF().Bz()(ix).real()
                                    -   netcurrent(2,ix) * Yin.EMF().By()(ix).real()); //+ Rie/ni;
        Hslope.vx(ix) -= netcurrent(0,ix) * (idx * d_vx[ix]);
    }

    Hslope.vx(szx - 1) += -idx * eventuallynetPovern[szx - 1] / Yin.HYDRO().mass()
                                    +   Zeovermi * eventuallychargeseparation[szx - 1] * (Yin.EMF().Ex()(szx - 1).real()
                                    +   netcurrent(1,szx - 1) * Yin.EMF().Bz()(szx - 1).real()
                                    -   netcurrent(2,szx - 1) * Yin.EMF().By()(szx - 1).real()); //+ Rie/ni
    Hslope.vx(szx - 1) -= netcurrent(0,szx - 1) * (idx * d_vx[szx - 1]);

    /// ----------------------------------------------------------------------
    /// Momentum equation - y direction
    /// ----------------------------------------------------------------------
    Hslope.vy(0) +=  Zeovermi *
                            (   netcurrent(2,0) * Yin.EMF().Bx()(0).real()
                            -   netcurrent(0,0) * Yin.EMF().Bz()(0).real()
                            +   eventuallychargeseparation[0] * Yin.EMF().Ey()(0).real());

    for (size_t ix(1);ix<szx-1;++ix)
    {
        Hslope.vy(ix) +=  Zeovermi *
                            (   netcurrent(2,ix) * Yin.EMF().Bx()(ix).real()
                            -   netcurrent(0,ix) * Yin.EMF().Bz()(ix).real()
                            +   eventuallychargeseparation[ix] * Yin.EMF().Ey()(ix).real());
    }

    Hslope.vy(szx - 1) +=  Zeovermi *
                            (   netcurrent(2,szx - 1) * Yin.EMF().Bx()(szx - 1).real()
                            -   netcurrent(0,szx - 1) * Yin.EMF().Bz()(szx - 1).real()
                            +   eventuallychargeseparation[szx - 1] * Yin.EMF().Ey()(szx - 1).real());


    /// ----------------------------------------------------------------------
    /// Momentum equation - z direction
    /// ----------------------------------------------------------------------
    Hslope.vz(0) +=  Zeovermi *
                            (   netcurrent(0,0) * Yin.EMF().By()(0).real()
                            -   netcurrent(1,0) * Yin.EMF().Bx()(0).real()
                            +   eventuallychargeseparation[0] * Yin.EMF().Ez()(0).real());

    for (size_t ix(1);ix<szx-1;++ix)
    {
        Hslope.vz(ix) +=  Zeovermi *
                            (   netcurrent(0,ix) * Yin.EMF().By()(ix).real()
                            -   netcurrent(1,ix) * Yin.EMF().Bx()(ix).real()
                            +   eventuallychargeseparation[ix] * Yin.EMF().Ez()(ix).real());
    }

    Hslope.vz(szx - 1) +=  Zeovermi *
                            (   netcurrent(0,szx - 1) * Yin.EMF().By()(szx - 1).real()
                            -   netcurrent(1,szx - 1) * Yin.EMF().Bx()(szx - 1).real()
                            +   eventuallychargeseparation[szx - 1] * Yin.EMF().Ez()(szx - 1).real());
}
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
void Fluid_Equation_1D::updateE(Hydro1D& HYDRO,
                                EMF1D& EMF){

    



    for (size_t ix(0);ix<szx;++ix){
        // EMF.Ex()(ix) -= 0.0;//dCdt[ix]; //+0.5*dC2xdx[ix]+(hydrocharge*electrondensity[ix]-hydrodensity[ix]); // Needs JxB in 2D
        // EMF.Ey()(ix) += -1.0*HYDRO.Z(ix)*HYDRO.charge()*HYDRO.vx(ix)*EMF.Bz()(ix);//+jx[ix]*Bz(ix)-jz[ix]*Bx(ix);
        // EMF.Ez()(ix) += HYDRO.Z(ix)*HYDRO.charge()*HYDRO.vx(ix)*EMF.By()(ix);//-jx[ix]*By(ix)-jy[ix]*Bx(ix);

        // EMF.Ex()(ix) += HYDRO.Z(ix)*HYDRO.charge()*HYDRO.vx(ix)*HYDRO.density(ix);
        // EMF.Ey()(ix) += HYDRO.Z(ix)*HYDRO.charge()*HYDRO.vy(ix)*HYDRO.density(ix);
        // EMF.Ez()(ix) += HYDRO.Z(ix)*HYDRO.charge()*HYDRO.vz(ix)*HYDRO.density(ix);
    }
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

// void Fluid_Equation_1D:: calcquantities(const State1D& Yin){

//    Yin.HYDRO().kineticspeciespressurearray()=0.0;

//    for (size_t s(0);s<Yin.Species();++s){

//        Yin.HYDRO().kineticspeciespressurearray() += Yin.DF(s).getpressure();
//        // Yin.HYDRO().jxarray() += (Yin.DF(s).getcurrent(0));
//        // Yin.HYDRO().jyarray() += (Yin.DF(s).getcurrent(1));
//        // Yin.HYDRO().jzarray() += (Yin.DF(s).getcurrent(2));
//        // std::cout << "kinetic pressure[" << ix << "]=" << Yin.HYDRO().kpressure(ix) << "\n";
//        // }
//    }

//    Yin.HYDRO().electrondensityarray() =  Yin.DF(0).getdensity();

// }

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

//******************************************************************************************
//------------------------------------------------------------------------------------------
Hydro_Advection_1D::Hydro_Advection_1D(size_t Nl, size_t Nm,
                         double pmin, double pmax, size_t Np,
                         double xmin, double xmax, size_t Nx)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
   : XX1(Nl+1,Nm+1), XX2(Nl+1,Nm+1), XX3(Nl+1,Nm+1), XX4(Nl+1,Nm+1),
     A1(Nl+1,Nm+1), A2(Nl+1,Nm+1), fd1(Nl+1,Nm+1), fd2(Nl+1,Nm+1),
     Hp0(Nl+1), H(Np,Nx), G(Np,Nx), TMP(Np,Nx),
     pr(Algorithms::MakeCAxis(static_cast<complex<double> >(pmin),
                             static_cast<complex<double> >(pmax),
                             Np)),
     invpr(pr), szx(Nx)
     {
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -
         complex<double> lc, mc;

//       Inverted momentum axis
         for (size_t i(0); i < pr.size(); ++i) {
             invpr[i] = 1.0/pr[i];
         }
         idp = (-1.0)/ (2.0*(pmax-pmin)/double(Np)); // -1/(2dp)
         idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx)); // -1/(2dx)

         Nbc = Input::List().BoundaryCells;

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Calculate the "A1, A2" parameters
         for (size_t l(0); l < Nl+1; ++l){
             for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 XX1(l,m) = (lc+mc    ) * (lc-mc    ) / (2.0*lc-1.0) / (2.0*lc+1.0) ;
                 XX2(l,m) = (lc-mc+1.0) * (lc+mc+1.0) / (2.0*lc+3.0) / (2.0*lc+1.0) ;

                 XX3(l,m) = (lc-mc+1.0) * (lc-mc+2.0) / (2.0*lc+3.0) / (2.0*lc+1.0) ;
                 XX4(l,m) = (lc+mc-1.0) * (lc+mc    ) / (2.0*lc-1.0) / (2.0*lc+1.0) ; //idx *(-1.0) * (lc-mc+1.0) / (2.0*lc+1.0);
             }
         }
         A2(0,0) = 1.0;
         //       -2*Dp*H at the 0 momentum cell
         Hp0[0] = 1.0 / pr[0];
         for (size_t l(1); l < Nl+1; ++l) {
             double ld(l);
             Hp0[l] = Hp0[l-1] * (pr[0]/pr[1]) * (2.0*ld+1.0)/(2.0*ld-1.0);
         }
         Hp0 *= (-2.0)*(pr[1]-pr[0]);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//   Advection in x
void Hydro_Advection_1D::operator()(const DistFunc1D& Din, const Hydro1D& hydro, DistFunc1D& Dh) {
//--------------------------------------------------------------

        valarray<complex<double> > vt(pr); vt *= 1.0/Din.mass();
        valarray<complex<double> > tempv(vt);


        valarray<complex<double> > Ux(0.0,(hydro.vxarray()).size());
        valarray<complex<double> > dUxdx(0.0,(hydro.vxarray()).size());


        Ux[0] = hydro.vx(0);
        dUxdx[0] = idx*(hydro.vx(1)-hydro.vx(0));

        Ux[szx-1] = hydro.vx(szx-1);
        dUxdx[szx-1] = idx*(hydro.vx(szx-1)-hydro.vx(szx-2));
        for (size_t i(1);i<szx-1;++i)
        {
            Ux[i] = hydro.vx(i);
            // dUxdx[i] = idx/12.0*(-hydro.vx(i+2)+8.0*(hydro.vx(i+1)-hydro.vx(i-1))+hydro.vx(i-2));
            dUxdx[i] = idx*0.5*(hydro.vx(i+1)-hydro.vx(i-1));
            // dUxdx[i] = idx/3.0*(-Cx[i+3]+6.0*Cx[i+2]-3.0*Cx[i+1]-2.0*Cx[i]);
        }


        // vt += fluidvx;
        size_t l0(Din.l0());
        size_t m0(Din.m0());

        SHarmonic1D tG(G), tH(H);
        
        // dCdx df/dv term
        TMP = Din(0,0);     //fd1 = fd1.Dx();
        MakeGH(TMP, 0);            MakeG00(TMP);     tH=H;tG=G;
        
        tempv = vt *  XX1(0,0);    tH.mpaxis(tempv); tH.mxaxis(dUxdx);   Dh(0,0) += tH; tH = H;  
        tempv = vt *  XX2(0,0);    tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(0,0) += tG; tG = G; 
        tempv = vt *  XX3(0,0);    tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(2,0) += tG;

        // Dh.checknan();

        TMP = Din(1,0);     //fd1 = fd1.Dx();
        MakeGH(TMP,1);             tH=H;tG=G;
        tempv = vt *  XX1(1,0);    tH.mpaxis(tempv); tH.mxaxis(dUxdx);   Dh(1,0) += tH; tH = H;  
        tempv = vt *  XX2(1,0);    tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(1,0) += tG; tG = G; 
        tempv = vt *  XX3(1,0);    tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(3,0) += tG; 
    
        // Dh.checknan();

        TMP = Din(1,1);    
        MakeGH(TMP,1);              tH=H;tG=G;
        tempv = vt *  XX1(1,1);     tH.mpaxis(tempv); tH.mxaxis(dUxdx);   Dh(1,1) += tH; tH = H;  
        tempv = vt *  XX2(1,1);     tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(1,1) += tG; tG = G; 
        tempv = vt *  XX3(1,1);     tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(3,1) += tG; tG = G; 
        // tempv = vt *  1.0;
        
        MakeGH(TMP,l0);             tH=H;tG=G;
        for (size_t m(0); m<((m0<l0-1)?(m0+1):(l0-1)); ++m){
            TMP = Din(l0,m);  
            tempv = vt *  XX1(l0,m);   tH.mpaxis(tempv); tH.mxaxis(dUxdx);   Dh(l0,m) += tH; tH = H; 
            tempv = vt *  XX2(l0,m);   tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(l0,m) += tG; tG = G; 
            tempv = vt *  XX4(l0,m);   tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(l0-2,m) += tG; tG = G; 
            // tempv = vt *  1.0);            
        }
        for (size_t m(l0-2); m<((m0<l0)?(m0):(l0+1)); ++m){
            // TMP = Din(l0,m);     //fd1 = fd1.Dx();
            tempv = vt *  XX1(l0,m);   tH.mpaxis(tempv); tH.mxaxis(dUxdx);   Dh(l0,m) += tH; tH = H; 
            tempv = vt *  XX2(l0,m);   tG.mpaxis(tempv); tG.mxaxis(dUxdx);   Dh(l0,m) += tG; tG = G; 
            // tempv = vt *  XX4(l0,m));      TMP =  H.mpaxis(tempv);            Dh(l0-2,m) += TMP;
            // tempv = vt *  1.0);            
        }
        

        MakeGH(TMP,l0-1);       tH=H;tG=G;
        for (size_t m(0); m<((m0<l0-2)?(m0+1):(l0-2)); ++m){
            // TMP = Din(l0-1,m);     //fd1 = fd1.Dx();
            tempv = vt *  XX1(l0-1,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l0-1,m) += tH; tH = H; 
            tempv = vt *  XX2(l0-1,m);  tG.mpaxis(tempv); tG.mxaxis(dUxdx);  Dh(l0-1,m) += tG; tG = G;         
            tempv = vt *  XX4(l0-1,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l0-3,m) += tH; tH = H; 
            // tempv = vt *  1.0;                        
        }

        for (size_t m(l0-3); m<((m0<l0)?(m0+1):(l0+1)); ++m){
            // TMP = Din(l0-1,m);     //fd1 = fd1.Dx();            
            tempv = vt *  XX1(l0-1,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l0-1,m) += tH; tH = H;             
            tempv = vt *  XX2(l0-1,m);  tG.mpaxis(tempv); tG.mxaxis(dUxdx);  Dh(l0-1,m) += tG; tG = G; 
            // tempv = vt *  1.0,m);                        
        }
    
        for (size_t l(2); l<l0-1; ++l){
            MakeGH(TMP,l);      tH=H;tG=G;
            for (size_t m(0); m < ((m0<l+1)?(m0+1):(l+1)); ++m){
                tempv = vt *  XX1(l,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l,m) += tH; tH = H; 
                tempv = vt *  XX2(l,m);  tG.mpaxis(tempv); tG.mxaxis(dUxdx);  Dh(l,m) += tG; tG = G; 
                // tempv = vt *  1.0;               
            }

            for (size_t m(0); m < ((m0<l-1)?(m0+1):(l-1)); ++m){
                tempv = vt *  XX4(l,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l-2,m) += tH; tH = H;
            }

            for (size_t m(0); m < ((m0<l+3)?(m0+1):(l+3)); ++m){                
                tempv = vt *  XX3(l,m);  tH.mpaxis(tempv); tH.mxaxis(dUxdx);  Dh(l+2,m) += tH; tH = H; 
                // tempv = vt *  1.0;               
            }
        }

        // C df/dr term
        for (size_t m(0); m < ((m0<l0)?(m0+1):m0); ++m){
            
            tH = Din(m,m);                             tH = tH.Dx();
            for (int ip(0); ip<tH.nump(); ++ip) tH(ip,Din(0,0).numx()-1)  = 0.0; 
            Dh(m,m) += tH.mxaxis(Ux);

            for (size_t l(m+1); l < l0; ++l) {
               tH = Din(l,m);                          tH = tH.Dx();  
               for (int ip(0); ip<tH.nump(); ++ip) tH(ip,Din(0,0).numx()-1)  = 0.0;  
               
               Dh(l,m) += tH.mxaxis(Ux);

            }
            
            tH = Din(l0,m);                            tH = tH.Dx();
            for (int ip(0); ip<tH.nump(); ++ip) tH(ip,Din(0,0).numx()-1)  = 0.0; 

            Dh(l0,m) += tH.mxaxis(Ux);
        
        
        }
}

//--------------------------------------------------------------
//--------------------------------------------------------------
//  Make derivatives 2*Dp*(l+1/l)*G and -2*Dp*H for a given f
    void Hydro_Advection_1D::MakeGH(SHarmonic1D& f, size_t el){
//--------------------------------------------------------------
        valarray<complex<double> > invpax(invpr); 
        complex<double> ld(el); 

        invpax *= (-2.0)*(ld+1.0) * (pr[1]-pr[0]);
       
        G = f;                   H = f;
        G = G.Dp();
                                 H  = H.mpaxis(invpax);
                                 H += G;
        G *= -(2.0*ld+1.0)/ld;
        G += H;

        for (size_t i(0); i < G.numx(); ++i) G(0,i) = 0.0;
        for (size_t i(0); i < H.numx(); ++i) H(0,i) = f(1,i) * Hp0[el];
    }
//--------------------------------------------------------------
//  Calculation of G00 = -2*Dp* df/dp(p0)
    void Hydro_Advection_1D::MakeG00(SHarmonic1D& f) {
//--------------------------------------------------------------
        G = f; G = G.Dp();                            

        complex<double> p0p1_sq( pr[0]*pr[0]/(pr[1]*pr[1]) ),
               inv_mp0p1_sq( 1.0/(1.0-p0p1_sq) ),
               g_r = -4.0*(pr[1]-pr[0]) * pr[0]/(pr[1]*pr[1]),
               f00;

        for (size_t i(0); i < f.numx(); ++i) {
             f00    = ( f(0,i) - f(1,i) * p0p1_sq) * inv_mp0p1_sq;  
             G(0,i) = ( f(1,i) - f00) * g_r;
        }
    }

//--------------------------------------------------------------