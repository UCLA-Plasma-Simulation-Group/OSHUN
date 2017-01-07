/*! \brief  Fluid Equation - Definitions
 * \author Archis Joglekar
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
    idx(numx/(xmax-xmin)),szx(numx),FEQ(xmin,xmax,numx)             {}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
void Hydro_Functor::operator()(const State1D& Yin, State1D& Yslope){

    Yslope = 0.0;
    

                                        FEQ.calcquantities(Yin);
    Yslope.HYDRO().densityarray()  =    FEQ.density(Yin.HYDRO());
    Yslope.HYDRO().velocityarray() =    FEQ.velocity(Yin);
                                        FEQ.updateE(Yin.HYDRO().velocityarray(),Yslope.HYDRO().velocityarray(),Yslope.EMF());

}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
void Hydro_Functor::operator()(const State1D& Yin, State1D& Yslope, size_t dir){}
void Hydro_Functor::operator()(const State1D& Y1in, const State1D& Y2in, State1D& Yslope){}
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
//  Fluid Equation Class Constructor
 
Fluid_Equation_1D::Fluid_Equation_1D(double xmin, double xmax, size_t Nx) 
  : idx(Nx/(xmax-xmin)), szx(Nx), dummy(0.0,Nx)
    {
        Nbc = Input::List().RKLevel;
    }

//--------------------------------------------------------------
//--------------------------------------------------------------

valarray<double> Fluid_Equation_1D::density(const Hydro1D& Hin){

    dummy[0]    = idx*(Hin.density(1)*Hin.velocity(1)-Hin.density(0)*Hin.velocity(0));            

    for (size_t ix(1);ix<szx-1;++ix)
    {
        dummy[ix]     = 0.5*idx*(Hin.density(ix+1)*Hin.velocity(ix+1)-Hin.density(ix-1)*Hin.velocity(ix-1));            
        // std::cout << "\n\n hydrodensity[" << ix << "] = " << Hin.density(ix) ;       
    }

    dummy[szx-1]     = idx*(Hin.density(szx-1)*Hin.velocity(szx-1)-Hin.density(szx-2)*Hin.velocity(szx-2));  

    return dummy; 
    
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
//
////--------------------------------------------------------------
//--------------------------------------------------------------

valarray<double>  Fluid_Equation_1D::velocity(const State1D& Yin){
    // State1D& Yin, double deltat){

    valarray<double> totalpressure(Yin.HYDRO().kineticspeciespressurearray());
    valarray<double> nh(Yin.HYDRO().densityarray());
    valarray<double> Th(Yin.HYDRO().temperaturearray());
    valarray<double> Ux(Yin.HYDRO().velocityarray());

    valarray<double> dUxdx(0.0,Ux.size());
    

    totalpressure                  *=  Yin.HYDRO().charge()/Yin.HYDRO().mass();

    dUxdx[0]                        = idx*(Ux[1]-Ux[0]);

    dummy[0]                        =   (      -idx*(
                                        (totalpressure[1]-totalpressure[0])
                                        +   (nh[1]*Th[1]-nh[0]*Th[0])/Yin.HYDRO().mass())   
            // +   Yin.HYDRO().jy(ix)*Yin.EMF().Bz()(ix)-Yin.HYDRO().jz(ix)*Yin.EMF().By()(ix)   
                                        +   (Yin.EMF().Ex()(0)).real()*(Yin.HYDRO().charge()*nh[0]-Yin.HYDRO().ne(0))/Yin.HYDRO().mass());


    dummy[0]                       -= Ux[0] * dUxdx[0];
    // Yin.HYDRO().velocity(0)        = Yin.HYDRO().velocity(0) + deltat * (dCdt[0]);



    for (size_t ix(1);ix<szx-1;++ix)
    {
        
        dUxdx[ix]                       = 0.5*idx*(Ux[ix+1]-Ux[ix-1]);

        dummy[ix]                       =   (      -0.5*idx*(
                                            (totalpressure[ix+1]-totalpressure[ix-1])
                                            +   (nh[ix+1]*Th[ix+1]-nh[ix-1]*Th[ix-1])/Yin.HYDRO().mass())   
            // +   Yin.HYDRO().jy(ix)*Yin.EMF().Bz()(ix)-Yin.HYDRO().jz(ix)*Yin.EMF().By()(ix)   
                                            +   (Yin.EMF().Ex()(ix)).real()*(Yin.HYDRO().charge()*nh[ix]-Yin.HYDRO().ne(ix))/Yin.HYDRO().mass());

        dummy[ix]                       -= Ux[ix] * dUxdx[ix];
        // std::cout << "\n\n hydrodensity[" << ix << "] = " << dummy[ix] ;
    
    }
    
    dUxdx[szx-1]                        = idx*(Ux[szx-1]-Ux[szx-2]);

    dummy[szx-1]                        =   (      -idx*(
                                            (totalpressure[szx-1]-totalpressure[szx-2])
                                        +   (nh[szx-1]*Th[szx-1]-nh[szx-2]*Th[szx-2])/Yin.HYDRO().mass())   
            // +   Yin.HYDRO().jy(ix)*Yin.EMF().Bz()(ix)-Yin.HYDRO().jz(ix)*Yin.EMF().By()(ix)   
                                        +   (Yin.EMF().Ex()(szx-1)).real()*(Yin.HYDRO().charge()*nh[szx-1]-Yin.HYDRO().ne(szx-1))/Yin.HYDRO().mass());


    dummy[szx-1]                       -= Ux[szx-1] * dUxdx[szx-1];
    
    return dummy;
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

void Fluid_Equation_1D::updateE(const valarray<double>& hydrovelocity,const valarray<double>& dCdt,
                                EMF1D& EMF){

    for (size_t ix(Nbc);ix<szx-Nbc;++ix){        
        EMF.Ex()(ix) -= dCdt[ix]; //+0.5*dC2xdx[ix]+(hydrocharge*electrondensity[ix]-hydrodensity[ix]); // Needs JxB in 2D
        EMF.Ey()(ix) += -1.0*hydrovelocity[ix]*EMF.Bz()(ix);//+jx[ix]*Bz(ix)-jz[ix]*Bx(ix);
        EMF.Ez()(ix) += hydrovelocity[ix]*EMF.By()(ix);//-jx[ix]*By(ix)-jy[ix]*Bx(ix);
    }
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

void Fluid_Equation_1D:: calcquantities(const State1D& Yin){
        
    Yin.HYDRO().kineticspeciespressurearray()=0.0;
   
    for (size_t s(0);s<Yin.Species();++s){

        Yin.HYDRO().kineticspeciespressurearray() += Yin.DF(s).getpressure();   
        Yin.HYDRO().jxarray() += (Yin.DF(s).getcurrent(0));
        Yin.HYDRO().jyarray() += (Yin.DF(s).getcurrent(1));
        Yin.HYDRO().jzarray() += (Yin.DF(s).getcurrent(2));         
        // std::cout << "kinetic pressure[" << ix << "]=" << Yin.HYDRO().kpressure(ix) << "\n";
        // } 
    }           
   
    Yin.HYDRO().electrondensityarray() =  Yin.DF(0).getdensity();

}   



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
     pr(Algorithms::MakeAxis(static_cast<complex<double> >(pmin),
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
         idp = (-1.0)/ (2.0*(pmax-pmin)/double(Np-1)); // -1/(2dp) 
         idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx)); // -1/(2dx) 

         Nbc = Input::List().RKLevel;

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
        

        valarray<complex<double> > Ux(0.0,(hydro.velocityarray()).size());
        valarray<complex<double> > dUxdx(0.0,(hydro.velocityarray()).size());


        Ux[0] = hydro.velocity(0);
        dUxdx[0] = idx*(hydro.velocity(1)-hydro.velocity(0));

        Ux[szx-1] = hydro.velocity(szx-1);
        dUxdx[szx-1] = idx*(hydro.velocity(szx-1)-hydro.velocity(szx-2));
        for (size_t i(1);i<szx-1;++i)
        {
            Ux[i] = hydro.velocity(i);
            // dUxdx[i] = idx/12.0*(-hydro.velocity(i+2)+8.0*(hydro.velocity(i+1)-hydro.velocity(i-1))+hydro.velocity(i-2));
            dUxdx[i] = idx*0.5*(hydro.velocity(i+1)-hydro.velocity(i-1));
            // dUxdx[i] = idx/3.0*(-Cx[i+3]+6.0*Cx[i+2]-3.0*Cx[i+1]-2.0*Cx[i]);
        }


        // vt += fluidvelocity;       
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
             G(0,i) = ( f(1,0) - f00) * g_r;
        }
    }

//--------------------------------------------------------------