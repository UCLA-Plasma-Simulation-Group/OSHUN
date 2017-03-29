/*! \brief Plasma Formulary and Units - Definitions
* \author PICKSC
 * \date   September 1, 2016
 * \file   formulary.cpp
 * 
 * Definitions for data structures that are related to
 *   simple physics concepts:
 *  1) units: a concrete class that combines a label and a
 *       numerical value
 *  2) Formulary: A class that enables conversion of units
 *       and implements expression from the NRL plasma
 *       formulary. 
 * 
 */
//-------------------------------------------------------------------

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"

//  Declerations
#include "input.h"
#include "formulary.h"


//-------------------------------------------------------------------
//-------------------------------------------------------------------
//  List all the units as in "label, conversion_factor" 
//-------------------------------------------------------------------
Formulary::Formulary() : n(Input::List().density_np),
                        wp(5.64*1.0e+4*sqrt(n)),
                        skindepth(cL/wp),
                        B0(-1.0*wp*me/qe),
                        T0(pow(Input::List().pth_ref,2.0)*keVnorm*1e3),
                        Zeta(Input::List().hydrocharge)
                        // ref_nuei(sqrt(2.0/pi)*LOGei(1.0,Input::List().pth_ref,Zeta)/exp(LOGei(1.0,Input::List().pth_ref,Zeta)) * wp)
                        {
    // Physically: label * d = const in any system

    // velocity - c
    D["Velocity"]     = units("c",1.0);    
    D["Velocity_cgs"] = units("cm/sec", c );    
    D["Velocity_si"]  = units("m/sec",  c * 0.01 );    

    D["v"]     = units("c",1.0);    
    D["v_cgs"] = units("cm/sec", c );    
    D["v_si"]  = units("m/sec",  c * 0.01 );    

    D["vx"]     = units("c",1.0);    
    D["vx_cgs"] = units("cm/sec", c );    
    D["vx_si"]  = units("m/sec",  c * 0.01 );    

    D["vy"]     = units("c",1.0);    
    D["vy_cgs"] = units("cm/sec", c );    
    D["vy_si"]  = units("m/sec",  c * 0.01 );    

    D["vz"]     = units("c",1.0);    
    D["vz_cgs"] = units("cm/sec", c );    
    D["vz_si"]  = units("m/sec",  c * 0.01 );    

    // time - 1/wp
    D["Time"]     = units("1/wp",1.0);    
    D["Time_cgs"] = units("sec", 1.0/wp);    
    D["Time_si"]  = units("sec", 1.0/wp);    
    D["Time_fs"]  = units("fsec", 1.0/wp*1.0e+15);    
    D["Time_ps"]  = units("psec", 1.0/wp*1.0e+12);    

    D["t"]     = units("1/wp",1.0);    
    D["t_cgs"] = units("sec", 1.0/wp);    
    D["t_si"]  = units("sec", 1.0/wp);    

    // space - c/wp
    D["Space"]     = units("c/wp",1.0);    
    D["Space_cgs"] = units("cm",       c/wp);    
    D["Space_si"]  = units("m", 0.01 * c/wp); // 1 cm = 0.01 m  

    D["x"]     = units("c/wp",1.0);    
    D["x_cgs"] = units("cm",       c/wp);    
    D["x_si"]  = units("m", 0.01 * c/wp); // 1 cm = 0.01 m  

    D["y"]     = units("c/wp",1.0);    
    D["y_cgs"] = units("cm",       c/wp);    
    D["y_si"]  = units("m", 0.01 * c/wp); // 1 cm = 0.01 m  

    D["z"]     = units("c/wp",1.0);    
    D["z_cgs"] = units("cm",       c/wp);    
    D["z_si"]  = units("m", 0.01 * c/wp); // 1 cm = 0.01 m  

    // density - n
    D["Density"]     = units("n",1.0);    
    D["Density_cgs"] = units("cm^-3",         n);    
    D["Density_si"]  = units("m^-3", 1.0e+6 * n);   // cm^-3 = 10^6 m^-3 

    D["n"]     = units("n",1.0);    
    D["n_cgs"] = units("cm^-3",         n);    
    D["n_si"]  = units("m^-3", 1.0e+6 * n);   // cm^-3 = 10^6 m^-3 

    // charge - e
    D["Charge"]     = units("e",1.0);    
    D["Charge_cgs"] = units("Fr",           e);    
    D["Charge_si"]  = units("C", (10.0/c) * e);  // 1 Fr = (10/c) * C  

    D["q"]     = units("e",1.0);    
    D["q_cgs"] = units("Fr",           e);    
    D["q_si"]  = units("C", (10.0/c) * e);  // 1 Fr = (10/c) * C  

    // momentum - mc
    D["Momentum"]     = units("mc",1.0);    
    D["Momentum_cgs"] = units("gr*cm/sec", m*c );    
    D["Momentum_si"]  = units("kg*m/sec",  m*c * 1.0e-5 );    

    D["p"]     = units("mc",1.0);    
    D["p_cgs"] = units("gr*cm/sec", m*c );    
    D["p_si"]  = units("kg*m/sec",  m*c * 1.0e-5 );    

    D["px"]     = units("mc",1.0);    
    D["px_cgs"] = units("gr*cm/sec", m*c );    
    D["px_si"]  = units("kg*m/sec",  m*c * 1.0e-5 );    

    D["py"]     = units("mc",1.0);    
    D["py_cgs"] = units("gr*cm/sec", m*c );    
    D["py_si"]  = units("kg*m/sec",  m*c * 1.0e-5 );    

    D["pz"]     = units("mc",1.0);    
    D["pz_cgs"] = units("gr*cm/sec", m*c );    
    D["pz_si"]  = units("kg*m/sec",  m*c * 1.0e-5 );    
//-----

    // energy - mc^2
    D["Energy"]     = units("mc^2",1.0);    
    D["Energy_cgs"] = units("erg",          m*pow(c,2)   );       
    D["Energy_si"]  = units("J",   1.0e-7 * m*pow(c,2)   );  // erg = 10^-7 J  
    D["Energy_eV"]  = units("eV",  1.0e-7 * m*pow(c,2)/((10.0/c) * e) );    

    D["T"]     = units("mc^2",1.0);    
    D["T_cgs"] = units("erg",          m*pow(c,2)   );       
    D["T_si"]  = units("J",   1.0e-7 * m*pow(c,2)   );  // erg = 10^-7 J  
    D["T_eV"]  = units("eV",  1.0e-7 * m*pow(c,2)/((10.0/c) * e) );    

    // current - ewp
    D["Current"]     = units("ewp",1.0);    
    D["Current_cgs"] = units("Fr/s",        e*wp );   
    D["Current_si"]  = units("A",  10.0/c * e*wp );  // Fr/s = 10/c A

    D["Jx"]     = units("ewp",1.0);    
    D["Jx_cgs"] = units("Fr/s",        e*wp );   
    D["Jx_si"]  = units("A",  10.0/c * e*wp );  // Fr/s = 10/c A

    D["Jy"]     = units("ewp",1.0);    
    D["Jy_cgs"] = units("Fr/s",        e*wp );   
    D["Jy_si"]  = units("A",  10.0/c * e*wp );  // Fr/s = 10/c A

    D["Jz"]     = units("ewp",1.0);    
    D["Jz_cgs"] = units("Fr/s",        e*wp );   
    D["Jz_si"]  = units("A",  10.0/c * e*wp );  // Fr/s = 10/c A

    // pressure - n*mc^2 
    D["Pressure"]      = units("n*mc^2",1.0);    
    D["Pressure_cgs"]  = units("Ba",             n*m*c*c ); 
    D["Pressure_si"]   = units("Pa",       0.1 * n*m*c*c ); // Ba = 0.1 Pa  
    D["Pressure_Mbar"] = units("Mbar", 1.0e-12 * n*m*c*c ); // Ba = 10^-12Mbar 

    D["P"]      = units("n*mc^2",1.0);    
    D["P_cgs"]  = units("Ba",             n*m*c*c ); 
    D["P_si"]   = units("Pa",       0.1 * n*m*c*c ); // Ba = 0.1 Pa  
    D["P_Mbar"] = units("Mbar", 1.0e-12 * n*m*c*c ); // Ba = 10^-12Mbar 
//-----

    // Electric field - mcwp/e
    D["Efield"]     = units("mcwp/e",1.0);    
    D["Efield_cgs"] = units("statV/cm",       m*c*wp/e);         
    D["Efield_si"]  = units("V/m", c*1.0e-6 * m*c*wp/e);  // statV/m = 10^-6*c * V/m     

    D["Ex"]     = units("mcwp/e",1.0);    
    D["Ex_cgs"] = units("statV/cm",       m*c*wp/e);         
    D["Ex_si"]  = units("V/m", c*1.0e-6 * m*c*wp/e);  // statV/m = 10^-6*c * V/m     

    D["Ey"]     = units("mcwp/e",1.0);    
    D["Ey_cgs"] = units("statV/cm",       m*c*wp/e);         
    D["Ey_si"]  = units("V/m", c*1.0e-6 * m*c*wp/e);  // statV/m = 10^-6*c * V/m     

    D["Ez"]     = units("mcwp/e",1.0);    
    D["Ez_cgs"] = units("statV/cm",       m*c*wp/e);         
    D["Ez_si"]  = units("V/m", c*1.0e-6 * m*c*wp/e);  // statV/m = 10^-6*c * V/m     

    // Magnetic field - mcwp/e 
    D["Bfield"]     = units("mcwp/e",1.0);    
    D["Bfield_cgs"] = units("Gauss",          m*c*wp/e);           
    D["Bfield_si"]  = units("Tesla", 1.0e-4 * m*c*wp/e);  // Gauss = 10^-4 Tesla           

    D["Bx"]     = units("mcwp/e",1.0);    
    D["Bx_cgs"] = units("Gauss",          m*c*wp/e);           
    D["Bx_si"]  = units("Tesla", 1.0e-4 * m*c*wp/e);  // Gauss = 10^-4 Tesla           

    D["By"]     = units("mcwp/e",1.0);    
    D["By_cgs"] = units("Gauss",          m*c*wp/e);           
    D["By_si"]  = units("Tesla", 1.0e-4 * m*c*wp/e);  // Gauss = 10^-4 Tesla           

    D["Bz"]     = units("mcwp/e",1.0);    
    D["Bz_cgs"] = units("Gauss",          m*c*wp/e);           
    D["Bz_si"]  = units("Tesla", 1.0e-4 * m*c*wp/e);  // Gauss = 10^-4 Tesla           

    // force - mcwp 
    D["Force"]     = units("mcwp",1.0);    
    D["Force_cgs"] = units("dyne",       m*c*wp); 
    D["Force_si"]  = units("N", 1.0e-5 * m*c*wp ); // dyne = 10^-5 N  

}

//-------------------------------------------------------------------
//   Calculate the Coulomb logarithm 
//-------------------------------------------------------------------
double Formulary::LOGee(double ne, double Te) {   
//   ne =  density/np, Te = energy/mc^2
//   Note: we assume nonrelativistic distribution functions
    if (ne < nmin) return 2.0; 

    // Te /= (3.0*ne);
    Te *= Units("Energy","eV").d; // Temperature in eV
    ne *= Units("Density","cgs").d;

    Te = log(Te); 
    ne = log(ne);

    return max(2.0,23.5 - 0.5*ne + 1.25*Te - sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0)));
}

double Formulary::LOGee(double ne, string un, double Te, string uT) {   
//   where "un" are the units of ne, and uT of Te
    Te /= Units("Energy",uT).d; 
    ne /= Units("Density",un).d;

    return LOGee(ne,Te);
}

double Formulary::LOGei(double ne, double Te, double Z) {   
//   ne =  density/np, Te = energy/mc^2
//   Note: we assume nonrelativistic distribution functions
    if (ne < nmin) return 2.0; 

    Te *= Units("Energy","eV").d; // Temperature in eV
    ne *= Units("Density","cgs").d;

    if ( Te > 10.0*Z*Z ) {
        return max(2.0, 24.0 - 0.5*log(ne) + log(Te));
    }
    return max(2.0, 23.0 - 0.5*log(ne) + 1.5*log(Te) - log(Z));
}

double Formulary::LOGei(double ne, string un, double Te, string uT, double Z) {   
//   where "un" are the units of ne, and uT of Te
    Te /= Units("Energy",uT).d; // Temperature in eV
    ne /= Units("Density",un).d;

    return LOGei(ne,Te,Z);
}

double Formulary::LOGii(double m1, double Z1, double n1, double T1, 
    double m2, double Z2, double n2, double T2) {   
//   ne =  density/np, Te = energy/mc^2
//   Note: we assume nonrelativistic distribution functions
    double mp=1;
    if (n1 < nmin || n2 < nmin) return 2.0; 

    T1 *= Units("Energy","eV").d; // Temperature in eV
    T2 *= Units("Energy","eV").d; // Temperature in eV
    n1 *= Units("Density","cgs").d;
    n2 *= Units("Density","cgs").d;

    return max(2.0, 23.0 - log(Z1*Z2*(m1/mp+m2/mp)/(m1/mp*T2+m2/mp*T1)*sqrt(n1*Z1*Z1/T1+n2*Z2*Z2/T2)));
}

//-------------------------------------------------------------------
//   Simple formulary expresions
//-------------------------------------------------------------------

//   Thermal velocity: NRL Formulary 2009, p.29
double Formulary::vth(double Te){
    Te *= Units("Energy","eV").d;
    return sqrt(Te)*4.19e+7/c;
}
double Formulary::vth(double Te, string uT){
    Te /= Units("Energy", uT).d;
    return vth(Te);
}
//   Electron-ion self-collision time: NRL Formulary 2009, p.37
double Formulary::Tau_i(double ne, double Te, double Zeta){
    return 2.09e+7*(pow(sqrt(Te/Zeta*Units("Energy","eV").d),3)/
            (ne*Units("Density","cgs").d*LOGei(ne,Te,Zeta))) /sqrt(me_over_mp);
}

//   Electron self-collision time: NRL Formulary 2009, p.37
double Formulary::Tau_e(double ne, double Te){
    return 3.44e+5*(pow(sqrt(Te*Units("Energy","eV").d),3)/
            (ne*Units("Density","cgs").d*LOGee(ne,Te))) ; //
            ///Units("Time","cgs").d;
}
double Formulary::Tau_e(double ne,  string un, double Te,  string uT){
    Te /= Units("Energy",uT).d; 
    ne /= Units("Density",un).d;
    return Tau_e(ne,Te);
}

//   Electron mean free path
double Formulary::MFP(double ne, double Te){
    return vth(Te)*Tau_e(ne,Te);
}
double Formulary::MFP(double ne,  string un, double Te,  string uT){
    Te /= Units("Energy",uT).d; 
    ne /= Units("Density",un).d;
    return MFP(ne,Te);
}
//-------------------------------------------------------------------

Formulary& formulary() {
    static Formulary f;
    return f;
}
//-------------------------------------------------------------------
//*******************************************************************
