/*! \brief Plasma Formulary and Units - Declarations
* \author PICKSC
 * \date   September 1, 2016
 * \file   formulary.h
 * 
 * Declarations for data structures that are related to
 *   simple physics concepts:
 *   1) units: a concrete class that combines a label and a
 *       numerical value
 *   2) Formulary: A class that enables conversion of units
 *       and implements expression from the NRL plasma
 *       formulary. 
 * 
 */
//-------------------------------------------------------------------


    #ifndef DECL_FORMULARY_H
    #define DECL_FORMULARY_H


//--------------------------------------------------------------
// Make Gaussian
template<class T> 
valarray<T> Gaussian(const size_t N, const T vmin, const T vmax, const T vth){

//  Make axis first
    valarray<T> G(N);
    for (size_t i(0); i < N; ++i) {
        G[i] = static_cast<T>(i);
    }
    G *= (vmax-vmin)/(static_cast<T>(N-1));
    G += vmin;

//  Make Gaussian
    T           C(pow( 1.0/ (sqrt(2.0*M_PI)*vth), 3));
    T           al( (-0.5) / (vth*vth));
    for (size_t i(0); i < N; ++i) {
        G[i]  = exp( al * G[i]*G[i] ); 
    }
    G *= C;

    return G;                              
}
//-------------------------------------------------------------------

//--------------------------------------------------------------
class units {
    public:
//      Contents
        string label;   // e.g. label = sec
        double d;  // e.g. x[label] = c * x[1/wp] => c = 1/wp 

//      Constructors
        units() : label("default"), d(0.0) {}
        units(string _x, double _d) : label(_x), d(_d){}

//      Copy constructor 
        units(const units& other) { 
           label = other.label; 
           d = other.d;
        } 
        ~units(){ }
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
class Formulary {
    public:
//    Construct an underlying dictionary for units
      Formulary();

//    Access to unit systems
      units Units(string key) { return D[key]; }
      units Units(string key1, string key2) { return D[key1+"_"+key2]; }
      string Label(string key) { return D[key].label; }
      string Label(string key1, string key2) { return D[key1+"_"+key2].label; }
      double Uconv(string key) { return D[key].d; }
      double Uconv(string key1, string key2) { return D[key1+"_"+key2].d; }

//    Coulomb logarithms
      double LOGee(double ne, double Te);
      double LOGee(double ne, string un, double Te, string uT);
      double LOGei(double ne, double Te, double Z);
      double LOGei(double ne, string un, double Te, string uT, double Z);
      double LOGii(double m1, double Z1, double n1, double T1, double m2, double Z2, double n2, double T2);   

//    Formulary functions
      double vth(double Te);
      double vth(double Te, string uT);
      double Tau_e(double ne, double Te);
      double Tau_e(double ne, string un, double Te, string uT);
      double Tau_i(double ne, double Te, double Zeta);
      double MFP(double ne, double Te);
      double MFP(double ne, string un, double Te, string uT);

//    Normalization density 
      //const double n = 1.0e+21;               // cm-3
      double n; // cm-3
      double wp;
      double skindepth;
      double B0;
      double T0;
      double Zeta;
      // const double ref_nuei;
      
//    Physical constants

      // double oneover_atomic_Z = 1.0/Input::List().Zeta;
      static constexpr double pi=3.141592653589793238;

      static constexpr double cL = 299792458;//3e8; // speed of light over v electron thermal
      static constexpr double eps0=8.854187817e-12;//8.85e-12;
      static constexpr double qe = -1.60217662e-19; //1.6e-19 // electron charge.
      static constexpr double me = 9.10938356e-31; //9.11e-31 //electron mass
      static constexpr double me_over_mp = 0.000544617024;
      static constexpr double keVnorm = 510.9989461; //to be multiplied by v^2/c^2

      static constexpr double c = 2.99792458*1.0e+10;    // cm/sec
      static constexpr double e = 4.80320425*1.0e-10;    // Franklin or statC 
      static constexpr double m = 9.10938215*1.0e-28;    // gr
       
    private:
      map<string,units> D;

      static constexpr double nmin = 1.0e-8;
};
//--------------------------------------------------------------

Formulary& formulary();
//--------------------------------------------------------------
//**************************************************************


    #endif
