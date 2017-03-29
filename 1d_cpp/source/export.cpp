/*! \brief  Export Files - Definitions
 * \author PICKSC
 *  \date   September 1, 2016
 *  \file   export.cpp
 * 
 */

//  Standard libraries
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cstring>

#include <math.h>
#include <map>

#include <sys/stat.h>
#include <sys/types.h>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"
#include "H5Cpp.h"
#include "exprtk.hpp"

//  Declarations
#include "input.h"
#include "state.h"
#include "formulary.h"
#include "setup.h"
#include "nmethods.h"
#include "parallel.h"
#include "nmethods.h"
#include "export.h"



//------------------------------------------------------------------------------
// Create a folder
//
// @param[in]  _name  The name of the folder
//
int Export_Files::Makefolder(string _name){

    mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    // char*   foldername = new char [_name.length() + 1];
    

    const char*   foldername = _name.data();
    

    // strcpy(foldername, _name.c_str());
    
    // vector<char> foldername(_name.c_str(), _name.c_str() + _name.size() + 1);

    int   status(mkdir(foldername,_permissions));

    // if (status != 0) cout << "Warning: Folder "<< _name<< " exists\n";

    // delete[] foldername;

    return status;
}
//--------------------------------------------------------------
//---------------------------------------------------------------
void Export_Files::Folders(){
//--------------------------------------------------------------
//   create the directory tree
//--------------------------------------------------------------

    if (Makefolder("RESTART") != 0) cout << "Warning: Folder 'RESTART' exists" << endl;

    if (Makefolder("OUTPUT") != 0) cout << "Warning: Folder 'OUTPUT' exists" << endl;

    if ( Input::List().o_EHist ) {
        if ( Makefolder("OUTPUT/NUM") != 0) cout << "Warning: Folder 'OUTPUT/NUM' exists" << endl;
    }

    if ( Input::List().o_Ex ||  Input::List().o_Ey || Input::List().o_Ez ||
         Input::List().o_Bx ||  Input::List().o_By || Input::List().o_Bz )  {
        if (Makefolder("OUTPUT/FLD") != 0)
            cout<<"Warning: Folder 'OUTPUT/FLD' exists" << endl;

        if (Input::List().o_Ex) {
            if (Makefolder("OUTPUT/FLD/Ex") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/EX' exists" << endl;
        }
        if (Input::List().o_Ey) {
            if (Makefolder("OUTPUT/FLD/Ey") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/EY' exists" << endl;
        }
        if (Input::List().o_Ez) {
            if (Makefolder("OUTPUT/FLD/Ez") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/EZ' exists" << endl;
        }
        if (Input::List().o_Bx) {
            if (Makefolder("OUTPUT/FLD/Bx") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/BX' exists" << endl;
        }
        if (Input::List().o_By) {
            if (Makefolder("OUTPUT/FLD/By") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/BY' exists" << endl;
        }
        if (Input::List().o_Bz) {
            if (Makefolder("OUTPUT/FLD/Bz") != 0)
                cout<<"Warning: Folder 'OUTPUT/FLD/BZ' exists" << endl;
        }
    }

    if ( Input::List().o_x1x2 ||  Input::List().o_pth ||
         Input::List().o_G  ||
         Input::List().o_Px || Input::List().o_PxPx ||
         Input::List().o_Py || Input::List().o_PxPy || Input::List().o_PyPy ||
         Input::List().o_Pz || Input::List().o_PxPz || Input::List().o_PyPz || Input::List().o_PzPz ||
         Input::List().o_Vx || Input::List().o_VxVx ||
         Input::List().o_Vy || Input::List().o_VxVy || Input::List().o_VyVy ||
         Input::List().o_VxVz || Input::List().o_VyVz || Input::List().o_VzVz ||
         Input::List().o_Vsq || Input::List().o_Temperature || Input::List().o_Pressure ||
         Input::List().o_Qx || Input::List().o_Qy || Input::List().o_Qz ||
         Input::List().o_p1x1 )  {

        if (Makefolder("OUTPUT/MOM") != 0)
            cout<<"Warning: Folder 'OUTPUT/MOM' exists" << endl;

//          Relativistic Energy - Momentum Tensor  
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_x1x2) {
            if (Makefolder("OUTPUT/MOM/n") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/N' exists" << endl;
        }
        if (Input::List().o_pth) {
            if (Makefolder("OUTPUT/MOM/Pth") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Pth' exists" << endl;
        }
        if (Input::List().o_G) {
            if (Makefolder("OUTPUT/MOM/Gam") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Gam' exists" << endl;
        }
        if (Input::List().o_Px) {
            if (Makefolder("OUTPUT/MOM/Px") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Px' exists" << endl;
        }
        if (Input::List().o_PxPx) {
            if (Makefolder("OUTPUT/MOM/PxPx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PxPx' exists" << endl;
        }
        if (Input::List().o_Py) {
            if (Makefolder("OUTPUT/MOM/Py") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Py' exists" << endl;
        }
        if (Input::List().o_PxPy) {
            if (Makefolder("OUTPUT/MOM/PxPy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PxPy' exists" << endl;
        }
        if (Input::List().o_PyPy) {
            if (Makefolder("OUTPUT/MOM/PyPy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PyPy' exists" << endl;
        }
        if (Input::List().o_Pz) {
            if (Makefolder("OUTPUT/MOM/Pz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Pz' exists" << endl;
        }
        if (Input::List().o_PxPz) {
            if (Makefolder("OUTPUT/MOM/PxPz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PxPz' exists" << endl;
        }
        if (Input::List().o_PyPz) {
            if (Makefolder("OUTPUT/MOM/PyPz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PyPz' exists" << endl;
        }
        if (Input::List().o_PzPz) {
            if (Makefolder("OUTPUT/MOM/PzPz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/PzPz' exists"  << endl;
        }

//          Current
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_Jx) {
            if (Makefolder("OUTPUT/MOM/Jx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Jx' exists" << endl;
        }
        if (Input::List().o_Jy) {
            if (Makefolder("OUTPUT/MOM/Jy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Jy' exists" << endl;
        }
        if (Input::List().o_Jz) {
            if (Makefolder("OUTPUT/MOM/Jz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Jz' exists" << endl;
        }


//          Nonrelativistic output
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_Vx) {
            if (Makefolder("OUTPUT/MOM/Vx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Vx' exists" << endl;
        }
        if (Input::List().o_VxVx) {
            if (Makefolder("OUTPUT/MOM/VxVx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VxVx' exists" << endl;
        }
        if (Input::List().o_Vy) {
            if (Makefolder("OUTPUT/MOM/Vy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Vy' exists" << endl;
        }
        if (Input::List().o_VxVy) {
            if (Makefolder("OUTPUT/MOM/VxVy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VxVy' exists" << endl;
        }
        if (Input::List().o_VyVy) {
            if (Makefolder("OUTPUT/MOM/VyVy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VyVy' exists" << endl;
        }
        if (Input::List().o_VxVz) {
            if (Makefolder("OUTPUT/MOM/VxVz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VxVz' exists" << endl;
        }
        if (Input::List().o_VyVz) {
            if (Makefolder("OUTPUT/MOM/VyVz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VyVz' exists" << endl;
        }
        if (Input::List().o_VzVz) {
            if (Makefolder("OUTPUT/MOM/VzVz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/VzVz' exists" << endl;
        }
        if (Input::List().o_Vsq) {
            if (Makefolder("OUTPUT/MOM/Vsq") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Vsq' exists" << endl;
        }
        if (Input::List().o_Temperature) {
//                 if (Makefolder("OUTPUT/MOM/T_eV") != 0) 
//                     cout<<"Warning: Folder 'OUTPUT/MOM/T_eV' exists" << endl;
            if (Makefolder("OUTPUT/MOM/T") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/T' exists" << endl;
        }
        if (Input::List().o_Pressure) {
            if (Makefolder("OUTPUT/MOM/P_Mbar") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/P_Mbar' exists" << endl;
        }
        if (Input::List().o_ND) {
            if (Makefolder("OUTPUT/MOM/ND") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/ND' exists" << endl;
        }
        if (Input::List().o_Nu) {
            if (Makefolder("OUTPUT/MOM/Nu") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Nu' exists" << endl;
        }
        if (Input::List().o_Qx) {
            if (Makefolder("OUTPUT/MOM/Qx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Qx' exists" << endl;
        }
        if (Input::List().o_Qy) {
            if (Makefolder("OUTPUT/MOM/Qy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Qy' exists" << endl;
        }
        if (Input::List().o_Qz) {
            if (Makefolder("OUTPUT/MOM/Qz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Qz' exists" << endl;
        }
        if (Input::List().o_vNx) {
            if (Makefolder("OUTPUT/MOM/vNx") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/vNx' exists" << endl;
        }
        if (Input::List().o_vNy) {
            if (Makefolder("OUTPUT/MOM/vNy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/vNy' exists" << endl;
        }
        if (Input::List().o_vNz) {
            if (Makefolder("OUTPUT/MOM/vNz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/vNz' exists" << endl;
        }
        if (Input::List().o_Ux) {
            if (Makefolder("OUTPUT/MOM/Ux") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Ux' exists" << endl;
        }
        if (Input::List().o_Uy) {
            if (Makefolder("OUTPUT/MOM/Uy") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Uy' exists" << endl;
        }
        if (Input::List().o_Uz) {
            if (Makefolder("OUTPUT/MOM/Uz") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Ux' exists" << endl;
        }
        if (Input::List().o_Z) {
            if (Makefolder("OUTPUT/MOM/Z") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Z' exists" << endl;
        }
        if (Input::List().o_ni) {
            if (Makefolder("OUTPUT/MOM/ni") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/ni' exists" << endl;
        }
        if (Input::List().o_Ti) {
            if (Makefolder("OUTPUT/MOM/Ti") != 0)
                cout<<"Warning: Folder 'OUTPUT/MOM/Ti' exists" << endl;
        }
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//            if (Input::List().o_p1x1) { 
//                if (Makefolder("OUTPUT/MOM/p1x1") != 0) 
//                    cout<<"Warning: Folder 'OUTPUT/MOM/p1x1' exists\n";
//            }
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

    }

    if ( Input::List().o_f0x1 || Input::List().o_p1x1 
        ||  Input::List().o_f10x1 ||  Input::List().o_f11x1 ||  Input::List().o_f20x1
        || Input::List().o_fl0x1) {
        if (Makefolder("OUTPUT/DISTR") != 0)
            cout<<"Warning: Folder 'OUTPUT/DISTR' exists" << endl;

    }


}
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
// List of the acceptable Tags 
// Constructor 
//------------------------------------------------------------------------------
// List of the acceptable Tags Constructor
//
// @param[in]  species  The species
//
Export_Files::DefaultTags::DefaultTags(size_t species){

//  Time 
    time.push_back( "Time_cgs");
    time.push_back( "Time_si" );
    time.push_back( "Time_fs" );
    time.push_back( "Time_ps" );
    time.push_back( "Time"    ); // Default tag

//  Space
    space.push_back( "Space_cgs");
    space.push_back( "Space_si" );
    space.push_back( "Space");

//  Fields
    fld.push_back( "Ex"     );
    fld.push_back( "Ex_cgs" );
    fld.push_back( "Ex_si"  );
    fld.push_back( "Ey"     );
    fld.push_back( "Ey_cgs" );
    fld.push_back( "Ey_si"  );
    fld.push_back( "Ez"     );
    fld.push_back( "Ez_cgs" );
    fld.push_back( "Ez_si"  );
    fld.push_back( "Bx"     );
    fld.push_back( "Bx_cgs" );
    fld.push_back( "Bx_si"  );
    fld.push_back( "By"     );
    fld.push_back( "By_cgs" );
    fld.push_back( "By_si"  );
    fld.push_back( "Bz"     );
    fld.push_back( "Bz_cgs" );
    fld.push_back( "Bz_si"  );

//  Moments
    mom.push_back( "P"      );
    mom.push_back( "P_cgs"  );
    mom.push_back( "P_si"   );
    mom.push_back( "P_Mbar" );
    mom.push_back( "T"      );
    mom.push_back( "T_cgs"  );
    mom.push_back( "T_si"   );
    mom.push_back( "T_eV" );
    mom.push_back( "n"      );
    mom.push_back( "n_cgs"  );
    mom.push_back( "n_si"   );
    mom.push_back( "Qx"     );
    mom.push_back( "Qy"     );
    mom.push_back( "Qz"     );
    mom.push_back( "vNx"    );
    mom.push_back( "vNy"    );
    mom.push_back( "vNz"    );
    mom.push_back( "Jx"     );
    mom.push_back( "Jx_cgs" );
    mom.push_back( "Jx_si"  );
    mom.push_back( "Jy"     );
    mom.push_back( "Jy_cgs" );
    mom.push_back( "Jy_si"  );
    mom.push_back( "Jz"     );
    mom.push_back( "Jz_cgs" );
    mom.push_back( "Jz_si"  );

    mom.push_back( "Ux"      );
    mom.push_back( "Uy"      );
    mom.push_back( "Uz"      );
    mom.push_back( "Z"      );
    mom.push_back( "ni"      );
    mom.push_back( "Ti"      );




//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px-x");
        pvsx.push_back( "py-x");
        pvsx.push_back( "pz-x");
    }

//  f-x
    for (size_t s(0); s < species; ++s) {
        fvsx.push_back( "f0-x");//+stringify(s) );
        fvsx.push_back( "f10-x");
        fvsx.push_back( "f11-x");
        fvsx.push_back( "f20-x");
        fvsx.push_back( "fl0-x");
    }

//  p-p-x
    for (size_t s(0); s < species; ++s) {
        pvspvsx.push_back( "pxpy-x");
        pvspvsx.push_back( "pypz-x");
        pvspvsx.push_back( "pxpz-x");
    }

}



//--------------------------------------------------------------


// Definition of the output axis
// Constructor 
Export_Files::oAxis::oAxis() : label(""), units(""), min(0.0), max(1.0), sz(3) {}
Export_Files::oAxis::oAxis( const float _m, const float _M,
                            const size_t _sz)
        : label(""), units(""), min(_m), max(_M), sz(_sz) {}
Export_Files::oAxis::oAxis(const string _l, const string _u, const float _m, const float _M,
                           const size_t _sz)
        : label(_l), units(_u), min(_m), max(_M), sz(_sz) {}
// Copy constructor 
Export_Files::oAxis::oAxis(const oAxis& other) {
    label = other.label;
    units = other.units;
    min   = other.min;
    max   = other.max;
    sz    = other.sz;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
// 1D header constructor
Export_Files::Header::Header(oAxis _x,
                             string _Ql, float _Qc,
                             string _tl, string _tu, float _tc,
                             string _oD)
        : title(_Ql), titleC(_Qc),
          time(_tl),  timeU(_tu),  timeC(_tc),
          oDir(_oD) {
    xyz_axis.push_back(_x);
}

// 2D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y,
                             string _Ql, float _Qc,
                             string _tl, string _tu, float _tc,
                             string _oD)
        : title(_Ql),  time(_tl), timeU(_tu),
          titleC(_Qc), timeC(_tc),
          oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
}

// 3D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z,
                             string _Ql, float _Qc,
                             string _tl, string _tu, float _tc,
                             string _oD)
        : title(_Ql), time(_tl), timeU(_tu),
          titleC(_Qc), timeC(_tc),
          oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
    xyz_axis.push_back(_z);
}

// xD header constructor
Export_Files::Header::Header(vector< oAxis > _xyz,
                             string _Ql, float _Qc,
                             string _tl, string _tu, float _tc,
                             string _oD)
        : xyz_axis(_xyz), title(_Ql), time(_tl), timeU(_tu),
          titleC(_Qc), timeC(_tc),
          oDir(_oD) {}

// number of header dimensions
size_t Export_Files::Header::dim() {
    return xyz_axis.size();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

valarray<float> Export_Files::Header::axis(const size_t i) {
//    return Algorithms::MakeAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz);
    return Algorithms::MakeCAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz);
}
string Export_Files::Header::label(const size_t i) {
    return xyz_axis[i].label;
}
string Export_Files::Header::units(const size_t i) {
    return xyz_axis[i].units;
}

string Export_Files::Header::Title_label() { return title; }
float Export_Files::Header::Title_conv()  { return titleC; }
string Export_Files::Header::Time_label()  { return time; }
float Export_Files::Header::Time_conv()   { return timeC; }
string Export_Files::Header::Directory()   { return oDir; }
//--------------------------------------------------------------


//--------------------------------------------------------------
// Constructor of the export facility for data structures
Export_Files::Xport::Xport(const Algorithms::AxisBundle<double>& _axis,
                           const vector< string > oTags,
                           string homedir){

    size_t species(_axis.pdim());
    DefaultTags dTags(species);

    vector< oAxis > xyz, pxyz, imre, pr;
    xyz.push_back( Export_Files::oAxis(_axis.xgmin(0), _axis.xgmax(0), _axis.Nxg(0)) );

    for (size_t s(0); s < species; ++s) {
        pr.push_back( Export_Files::oAxis(_axis.pmin(s), _axis.pmax(s), _axis.Np(s)));
        pxyz.push_back( Export_Files::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npx(s)) );
        pxyz.push_back( Export_Files::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npy(s)) );
        pxyz.push_back( Export_Files::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npz(s)) );
    }
    imre.push_back( Export_Files::oAxis(0,1,2) );
//  Time 
    size_t tloc(0);                  // Find the location of the right tag
    while ( ( tloc < dTags.time.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.time[tloc]) == oTags.end() ) ) {
        ++tloc;
    }
    string tlabel = "t[" +formulary().Label(dTags.time[tloc])+"]";
    string tunits = formulary().Label(dTags.time[tloc]);
    float  tconv  =  formulary().Uconv(dTags.time[tloc]);

//  xyz Axis
    size_t xloc(0);                  // Find the location of the right tag
    while ( ( xloc < dTags.space.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.space[xloc]) == oTags.end() ) ) {
        ++xloc;
    }
    xyz[0].label = "x["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 1 ) xyz[1].label = "y["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 2 ) xyz[2].label = "z["+ formulary().Label(dTags.space[xloc]) +"]";
    xyz[0].units = formulary().Label(dTags.space[xloc]);
    if ( xyz.size() > 1 ) xyz[1].units = formulary().Label(dTags.space[xloc]);
    if ( xyz.size() > 2 ) xyz[2].units = formulary().Label(dTags.space[xloc]);


//  pxyz Axis
    for (size_t i(0); i < species; ++i) {
        pr[i].label = "p";
        pr[i].units = "m_e c";

        pxyz[i].label = "px[mc]";
        pxyz[i].units = "m_e c";
        // if ( pxyz.size()/species > 1 ) 
        pxyz[1+i].label = "py[mc]";
        pxyz[1+i].units = "m_e c";
        // if ( pxyz.size()/species > 2 )
        pxyz[2+i].label = "pz[mc]";
        pxyz[2+i].units = "m_e c";

        // if ( pxyz.size()/species > 1 ) 
        // if ( pxyz.size()/species > 2 ) pxyz[2*species+i].units = "m_e c";
    }

//  Tags for Fields -->
    for (size_t i(0); i < dTags.fld.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fld[i]) != oTags.end() ) {

            string nounits = dTags.fld[i].substr(0, dTags.fld[i].find("_"));
            string folder = homedir + "OUTPUT/FLD/" + nounits + "/";
            //         Generate a header file for this tag
            Hdr[dTags.fld[i]] = Header(xyz,
                                       nounits+"["+formulary().Label(dTags.fld[i])+"]",
                                       formulary().Uconv(dTags.fld[i]),
                                       tlabel, tunits, tconv, folder);
        }
    } // <--

//  Tags for Moments -->
    for (size_t i(0); i < dTags.mom.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.mom[i]) != oTags.end() ) {

            string nounits = dTags.mom[i].substr(0, dTags.mom[i].find("_"));
            string folder = homedir + "OUTPUT/MOM/" + nounits + "/";
            //         Generate a header file for this tag
            Hdr[dTags.mom[i]] = Header(xyz,
                                       nounits+"["+formulary().Label(dTags.mom[i])+"]",
                                       formulary().Uconv(dTags.mom[i]),
                                       tlabel, tunits, tconv, folder);
        }
    } // <--

//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx.size(); ++i) {

        // for (size_t k(0); k < oTags.size(); ++k)
        //         std::cout << " \n \n k = " << oTags[k];

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.pvsx[i]) != oTags.end() ) {

            string folder = homedir + "OUTPUT/DISTR/" + dTags.pvsx[i] + "/";
            Makefolder(folder);


//          Generate a header file for this tag
//          For each 9 you have a different species 
            Hdr[dTags.pvsx[i]] = Header( pr[i/3], xyz[0],
                                         "f"+stringify(i/3), 1.0, tlabel, tunits, tconv, folder);
            // std::cout << " \n k = " << dTags.pvsx[i] << "\n";
            // std::cout << " \n \n 11 ";
        }
    } //<--

//  Tags for f-x -->
    for (size_t i(0); i < dTags.fvsx.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fvsx[i]) != oTags.end() ) {

            string folder = homedir + "OUTPUT/DISTR/" + dTags.fvsx[i] + "/";
            Makefolder(folder);

//          Generate a header file for this tag
            Hdr[dTags.fvsx[i]] = Header( pr[i/5], xyz[0], imre[0],
                                         "f", 1.0, tlabel, tunits, tconv, folder);

        }
    } //<--

    //  Tags for p-p-x -->
    for (size_t i(0); i < dTags.pvspvsx.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.pvspvsx[i]) != oTags.end() ) {

            string folder = homedir + "OUTPUT/DISTR/" + dTags.pvspvsx[i] + "/";
            Makefolder(folder);

//          Generate a header file for this tag
            Hdr[dTags.pvspvsx[i]] = Header( pxyz[i/1+species],pxyz[(i+1)/1+species], xyz[0],
                                            "f"+stringify(i/1), 1.0, tlabel, tunits, tconv, folder);

            // std::cout << " \n \n 11 ";
        }
    } //<--
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//--------------------------------------------------------------
//  Adjust H5 filenames with zeros to reach some prescribed length. 
//  Add the H5 filename extension.  
string Export_Files::Xport::oH5Fextension(size_t step, int species){

    stringstream sFilename;

    if(species >= 0) sFilename << "_s" << species;

    sFilename << "_";

    // Number of zeros to add to the filename
    int Nzeros(ofconventions::ofile_digits - stringify(step).length());
    while (Nzeros-- > 0) {
        sFilename << "0";
    }

    sFilename << step << ofconventions::h5file_extension;

    return sFilename.str();
}
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
Export_Files::Restart_Facility::Restart_Facility(const int rank, string homedir) {
    hdir = homedir;

    if (!rank) Makefolder(hdir+"RESTART/");

//        if (Makefolder(hdir+"RESTART/") != 0) cout<<"Warning: Folder "<< hdir+"RESTART/"<<" exists\n";
}

//--------------------------------------------------------------
//  Read restart file
void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, State1D& Y) {

//      Generate filename 
    string   filename(hdir+"RESTART/re_1D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ifstream  fin(filename.c_str(), ios::binary);

    if (fin)
    {
    //      Read distribution functions
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
                for(size_t i(0); i < (Y.DF(s))(nh).dim(); ++i) {
                    fin.read((char *)&(Y.DF(s))(nh)(i), sizeof((Y.DF(s))(nh)(i)));
                }
            }
        }
    }
    else {
        
        if (!rank) std::cout << "\n\n ERROR :: No files to read! \n\n";
        exit(1);

    }

//      Read Ex
    for(size_t i(0); i < Y.EMF().Ex().numx(); ++i){
        fin.read((char *)(&Y.EMF().Ex()(i)), sizeof(Y.EMF().Ex()(i)));
    }

    fin.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Write restart file
void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, State1D& Y) {

//      Generate filename 
    string   filename(hdir+"RESTART/re_1D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ofstream  fout(filename.c_str(), ios::binary);

//      Write distribution functions
    for(size_t s(0); s < Y.Species(); ++s) {
        for (size_t s(0); s < Y.Species(); ++s) {
            for (size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
                for (size_t i(0); i < (Y.DF(s))(nh).dim(); ++i) {
                    fout.write((char *) &(Y.DF(s))(nh)(i), sizeof((Y.DF(s))(nh)(i)));
                }
            }
        }
    }

//      Write Ex
    for(size_t i(0); i < Y.EMF().Ex().numx(); ++i) {
        fout.write((char *)(&Y.EMF().Ex()(i)), sizeof(Y.EMF().Ex()(i)));
    }

    fout.flush();
    fout.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension. 
string Export_Files::Restart_Facility::rFextension(const int rank, const size_t rstep){
    stringstream sFilename;

    // Number of zeros to add to the filename
    int Nzeros(ofconventions::rank_digits - stringify(rank).length());
    while (Nzeros-- > 0) {
        sFilename << "0";
    }

    sFilename << rank << "_";

    // Number of zeros to add to the filename
    Nzeros = ofconventions::rfile_digits - stringify(rstep).length();
    while (Nzeros-- > 0) {
        sFilename << "0";
    }

    sFilename << rstep << ofconventions::rfile_extension;

    return sFilename.str();
}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Output_Data::PLegendre1D::PLegendre1D( size_t Nl, size_t Nm, size_t Np,
                                       float pmin, float pmax, size_t Npx ) {

    size_t sz(0);

    valarray<float> p( Algorithms::MakeAxis(pmin, pmax, Np)),
            px(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npx));

//  Generate the structure to save the polynomials
    sz = ((Nm+1)*(2*Nl-Nm+2))/2;
//  
    plegendre = new vector< Array2D<float>>(sz,Array2D<float>(Npx,Np))  ;
    //    for (size_t m(0); m < Nm+1; ++m) {
    //     for (size_t l(m); l < Nl+1; ++l) { 
    //         (*plegendre).push_back( Array2D<float>(Npx,Np) );
    //         ++sz;
    //  }
    // }
    size_t k(0);
//  Generate the polynomial values for each cos(theta) = px/p 
    for (size_t j(0); j < p.size(); ++j) {
        float invp(1.0/p[j]);
        for (size_t i(0); i < px.size(); ++i) {

//          For given px/p generate all the polynomial values ...
//             valarray<float> vL( Algorithms::Legendre( px[i]*invp, Nl, Nm ) );
            Array2D<float> vL( Algorithms::Legendre( px[i]*invp, Nl, Nm ) );
            // std::cout << "\n x = " << px[i]*invp << ", sqrt(1-x^2) = " << sqrt(1.0-px[i]*invp*px[i]*invp) ;
//          ... and save them
            // std::cout << "dim = " << vL.dim() << "\n";
            // std::cout << "dim1:" << vL.dim1() << "\n";
            // std::cout << "dim2:" << vL.dim1() << "\n\n";
            for (size_t l(0); l < Nl+1; ++l){
                for (size_t m=0; m<((Nm<l)?Nm:l)+1; ++m){
                    // for (size_t k(0); k < sz; ++k) {
                    k = ((l < Nm+1)?((l*(l+1))/2+m):(l*(Nm+1)-(Nm*(Nm+1))/2 + m));

                    (*plegendre)[k](i,j) = vL(l,m);

                    std::cout << "\n LP(" << l << "," << m << "," << k << ") = " <<  vL(l,m);

                }
                // (*plegendre)[k](i,j) = vL(k);
            }

            // exit(1);
        }
    }

}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre1D::PLegendre1D( const PLegendre1D& other ) {

//  Generate the structure to save the polynomials
    plegendre = new vector< Array2D<float> > ;
    for (size_t i(0); i < other.dim(); ++i) {
        (*plegendre).push_back( other(i) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre1D::~PLegendre1D(){
    delete plegendre;
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Output_Data::PLegendre2D::PLegendre2D( size_t Nl, size_t Nm, size_t Np,
                                       float pmin, float pmax, size_t Npx, size_t Npy ) {

    size_t sz(0);

    valarray<float> p( Algorithms::MakeAxis(pmin, pmax, Np)),
            px(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npx)),
            py(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npy));

//  Generate the structure to save the polynomials
    sz = ((Nm+1)*(2*Nl-Nm+2))/2;
//  
    plegendre = new vector< Array2D<float>>(sz,Array2D<float>(Npx,Npy))  ;
    //    for (size_t m(0); m < Nm+1; ++m) {
    //     for (size_t l(m); l < Nl+1; ++l) { 
    //         (*plegendre).push_back( Array2D<float>(Npx,Np) );
    //         ++sz;
    //  }
    // }
    size_t k(0);
    // float invpmax(1.0/pmax);
//  Generate the polynomial values for each cos(theta) = px/p 
    for (size_t j(0); j < py.size(); ++j) {

        for (size_t i(0); i < px.size(); ++i) {
            float invp(sqrt(py[j]*py[j]+px[i]+px[i]));
            if (invp > 0.0 ) invp=1.0/invp;
            else if (invp > pmax) invp *= 2;
            else invp = 0.0;
//          For given px/p generate all the polynomial values ...
//             valarray<float> vL( Algorithms::Legendre( px[i]*invp, Nl, Nm ) );

            Array2D<float> vL( Algorithms::Legendre( px[i]*invp, Nl, Nm ) );
            // std::cout << "\n x = " << px[i]*invp << ", sqrt(1-x^2) = " << sqrt(1.0-px[i]*invp*px[i]*invp) ;
            //          ... and save them
            // std::cout << "dim = " << vL.dim() << "\n";
            // std::cout << "dim1:" << vL.dim1() << "\n";
            // std::cout << "dim2:" << vL.dim1() << "\n\n";
            for (size_t l(0); l < Nl+1; ++l){
                for (size_t m(0); m<((Nm<l)?Nm:l)+1; ++m){
                    // for (size_t k(0); k < sz; ++k) {
                    k = ((l < Nm+1)?((l*(l+1))/2+m):(l*(Nm+1)-(Nm*(Nm+1))/2 + m));

                    (*plegendre)[k](i,j) = vL(l,m);
                    // if (invp == 0.0) (*plegendre)[k](i,j) = 0.0;

                    // std::cout << "\n LP(" << l << "," << m << "," << k << ") = " <<  vL(l,m);

                }
                // (*plegendre)[k](i,j) = vL(k);
            }
        }
        // exit(1);
    }

}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre2D::PLegendre2D( const PLegendre2D& other ) {

//  Generate the structure to save the polynomials
    plegendre = new vector< Array2D<float> > ;
    for (size_t i(0); i < other.dim(); ++i) {
        (*plegendre).push_back( other(i) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre2D::~PLegendre2D(){
    delete plegendre;
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Output_Data::p1x1_1D::p1x1_1D( const Grid_Info& _G) {

//  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pmax.push_back( static_cast<float>(_G.axis.pmax(s)) );
        pmin.push_back( static_cast<float>(-1.0*_G.axis.pmax(s)) );
//        Pl.push_back( PLegendre1D( _G.l0[s], _G.m0[s], _G.axis.Np(s), static_cast<float>(_G.axis.pmin(s)),
//                                   static_cast<float>(_G.axis.pmax(s)), _G.axis.Npx(s) ) );
//
//        polarR.push_back( Array2D<float>( _G.axis.Npx(s), _G.axis.Np(s) ) );

        p1x1.push_back( valarray<float>( _G.axis.Npx(s)) );

    }

/*//  Calculate a 2D array for the polar radius for each px and pr
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        valarray<float> p( Algorithms::MakeAxis( pmin[s], pmax[s], _G.axis.Np(s) )),
                px(Algorithms::MakeAxis( float(-1.0)* pmax[s], pmax[s], _G.axis.Npx(s) ));
//      --->
        for (size_t i(0); i < _G.axis.Npx(s); ++i){
            for (size_t j(0); j < _G.axis.Np(s); ++j){
                float polarR_sq  = p[j] * p[j] - px[i] * px[i];
                if (polarR_sq < 0.0) {
                    polarR[s](i,j) = -1.0;
                }
                else {
                    polarR[s](i,j)  = sqrt(abs(polarR_sq));
                }
                // std::cout << "\n polarR[s](" << px[i] << ',' << p[j] << ") = " << polarR[s](i,j) << "\n";
            }
        }
//      <---
    }*/
//    p1x1 = new valarray<float>(Npx);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p1x1_1D::p1x1_1D( const p1x1_1D& other) {

    for (size_t s(0); s < other.Species(); ++s) {
        pmin.push_back( other.Pmin(s) );
        pmax.push_back( other.Pmax(s) );
//        Pl.push_back( other.PL(s) );
//        polarR.push_back( other.PolarR(s) );
        p1x1.push_back( other.p1_x1(s) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p1x1_1D::~p1x1_1D(){
    // delete p1x1;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location x0 
//  into a cartesian grid.
valarray<float>  Output_Data::p1x1_1D::operator()(DistFunc1D& df, size_t x0, size_t s) {

    p1x1[s] = 0.0; // this is a valarray<float>
    vector<complex<double> > dummyvec;

/*    float p_im1(0.0), p_ip1(0.0);
    float integrant_low(0.0),
            integrant_high(0.0);

//  Calculate the integral for each harmonic separately
//     for (size_t l(0); l < Nl(s); ++l) { 
    for (size_t i(0); i < (Pl[s]).dim(); ++i) {

        valarray<float> InPx( Npx(s) );
        for (size_t ipx(0); ipx < Npx(s)-1; ++ipx) { // at each location in px

// //          At each location in |p|
            size_t ip(0);
            while (polarR[s](ipx,ip) < 0 ) { ++ip; }

            p_ip1 = polarR[s](ipx,ip);
            integrant_high =(float( (df(i))(ip,x0).real() ))*Pl[s](i)(ipx, ip) ;



            InPx[ipx] += p_ip1 * (0.5*p_ip1) * integrant_high;
            ++ip;
            while ( (ip < Np(s) ) && (polarR[s](ipx,ip) > 0)  )   {
                p_im1 = p_ip1;
                integrant_low = integrant_high;

                p_ip1 = polarR[s](ipx,ip);
                integrant_high =(float( (df(i))(ip,x0).real() )) * Pl[s](i)(ipx, ip);

                InPx[ipx] += p_im1 * (0.5*(p_ip1-p_im1)) * integrant_low;
                InPx[ipx] += p_ip1 * (0.5*(p_ip1-p_im1)) * integrant_high;

                ++ip;
            }
            // if (i==3) std::cout << "\n p1x1[s](" << ipx << ',' << ip << ") = " << InPx[ipx] << "\n";
            (p1x1[s])[ipx] = integrant_high;
        }
    }*/

    int sgn(1);
    size_t midpoint,ipxp,ipxm;

    for (size_t i(0); i < df.dim(); ++i)
    {
        dummyvec = df(i).xVec(x0);

        midpoint = df(0,0).nump();
        p1x1[s][midpoint] = (float) ((dummyvec[0]).real());
        for (size_t ip(0); ip < df(0,0).nump(); ++ip)
        {
            ipxm = midpoint - 1 - ip;
            ipxp = midpoint + 1 + ip;

            (p1x1[s])[ipxm] += (float) (sgn*(dummyvec[ip]).real());
            (p1x1[s])[ipxp] += (float) ((dummyvec[ip]).real());
        }
        sgn *= -1;
    }

//    p1x1[s] *= 2.0 * M_PI;
    // exit(1);
    return p1x1[s];
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Output_Data::fx1_1D::fx1_1D(const Grid_Info& _G) {

//  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pmax.push_back( static_cast<float>(_G.axis.pmax(s)) );
        pmin.push_back( static_cast<float>(_G.axis.pmin(s)) );
        nump.push_back( static_cast<float>(_G.axis.Np(s)));
        // f0x1.push_back( valarray<float>( _G.axis.Npx(s)) );
    }
//    f0x1 = new valarray<float>(Npx);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::fx1_1D::fx1_1D( const fx1_1D& other) {

    for (size_t s(0); s < other.Species(); ++s) {
        pmin.push_back( other.Pmin(s) );
        pmax.push_back( other.Pmax(s) );
        nump.push_back( other.Np(s));
        // fx1.push_back( other.f_x1(s) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::fx1_1D::~fx1_1D(){

}
//-------------------------------------------------------------
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location x0 
//  into a cartesian grid.
Array2D<float>  Output_Data::fx1_1D::operator()(DistFunc1D& df, size_t l, size_t m, size_t x0, size_t s) {

    // valarray<float> fx1(0.0,Npx(s)); // this is a valarray<float>
    Array2D<float> fx1(Np(s),2); // this is a valarray<float>// 

    for (size_t ip(0); ip < Np(s); ++ip) {
        //  At each location in |p|
        // fx1[0.5*(Npx(s)-1)+ip] = (float( (df(whichdist))(ip,x0).real() ));
        // fx1[0.5*(Npx(s)-1)-ip] = (float( (df(whichdist))(ip,x0).real() ));

        // fx1(0.5*(Npx(s)-1)+1+ip,1) = (float( (df(whichdist))(ip,x0).real() ));
        // fx1(0.5*(Npx(s)-1)-1-ip,1) = (float( (df(whichdist))(ip,x0).real() ));

        // fx1(0.5*(Npx(s)-1)+1+ip,2) = (float( (df(whichdist))(ip,x0).imag() ));
        // fx1(0.5*(Npx(s)-1)-1-ip,2) = (float( (df(whichdist))(ip,x0).imag() ));
        // 
        fx1(ip,0) = (float( (df(l,m))(ip,x0).real() ));
        fx1(ip,1) = (float( (df(l,m))(ip,x0).imag() ));

        // std::cout << "fx1(" << ip << ") = " << (df(whichdist))(ip,x0).real() << "\n";
        // << fx1(0.5*(Npx(s)-1)+1+ip,1) << "," << fx1(0.5*(Npx(s)-1)+1+ip,2) << "\n";
    }
    // fx1[Npx(s)] = (float( (df(whichdist))(0,x0).real() ));
    // fx1(,1) = (float( (df(whichdist))(0,x0).real() ));
    // fx1(0.5*(Npx(s)-1),2) = (float( (df(whichdist))(0,x0).imag() ));
    return fx1;
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Output_Data::p2p1x1_1D::p2p1x1_1D( const Grid_Info& _G) {

//  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pmax.push_back( static_cast<float>(_G.axis.pmax(s)) );
        pmin.push_back( static_cast<float>(_G.axis.pmin(s)) );
        numl.push_back( _G.l0[s]);
        numm.push_back( _G.m0[s]);

        Pl.push_back( PLegendre2D( _G.l0[s], _G.m0[s], _G.axis.Np(s), static_cast<float>(_G.axis.pmin(s)),
                                   static_cast<float>(_G.axis.pmax(s)), _G.axis.Npx(s), _G.axis.Npy(s) ) );
        pr.push_back( Array2D<float>( _G.axis.Npx(s),_G.axis.Npy(s)));


        nextpcell.push_back( Array2D<size_t>( _G.axis.Npx(s),_G.axis.Npy(s)));
        distancetothatcell.push_back( Array2D<float>( _G.axis.Npx(s),_G.axis.Npy(s)));


        p2p1x1.push_back( Array2D<float>( _G.axis.Npx(s),_G.axis.Npy(s)) );
    }

    size_t k(0);
    //  Calculate a 2D array for the pr for each px and py   
    for (size_t s(0); s < _G.axis.pdim(); ++s) {

        valarray<float> p( Algorithms::MakeAxis( pmin[s], pmax[s], _G.axis.Np(s) )),
                px(Algorithms::MakeAxis( float(-1.0)* pmax[s], pmax[s], _G.axis.Npx(s) )),
                py(Algorithms::MakeAxis( float(-1.0)* pmax[s], pmax[s], _G.axis.Npy(s) ));
//      --->
        for (size_t i(0); i < _G.axis.Npx(s); ++i){
            for (size_t j(0); j < _G.axis.Npy(s); ++j){
                pr[s](i,j) = sqrt(px[i]*px[i]+py[j]*py[j]);

                while ((pr[s](i,j) > p[k]) && (k < _G.axis.Np(s))) ++k;
                if (k == 0) ++k;
                nextpcell[s](i,j) = k-1;
                distancetothatcell[s](i,j) = (p[k] - pr[s](i,j))/(p[1]-p[0]);
                // std::cout << "\n\n (px,py,pr,p[k],nextpcell(i,j)) = " << px[i]
                // << ", " << py[j] 
                // << ", " << pr[s](i,j) 
                // << ", " << p[k] 
                // << ", " << nextpcell[s](i,j);

                k = 0;

            }
        }
//      <---
    }


//    f0x1 = new valarray<float>(Npx);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p2p1x1_1D::p2p1x1_1D( const p2p1x1_1D& other) {

    for (size_t s(0); s < other.Species(); ++s) {
        pmin.push_back( other.Pmin(s) );
        pmax.push_back( other.Pmax(s) );
        numl.push_back( other.Nl(s));
        numm.push_back( other.Nm(s));
        Pl.push_back(other.PL(s));
        pr.push_back(other.prad(s));
        nextpcell.push_back(other.npcell(s));
        distancetothatcell.push_back(other.dcell(s));
        p2p1x1.push_back( other.p2p1_x1(s) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p2p1x1_1D::~p2p1x1_1D(){

}
//-------------------------------------------------------------
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location x0 
//  into a cartesian grid.
Array2D<float>  Output_Data::p2p1x1_1D::operator()(DistFunc1D& df, size_t x0, size_t s) {

    p2p1x1[s] = 0.0; // this is a Array2D<float>
    // std::cout << "\n\n x0 = " << x0 << " \n";
    size_t i(0), temploc(0);
    float tempdist(0.0);

    // valarray<float> p( Algorithms::MakeAxis(pmin, pmax, Np));

    for (size_t l(0); l < Nl(s)+1; ++l){
        for (size_t m(0); m<((Nm(s)<l)?Nm(s):l)+1; ++m){

            i = ((l < Nm(s)+1)?((l*(l+1))/2+m):(l*(Nm(s)+1)-(Nm(s)*(Nm(s)+1))/2 + m));

            for (size_t ipy(0); ipy < Npy(s); ++ipy) {
                for (size_t ipx(0); ipx < Npx(s); ++ipx) {
                    //  At each location in px,py
                    if (nextpcell[s](ipx,ipy) == 0)
                    {
                        temploc = 1;
                        tempdist = 1.0;
                    }
                    else
                    {
                        temploc = nextpcell[s](ipx,ipy);
                        tempdist = distancetothatcell[s](ipx,ipy);
                    }

                    p2p1x1[s](ipx,ipy) += Pl[s](i)(ipx, ipy) *
                                          ((float((df(i))(temploc  ,x0).real())*(1.0-tempdist))+
                                           (float((df(i))(temploc-1,x0).real())*(    tempdist)));
                }
            }
        }
    }
    return p2p1x1[s];
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor_1D::operator()(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                                     const Parallel_Environment_1D& PE) {

    if (Input::List().o_Ex) {
        Ex( Y, grid, tout, PE );
    }

    if (Input::List().o_Ey) {
        Ey( Y, grid, tout, PE );
    }

    if (Input::List().o_Ez) {
        Ez( Y, grid, tout, PE );
    }

    if (Input::List().o_Bx) {
        Bx( Y, grid, tout, PE );
    }

    if (Input::List().o_By) {
        By( Y, grid, tout, PE );
    }

    if (Input::List().o_Bz) {
        Bz( Y, grid, tout, PE );
    }

    if (Input::List().o_x1x2) {
        n( Y, grid, tout, PE );
    }

    if (Input::List().o_Temperature) {
        T( Y, grid, tout, PE );
    }

    if (Input::List().o_Jx) {
        Jx( Y, grid, tout, PE );
    }

    if (Input::List().o_Jy) {
        Jy( Y, grid, tout, PE );
    }

    if (Input::List().o_Jz) {
        Jz( Y, grid, tout, PE );
    }

    if (Input::List().o_Qx) {
        Qx( Y, grid, tout, PE );
    }

    if (Input::List().o_Qy) {
        Qy( Y, grid, tout, PE );
    }

    if (Input::List().o_Qz) {
        Qz( Y, grid, tout, PE );
    }
    if (Input::List().o_vNx) {
        vNx( Y, grid, tout, PE );
    }
    if (Input::List().o_vNy) {
        vNy( Y, grid, tout, PE );
    }
    if (Input::List().o_vNz) {
        vNz( Y, grid, tout, PE );
    }

//    if (Input::List().o_p1x1){
//        px( Y, grid, tout, PE );
//    }
//    if (Input::List().o_f0x1){
//        f0( Y, grid, tout, PE );
//    }
//    if (Input::List().o_f10x1){
//        f10( Y, grid, tout, PE );
//    }
//    if (Input::List().o_f11x1){
//        f11( Y, grid, tout, PE );
//    }
//    if (Input::List().o_p2p1x1){
//        pxpy( Y, grid, tout, PE );
//    }

    if (Input::List().o_Ux) {
        Ux( Y, grid, tout, PE );
    }

    if (Input::List().o_Uy) {
        Uy( Y, grid, tout, PE );
    }

    if (Input::List().o_Uz) {
        Uz( Y, grid, tout, PE );
    }

    if (Input::List().o_Z) {
        Z( Y, grid, tout, PE );
    }

    if (Input::List().o_ni) {
        ni( Y, grid, tout, PE );
    }

    if (Input::List().o_Ux) {
        Ti( Y, grid, tout, PE );
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor_1D::distdump(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                                   const Parallel_Environment_1D& PE) {

    if (Input::List().o_p1x1){
        px( Y, grid, tout, PE );
    }
    if (Input::List().o_f0x1){
        f0( Y, grid, tout, PE );
    }
    if (Input::List().o_f10x1){
        f10( Y, grid, tout, PE );
    }
    if (Input::List().o_f11x1){
        f11( Y, grid, tout, PE );
    }
    if (Input::List().o_f20x1){
        f20( Y, grid, tout, PE );
    }
    if (Input::List().o_fl0x1){
        fl0( Y, grid, tout, PE );
    }
//    if (Input::List().o_p2p1x1){
//        pxpy( Y, grid, tout, PE );
//    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//     void Parallel_Output::parallel_Ex(size_t step, State1D& Y) {
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Ex(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Exbuf = new float[msg_sz];
    valarray<float> ExGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Exbuf[i] = static_cast<float>( Y.EMF().Ex()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Exbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ExGlobal[i] = Exbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Exbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    ExGlobal[i + outNxLocal*rr] = Exbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            ExGlobal[i] = Exbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ex", ExGlobal, tout);

    delete[] Exbuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Ey(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Eybuf = new float[msg_sz];
    valarray<float> EyGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Eybuf[i] = static_cast<float>( Y.EMF().Ey()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Eybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EyGlobal[i] = Eybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Eybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    EyGlobal[i + outNxLocal*rr] = Eybuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            EyGlobal[i] = Eybuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ey", EyGlobal, tout);

    delete[] Eybuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Ez(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Ezbuf = new float[msg_sz];
    valarray<float> EzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Ezbuf[i] = static_cast<float>( Y.EMF().Ez()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Ezbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EzGlobal[i] = Ezbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Ezbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    EzGlobal[i + outNxLocal*rr] = Ezbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            EzGlobal[i] = Ezbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ez", EzGlobal, tout);

    delete[] Ezbuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Bx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Bxbuf = new float[msg_sz];
    valarray<float> BxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bxbuf[i] = static_cast<float>( Y.EMF().Bx()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BxGlobal[i] = Bxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Bxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    BxGlobal[i + outNxLocal*rr] = Bxbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            BxGlobal[i] = Bxbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Bx", BxGlobal, tout);

    delete[] Bxbuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::By(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Bybuf = new float[msg_sz];
    valarray<float> ByGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bybuf[i] = static_cast<float>( Y.EMF().By()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ByGlobal[i] = Bybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Bybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    ByGlobal[i + outNxLocal*rr] = Bybuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            ByGlobal[i] = Bybuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("By", ByGlobal, tout);

    delete[] Bybuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Bz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Bzbuf = new float[msg_sz];
    valarray<float> BzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bzbuf[i] = static_cast<float>( Y.EMF().Bz()(Nbc+i).real() );
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BzGlobal[i] = Bzbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Bzbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    BzGlobal[i + outNxLocal*rr] = Bzbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            BzGlobal[i] = Bzbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Bz", BzGlobal, tout);

    delete[] Bzbuf;

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::px(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));



    for(int s(0); s < Y.Species(); ++s) {
        size_t Npx(2*f_x.Np(s)+1);
        int msg_sz(outNxLocal*Npx);
        Array2D<float> p1x1Global(Npx,outNxGlobal); //, yglob_axis.dim());
        float* pxbuf = new float[Npx*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) {

            valarray<float> data1D = px_x( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npx; ++j) {
                pxbuf[j+i*Npx]=data1D[j];
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    for (size_t j(0); j < Npx; ++j) {
                        p1x1Global(j,i) = pxbuf[j+i*Npx];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(pxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < Npx; ++j) {
                            p1x1Global(j,i + outNxLocal*rr) = pxbuf[j+i*Npx];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                for (size_t j(0); j < Npx; ++j) {
                    p1x1Global(j,i) = pxbuf[j+i*Npx];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("px-x", p1x1Global, tout, s);

        delete[] pxbuf;
    }

}
//--------------------------------------------------------------
////--------------------------------------------------------------
//void Output_Data::Output_Preprocessor_1D::pxpy(const State1D& Y, const Grid_Info& grid, const size_t tout,
//                                               const Parallel_Environment_1D& PE) {
//
//    size_t Nbc = Input::List().BoundaryCells;
//    MPI_Status status;
//    size_t st(0), bi(0);
//    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
//    size_t outNxGlobal(grid.axis.Nxg(0));
//    size_t ind(0);
//    for(int s(0); s < Y.Species(); ++s) {
//
//        int msg_sz(outNxLocal*pxpy_x.Npx(s)*pxpy_x.Npy(s));
//        Array3D<float> pxpyGlobal(pxpy_x.Npx(s),pxpy_x.Npy(s),outNxGlobal); //, yglob_axis.dim());
//        float* pbuf = new float[pxpy_x.Npx(s)*pxpy_x.Npy(s)*outNxGlobal];
//
//        for (size_t i(0); i < outNxLocal; ++i) {
//
//            Array2D<float> data2D = pxpy_x( Y.DF(s), i+Nbc, s);
//
//            for (size_t j(0); j < pxpy_x.Npx(s); ++j) {
//                for (size_t k(0); k < pxpy_x.Npy(s); ++k) {
//                    pbuf[ind+i*(pxpy_x.Npx(s)*pxpy_x.Npy(s))]=data2D(j,k);
//                    ++ind;
//                }
//            }
//            ind = 0;
//        }
//
//        if (PE.NODES() > 1) {
//            if (PE.RANK()!=0) {
//                MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
//            }
//            else {
//                // Fill data for rank = 0
//                for(size_t i(0); i < outNxLocal; i++) {
//                    ind = 0;
//                    for (size_t j(0); j < pxpy_x.Npx(s); ++j) {
//                        for (size_t k(0); k < pxpy_x.Npy(s); ++k) {
//                            pxpyGlobal(k,j,i) = pbuf[ind+i*(pxpy_x.Npx(s)*pxpy_x.Npy(s))];
//                            ++ind;
//                        }
//                    }
//                }
//                // Fill data for rank > 0
//                for (int rr = 1; rr < PE.NODES(); ++rr){
//                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
//                    for(size_t i(0); i < outNxLocal; i++) {
//                        ind = 0;
//                        for (size_t j(0); j < pxpy_x.Npx(s); ++j) {
//                            for (size_t k(0); k < pxpy_x.Npy(s); ++k) {
//                                pxpyGlobal(k,j,i + outNxLocal*rr) = pbuf[ind+i*(pxpy_x.Npx(s)*pxpy_x.Npy(s))];
//                                // p1x1Global(j,i + outNxLocal*rr) = pxbuf[j+i*f_x.Np(s)];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        else {
//            for(size_t i(0); i < outNxLocal; i++) {
//                ind = 0;
//                for (size_t j(0); j < pxpy_x.Npx(s); ++j) {
//                    for (size_t k(0); k < pxpy_x.Npy(s); ++k) {
//                        pxpyGlobal(k,j,i) = pbuf[ind+i*(pxpy_x.Npx(s)*pxpy_x.Npy(s))];
//                        ++ind;
//                    }
//                }
//            }
//        }
//
//        if (PE.RANK() == 0) expo.Export_h5("pxpy-x", pxpyGlobal, tout, s);
//
//        delete[] pbuf;
//    }
//
//}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::f0(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float* f0xbuf = new float[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s),0,0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(j,i + outNxLocal*rr,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(j,i + outNxLocal*rr,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f0-x", f0x1Global, tout, s);

        delete[] f0xbuf;
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------  
//--------------------------------------------------------------
// void Output_Data::Output_Preprocessor_1D::f0(const State1D& Y, const Grid_Info& grid, const size_t tout,
//                                                             const Parallel_Environment_1D& PE) { 
//         // std::cout << "0 \n";
//         size_t Nbc = Input::List().BoundaryCells;
//         MPI_Status status;

//         size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
//         size_t outNxGlobal(grid.axis.Nxg(0));

//         for(int s(0); s < Y.Species(); ++s) {     

//             int msg_sz(outNxLocal*f_x.Npx(s)); 
//             Array2D<float> f0x1Global(f_x.Npx(s),outNxGlobal); //, yglob_axis.dim());
//             float* f0xbuf = new float[f_x.Npx(s)*outNxLocal];

//             for (size_t i(0); i < outNxLocal; ++i) {

//                 valarray<float> data1D = f_x( Y.DF(s),0, i+Nbc, s);

//                 for (size_t j(0); j < f_x.Npx(s); ++j) {
//                     f0xbuf[j+i*f_x.Npx(s)]=data1D[j];
//                 }
//             }

//             if (PE.NODES() > 1) { 
//                 if (PE.RANK()!=0) {
//                    MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
//                 } 
//                 else {            
//                     // Fill data for rank = 0 
//                     for(size_t i(0); i < outNxLocal; i++) {
//                         for (size_t j(0); j < f_x.Npx(s); ++j) {
//                             f0x1Global(j,i) = f0xbuf[j+i*f_x.Npx(s)];
//                         }
//                     }
//                     // Fill data for rank > 0 
//                     for (int rr = 1; rr < PE.NODES(); ++rr){
//                         MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
//                         for(size_t i(0); i < outNxLocal; i++) {
//                             for (size_t j(0); j < f_x.Npx(s); ++j) {
//                                 f0x1Global(j,i + outNxLocal*rr) = f0xbuf[j+i*f_x.Npx(s)];
//                             }
//                         }
//                     }
//                 }
//             }
//             else { 
//                 for(size_t i(0); i < outNxGlobal; i++) {
//                     for (size_t j(0); j < f_x.Npx(s); ++j) {
//                         f0x1Global(j,i) = f0xbuf[j+i*f_x.Npx(s)];
//                     }
//                 }
//             }

//             if (PE.RANK() == 0) expo.Export_h5("f0-x", f0x1Global, tout, s);

//             delete[] f0xbuf;
//         }

//     }
// //--------------------------------------------------------------
// //--------------------------------------------------------------  
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::f10(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float* f0xbuf = new float[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s), 1, 0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(j,i + outNxLocal*rr,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(j,i + outNxLocal*rr,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f10-x", f0x1Global, tout, s);

        delete[] f0xbuf;
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::f11(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float* f0xbuf = new float[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s),1,1, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(j,i + outNxLocal*rr,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(j,i + outNxLocal*rr,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f11-x", f0x1Global, tout, s);

        delete[] f0xbuf;
    }


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::f20(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float* f0xbuf = new float[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s), 2, 0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(j,i + outNxLocal*rr,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(j,i + outNxLocal*rr,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f20-x", f0x1Global, tout, s);

        delete[] f0xbuf;
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::fl0(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float* f0xbuf = new float[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s), Y.DF(s).l0()-1,0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(j,i + outNxLocal*rr,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(j,i + outNxLocal*rr,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(j,i,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(j,i,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("fl0-x", f0x1Global, tout, s);

        delete[] f0xbuf;
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::n(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                            const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* nbuf = new float[msg_sz];
    valarray<float> nGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        for(size_t i(0); i < msg_sz; ++i) {

            // std::cout << "f00n[" << i << "]=" << (Y.SH(s,0,0)).xVec(i+Nbc)[0] << "\n";

            nbuf[i] = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(nbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    nGlobal[i] = nbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(nbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        nGlobal[i + outNxLocal*rr] = nbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                nGlobal[i] = nbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("n", nGlobal, tout, s);

    }

    delete[] nbuf;

}
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor_1D::T(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                            const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    // float* nbuf = new float[msg_sz];
    float* tbuf = new float[msg_sz];
    valarray<float> tGlobal(outNxGlobal); //, yglob_axis.dim());

    double convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        // float deltapra = pra[1] - pra[0];

        // need nbuf for density normalization
        for(size_t i(0); i < msg_sz; ++i) {
            tbuf[i] = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
            tbuf[i] /= 3.0*4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);

            // tbuf[i] *= convert_factor/Y.DF(s).mass();

            tbuf[i] *= 1.0/Y.DF(s).mass();

        }

//          for (int ix(0); ix < outNxLocal; ++ix){
//              double accumulator = 0.0;
//              for(int ip(0); ip < pra.size(); ++ip) {
//                  double p_squared = pra[ip]*pra[ip];
//                  double val = Y.SH(s,0,0)(ip,ix+Nbc).real();
//                  // pr/dx() can be moved outside loop if anyone cares.
//                  val *= (deltapra) *p_squared*p_squared;
//                  accumulator += val;
//              }
//              tbuf[ix] = accumulator; 
//              // do the normalization
//              tbuf[ix] *= ( 4.0*M_PI*2.0/3.0 );
//              // do the density part of the normalization...
//              tbuf[ix] /= nbuf[ix]; 

// //               double pth2 = ( TemperatureFromPthSq[ix] );
//              tbuf[ix] *= convert_factor;//pth2 / (1.0 + pth2) * convert_factor;
//          }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(tbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    tGlobal[i] = tbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(tbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        tGlobal[i + outNxLocal*rr] = tbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                tGlobal[i] = tbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("T", tGlobal, tout, s);

    }

    // delete[] nbuf;
    delete[] tbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Jx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Jxbuf = new float[msg_sz];
    valarray<float> JxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jxbuf[i] = Y.DF(s).q()/Y.DF(s).mass()*4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JxGlobal[i] = Jxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Jxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        JxGlobal[i + outNxLocal*rr] = Jxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                JxGlobal[i] = Jxbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jx", JxGlobal, tout, s);

    }

    delete[] Jxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Jy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Jybuf = new float[msg_sz];
    valarray<float> JyGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jybuf[i] = Y.DF(s).q()/Y.DF(s).mass()*8.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JyGlobal[i] = Jybuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Jybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        JyGlobal[i + outNxLocal*rr] = Jybuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                JyGlobal[i] = Jybuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jy", JyGlobal, tout, s);

    }

    delete[] Jybuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Jz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Jzbuf = new float[msg_sz];
    valarray<float> JzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        for(size_t i(0); i < msg_sz; ++i) {
            Jzbuf[i] = Y.DF(s).q()/Y.DF(s).mass()*-8.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JzGlobal[i] = Jzbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Jzbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        JzGlobal[i + outNxLocal*rr] = Jzbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                JzGlobal[i] = Jzbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jz", JzGlobal, tout, s);

    }

    delete[] Jzbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Qx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Vxbuf;
    double Vybuf;
    double Vzbuf;
    double VxVxbuf;
    double VyVybuf;
    double VzVzbuf;
    double VxVybuf;
    double VxVzbuf;
    double VyVzbuf;
    double nbuf;
    double Ubuf;

    float* Qxbuf = new float[msg_sz];
    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        for(size_t i(0); i < msg_sz; ++i) {
//                nbuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
//                Ubuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
//
//                Vxbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
//                Vybuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vzbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vxbuf /= nbuf;
//                Vybuf /= nbuf;
//                Vzbuf /= nbuf;
//
//                VxVxbuf = 4.0/3.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,0)).xVec(i+Nbc) ), pra, 4);
//                VyVybuf = 4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//                VzVzbuf = -1.0*VyVybuf;
//
//                VyVybuf -= VxVxbuf/2.0;
//                VzVzbuf -= VxVxbuf/2.0;
//
//                VxVxbuf += Ubuf/2.0;
//                VyVybuf += Ubuf/2.0;
//                VzVzbuf += Ubuf/2.0;
//
//                VxVybuf = 4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VxVzbuf = -4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VyVzbuf = -4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//
//                VxVxbuf /= nbuf;
//                VyVybuf /= nbuf;
//                VzVzbuf /= nbuf;
//                VxVybuf /= nbuf;
//                VxVzbuf /= nbuf;
//                VyVzbuf /= nbuf;

            Qxbuf[i] = 4.0*M_PI/3.0*Algorithms::moment(  vfloat( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        QxGlobal[i + outNxLocal*rr] = Qxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                QxGlobal[i] = Qxbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qx", QxGlobal, tout, s);

    }

    delete[] Qxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Qy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Vxbuf;
    double Vybuf;
    double Vzbuf;
    double VxVxbuf;
    double VyVybuf;
    double VzVzbuf;
    double VxVybuf;
    double VxVzbuf;
    double VyVzbuf;
    double nbuf;
    double Ubuf;

    float* Qxbuf = new float[msg_sz];
    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        for(size_t i(0); i < msg_sz; ++i) {
//                nbuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
//                Ubuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
//
//                Vxbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
//                Vybuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vzbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vxbuf /= nbuf;
//                Vybuf /= nbuf;
//                Vzbuf /= nbuf;
//
//                VxVxbuf = 4.0/3.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,0)).xVec(i+Nbc) ), pra, 4);
//                VyVybuf = 4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//                VzVzbuf = -1.0*VyVybuf;
//
//                VyVybuf -= VxVxbuf/2.0;
//                VzVzbuf -= VxVxbuf/2.0;
//
//                VxVxbuf += Ubuf/2.0;
//                VyVybuf += Ubuf/2.0;
//                VzVzbuf += Ubuf/2.0;
//
//                VxVybuf = 4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VxVzbuf = -4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VyVzbuf = -4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//
//                VxVxbuf /= nbuf;
//                VyVybuf /= nbuf;
//                VzVzbuf /= nbuf;
//                VxVybuf /= nbuf;
//                VxVzbuf /= nbuf;
//                VyVzbuf /= nbuf;

            Qxbuf[i] = 8.0*M_PI/3.0*Algorithms::moment(  vfloat( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        QxGlobal[i + outNxLocal*rr] = Qxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                QxGlobal[i] = Qxbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qy", QxGlobal, tout, s);

    }

    delete[] Qxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Qz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Vxbuf;
    double Vybuf;
    double Vzbuf;
    double VxVxbuf;
    double VyVybuf;
    double VzVzbuf;
    double VxVybuf;
    double VxVzbuf;
    double VyVzbuf;
    double nbuf;
    double Ubuf;

    float* Qxbuf = new float[msg_sz];
    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis( float(0.0),f_x.Pmax(s), f_x.Np(s) ) );
        for(size_t i(0); i < msg_sz; ++i) {
//                nbuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
//                Ubuf = 4.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
//
//                Vxbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
//                Vybuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vzbuf  = 4.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
//                Vxbuf /= nbuf;
//                Vybuf /= nbuf;
//                Vzbuf /= nbuf;
//
//                VxVxbuf = 4.0/3.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,0)).xVec(i+Nbc) ), pra, 4);
//                VyVybuf = 4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//                VzVzbuf = -1.0*VyVybuf;
//
//                VyVybuf -= VxVxbuf/2.0;
//                VzVzbuf -= VxVxbuf/2.0;
//
//                VxVxbuf += Ubuf/2.0;
//                VyVybuf += Ubuf/2.0;
//                VzVzbuf += Ubuf/2.0;
//
//                VxVybuf = 4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VxVzbuf = -4.0*2.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,1)).xVec(i+Nbc) ), pra, 4);
//                VyVzbuf = -4.0*4.0/5.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,2,2)).xVec(i+Nbc) ), pra, 4);
//
//                VxVxbuf /= nbuf;
//                VyVybuf /= nbuf;
//                VzVzbuf /= nbuf;
//                VxVybuf /= nbuf;
//                VxVzbuf /= nbuf;
//                VyVzbuf /= nbuf;

            Qxbuf[i] = -8.0*M_PI/3.0*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        QxGlobal[i + outNxLocal*rr] = Qxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                QxGlobal[i] = Qxbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qz", QxGlobal, tout, s);

    }

    delete[] Qxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::vNx(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* vNxbuf = new float[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<float> pra( Algorithms::MakeCAxis( float(0.0),f_x.Pmax(s), f_x.Np(s) ) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = (float) (1.0 / 6.0 * (Algorithms::moment(vfloat((Y.SH(s, 1, 0)).xVec(i + Nbc) ), pra, 6)
                                              / Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5)));
        }

        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        vNxGlobal[i + outNxLocal*rr] = vNxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                vNxGlobal[i] = vNxbuf[i];
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("vNx", vNxGlobal, tout, 0);

    }

    delete[] vNxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::vNy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* vNxbuf = new float[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal);

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( Algorithms::MakeCAxis( float(0.0), f_x.Pmax(s), f_x.Np(s) ) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = (float) (2.0 / 6.0 * (Algorithms::moment(vfloat((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
                                              / Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5)));
        }
        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        vNxGlobal[i + outNxLocal*rr] = vNxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                vNxGlobal[i] = vNxbuf[i];
            }
        }
        if (PE.RANK() == 0) expo.Export_h5("vNy", vNxGlobal, tout, 0);
    }
    delete[] vNxbuf;
}
// --------------------------------------------------------------
// --------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::vNz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    float* vNxbuf = new float[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal);

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<float> pra( Algorithms::MakeCAxis( float(0.0), f_x.Pmax(s), f_x.Np(s) ) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = (float) (-2.0 / 6.0 * (Algorithms::moment(vfloat_complex((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
                                               / Algorithms::moment(  vfloat( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5)));
        }
        if (PE.NODES() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.NODES(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        vNxGlobal[i + outNxLocal*rr] = vNxbuf[i];
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                vNxGlobal[i] = vNxbuf[i];
            }
        }
        if (PE.RANK() == 0) expo.Export_h5("vNz", vNxGlobal, tout, 0);
    }
    delete[] vNxbuf;
}
// --------------------------------------------------------------
// --------------------------------------------------------------


//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Ux(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Uxbuf = new float[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vx(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    UxGlobal[i + outNxLocal*rr] = Uxbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            UxGlobal[i] = Uxbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ux", UxGlobal, tout, 0);



    delete[] Uxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Uy(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Uxbuf = new float[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vy(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    UxGlobal[i + outNxLocal*rr] = Uxbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            UxGlobal[i] = Uxbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Uy", UxGlobal, tout, 0);



    delete[] Uxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Uz(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Uxbuf = new float[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vz(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    UxGlobal[i + outNxLocal*rr] = Uxbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            UxGlobal[i] = Uxbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Uz", UxGlobal, tout, 0);



    delete[] Uxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Z(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                            const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Uxbuf = new float[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().Z(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    UxGlobal[i + outNxLocal*rr] = Uxbuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            UxGlobal[i] = Uxbuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Z", UxGlobal, tout, 0);



    delete[] Uxbuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::ni(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* nibuf = new float[msg_sz];
    valarray<float> niGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        nibuf[i] = static_cast<float>(Y.HYDRO().density(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(nibuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                niGlobal[i] = nibuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(nibuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    niGlobal[i + outNxLocal*rr] = nibuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            niGlobal[i] = nibuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("ni", niGlobal, tout, 0);



    delete[] nibuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor_1D::Ti(const State1D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float* Thydrobuf = new float[msg_sz];
    valarray<float> ThydroGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Thydrobuf[i] = static_cast<float>(511000.0/3.0*Y.HYDRO().temperature(i+Nbc));
    }

    if (PE.NODES() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Thydrobuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ThydroGlobal[i] = Thydrobuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.NODES(); ++rr){
                MPI_Recv(Thydrobuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < outNxLocal; i++) {
                    ThydroGlobal[i + outNxLocal*rr] = Thydrobuf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < outNxGlobal; i++) {
            ThydroGlobal[i] = Thydrobuf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ti", ThydroGlobal, tout, 0);



    delete[] Thydrobuf;

}
//--------------------------------------------------------------
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
                                     std::valarray<float> data,
                                     const size_t& step,
                                     const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

//  Check Header file correctness
    if (Hdr[tag].dim() != 1) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 1D structure\n";
        exit(1);
    }

//  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    H5::H5File file = hmake_file(filename);

    hsize_t dimsf[1] = {data.size()};              // dataset dimensions
    H5::DataSpace dataspace( 1, dimsf );
    H5::DataSet dataset = file.createDataSet( tag, H5::PredType::NATIVE_FLOAT, dataspace );

    float xmin = {  static_cast<float>( Hdr[tag].min(0)) };
    float xmax = {  static_cast<float>( Hdr[tag].max(0)) };

    hinit_attr(file, tag, step, xmax, xmin);

    H5::Group group = file.createGroup("/AXIS");

    std::size_t pos = Hdr[tag].label(0).find("[");

    string axismainname = "AXIS1";
    float axisrange[2];
    axisrange[0] = xmin;
    axisrange[1] = xmax;
    string axislongname = Hdr[tag].label(0).substr(0,pos);
    string axisname = Hdr[tag].label(0).substr(0,pos);
    string axistype = "linear";
    string axisunits = Hdr[tag].units(0);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    group.close();

    dataset.write( &(data[0]) , H5::PredType::NATIVE_FLOAT );
    hfile_add_attr_todataset(dataset, "LONG_NAME", tag);
    hfile_add_attr_todataset(dataset, "UNITS", formulary().Label(tag));

    dataset.close();

    hclose_file(file);

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
                                     Array2D<float> data,
                                     const size_t& step,
                                     const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

//  Check Header file correctness
    if (Hdr[tag].dim() != 2) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 2D structure\n";
        exit(1);
    }

//  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    H5::H5File file = hmake_file(filename);

    hsize_t dimsf[2] = { data.dim2(), data.dim1() };              // dataset dimensions
    H5::DataSpace dataspace( 2, dimsf );
    H5::DataSet dataset = file.createDataSet( tag, H5::PredType::NATIVE_FLOAT, dataspace );

    float xmin[2] = {  static_cast<float>( Hdr[tag].min(0)),
                       static_cast<float>( Hdr[tag].min(1)) };

    float xmax[2] = {  static_cast<float>( Hdr[tag].max(0)),
                       static_cast<float>( Hdr[tag].max(1)) };

    hinit_attr2(file, tag, step, xmax, xmin);

    H5::Group group = file.createGroup("/AXIS");

    std::size_t pos = Hdr[tag].label(1).find("[");

    string axismainname = "AXIS1";
    float axisrange[2];
    axisrange[0] = xmin[1];
    axisrange[1] = xmax[1];
    string axislongname = Hdr[tag].label(1).substr(0,pos);
    string axisname = Hdr[tag].label(1).substr(0,pos);
    string axistype = "linear";
    string axisunits = Hdr[tag].units(1);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    pos = Hdr[tag].label(0).find("[");

    axismainname = "AXIS2";
    axisrange[0] = xmin[0];
    axisrange[1] = xmax[0];
    axislongname = Hdr[tag].label(0).substr(0,pos);
    axisname = Hdr[tag].label(0).substr(0,pos);
    axisunits = Hdr[tag].units(0);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    group.close();

    dataset.write( &(data.array())[0] , H5::PredType::NATIVE_FLOAT );
    hfile_add_attr_todataset(dataset, "LONG_NAME", tag);
    hfile_add_attr_todataset(dataset, "UNITS", formulary().Label(tag));

    dataset.close();

    hclose_file(file);

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
                                     Array3D<float> data,
                                     const size_t& step,
                                     const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

//  Check Header file correctness
    if (Hdr[tag].dim() != 3) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n";
        exit(1);
    }

//  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    H5::H5File file = hmake_file(filename);

    hsize_t dimsf[3] = { data.dim3(), data.dim2(), data.dim1() };              // dataset dimensions
    H5::DataSpace dataspace( 3, dimsf );
    H5::DataSet dataset = file.createDataSet( tag, H5::PredType::NATIVE_FLOAT, dataspace );

    float xmin[3] = {  static_cast<float>( Hdr[tag].min(0)),
                       static_cast<float>( Hdr[tag].min(1)),
                       static_cast<float>( Hdr[tag].min(2))};
    float xmax[3] = {  static_cast<float>( Hdr[tag].max(0)),
                       static_cast<float>( Hdr[tag].max(1)),
                       static_cast<float>( Hdr[tag].max(2))};

    hinit_attr2(file, tag, step, xmax, xmin);

    H5::Group group = file.createGroup("/AXIS");

    std::size_t pos = Hdr[tag].label(1).find("[");

    string axismainname = "AXIS1";
    float axisrange[2];
    axisrange[0] = xmin[1];
    axisrange[1] = xmax[1];
    string axislongname = Hdr[tag].label(1).substr(0,pos);
    string axisname = Hdr[tag].label(1).substr(0,pos);
    string axistype = "linear";
    string axisunits = Hdr[tag].units(1);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    pos = Hdr[tag].label(0).find("[");

    axismainname = "AXIS2";
    axisrange[0] = xmin[0];
    axisrange[1] = xmax[0];
    axislongname = Hdr[tag].label(0).substr(0,pos);
    axisname = Hdr[tag].label(0).substr(0,pos);
    axisunits = Hdr[tag].units(0);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    pos = Hdr[tag].label(2).find("[");

    axismainname = "AXIS3";
    axisrange[0] = xmin[2];
    axisrange[1] = xmax[2];
    axislongname = Hdr[tag].label(2).substr(0,pos);
    axisname = Hdr[tag].label(2).substr(0,pos);
    axisunits = Hdr[tag].units(2);

    haxis(group, axismainname, axisrange, axislongname, axisname, axistype, axisunits);

    group.close();

    dataset.write( &(data.array())[0] , H5::PredType::NATIVE_FLOAT );
    hfile_add_attr_todataset(dataset, "LONG_NAME", tag);
    hfile_add_attr_todataset(dataset, "UNITS", formulary().Label(tag));

    dataset.close();

    hclose_file(file);

}
//--------------------------------------------------------------
//
//--------------------------------------------------------------
H5::H5File Export_Files::Xport:: hmake_file(string ofilename) {
//--------------------------------------------------------------
//    Create HDF5 file
//--------------------------------------------------------------

    const H5std_string  FILE_NAME( ofilename );
    H5::H5File file( FILE_NAME, H5F_ACC_TRUNC );
    return file;

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hclose_file(H5::H5File &file) {
//--------------------------------------------------------------
//    Close HDF5 file
//--------------------------------------------------------------

    file.close();

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hinit_attr(H5::H5File &hfilehandle, const std::string tag,
                                      size_t step,
                                      float xmax, float xmin) {
//--------------------------------------------------------------
//    Add initial attributes:
//    time step, iteration number, name of diagnostic,
//    physical time, time units, type, and max and min of axes ranges
//--------------------------------------------------------------

    double dt_out(Input::List().t_stop / Input::List().n_outsteps);

    float dt = { static_cast<float>(
                         dt_out/static_cast<double>(
                                 size_t(static_cast<int>(dt_out
                                                         /Input::List().clf_dp))+1
                         ))};
    
    hfile_add_attr(hfilehandle, "DT", dt);

    float time = { static_cast<float>(dt_out*step) };
    hfile_add_attr(hfilehandle, "TIME", time);

    hfile_add_attr(hfilehandle, "NAME", tag);

    // float time = { static_cast<float>(dt_out*step*dt) };
    // hfile_add_attr(hfilehandle, "TIME", time);

    string timeunits = "1/\\omega_p";
    hfile_add_attr(hfilehandle, "TIME UNITS", timeunits);

    string typegrid = "grid";
    hfile_add_attr(hfilehandle, "TYPE", typegrid);

    hfile_add_attr(hfilehandle, "XMAX", xmax);

    hfile_add_attr(hfilehandle, "XMIN", xmin);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Export_Files::Xport:: hinit_attr2(H5::H5File &hfilehandle, const std::string tag,
                                       size_t step,
                                       float xmax[], float xmin[]) {
//--------------------------------------------------------------
//    Add initial attributes:
//    time step, iteration number, name of diagnostic,
//    physical time, time units, type, and max and min of axes ranges
//--------------------------------------------------------------

    double dt_out(Input::List().t_stop / Input::List().n_outsteps);

    float dt = { static_cast<float>(
                         dt_out/static_cast<double>(
                                 size_t(static_cast<int>(dt_out
                                                         /Input::List().clf_dp))+1
                         ))};
    hfile_add_attr(hfilehandle, "DT", dt);

    float time = { static_cast<float>(dt_out*step) };
    hfile_add_attr(hfilehandle, "TIME", time);

    hfile_add_attr(hfilehandle, "NAME", tag);

    string timeunits = "1/\\omega_p";
    hfile_add_attr(hfilehandle, "TIME UNITS", timeunits);

    string typegrid = "grid";
    hfile_add_attr(hfilehandle, "TYPE", typegrid);

    hfile_add_attr2(hfilehandle, "XMAX", xmax);

    hfile_add_attr2(hfilehandle, "XMIN", xmin);

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: haxis(H5::Group &hgrouphandle, string axismainname, float axisrange[2],
                                 string axislongname, string axisname,
                                 string axistype, string axisunits) {
//--------------------------------------------------------------
//   Add axis attributes
//--------------------------------------------------------------

    hsize_t dims2[1] = { 2 };
    H5::DataSpace dataspaceg1( 1, dims2 );
    H5::DataSet datasetg1 = hgrouphandle.createDataSet(axismainname, H5::PredType::NATIVE_FLOAT,
                                                       dataspaceg1);
    datasetg1.write( axisrange, H5::PredType::NATIVE_FLOAT );
    hfile_add_attr_todataset(datasetg1, "LONG_NAME", axislongname);
    hfile_add_attr_todataset(datasetg1, "NAME", axisname);
    hfile_add_attr_todataset(datasetg1, "TYPE", axistype);
    hfile_add_attr_todataset(datasetg1, "UNITS", axisunits);

    dataspaceg1.close();
    datasetg1.close();

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr(H5::H5File &hfilehandle, string attrname,
                                          int attrdata) {
//--------------------------------------------------------------
//   Add 1D integer attribute 
//--------------------------------------------------------------

    hsize_t dims1[1] = { 1 };
    H5::DataSpace attr_dataspace_1 = H5::DataSpace (1, dims1 );
    H5::Attribute attribute_int = hfilehandle.createAttribute( attrname, H5::PredType::NATIVE_INT,
                                                               attr_dataspace_1);
    attribute_int.write( H5::PredType::NATIVE_INT, &attrdata);
    attr_dataspace_1.close();
    attribute_int.close();

}

//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr2(H5::H5File &hfilehandle, string attrname,
                                           int attrdata[2]) {
//--------------------------------------------------------------
//   Add 2D integer attribute 
//--------------------------------------------------------------

    hsize_t dims1[1] = { 2 };
    H5::DataSpace attr_dataspace_2 = H5::DataSpace (1, dims1 );
    H5::Attribute attribute_int = hfilehandle.createAttribute( attrname, H5::PredType::NATIVE_INT,
                                                               attr_dataspace_2);
    attribute_int.write( H5::PredType::NATIVE_INT, attrdata);
    attr_dataspace_2.close();
    attribute_int.close();

}

//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr(H5::H5File &hfilehandle, string attrname,
                                          float attrdata) {
//--------------------------------------------------------------
//   Add 1D float attribute 
//--------------------------------------------------------------

    hsize_t dims1[1] = { 1 };
    H5::DataSpace attr_dataspace_1 = H5::DataSpace (1, dims1 );
    H5::Attribute attribute_float = hfilehandle.createAttribute( attrname, H5::PredType::NATIVE_FLOAT,
                                                                 attr_dataspace_1);
    attribute_float.write( H5::PredType::NATIVE_FLOAT, &attrdata);
    attr_dataspace_1.close();
    attribute_float.close();

}

//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr2(H5::H5File &hfilehandle, string attrname,
                                           float attrdata[2]) {
//--------------------------------------------------------------
//   Add 2D float attribute 
//--------------------------------------------------------------

    hsize_t dims1[1] = { 2 };
    H5::DataSpace attr_dataspace_2 = H5::DataSpace (1, dims1 );
    H5::Attribute attribute_float = hfilehandle.createAttribute( attrname, H5::PredType::NATIVE_FLOAT,
                                                                 attr_dataspace_2);
    attribute_float.write( H5::PredType::NATIVE_FLOAT, attrdata);
    attr_dataspace_2.close();
    attribute_float.close();

}

//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr(H5::H5File &hfilehandle, string attrname,
                                          string attrdata) {
//--------------------------------------------------------------
//   Add string attribute 
//--------------------------------------------------------------

    hsize_t dims1[1] = { 1 };
    H5::DataSpace attr_dataspace_1 = H5::DataSpace (1, dims1 );
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    H5std_string strwritebuf ( attrdata );
    H5::Attribute attribute_string = hfilehandle.createAttribute( attrname, strdatatype, attr_dataspace_1);
    attribute_string.write(strdatatype, strwritebuf);
    attr_dataspace_1.close();
    attribute_string.close();

}

//--------------------------------------------------------------


//--------------------------------------------------------------
void Export_Files::Xport:: hfile_add_attr_todataset(H5::DataSet &hdatasethandle, string attrname,
                                                    string attrdata) {
//--------------------------------------------------------------
//   Add string attribute to a specific dataset
//--------------------------------------------------------------

    hsize_t dims1[1] = { 1 };
    H5::DataSpace attr_dataspace_1 = H5::DataSpace (1, dims1 );
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    H5std_string strwritebuf ( attrdata );
    H5::Attribute attribute_string = hdatasethandle.createAttribute( attrname, strdatatype, attr_dataspace_1);
    attribute_string.write(strdatatype, strwritebuf);
    attr_dataspace_1.close();
    attribute_string.close();

}

//-------------------------------------------------------------



//**************************************************************
//**************************************************************
