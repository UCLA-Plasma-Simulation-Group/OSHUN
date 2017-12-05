/*! \brief  Export Files - Definitions
 * \author  PICKSC
 *  \date   2017
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

//  Declarations
#include "external/spline.h"
#include "state.h"
#include "formulary.h"
#include "setup.h"
#include "parallel.h"
#include "nmethods.h"
#include "input.h"
#include "export.h"

//------------------------------------------------------------------------------
// Create a folder
//
// @param[in]  _name  The name of the folder
//
int Export_Files::Makefolder(string _name){

    mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    const char*   foldername = _name.data();
    int   status(mkdir(foldername,_permissions));
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

        if (Input::List().particlepusher) {
            if (Makefolder("OUTPUT/PARTICLES") != 0)
                cout<<"Warning: Folder 'OUTPUT/PARTICLES' exists" << endl;
            if (Makefolder("OUTPUT/PARTICLES/prtx") != 0)
                cout<<"Warning: Folder 'OUTPUT/PARTICLES/prtx' exists" << endl;
            if (Makefolder("OUTPUT/PARTICLES/prtpx") != 0)
                cout<<"Warning: Folder 'OUTPUT/PARTICLES/prtpx' exists" << endl;
            if (Makefolder("OUTPUT/PARTICLES/prtpy") != 0)
                cout<<"Warning: Folder 'OUTPUT/PARTICLES/prtpy' exists" << endl;
            if (Makefolder("OUTPUT/PARTICLES/prtpz") != 0)
                cout<<"Warning: Folder 'OUTPUT/PARTICLES/prtpz' exists" << endl;
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

    if (  Input::List().o_p1x1 || Input::List().o_p2x1 || Input::List().o_p3x1 ||
      Input::List().o_p1p3x1 || Input::List().o_p1p2x1 || Input::List().o_p2p3x1 ||
      Input::List().o_f0x1 ||  Input::List().o_f10x1 ||  Input::List().o_f11x1 
      ||  Input::List().o_f20x1 || Input::List().o_fl0x1) {
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

    fld2d.push_back( "Ex-2d"     );
    fld2d.push_back( "Ey-2d"     );
    fld2d.push_back( "Ez-2d"     );
    fld2d.push_back( "Bx-2d"     );
    fld2d.push_back( "By-2d"     );
    fld2d.push_back( "Bz-2d"     );

//  Moments
    mom2d.push_back( "P-2d"      );
    mom2d.push_back( "T-2d"      );
    mom2d.push_back( "n-2d"      );
    mom2d.push_back( "Qx-2d"     );
    mom2d.push_back( "Qy-2d"     );
    mom2d.push_back( "Qz-2d"     );
    mom2d.push_back( "vNx-2d"    );
    mom2d.push_back( "vNy-2d"    );
    mom2d.push_back( "vNz-2d"    );
    mom2d.push_back( "Jx-2d"     );
    mom2d.push_back( "Jy-2d"     );
    mom2d.push_back( "Jz-2d"     );
    

    // mom2d.push_back( "Ux"      );
    // mom2d.push_back( "Uy"      );
    // mom2d.push_back( "Uz"      );
    // mom2d.push_back( "Z"      );
    // mom2d.push_back( "ni"      );
    // mom2d.push_back( "Ti"      );

    part.push_back("prtx");
    part.push_back("prtpx");
    part.push_back("prtpy");
    part.push_back("prtpz");

//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px-x");
        pvsx.push_back( "py-x");
        pvsx.push_back( "pz-x");
    }

    for (size_t s(0); s < species; ++s) {
        pvsx2d.push_back( "px-xy");
        pvsx2d.push_back( "py-xy");
        pvsx2d.push_back( "pz-xy");
    }

//  f-x
    for (size_t s(0); s < species; ++s) {
        fvsx.push_back( "f0-x");//+stringify(s) );
        fvsx.push_back( "f10-x");
        fvsx.push_back( "f11-x");
        fvsx.push_back( "f20-x");
        fvsx.push_back( "fl0-x");
    }

    for (size_t s(0); s < species; ++s) {
        fvsx2d.push_back( "f0-xy");//+stringify(s) );
        fvsx2d.push_back( "f10-xy");
        fvsx2d.push_back( "f11-xy");
        fvsx2d.push_back( "f20-xy");
        fvsx2d.push_back( "fl0-xy");
    }

//  p-p-x
    for (size_t s(0); s < species; ++s) {
        pvspvsx.push_back( "pxpy-x");
        pvspvsx.push_back( "pypz-x");
        pvspvsx.push_back( "pxpz-x");
    }

    for (size_t s(0); s < species; ++s) {
        pvspvsx2d.push_back( "pxpy-xy");
        pvspvsx2d.push_back( "pypz-xy");
        pvspvsx2d.push_back( "pxpz-xy");
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

// 3D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z, oAxis _imre,
 string _Ql, float _Qc,
 string _tl, string _tu, float _tc,
 string _oD)
: title(_Ql), time(_tl), timeU(_tu),
titleC(_Qc), timeC(_tc),
oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
    xyz_axis.push_back(_z);
    xyz_axis.push_back(_imre);
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
   string homedir)

{
    size_t species(_axis.pdim());
    DefaultTags dTags(species);

    vector< oAxis > xyz, pxyz, imre, pr, parts;
    for (size_t i(0); i < _axis.xdim(); ++i) {
        xyz.push_back( Export_Files::oAxis(_axis.xgmin(i), _axis.xgmax(i), _axis.Nxg(i)) );
    }
    parts.push_back( Export_Files::oAxis(0.5, Input::List().numparticles+0.5, Input::List().numparticles));

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
        ( find(oTags.begin(),oTags.end(), dTags.time[tloc]) == oTags.end() ) ) 
    {
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

    parts[0].label = "particle index";
    parts[0].units = "#";

    //  pxyz Axis
    for (size_t i(0); i < species; ++i) 
    {
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
    for (size_t i(0); i < dTags.fld.size(); ++i) 
    {
        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fld[i]) != oTags.end() ) 
        {
            string nounits = dTags.fld[i].substr(0, dTags.fld[i].find("_"));
            string folder = homedir + "OUTPUT/FLD/" + nounits + "/";
                //         Generate a header file for this tag
            Hdr[dTags.fld[i]] = Header(xyz[0],
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
            Hdr[dTags.mom[i]] = Header(xyz[0],
               nounits+"["+formulary().Label(dTags.mom[i])+"]",
               formulary().Uconv(dTags.mom[i]),
               tlabel, tunits, tconv, folder);
        }
    } // <--

//  Tags for Particles -->
    for (size_t i(0); i < dTags.part.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.part[i]) != oTags.end() ) {

            string nounits = dTags.part[i].substr(0, dTags.part[i].find("_"));
            string folder = homedir + "OUTPUT/PARTICLES/" + nounits + "/";
            //         Generate a header file for this tag
            

            Hdr[dTags.part[i]] = Header(parts,
             "particles", //nounits+"["+formulary().Label(dTags.mom[i])+"]",
             1.0,//formulary().Uconv(dTags.mom[i]),
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
      //
      //
    //  Tags for Fields -->
    for (size_t i(0); i < dTags.fld2d.size(); ++i) 
    {
        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fld2d[i]) != oTags.end() ) 
        {
            string nounits = dTags.fld2d[i].substr(0, dTags.fld2d[i].find("-"));
            string folder = homedir + "OUTPUT/FLD/" + nounits + "/";
                //         Generate a header file for this tag
            Hdr[dTags.fld2d[i]] = Header(xyz[0],xyz[1],
               nounits+"["+formulary().Label(dTags.fld2d[i])+"]",
               formulary().Uconv(dTags.fld2d[i]),
               tlabel, tunits, tconv, folder);
        }
    } // <--

//  Tags for Moments -->
    for (size_t i(0); i < dTags.mom2d.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.mom2d[i]) != oTags.end() ) {

            string nounits = dTags.mom2d[i].substr(0, dTags.mom2d[i].find("-"));
            string folder = homedir + "OUTPUT/MOM/" + nounits + "/";
            //         Generate a header file for this tag
            Hdr[dTags.mom2d[i]] = Header(xyz[0],xyz[1],
               nounits+"["+formulary().Label(dTags.mom2d[i])+"]",
               formulary().Uconv(dTags.mom2d[i]),
               tlabel, tunits, tconv, folder);
        }
    } // <--

// //  Tags for Particles -->
//     for (size_t i(0); i < dTags.part.size(); ++i) {

//         //     If this tag is an output tag
//         if ( find(oTags.begin(),oTags.end(), dTags.part[i]) != oTags.end() ) {

//             string nounits = dTags.part[i].substr(0, dTags.part[i].find("_"));
//             string folder = homedir + "OUTPUT/PARTICLES/" + nounits + "/";
//             //         Generate a header file for this tag
            

//             Hdr[dTags.part[i]] = Header(parts,
//              "particles", //nounits+"["+formulary().Label(dTags.mom[i])+"]",
//              1.0,//formulary().Uconv(dTags.mom[i]),
//              tlabel, tunits, tconv, folder);
//         }
//     } // <--


//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx2d.size(); ++i) {

        // for (size_t k(0); k < oTags.size(); ++k)
        //         std::cout << " \n \n k = " << oTags[k];

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.pvsx2d[i]) != oTags.end() ) {

            string folder = homedir + "OUTPUT/DISTR/" + dTags.pvsx2d[i] + "/";
            Makefolder(folder);


//          Generate a header file for this tag
//          For each 9 you have a different species 
            Hdr[dTags.pvsx2d[i]] = Header( pr[i/3], xyz[0], xyz[1],
             "f"+stringify(i/3), 1.0, tlabel, tunits, tconv, folder);
            // std::cout << " \n k = " << dTags.pvsx2d[i] << "\n";
            // std::cout << " \n \n 11 ";
        }
    } //<--

//  Tags for f-x -->
    for (size_t i(0); i < dTags.fvsx2d.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fvsx2d[i]) != oTags.end() ) {

            string folder = homedir + "OUTPUT/DISTR/" + dTags.fvsx2d[i] + "/";
            Makefolder(folder);

//          Generate a header file for this tag
            Hdr[dTags.fvsx2d[i]] = Header( pr[i/5], xyz[0], xyz[1], imre[0],
             "f", 1.0, tlabel, tunits, tconv, folder);

        }
    } //<--

    //  Tags for p-p-x -->
//     for (size_t i(0); i < dTags.pvspvsx.size(); ++i) {

//         //     If this tag is an output tag
//         if ( find(oTags.begin(),oTags.end(), dTags.pvspvsx[i]) != oTags.end() ) {

//             string folder = homedir + "OUTPUT/DISTR/" + dTags.pvspvsx[i] + "/";
//             Makefolder(folder);

// //          Generate a header file for this tag
//             Hdr[dTags.pvspvsx[i]] = Header( pxyz[i/1+species],pxyz[(i+1)/1+species], xyz[0],
//                 "f"+stringify(i/1), 1.0, tlabel, tunits, tconv, folder);

//             // std::cout << " \n \n 11 ";
//         }
//     } //<--
    

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

//      Read fields
    for(size_t ifields(0); ifields < Y.Fields(); ++ifields){
        for(size_t i(0); i < Y.EMF().Ex().numx(); ++i) {
            fin.read((char *)(&(Y.FLD(ifields))(i)), sizeof(Y.FLD(ifields)(i)));
        }
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
    for (size_t s(0); s < Y.Species(); ++s) {
        for (size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
            for (size_t i(0); i < (Y.DF(s))(nh).dim(); ++i) {
                fout.write((char *) &(Y.DF(s))(nh)(i), sizeof((Y.DF(s))(nh)(i)));
            }
        }
    }

//      Write fields
    for(size_t ifields(0); ifields < Y.Fields(); ++ifields){
        for(size_t i(0); i < Y.EMF().Ex().numx(); ++i) {
            fout.write((char *)(&(Y.FLD(ifields))(i)), sizeof(Y.FLD(ifields)(i)));
        }
    }


    fout.flush();
    fout.close();
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Read restart file
void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, State2D& Y) {

//      Generate filename 
    string   filename(hdir+"RESTART/re_2D_");
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

//      Read fields
    for(size_t ifields(0); ifields < Y.Fields(); ++ifields){
        for(size_t ix(0); ix < Y.EMF().Ex().numx(); ++ix) 
        {
            for(size_t iy(0); iy < Y.EMF().Ex().numy(); ++iy) 
            {
                fin.read((char *)(&(Y.FLD(ifields))(ix,iy)), sizeof(Y.FLD(ifields)(ix,iy)));
            }
        }
    }

    fin.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Write restart file
void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, State2D& Y) {

//      Generate filename 
    string   filename(hdir+"RESTART/re_2D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ofstream  fout(filename.c_str(), ios::binary);

//      Write distribution functions
    
    for (size_t s(0); s < Y.Species(); ++s) {
        for (size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
            for (size_t i(0); i < (Y.DF(s))(nh).dim(); ++i) {
                fout.write((char *) &(Y.DF(s))(nh)(i), sizeof((Y.DF(s))(nh)(i)));
            }
        }
    }
    

//      Write fields
    for(size_t ifields(0); ifields < Y.Fields(); ++ifields){
        for(size_t ix(0); ix < Y.EMF().Ex().numx(); ++ix) 
        {
            for(size_t iy(0); iy < Y.EMF().Ex().numy(); ++iy) 
            {
                fout.write((char *)(&(Y.FLD(ifields))(ix,iy)), sizeof(Y.FLD(ifields)(ix,iy)));
            }
        }
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
Output_Data::PLegendre2D::PLegendre2D( size_t Nl, size_t Nm,
    double pmax, valarray<double> px, valarray<double> py){
    // size_t Npx, size_t Npy ) {

    size_t sz(0);


    // valarray<float> //p( Algorithms::MakeAxis(pmin, pmax, Np)),
    // px(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npx)),
    // py(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npy));

//  Generate the structure to save the polynomials
    sz = ((Nm+1)*(2*Nl-Nm+2))/2;

    // plegendre = new vector< Array2D<double> >(sz,Array2D<double>(px.size(),py.size()))  ;
    for (size_t i(0);i < sz; ++i)
    {
        plegendre.push_back(Array2D<double>(px.size(),py.size()));
    }
    // 
    // vector< Array2D<double> > plegendre;

    
//  Generate the polynomial values for each cos(theta) = px/p 
    if (Nl > 1)
    {
        for (size_t ipy(0); ipy < py.size(); ++ipy) {

            for (size_t ipx(0); ipx < px.size(); ++ipx) {

                double invp(sqrt(py[ipy] * py[ipy] + px[ipx] * px[ipx]));

                if (invp > pmax) invp = 0;
                if (invp > 0.0) invp = 1.0 / invp;

    //          For given px/p generate all the polynomial values ...
                Array2D<double> vL(Algorithms::Legendre(static_cast<double>(px[ipx]) * invp, Nl, Nm));
                //          ... and save them
                size_t k(0);
                for (size_t l(0); l < Nl + 1; ++l) {
                    for (size_t m(0); m < ((Nm < l) ? Nm : l) + 1; ++m) {
                        // for (size_t k(0); k < sz; ++k) {
                        // k = ((l < Nm + 1) ? ((l * (l + 1)) / 2 + m) : (l * (Nm + 1) - (Nm * (Nm + 1)) / 2 + m));

                        // if (isnan(vL(l,m)))
                        // {
                        //     std::cout << "\n\n is nan \n ipx,ipy = " << ipx << " , " << ipy;
                        //     std::cout << "\n invp = " << invp;
                        //     std::cout << "\n px[ipx] = " << px[ipx];
                        //     exit(1);
                        // }

                        // plegendre.push_back
                        // (*plegendre)[k](ipx, ipy) = vL(l, m);
                        // std::cout << "\n (l,m) = " << l << "," << m << " \n ";
                        
                        (plegendre)[k](ipx, ipy) = vL(l, m);

                        // std::cout << "\n ka = " << k << " \n ";
                        ++k;
                    }
                }
            }
        }
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre2D::PLegendre2D( const PLegendre2D& other ) {

//  Generate the structure to save the polynomials
    // plegendre = new vector< Array2D<double> > ;
    // vector< Array2D<double> > plegendre;
    for (size_t i(0); i < other.dim(); ++i) {
        // (*plegendre).push_back( other(i) );
        (plegendre).push_back( other(i) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre2D::~PLegendre2D(){
    // delete plegendre;
}
//--------------------------------------------------------------

//**************************************************************
//--------------------------------------------------------------
Output_Data::harmonicvsposition::harmonicvsposition(const Grid_Info& _G) {

//  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pmax.push_back( static_cast<double>(_G.axis.pmax(s)) );
        pmin.push_back( static_cast<double>(_G.axis.pmin(s)) );
        nump.push_back( static_cast<double>(_G.axis.Np(s))   );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::harmonicvsposition::harmonicvsposition( const harmonicvsposition& other) {

    for (size_t s(0); s < other.Species(); ++s) {
        pmin.push_back( other.Pmin(s) );
        pmax.push_back( other.Pmax(s) );
        nump.push_back( other.Np(s)   );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::harmonicvsposition::~harmonicvsposition(){

}
//-------------------------------------------------------------
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location x0 
//  into a cartesian grid.
Array2D<float>  Output_Data::harmonicvsposition::operator()(DistFunc1D& df, size_t l, size_t m, size_t x0, size_t s) {
    
    Array2D<float> fx1(Np(s),2); 

    for (size_t ip(0); ip < Np(s); ++ip) {
        //  At each location in |p|
        fx1(ip,0) = static_cast<float>( (df(l,m))(ip,x0).real() );
        fx1(ip,1) = static_cast<float>( (df(l,m))(ip,x0).imag() );
    }

    return fx1;
}

Array2D<float>  Output_Data::harmonicvsposition::operator()(DistFunc2D& df, size_t l, size_t m, size_t x0, size_t y0, size_t s) {
    
    Array2D<float> fx1(Np(s),2); 

    for (size_t ip(0); ip < Np(s); ++ip) {
        //  At each location in |p|
        fx1(ip,0) = static_cast<float>( (df(l,m))(ip,x0,y0).real() );
        fx1(ip,1) = static_cast<float>( (df(l,m))(ip,x0,y0).imag() );
    }

    return fx1;
}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
Output_Data::fulldistvsposition::fulldistvsposition( const Grid_Info& _G): grid(_G) 
{
    
    //  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pvec.push_back(valtovec(_G.axis.p(s)));
        
        pxvec.push_back(_G.axis.px(s));
        pyvec.push_back(_G.axis.py(s));
        pzvec.push_back(_G.axis.pz(s));

        PL2D.push_back( PLegendre2D( _G.l0[s], _G.m0[s],
          (_G.axis.pmax(s)), _G.axis.px(s), _G.axis.py(s) ) );  
    
        pout1D_p1.push_back( valarray<float>( _G.axis.Npx(s))); 
        pout1D_p2.push_back( valarray<float>( _G.axis.Npy(s))); 
        pout1D_p3.push_back( valarray<float>( _G.axis.Npz(s)));         
        pout2D_p1p2.push_back( Array2D<float>( _G.axis.Npx(s),_G.axis.Npy(s))); 
        pout2D_p1p3.push_back( Array2D<float>( _G.axis.Npx(s),_G.axis.Npz(s))); 
        pout2D_p2p3.push_back( Array2D<float>( _G.axis.Npy(s),_G.axis.Npz(s))); 

        pout3D.push_back( Array3D<float>( _G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
    }
    
    

    double py_sq, pz_sqppy_sq;
    size_t pindlt;

    for (size_t s(0); s < _G.axis.pdim(); ++s) 
    {
        pradius.push_back(Array3D<double>(_G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
        // pind_disttolowerbound.push_back(Array3D<float>(_G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
        // pind_oneminusdisttolowerbound.push_back(Array3D<float>(_G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
        // pind_lowerbound.push_back(Array3D<size_t>(_G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
        phi.push_back(Array2D<double>(_G.axis.Npy(s),_G.axis.Npz(s)));  phi[s] = 0.;
        
        for (size_t ipy(0); ipy < _G.axis.Npy(s); ++ipy)
        {
            py_sq = (_G.axis.py(s))[ipy]*(_G.axis.py(s))[ipy];
            for (size_t ipz(0); ipz < _G.axis.Npz(s); ++ipz)
            {
                pz_sqppy_sq = py_sq + (_G.axis.pz(s))[ipz]*(_G.axis.pz(s))[ipz];
                for (size_t ipx(0); ipx < _G.axis.Npx(s); ++ipx)
                {

                    pradius[s](ipx,ipy,ipz) = sqrt((_G.axis.px(s))[ipx]*(_G.axis.px(s))[ipx]+pz_sqppy_sq);

                    if (pradius[s](ipx,ipy,ipz) < _G.axis.pmax(s))
                    {
                        // pindlt = floor(pradius[s](ipx,ipy,ipz)/_G.axis.dp(s)[0]+0.5)-1;
                        // pind_disttolowerbound[s](ipx,ipy,ipz) = abs(pradius[s](ipx,ipy,ipz)-(_G.axis.p(s))[pindlt])/_G.axis.dp(s)[0];
                        // pind_oneminusdisttolowerbound[s](ipx,ipy,ipz) = 1.0 - pind_disttolowerbound[s](ipx,ipy,ipz);
                        // pind_lowerbound[s](ipx,ipy,ipz) = pindlt;
                        phi[s](ipy,ipz) = (atan2((_G.axis.pz(s))[ipz],(_G.axis.py(s))[ipy]) + M_PI);
                    }
                    else
                    {
                        // pind_disttolowerbound[s](ipx,ipy,ipz) = 0.;
                        // pind_lowerbound[s](ipx,ipy,ipz) = 0;
                        // pind_oneminusdisttolowerbound[s](ipx,ipy,ipz) = 0;
                    }
                }
            }       
        }


        dpx.push_back(_G.axis.dpx(s));
        dpy.push_back(_G.axis.dpy(s));
        dpz.push_back(_G.axis.dpz(s));


    }
        
}
//-------------------------------------------------------------
/**
 * @brief      Destroys the object.
 */
//-------------------------------------------------------------
Output_Data::fulldistvsposition::~fulldistvsposition(){

}

//-------------------------------------------------------------
/**
 * @brief      Creates p1 from a 3D grid
 */
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p1(DistFunc1D& df, size_t x0, size_t s) {

    pout1D_p1[s] = 0.0; 


    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p1[s][ipx] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]*dpz[s][ipz]);
                    }
                    else pout1D_p1[s][ipx] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p1[s][ipx] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re*dpy[s][ipy]*dpz[s][ipz]);
                        }
                        else pout1D_p1[s][ipx] += 0.;

                    }
                }
            }
        }
    }
    return pout1D_p1[s];    
}
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p1(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout1D_p1[s] = 0.0; 


    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p1[s][ipx] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]*dpz[s][ipz]);
                    }
                    else pout1D_p1[s][ipx] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p1[s][ipx] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re*dpy[s][ipy]*dpz[s][ipz]);
                        }
                        else pout1D_p1[s][ipx] += 0.;

                    }
                }
            }
        }
    }
    return pout1D_p1[s];    
}

//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p2(DistFunc1D& df, size_t x0, size_t s) {

    pout1D_p2[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {
        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));
     
        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p2[s][ipy] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpz[s][ipz]);
                    }
                    else pout1D_p2[s][ipy] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;


            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p2[s][ipy] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpz[s][ipz]);
                        }
                        else pout1D_p2[s][ipy] += 0.;

                    }
                }
            }
        }
    }

    return pout1D_p2[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p2(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout1D_p2[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {
        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));
     
        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p2[s][ipy] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpz[s][ipz]);
                    }
                    else pout1D_p2[s][ipy] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;


            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p2[s][ipy] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpz[s][ipz]);
                        }
                        else pout1D_p2[s][ipy] += 0.;

                    }
                }
            }
        }
    }

    return pout1D_p2[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p3(DistFunc1D& df, size_t x0, size_t s) {

    pout1D_p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p3[s][ipz] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpy[s][ipy]);
                    }
                    else pout1D_p3[s][ipz] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p3[s][ipz] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpy[s][ipy]);
                        }
                        else pout1D_p3[s][ipz] += 0.;

                    }
                }
            }
        }
    }

    return pout1D_p3[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
valarray<float>  Output_Data::fulldistvsposition::p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout1D_p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);
                
        double YSH_re;

        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout1D_p3[s][ipz] += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpy[s][ipy]);
                    }
                    else pout1D_p3[s][ipz] += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {
            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout1D_p3[s][ipz] += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpy[s][ipy]);
                        }
                        else pout1D_p3[s][ipz] += 0.;

                    }
                }
            }
        }
    }

    return pout1D_p3[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p1p2(DistFunc1D& df, size_t x0, size_t s) {

    pout2D_p1p2[s] = 0.0; 
    size_t im(0);

    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;                   
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {

                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p1p2[s](ipx,ipy) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpz[s][ipz]);
                    }
                    else pout2D_p1p2[s](ipx,ipy) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;


            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p1p2[s](ipx,ipy) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpz[s][ipz]);
                        }
                        else pout2D_p1p2[s](ipx,ipy) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p1p2[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p2 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p1p2(DistFunc2D& df, size_t x0,  size_t y0, size_t s) {

    pout2D_p1p2[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {

                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p1p2[s](ipx,ipy) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpz[s][ipz]);
                    }
                    else pout2D_p1p2[s](ipx,ipy) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;


            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p1p2[s](ipx,ipy) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpz[s][ipz]);
                        }
                        else pout2D_p1p2[s](ipx,ipy) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p1p2[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p2p3 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p2p3(DistFunc1D& df, size_t x0, size_t s) {

    pout2D_p2p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p2p3[s](ipy,ipz) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]);
                    }
                    else pout2D_p2p3[s](ipy,ipz) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);

                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p2p3[s](ipy,ipz) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpx[s][ipx]);
                        }
                        else pout2D_p2p3[s](ipy,ipz) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p2p3[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p2p3 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p2p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout2D_p2p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p2p3[s](ipy,ipz) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]);
                    }
                    else pout2D_p2p3[s](ipy,ipz) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);

                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p2p3[s](ipy,ipz) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpx[s][ipx]);
                        }
                        else pout2D_p2p3[s](ipy,ipz) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p2p3[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p3 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p1p3(DistFunc1D& df, size_t x0, size_t s) {

    pout2D_p1p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p1p3[s](ipx,ipz) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]);
                    }
                    else pout2D_p1p3[s](ipx,ipz) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p1p3[s](ipx,ipz) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpy[s][ipy]);
                        }
                        else pout2D_p1p3[s](ipx,ipz) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p1p3[s];    
}
//-------------------------------------------------------------
/**
 * @brief      Creates p1p3 from a 3D grid
 */
//-------------------------------------------------------------
Array2D<float>  Output_Data::fulldistvsposition::p1p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout2D_p1p3[s] = 0.0; 
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout2D_p1p3[s](ipx,ipz) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]);
                    }
                    else pout2D_p1p3[s](ipx,ipz) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0,y0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0,y0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout2D_p1p3[s](ipx,ipz) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpy[s][ipy]);
                        }
                        else pout2D_p1p3[s](ipx,ipz) += 0.;

                    }
                }
            }
        }
    }

    return pout2D_p1p3[s];    
}
//--------------------------------------------------------------
Array3D<float>  Output_Data::fulldistvsposition::p1p2p3(DistFunc1D& df, size_t x0, size_t s){
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location (x0,y0) 
//  into a cartesian grid.
//--------------------------------------------------------------
    pout3D[s] = 0.0;
    size_t im(0);
    for(size_t il = 0; il < grid.l0[s]+1; ++il)
    {

        size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));


        vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
                 
        tk::spline splSH;
        splSH.set_points(pvec[s],shdata_real);

        double YSH_re;
                        
        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
        {
            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
            {
                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                {
                    if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                    {                       
                        YSH_re = (splSH(pradius[s](ipx,ipy,ipz)));
                        pout3D[s](ipx,ipy,ipz) += static_cast<float>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy));
                    }
                    else pout3D[s](ipx,ipy,ipz) += 0.;
                }
            }
        }
    }

    for (im = 1; im < grid.m0[s]+1; ++im)
    {
        for (size_t il(im); il < grid.l0[s]+1; ++il)
        {

            vector<double> shdata_real( vdouble_real(   df(il,im).xVec(x0)));
            vector<double> shdata_imag( vdouble_imag(   df(il,im).xVec(x0)));
                 
            tk::spline splSH_r, splSH_i;
            splSH_r.set_points(pvec[s],shdata_real);
            splSH_i.set_points(pvec[s],shdata_imag);

            size_t i_dist = ((il < grid.m0[s]+1)?((il*(il+1))/2+im):(il*(grid.m0[s]+1)-(grid.m0[s]*(grid.m0[s]+1))/2 + im));

            double YSH_re, YSH_im;

            double mphi,calcos,calsin;


            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        mphi = im*phi[s](ipy,ipz);
                        calcos = cos(mphi);
                        calsin = sin(mphi);


                        if (pradius[s](ipx,ipy,ipz) <= grid.axis.pmax(s))
                        {                       
                            YSH_re = (splSH_r(pradius[s](ipx,ipy,ipz)));
                            YSH_im = (splSH_i(pradius[s](ipx,ipy,ipz)));

                            YSH_re *= calcos; 
                            YSH_im *= calsin; 
                            YSH_re -= YSH_im;

                            pout3D[s](ipx,ipy,ipz) += static_cast<float>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re);//*dpy[s]*dpz[s];
                        }
                        else pout3D[s](ipx,ipy,ipz) += 0.;

                    }
                }
            }
        }
    }
    return pout3D[s];
}
//**************************************************************
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::operator()(const State1D& Y, const Grid_Info& grid, const size_t tout,
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

    if (Input::List().particlepusher) {
        particles_x( Y, grid, tout, PE );
        particles_px( Y, grid, tout, PE );
        particles_py( Y, grid, tout, PE );
        particles_pz( Y, grid, tout, PE );
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::distdump(const State1D& Y, const Grid_Info& grid, const size_t tout,
   const Parallel_Environment_1D& PE) 
{

    if (Input::List().o_p1x1){
        px( Y, grid, tout, PE );
    }
    if (Input::List().o_p2x1){
        py( Y, grid, tout, PE );
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
    if (Input::List().o_f20x1)
    {
        f20( Y, grid, tout, PE );
    }
    if (Input::List().o_fl0x1){
        fl0( Y, grid, tout, PE );
    }

}

void Output_Data::Output_Preprocessor::bigdistdump(const State1D& Y, const Grid_Info& grid, const size_t tout,
   const Parallel_Environment_1D& PE) 
{

    if (Input::List().o_p1p2x1)
    {
        pxpy( Y, grid, tout, PE );
    }
    if (Input::List().o_p1p3x1)
    {
        pxpz( Y, grid, tout, PE );
    }
    if (Input::List().o_p2p3x1)
    {
        pypz( Y, grid, tout, PE );
    }
    if (Input::List().o_p1p2p3x1)
    {
        // pxpypz( Y, grid, tout, PE );
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::operator()(const State2D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_2D& PE) {

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

    // if (Input::List().o_Ux) {
    //     Ux( Y, grid, tout, PE );
    // }

    // if (Input::List().o_Uy) {
    //     Uy( Y, grid, tout, PE );
    // }

    // if (Input::List().o_Uz) {
    //     Uz( Y, grid, tout, PE );
    // }

    // if (Input::List().o_Z) {
    //     Z( Y, grid, tout, PE );
    // }

    // if (Input::List().o_ni) {
    //     ni( Y, grid, tout, PE );
    // }

    // if (Input::List().o_Ux) {
    //     Ti( Y, grid, tout, PE );
    // }

    // if (Input::List().particlepusher) {
    //     particles_x( Y, grid, tout, PE );
    //     particles_px( Y, grid, tout, PE );
    //     particles_py( Y, grid, tout, PE );
    //     particles_pz( Y, grid, tout, PE );
    // }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::distdump(const State2D& Y, const Grid_Info& grid, const size_t tout,
   const Parallel_Environment_2D& PE) 
{
    if (Input::List().o_p1x1){
        px( Y, grid, tout, PE );
    }
    if (Input::List().o_p2x1){
        py( Y, grid, tout, PE );
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
    if (Input::List().o_f20x1)
    {
        f20( Y, grid, tout, PE );
    }
    if (Input::List().o_fl0x1){
        fl0( Y, grid, tout, PE );
    }
}

void Output_Data::Output_Preprocessor::bigdistdump(const State2D& Y, const Grid_Info& grid, const size_t tout,
   const Parallel_Environment_2D& PE) 
{
    if (Input::List().o_p1p2x1)
    {
        pxpy( Y, grid, tout, PE );
    }
    if (Input::List().o_p1p3x1)
    {
        pxpz( Y, grid, tout, PE );
    }
    if (Input::List().o_p2p3x1)
    {
        pypz( Y, grid, tout, PE );
    }
    if (Input::List().o_p1p2p3x1)
    {
        // pxpypz( Y, grid, tout, PE );
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ex(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    // float* Exbuf[msg_sz];
    float Exbuf[msg_sz];
    valarray<float> ExGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Exbuf[i] = static_cast<float>( Y.EMF().Ex()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Exbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ExGlobal[i] = Exbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ey(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Eybuf[msg_sz];
    valarray<float> EyGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Eybuf[i] = static_cast<float>( Y.EMF().Ey()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Eybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EyGlobal[i] = Eybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ez(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Ezbuf[msg_sz];
    valarray<float> EzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Ezbuf[i] = static_cast<float>( Y.EMF().Ez()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Ezbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EzGlobal[i] = Ezbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bx(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Bxbuf[msg_sz];
    valarray<float> BxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bxbuf[i] = static_cast<float>( Y.EMF().Bx()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BxGlobal[i] = Bxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::By(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Bybuf[msg_sz];
    valarray<float> ByGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bybuf[i] = static_cast<float>( Y.EMF().By()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ByGlobal[i] = Bybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bz(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Bzbuf[msg_sz];
    valarray<float> BzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Bzbuf[i] = static_cast<float>( Y.EMF().Bz()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BzGlobal[i] = Bzbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ex(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float Exbuf[msg_sz];
    Array2D<float> ExGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Exbuf[i] = static_cast<float>( ( Y.EMF().Ex()(Nbc+ix,Nbc+iy).real() ));
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Exbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    // for(size_t i(0); i < outNxLocal; i++) {
                    ExGlobal(ix,iy) = Exbuf[i]; ++i;
                }
            }


            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();
                i = 0;

                MPI_Recv(Exbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        ExGlobal(ix + outNxLocal*rankx,iy + outNyLocal*ranky) = Exbuf[i]; ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                ExGlobal(ix,iy) = Exbuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ex-2d", ExGlobal, tout);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ey(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float Eybuf[msg_sz];
    Array2D<float> EyGlobal(outNxGlobal,outNyGlobal);

    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Eybuf[i] = static_cast<float>( Y.EMF().Ey()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Eybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    // for(size_t i(0); i < outNxLocal; i++) {
                    EyGlobal(ix,iy) = Eybuf[i]; ++i;
                }
            }


            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();
                i = 0;

                MPI_Recv(Eybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        EyGlobal(ix + outNxLocal*rankx,iy + outNyLocal*ranky) = Eybuf[i]; ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                EyGlobal(ix,iy) = Eybuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ey-2d", EyGlobal, tout);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ez(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float Ezbuf[msg_sz];
    Array2D<float> EzGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Ezbuf[i] = static_cast<float>( Y.EMF().Ez()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Ezbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    // for(size_t i(0); i < outNxLocal; i++) {
                    EzGlobal(ix,iy) = Ezbuf[i]; ++i;
                }
            }


            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();
                i = 0;

                MPI_Recv(Ezbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        EzGlobal(ix + outNxLocal*rankx,iy + outNyLocal*ranky) = Ezbuf[i]; ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                EzGlobal(ix,iy) = Ezbuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Ez-2d", EzGlobal, tout);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bx(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float Bxbuf[msg_sz];
    Array2D<float> BxGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bxbuf[i] = static_cast<float>( Y.EMF().Bx()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    // for(size_t i(0); i < outNxLocal; i++) {
                    BxGlobal(ix,iy) = Bxbuf[i]; ++i;
                }
            }


            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();
                i = 0;

                MPI_Recv(Bxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        BxGlobal(ix + outNxLocal*rankx,iy + outNyLocal*ranky) = Bxbuf[i]; ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                BxGlobal(ix,iy) = Bxbuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Bx-2d", BxGlobal, tout);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::By(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float Bybuf[msg_sz];
    Array2D<float> ByGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bybuf[i] = static_cast<float>( Y.EMF().By()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    // for(size_t i(0); i < outNxLocal; i++) {
                    ByGlobal(ix,iy) = Bybuf[i]; ++i;
                }
            }


            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();
                i = 0;

                MPI_Recv(Bybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        ByGlobal(ix + outNxLocal*rankx,iy + outNyLocal*ranky) = Bybuf[i]; ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                ByGlobal(ix,iy) = Bybuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("By-2d", ByGlobal, tout);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bz(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal);
    float Bzbuf[msg_sz];
    Array2D<float> BzGlobal(outNxGlobal,outNyGlobal);
    size_t i(0);
    int rankx,ranky;

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bzbuf[i] = static_cast<float>( Y.EMF().Bz()(Nbc + ix, Nbc + iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            i = 0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    BzGlobal(ix,iy) = Bzbuf[i]; ++i;
                }
            }

            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){

                rankx = rr % PE.MPI_X();
                ranky = rr / PE.MPI_X();

                i = 0;

                MPI_Recv(Bzbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        BzGlobal(ix + outNxLocal * rankx , iy + outNyLocal * ranky) = Bzbuf[i];
                        ++i;
                    }
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                BzGlobal(ix,iy) = Bzbuf[i];++i;
            }
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("Bz-2d", BzGlobal, tout);

}
//--------------------------------------------------------------   
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::px(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npx(grid.axis.Npx(s));
        int msg_sz(outNxLocal*Npx);
        Array2D<float> p1x1Global(Npx,outNxGlobal); //, yglob_axis.dim());
        float pxbuf[Npx*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) {

            // valarray<float> data1D = px_x( Y.DF(s), i+Nbc, s);
            valarray<float> data1D = p_x.p1( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npx; ++j) {
                pxbuf[j+i*Npx]=data1D[j];
            }
        }

        if (PE.MPI_Processes() > 1) {
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
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::py(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));



    for(int s(0); s < Y.Species(); ++s) {
        size_t Npy(grid.axis.Npy(s));
        int msg_sz(outNxLocal*Npy);
        Array2D<float> p2x1Global(Npy,outNxGlobal); //, yglob_axis.dim());
        float pybuf[Npy*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) {

            // valarray<float> data1D = py_x( Y.DF(s), i+Nbc, s);
            valarray<float> data1D = p_x.p2( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npy; ++j) {
                pybuf[j+i*Npy]=data1D[j];
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    for (size_t j(0); j < Npy; ++j) {
                        p2x1Global(j,i) = pybuf[j+i*Npy];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < Npy; ++j) {
                            p2x1Global(j,i + outNxLocal*rr) = pybuf[j+i*Npy];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                for (size_t j(0); j < Npy; ++j) {
                    p2x1Global(j,i) = pybuf[j+i*Npy];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("py-x", p2x1Global, tout, s);

    }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pz(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));



    for(int s(0); s < Y.Species(); ++s) {
        size_t Npz(grid.axis.Npz(s));
        int msg_sz(outNxLocal*Npz);
        Array2D<float> p3x1Global(Npz,outNxGlobal); //, yglob_axis.dim());
        float pzbuf[Npz*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) {

            // valarray<float> data1D = py_x( Y.DF(s), i+Nbc, s);
            valarray<float> data1D = p_x.p3( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npz; ++j) {
                pzbuf[j+i*Npz]=data1D[j];
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    for (size_t j(0); j < Npz; ++j) {
                        p3x1Global(j,i) = pzbuf[j+i*Npz];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pzbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < Npz; ++j) {
                            p3x1Global(j,i + outNxLocal*rr) = pzbuf[j+i*Npz];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                for (size_t j(0); j < Npz; ++j) {
                    p3x1Global(j,i) = pzbuf[j+i*Npz];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("pz-x", p3x1Global, tout, s);

    }

}
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::px(const State2D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npx(grid.axis.Npx(s));
        int msg_sz(outNxLocal*outNyLocal*Npx);
        Array3D<float> p1x1Global(Npx,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pxbuf[Npx*outNxLocal*outNyLocal];
        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {

                valarray<float> data1D = p_x.p1( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npx; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                size_t counter(0);
                // Fill data for rank = 0
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < Npx; ++j) 
                        {
                            p1x1Global(j,ix,iy) = pxbuf[counter];
                            ++counter;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr)
                {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(pxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    counter = 0;    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy) 
                        {
                            for (size_t j(0); j < Npx; ++j) 
                            {
                                p1x1Global(j,ix + outNxLocal*rankx, iy + outNyLocal*ranky) = pxbuf[counter];
                                ++counter;
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            size_t counter(0);            
            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < Npx; ++j) 
                    {
                        p1x1Global(j,ix,iy) = pxbuf[counter];
                        ++counter;
                    }
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("px-2d", p1x1Global, tout, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::py(const State2D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npy(grid.axis.Npy(s));
        int msg_sz(outNxLocal*outNyLocal*Npy);
        Array3D<float> p1x1Global(Npy,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pxbuf[Npy*outNxLocal*outNyLocal];
        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {

                valarray<float> data1D = p_x.p2( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npy; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                size_t counter(0);
                // Fill data for rank = 0
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < Npy; ++j) 
                        {
                            p1x1Global(j,ix,iy) = pxbuf[counter];
                            ++counter;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr)
                {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(pxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    counter = 0;    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy) 
                        {
                            for (size_t j(0); j < Npy; ++j) 
                            {
                                p1x1Global(j,ix + outNxLocal*rankx, iy + outNyLocal*ranky) = pxbuf[counter];
                                ++counter;
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            size_t counter(0);            
            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < Npy; ++j) 
                    {
                        p1x1Global(j,ix,iy) = pxbuf[counter];
                        ++counter;
                    }
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("py-2d", p1x1Global, tout, s);

    }

}
//--------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::pz(const State2D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npz(grid.axis.Npz(s));
        int msg_sz(outNxLocal*outNyLocal*Npz);
        Array3D<float> p1x1Global(Npz,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pxbuf[Npz*outNxLocal*outNyLocal];
        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {

                valarray<float> data1D = p_x.p3( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npz; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                size_t counter(0);
                // Fill data for rank = 0
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < Npz; ++j) 
                        {
                            p1x1Global(j,ix,iy) = pxbuf[counter];
                            ++counter;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr)
                {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(pxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    counter = 0;    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy) 
                        {
                            for (size_t j(0); j < Npz; ++j) 
                            {
                                p1x1Global(j,ix + outNxLocal*rankx, iy + outNyLocal*ranky) = pxbuf[counter];
                                ++counter;
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            size_t counter(0);            
            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < Npz; ++j) 
                    {
                        p1x1Global(j,ix,iy) = pxbuf[counter];
                        ++counter;
                    }
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("pz-2d", p1x1Global, tout, s);

    }

}
//--------------------------------------------------------------
/**
 * @brief      { function_description }
 *
 * @param[in]  Y     { parameter_description }
 * @param[in]  grid  The grid
 * @param[in]  tout  The tout
 * @param[in]  PE    { parameter_description }
 */
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpy(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) 
{

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);
    
    for(int s(0); s < Y.Species(); ++s) {
        
        int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npy(s));        
        Array3D<float> pxpyGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal); //, yglob_axis.dim());        
        float pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal];        
        ind = 0;    

        for (size_t i(0); i < outNxLocal; ++i) 
        {
            Array2D<float> data2D = p_x.p1p2( Y.DF(s), i+Nbc, s);
            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
            {
                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                {
                    pbuf[ind]=data2D(j,k);
                    ++ind;
                }
            }
           
        }

        if (PE.MPI_Processes() > 1) {
           if (PE.RANK()!=0) {
               MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
           }
           else {
               // Fill data for rank = 0
               ind = 0;
               for(size_t i(0); i < outNxLocal; i++) 
               {
                    for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                        {
                            pxpyGlobal(k,j,i) = pbuf[ind];
                            ++ind;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    ind = 0;
                    for(size_t i(0); i < outNxLocal; i++) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                            {
                                pxpyGlobal(k,j,i + outNxLocal*rr) = pbuf[ind];
                                ++ind;
                            }
                        }
                    }
                }
            }
        }
       else {
           ind = 0;
           for(size_t i(0); i < outNxLocal; i++) {
               for (size_t j(0); j < grid.axis.Npx(s); ++j) {
                   for (size_t k(0); k < grid.axis.Npy(s); ++k) {
                       pxpyGlobal(k,j,i) = pbuf[ind];
                       ++ind;
                   }
               }
           }
       }

       if (PE.RANK() == 0) expo.Export_h5("pxpy-x", pxpyGlobal, tout, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpy(const State2D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;
    size_t counter;


    for(int s(0); s < Y.Species(); ++s) 
    {

        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<float> pxpyGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<float> data2D = p_x.p1p2( Y.DF(s), ix+Nbc, iy+Nbc, s);

                for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                {
                    for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                    {
                        pbuf[counter]=data2D(j,k);
                        ++counter;
                    }
                }
            }
        }

        if (PE.MPI_Processes() > 1) 
        {
            if (PE.RANK()!=0) {
                MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                // Fill data for rank = 0
                counter = 0;
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                            {
                                pxpyGlobal(k,j,ix,iy) = pbuf[counter];
                                ++counter;
                            }
                        }
                   }
               }
               // Fill data for rank > 0
               for (int rr(1); rr < PE.MPI_Processes(); ++rr)
               {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    counter = 0;

                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    pxpyGlobal(k,j,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
                                    ++counter;
                                }
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            counter = 0;

            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                        {
                            pxpyGlobal(k,j,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pxpy-2d", pxpyGlobal, tout, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
/**
 * @brief      { function_description }
 *
 * @param[in]  Y     { parameter_description }
 * @param[in]  grid  The grid
 * @param[in]  tout  The tout
 * @param[in]  PE    { parameter_description }
 */
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pypz(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(outNxLocal*grid.axis.Npy(s)*grid.axis.Npz(s));
        Array3D<float> pypzGlobal(grid.axis.Npy(s),grid.axis.Npz(s),outNxGlobal); //, yglob_axis.dim());
        float pbuf[grid.axis.Npy(s)*grid.axis.Npz(s)*outNxLocal];

        ind = 0;

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = p_x.p2p3( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < grid.axis.Npy(s); ++j) 
            {
               for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                {
                   pbuf[ind]=data2D(j,k);
                   ++ind;
                }
           }
       }

       if (PE.MPI_Processes() > 1) {
           if (PE.RANK()!=0) {
               MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
           }
           else 
           {
               // Fill data for rank = 0
                ind = 0;
                for(size_t i(0); i < outNxLocal; i++) 
                {
                    for (size_t j(0); j < grid.axis.Npy(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                        {
                           pypzGlobal(k,j,i) = pbuf[ind];
                           ++ind;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    ind = 0;
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < grid.axis.Npy(s); ++j) {
                            for (size_t k(0); k < grid.axis.Npz(s); ++k) {
                                pypzGlobal(k,j,i + outNxLocal*rr) = pbuf[ind];
                                ++ind;
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            ind = 0;
            for(size_t i(0); i < outNxLocal; i++) 
            {
                for (size_t j(0); j < grid.axis.Npy(s); ++j) 
                {
                    for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                    {
                        pypzGlobal(k,j,i) = pbuf[ind];
                        ++ind;
                    }
                }
            }
        }

       if (PE.RANK() == 0) expo.Export_h5("pypz-x", pypzGlobal, tout, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pypz(const State2D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;
    size_t counter;


    for(int s(0); s < Y.Species(); ++s) 
    {

        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<float> dataGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<float> data2D = p_x.p2p3( Y.DF(s), ix+Nbc, iy+Nbc, s);

                for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                {
                    for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                    {
                        pbuf[counter]=data2D(j,k);
                        ++counter;
                    }
                }
            }
        }

        if (PE.MPI_Processes() > 1) 
        {
            if (PE.RANK()!=0) {
                MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                // Fill data for rank = 0
                counter = 0;
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                            {
                                dataGlobal(k,j,ix,iy) = pbuf[counter];
                                ++counter;
                            }
                        }
                   }
               }
               // Fill data for rank > 0
               for (int rr(1); rr < PE.MPI_Processes(); ++rr)
               {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    counter = 0;

                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    dataGlobal(k,j,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
                                    ++counter;
                                }
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            counter = 0;

            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                        {
                            dataGlobal(k,j,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pypz-2d", dataGlobal, tout, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
/**
 * @brief      { function_description }
 *
 * @param[in]  Y     { parameter_description }
 * @param[in]  grid  The grid
 * @param[in]  tout  The tout
 * @param[in]  PE    { parameter_description }
 */
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpz(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

   size_t Nbc = Input::List().BoundaryCells;
   MPI_Status status;

   size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
   size_t outNxGlobal(grid.axis.Nxg(0));
   size_t ind(0);

   for(int s(0); s < Y.Species(); ++s) {

       int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npz(s));
       Array3D<float> pxpzGlobal(grid.axis.Npx(s),grid.axis.Npz(s),outNxGlobal); //, yglob_axis.dim());
       float pbuf[grid.axis.Npx(s)*grid.axis.Npz(s)*outNxLocal];

       for (size_t i(0); i < outNxLocal; ++i) {

           Array2D<float> data2D = p_x.p1p2( Y.DF(s), i+Nbc, s);
           ind = 0;

           for (size_t j(0); j < grid.axis.Npx(s); ++j) {
               for (size_t k(0); k < grid.axis.Npz(s); ++k) {
                   pbuf[ind]=data2D(j,k);
                   ++ind;
               }
           }
           
       }

       if (PE.MPI_Processes() > 1) {
           if (PE.RANK()!=0) {
               MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
           }
           else 
           {
               ind = 0;
               // Fill data for rank = 0
               for(size_t i(0); i < outNxLocal; i++) 
               {
                    for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                        {
                            pxpzGlobal(k,j,i) = pbuf[ind];
                            ++ind;
                        }
                    }
                }
               // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    ind = 0;
                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                            {
                                pxpzGlobal(k,j,i + outNxLocal*rr) = pbuf[ind];
                                ++ind;
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            ind = 0;
            for(size_t i(0); i < outNxLocal; i++) 
            {
                for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                {
                    for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                    {
                        pxpzGlobal(k,j,i) = pbuf[ind];
                        ++ind;
                    }
                }
            }
        }

       if (PE.RANK() == 0) expo.Export_h5("pxpz-x", pxpzGlobal, tout, s);

   }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpz(const State2D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx, ranky;
    size_t counter;


    for(int s(0); s < Y.Species(); ++s) 
    {

        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<float> dataGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        float pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<float> data2D = p_x.p1p3( Y.DF(s), ix+Nbc, iy+Nbc, s);

                for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                {
                    for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                    {
                        pbuf[counter]=data2D(j,k);
                        ++counter;
                    }
                }
            }
        }

        if (PE.MPI_Processes() > 1) 
        {
            if (PE.RANK()!=0) {
                MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                // Fill data for rank = 0
                counter = 0;
                for(size_t ix(0); ix < outNxLocal; ++ix) 
                {
                    for(size_t iy(0); iy < outNyLocal; ++iy) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                            {
                                dataGlobal(k,j,ix,iy) = pbuf[counter];
                                ++counter;
                            }
                        }
                   }
               }
               // Fill data for rank > 0
               for (int rr(1); rr < PE.MPI_Processes(); ++rr)
               {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    counter = 0;

                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    dataGlobal(k,j,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
                                    ++counter;
                                }
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            counter = 0;

            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                    {
                        for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                        {
                            dataGlobal(k,j,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pxpz-2d", dataGlobal, tout, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
/**
 * @brief      { function_description }
 *
 * @param[in]  Y     { parameter_description }
 * @param[in]  grid  The grid
 * @param[in]  tout  The tout
 * @param[in]  PE    { parameter_description }
 */
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpypz(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npy(s)*grid.axis.Npz(s));
        Array4D<float> pxpypzGlobal(grid.axis.Npx(s),grid.axis.Npy(s),grid.axis.Npz(s),outNxGlobal); //, yglob_axis.dim());
        float pbuf[msg_sz];

        ind = 0;

        for (size_t i(0); i < outNxLocal; ++i) {

            Array3D<float> data3D = p_x.p1p2p3( Y.DF(s), i+Nbc, s);

            for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
            {
                for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                {
                    for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                    {
                        pbuf[ind]=data3D(ipx,ipy,ipz);
                        ++ind;
                    }
                }
            }
        }

        if (PE.MPI_Processes() > 1) 
        {
            if (PE.RANK()!=0) 
            {
                MPI_Send(pbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else 
            {
                // Fill data for rank = 0
                ind = 0;
                for(size_t i(0); i < outNxLocal; i++) 
                {
                    
                    for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
                    {
                        for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                        {
                            for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                            {
                                pxpypzGlobal(ipx,ipy,ipz,i) = pbuf[ind];
                                ++ind;
                            }
                        }
                    }
                }
               // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr)
                {
                    MPI_Recv(pbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    ind = 0;
                    for(size_t i(0); i < outNxLocal; i++) 
                    {
                        for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
                        {
                            for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                            {
                                for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                                {
                                    pxpypzGlobal(ipx,ipy,ipz,i + outNxLocal*rr) = pbuf[ind];
                                    ++ind;
                                }
                            }
                        }
                    }
                }
            }
        }
        else 
        {
            ind = 0;
            for(size_t i(0); i < outNxLocal; i++) 
            {                
                for (size_t ipx(0); ipx < grid.axis.Npx(s); ++ipx) 
                {
                    for (size_t ipy(0); ipy < grid.axis.Npy(s); ++ipy) 
                    {
                        for (size_t ipz(0); ipz < grid.axis.Npz(s); ++ipz) 
                        {
                            pxpypzGlobal(ipx,ipy,ipz,i) = pbuf[ind];
                            ++ind;
                        }
                    }
                }
            }
        }

        // if (PE.RANK() == 0) expo.Export_h5("pxpypz-x", pxpypzGlobal, tout, s);
    }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f0(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s),0,0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
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
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f10(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s), 1, 0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
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
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f11(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s),1,1, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
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
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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
    }


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f20(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) 
    {
        if (Y.DF(s).l0() > 1)
        {
            int msg_sz(2*outNxLocal*f_x.Np(s));
            Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
            float f0xbuf[msg_sz];

            for (size_t i(0); i < outNxLocal; ++i) {

                Array2D<float> data2D = f_x( Y.DF(s), 2, 0, i+Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                    f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                    f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
                }
            }

            if (PE.MPI_Processes() > 1) {
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
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

        }
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::fl0(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<float> f0x1Global(f_x.Np(s),outNxGlobal,2); //, yglob_axis.dim());
        float f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<float> data2D = f_x( Y.DF(s), Y.DF(s).l0(),0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
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
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f0(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;

//    float buf[msg_sz];
//    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<float> global(f_x.Np(s),outNxGlobal,outNyGlobal,2); //, yglob_axis.dim());
        float buf[msg_sz];
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<float> data2D = f_x(Y.DF(s), 0, 0, ix + Nbc, iy + Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                    buf[i] = data2D(j, 0);
                    buf[i + 1] = data2D(j, 1);
                    i+=2;
                }
            }

        }

            if (PE.MPI_Processes() > 1) {
                if (PE.RANK()!=0) {
                    MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(j, ix, iy, 0) = buf[i];
                                global(j, ix, iy, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 0) = buf[i];
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                    }
                }
            }
            else {
                i=0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            global(j, ix, iy, 0) = buf[i];
                            global(j, ix, iy, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f0-xy", global, tout, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
    void Output_Data::Output_Preprocessor::f10(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                                 const Parallel_Environment_2D& PE) {
        size_t Nbc(Input::List().BoundaryCells);
        MPI_Status status;
        size_t i(0);
        size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
        size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
        size_t outNxGlobal(grid.axis.Nxg(0));
        size_t outNyGlobal(grid.axis.Nxg(1));
        int rankx, ranky;

//    float buf[msg_sz];
//    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

        for(int s(0); s < Y.Species(); ++s) {
            int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
            Array4D<float> global(f_x.Np(s),outNxGlobal,outNyGlobal,2); //, yglob_axis.dim());
            float buf[msg_sz];
            i=0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    Array2D<float> data2D = f_x(Y.DF(s), 1, 0, ix + Nbc, iy + Nbc, s);
                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                        buf[i] = data2D(j, 0);
                        buf[i + 1] = data2D(j, 1);
                        i+=2;
                    }
                }
            }

                if (PE.MPI_Processes() > 1) {
                    if (PE.RANK()!=0) {
                        MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
                    }
                    else {
                        // Fill data for rank = 0
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(j, ix, iy, 0) = buf[i];
                                    global(j, ix, iy, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                        // Fill data for rank > 0
                        for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                            MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                            rankx = rr % PE.MPI_X();
                            ranky = rr / PE.MPI_X();
                            i=0;
                            for (size_t ix(0); ix < outNxLocal; ++ix) {
                                for (size_t iy(0); iy < outNyLocal; ++iy) {
                                    for (size_t j(0); j < f_x.Np(s); ++j) {
                                        global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 0) = buf[i];
                                        global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 1) = buf[i+1];
                                        i+=2;
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(j, ix, iy, 0) = buf[i];
                                global(j, ix, iy, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                }

                if (PE.RANK() == 0) expo.Export_h5("f10-xy", global, tout, s);

            }

        }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f11(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;

//    float buf[msg_sz];
//    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<float> global(f_x.Np(s),outNxGlobal,outNyGlobal,2); //, yglob_axis.dim());
        float buf[msg_sz];
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<float> data2D = f_x(Y.DF(s), 1, 1, ix + Nbc, iy + Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                    buf[i] = data2D(j, 0);
                    buf[i + 1] = data2D(j, 1);
                    i+=2;
                }
            }

        }

            if (PE.MPI_Processes() > 1) {
                if (PE.RANK()!=0) {
                    MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(j, ix, iy, 0) = buf[i];
                                global(j, ix, iy, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 0) = buf[i];
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                    }
                }
            }
            else {
                i=0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            global(j, ix, iy, 0) = buf[i];
                            global(j, ix, iy, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f11-xy", global, tout, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f20(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;

//    float buf[msg_sz];
//    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<float> global(f_x.Np(s),outNxGlobal,outNyGlobal,2); //, yglob_axis.dim());
        float buf[msg_sz];
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<float> data2D = f_x(Y.DF(s), 2, 0, ix + Nbc, iy + Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                    buf[i] = data2D(j, 0);
                    buf[i + 1] = data2D(j, 1);
                    i+=2;
                }
            }

        }

            if (PE.MPI_Processes() > 1) {
                if (PE.RANK()!=0) {
                    MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(j, ix, iy, 0) = buf[i];
                                global(j, ix, iy, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 0) = buf[i];
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                    }
                }
            }
            else {
                i=0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            global(j, ix, iy, 0) = buf[i];
                            global(j, ix, iy, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f20-xy", global, tout, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::fl0(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;

//    float buf[msg_sz];
//    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<float> global(f_x.Np(s),outNxGlobal,outNyGlobal,2); //, yglob_axis.dim());
        float buf[msg_sz];
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<float> data2D = f_x(Y.DF(s), Y.DF(s).l0(), 0, ix + Nbc, iy + Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                    buf[i] = data2D(j, 0);
                    buf[i + 1] = data2D(j, 1);
                    i+=2;
                }
            }

        }

            if (PE.MPI_Processes() > 1) {
                if (PE.RANK()!=0) {
                    MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(j, ix, iy, 0) = buf[i];
                                global(j, ix, iy, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 0) = buf[i];
                                    global(j, ix + outNxLocal * rankx, iy + outNyLocal * ranky, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                    }
                }
            }
            else {
                i=0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            global(j, ix, iy, 0) = buf[i];
                            global(j, ix, iy, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("fl0-xy", global, tout, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::n(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float nbuf[msg_sz];
    valarray<float> nGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        // // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( vfloat(grid.axis.p(s)) );
        valarray<float> pra( vfloat(grid.axis.p(s)));
        for(size_t i(0); i < msg_sz; ++i) {

            // std::cout << "f00n[" << i << "]=" << (Y.SH(s,0,0)).xVec(i+Nbc)[0] << "\n";

            nbuf[i] = 4.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(nbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    nGlobal[i] = nbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_x(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    float buf[msg_sz];
    valarray<float> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().x(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < msg_sz; i++) {
                    pGlobal[i] += buf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < msg_sz; i++) {
            pGlobal[i] = buf[i];
        }
    }
    

    if (PE.RANK() == 0) expo.Export_h5("prtx", pGlobal, tout);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_px(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    float buf[msg_sz];
    valarray<float> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().px(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < msg_sz; i++) {
                    pGlobal[i] += buf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < msg_sz; i++) {
            pGlobal[i] = buf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("prtpx", pGlobal, tout);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_py(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    float buf[msg_sz];
    valarray<float> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().py(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < msg_sz; i++) {
                    pGlobal[i] += buf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < msg_sz; i++) {
            pGlobal[i] = buf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("prtpy", pGlobal, tout);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_pz(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    float buf[msg_sz];
    valarray<float> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().pz(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                for(size_t i(0); i < msg_sz; i++) {
                    pGlobal[i] += buf[i];
                }
            }
        }
    }
        // Fill data for Nodes = 0
    else {
        for(size_t i(0); i < msg_sz; i++) {
            pGlobal[i] = buf[i];
        }
    }

    if (PE.RANK() == 0) expo.Export_h5("prtpz", pGlobal, tout);

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::T(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    // float nbuf[msg_sz];
    float tbuf[msg_sz];
    valarray<float> tGlobal(outNxGlobal); //, yglob_axis.dim());

    double convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );
        // float deltapra = pra[1] - pra[0];

        // need nbuf for density normalization
        for(size_t i(0); i < msg_sz; ++i) {
            tbuf[i] = 4.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
            tbuf[i] /= 3.0*4.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);

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

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(tbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    tGlobal[i] = tbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jx(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Jxbuf[msg_sz];
    valarray<float> JxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jxbuf[i] = Y.DF(s).q()*4.0/3.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JxGlobal[i] = Jxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jy(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Jybuf[msg_sz];
    valarray<float> JyGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jybuf[i] = Y.DF(s).q()*8.0/3.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jybuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JyGlobal[i] = Jybuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jz(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Jzbuf[msg_sz];
    valarray<float> JzGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );
        for(size_t i(0); i < msg_sz; ++i) {
            Jzbuf[i] = Y.DF(s).q()*-8.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jzbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JzGlobal[i] = Jzbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qx(const State1D& Y, const Grid_Info& grid, const size_t tout,
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

    float Qxbuf[msg_sz];

    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );
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

            Qxbuf[i] = 4.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_real( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qy(const State1D& Y, const Grid_Info& grid, const size_t tout,
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

    float Qxbuf[msg_sz];
    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        // valarray<float> pra( Algorithms::MakeCAxis(float(0.0), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );
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

            Qxbuf[i] = 8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_real( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qz(const State1D& Y, const Grid_Info& grid, const size_t tout,
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

    float Qxbuf[msg_sz];
    valarray<float> QxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );
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

            Qxbuf[i] = -8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
//                Qxbuf[i] -= Vxbuf *(VxVxbuf +VyVybuf +VzVzbuf )*nbuf ;
//                Qxbuf[i] -= 2.0 *( VxVxbuf *Vxbuf  + VxVybuf *Vybuf + VxVzbuf*Vzbuf)*nbuf ;
//                Qxbuf[i] += 2.0 *( Vxbuf *Vxbuf  + Vybuf *Vybuf +Vzbuf*Vzbuf ) * Vxbuf  * nbuf ;
            Qxbuf[i] *= 0.5;
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNx(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float vNxbuf[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<float> pra( vfloat(grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<float>( (1.0 / 6.0 * (Algorithms::moment(vfloat_real((Y.SH(s, 1, 0)).xVec(i + Nbc) ), pra, 6)
              / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))) );
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNy(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float vNxbuf[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal);

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<float> pra( vfloat(grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<float>( (2.0 / 6.0 * (Algorithms::moment(vfloat_real((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
              / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))) );
        }
        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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
}
// --------------------------------------------------------------
// --------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNz(const State1D& Y, const Grid_Info& grid, const size_t tout,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    float vNxbuf[msg_sz];
    valarray<float> vNxGlobal(outNxGlobal);

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<float> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<float> pra( vfloat(grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<float>( (-2.0 / 6.0 * (Algorithms::moment(vfloat_complex((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
               / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))));
        }
        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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
    
}
// --------------------------------------------------------------
// --------------------------------------------------------------
void Output_Data::Output_Preprocessor::n(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                            const Parallel_Environment_2D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    float nbuf[msg_sz];
    Array2D<float> nGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    int rankx,ranky;

    for(int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );

        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                // std::cout << "f00n[" << i << "]=" << (Y.SH(s,0,0)).xVec(i+Nbc)[0] << "\n";
                nbuf[i] = 4.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 2);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(nbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                i=0;
                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        nGlobal(ix,iy) = nbuf[i]; ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    i=0;
                    MPI_Recv(nbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t ix(0); ix < outNxLocal; ++ix) {
                        for(size_t iy(0); iy < outNyLocal; ++iy) {
                            nGlobal(ix + outNxLocal*rankx, iy + outNyLocal*ranky) = nbuf[i]; 
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i=0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    nGlobal(ix,iy) = nbuf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("n-2d", nGlobal, tout, s);

    }

    

}
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::T(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                            const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    int outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    int outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal*outNyLocal);

    float tbuf[msg_sz];
    Array2D<float> tGlobal(outNxGlobal,outNyGlobal);

    double convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);

    for(size_t s(0); s < Y.Species(); ++s) 
    {
        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                tbuf[i] = Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 4);
                tbuf[i] /= 3.0*Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 2);

                tbuf[i] *= 1.0/Y.DF(s).mass();

                // std::cout << "T1[" << ix << "," << iy << "] = " << tbuf[i] << "\n";
                
                ++i;
     
     
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(tbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                i=0;
                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        tGlobal(ix,iy) = tbuf[i];++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(tbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i=0;
                    for(size_t ix(0); ix < outNxLocal; ++ix) {
                        for(size_t iy(0); iy < outNyLocal; ++iy) {
                            tGlobal(ix + outNxLocal*rankx, iy + outNyLocal*ranky) = tbuf[i]; 
                            ++i;
                        }
                    }
                }
            }
        }
        // Fill data for Nodes = 0
        else 
        {
            i=0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    tGlobal(ix,iy) = tbuf[i]; 
                    // std::cout << "T2[" << ix << "," << iy << "] = " << tbuf[i] << "\n";
                    ++i;
                    
                }
            }
            // exit(1);
        }

        if (PE.RANK() == 0) expo.Export_h5("T-2d", tGlobal, tout, s);

    }

    

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jx(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal*outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) 
    {
        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<float>(Y.DF(s).q()*4.0/3.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,1,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                i=0;
                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix,iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i=0;
                    for(size_t ix(0); ix < outNxLocal; ++ix) {
                        for(size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal*rankx, iy + outNyLocal*ranky) = buf[i]; 
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else 
        {
            i=0;
            for(size_t ix(0); ix < outNxLocal; ++ix) 
            {
                for(size_t iy(0); iy < outNyLocal; ++iy) 
                {
                    global(ix,iy) = buf[i]; ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jx-2d", global, tout, s);

    }

    

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jy(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal*outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal,outNyGlobal);

    for(int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<float>(Y.DF(s).q()*8.0/3.0*M_PI*Algorithms::moment(  vfloat_real( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                i=0;
                for(size_t ix(0); ix < outNxLocal; ++ix) {
                    for(size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix,iy) = buf[i];++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i=0;
                    for(size_t ix(0); ix < outNxLocal; ++ix) {
                        for(size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal*rankx, iy + outNyLocal*ranky) = buf[i]; ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i=0;
            for(size_t ix(0); ix < outNxLocal; ++ix) {
                for(size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix,iy) = buf[i]; ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jy-2d", global, tout, s);

    }

    

}
// //--------------------------------------------------------------
// //--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jz(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<float>(Y.DF(s).q()*-8.0/3.0*M_PI*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Jz-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qx(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*4.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_real( (Y.SH(s,1,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qx-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qy(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_real( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qy-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qz(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*-8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vfloat_complex( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("Qz-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNx(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<float>( (1.0 / 6.0 * (Algorithms::moment(vfloat_real((Y.SH(s, 1, 0)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
              / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("vNx-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNy(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<float>( (2.0 / 6.0 * (Algorithms::moment(vfloat_real((Y.SH(s, 1, 1)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
              / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("vNy-2d", global, tout, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNz(const State2D& Y, const Grid_Info& grid, const size_t tout,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx, ranky;
    int msg_sz(outNxLocal * outNyLocal);

    float buf[msg_sz];
    Array2D<float> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<float> pra( vfloat(grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for (size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<float>( (-2.0 / 6.0 * (Algorithms::moment(vfloat_complex((Y.SH(s, 1, 1)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
                        / Algorithms::moment(  vfloat_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
            } else {
                // Fill data for rank = 0
                i = 0;
                for (size_t ix(0); ix < outNxLocal; ++ix) {
                    for (size_t iy(0); iy < outNyLocal; ++iy) {
                        global(ix, iy) = buf[i];
                        ++i;
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr) {
                    rankx = rr % PE.MPI_X();
                    ranky = rr / PE.MPI_X();

                    MPI_Recv(buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);

                    i = 0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            global(ix + outNxLocal * rankx, iy + outNyLocal * ranky) = buf[i];
                            ++i;
                        }
                    }
                }
            }
        }
            // Fill data for Nodes = 0
        else {
            i = 0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    global(ix, iy) = buf[i];
                    ++i;
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("vNz-2d", global, tout, s);

    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------

//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ux(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Uxbuf[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vx(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Uy(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Uxbuf[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vy(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Uz(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Uxbuf[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().vz(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Z(const State1D& Y, const Grid_Info& grid, const size_t tout,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Uxbuf[msg_sz];
    valarray<float> UxGlobal(outNxGlobal); //, yglob_axis.dim());


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<float>(Y.HYDRO().Z(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::ni(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float nibuf[msg_sz];
    valarray<float> niGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        nibuf[i] = static_cast<float>(Y.HYDRO().density(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(nibuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                niGlobal[i] = nibuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ti(const State1D& Y, const Grid_Info& grid, const size_t tout,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    // size_t st(0), bi(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    float Thydrobuf[msg_sz];
    valarray<float> ThydroGlobal(outNxGlobal); //, yglob_axis.dim());

    for(size_t i(0); i < msg_sz; ++i) {
        Thydrobuf[i] = static_cast<float>(511000.0/3.0*Y.HYDRO().temperature(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Thydrobuf, msg_sz, MPI_FLOAT, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ThydroGlobal[i] = Thydrobuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
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
//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
                                     Array4D<float> data,
                                     const size_t& step,
                                     const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

//  Check Header file correctness
    if (Hdr[tag].dim() != 4) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 4D structure\n";
        exit(1);
    }

//  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    H5::H5File file = hmake_file(filename);

    hsize_t dimsf[4] = { data.dim4(), data.dim3(), data.dim2(), data.dim1() };              // dataset dimensions
    H5::DataSpace dataspace( 4, dimsf );
    H5::DataSet dataset = file.createDataSet( tag, H5::PredType::NATIVE_FLOAT, dataspace );

    float xmin[4] = {   static_cast<float>( Hdr[tag].axis(0)[ 0 ] ),
                        static_cast<float>( Hdr[tag].axis(1)[ 0 ] ),
                        static_cast<float>( Hdr[tag].axis(2)[ 0 ] ),
                        static_cast<float>( Hdr[tag].axis(3)[ 0 ] )};
    float xmax[4] = {   static_cast<float>( Hdr[tag].axis(0)[ (Hdr[tag].axis(0).size()-1) ] ),
                        static_cast<float>( Hdr[tag].axis(1)[ (Hdr[tag].axis(1).size()-1) ] ),
                        static_cast<float>( Hdr[tag].axis(2)[ (Hdr[tag].axis(1).size()-1) ] ),
                        static_cast<float>( Hdr[tag].axis(3)[ (Hdr[tag].axis(1).size()-1) ] )};

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

    pos = Hdr[tag].label(3).find("[");

    axismainname = "AXIS4";
    axisrange[0] = xmin[3];
    axisrange[1] = xmax[3];
    axislongname = Hdr[tag].label(3).substr(0,pos);
    axisname = Hdr[tag].label(3).substr(0,pos);
    axisunits = Hdr[tag].units(3);

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
                 /Input::List().dt))+1
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
                 /Input::List().dt))+1
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
