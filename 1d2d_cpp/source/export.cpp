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

#include "external/highfive/H5Attribute.hpp"
#include "external/highfive/H5File.hpp"
#include "external/highfive/H5Group.hpp"
#include "external/highfive/H5DataSet.hpp"
#include "external/highfive/H5DataSpace.hpp"

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

    if (Makefolder("restart") != 0) cout << "Warning: Folder 'restart' exists" << endl;

    if (Makefolder("output") != 0) cout << "Warning: Folder 'output' exists" << endl;

    if ( Input::List().o_EHist ) {
        if ( Makefolder("output/NUM") != 0) cout << "Warning: Folder 'output/NUM' exists" << endl;
    }

    if ( Input::List().o_Ex ||  Input::List().o_Ey || Input::List().o_Ez ||
     Input::List().o_Bx ||  Input::List().o_By || Input::List().o_Bz )  {
        if (Makefolder("output/fields") != 0)
            cout<<"Warning: Folder 'output/fields' exists" << endl;

        if (Input::List().o_Ex) {
            if (Makefolder("output/fields/Ex") != 0)
                cout<<"Warning: Folder 'output/fields/EX' exists" << endl;
        }
        if (Input::List().o_Ey) {
            if (Makefolder("output/fields/Ey") != 0)
                cout<<"Warning: Folder 'output/fields/EY' exists" << endl;
        }
        if (Input::List().o_Ez) {
            if (Makefolder("output/fields/Ez") != 0)
                cout<<"Warning: Folder 'output/fields/EZ' exists" << endl;
        }
        if (Input::List().o_Bx) {
            if (Makefolder("output/fields/Bx") != 0)
                cout<<"Warning: Folder 'output/fields/BX' exists" << endl;
        }
        if (Input::List().o_By) {
            if (Makefolder("output/fields/By") != 0)
                cout<<"Warning: Folder 'output/fields/BY' exists" << endl;
        }
        if (Input::List().o_Bz) {
            if (Makefolder("output/fields/Bz") != 0)
                cout<<"Warning: Folder 'output/fields/BZ' exists" << endl;
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

        if (Makefolder("output/moments") != 0)
            cout<<"Warning: Folder 'output/moments' exists" << endl;

//          Relativistic Energy - Momentum Tensor  
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_x1x2) {
            if (Makefolder("output/moments/n") != 0)
                cout<<"Warning: Folder 'output/moments/N' exists" << endl;
        }
        if (Input::List().o_pth) {
            if (Makefolder("output/moments/Pth") != 0)
                cout<<"Warning: Folder 'output/moments/Pth' exists" << endl;
        }
        if (Input::List().o_G) {
            if (Makefolder("output/moments/Gam") != 0)
                cout<<"Warning: Folder 'output/moments/Gam' exists" << endl;
        }
        if (Input::List().o_Px) {
            if (Makefolder("output/moments/Px") != 0)
                cout<<"Warning: Folder 'output/moments/Px' exists" << endl;
        }
        if (Input::List().o_PxPx) {
            if (Makefolder("output/moments/PxPx") != 0)
                cout<<"Warning: Folder 'output/moments/PxPx' exists" << endl;
        }
        if (Input::List().o_Py) {
            if (Makefolder("output/moments/Py") != 0)
                cout<<"Warning: Folder 'output/moments/Py' exists" << endl;
        }
        if (Input::List().o_PxPy) {
            if (Makefolder("output/moments/PxPy") != 0)
                cout<<"Warning: Folder 'output/moments/PxPy' exists" << endl;
        }
        if (Input::List().o_PyPy) {
            if (Makefolder("output/moments/PyPy") != 0)
                cout<<"Warning: Folder 'output/moments/PyPy' exists" << endl;
        }
        if (Input::List().o_Pz) {
            if (Makefolder("output/moments/Pz") != 0)
                cout<<"Warning: Folder 'output/moments/Pz' exists" << endl;
        }
        if (Input::List().o_PxPz) {
            if (Makefolder("output/moments/PxPz") != 0)
                cout<<"Warning: Folder 'output/moments/PxPz' exists" << endl;
        }
        if (Input::List().o_PyPz) {
            if (Makefolder("output/moments/PyPz") != 0)
                cout<<"Warning: Folder 'output/moments/PyPz' exists" << endl;
        }
        if (Input::List().o_PzPz) {
            if (Makefolder("output/moments/PzPz") != 0)
                cout<<"Warning: Folder 'output/moments/PzPz' exists"  << endl;
        }

//          Current
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_Jx) {
            if (Makefolder("output/moments/Jx") != 0)
                cout<<"Warning: Folder 'output/moments/Jx' exists" << endl;
        }
        if (Input::List().o_Jy) {
            if (Makefolder("output/moments/Jy") != 0)
                cout<<"Warning: Folder 'output/moments/Jy' exists" << endl;
        }
        if (Input::List().o_Jz) {
            if (Makefolder("output/moments/Jz") != 0)
                cout<<"Warning: Folder 'output/moments/Jz' exists" << endl;
        }


//          Nonrelativistic output
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
        if (Input::List().o_Vx) {
            if (Makefolder("output/moments/Vx") != 0)
                cout<<"Warning: Folder 'output/moments/Vx' exists" << endl;
        }
        if (Input::List().o_VxVx) {
            if (Makefolder("output/moments/VxVx") != 0)
                cout<<"Warning: Folder 'output/moments/VxVx' exists" << endl;
        }
        if (Input::List().o_Vy) {
            if (Makefolder("output/moments/Vy") != 0)
                cout<<"Warning: Folder 'output/moments/Vy' exists" << endl;
        }
        if (Input::List().o_VxVy) {
            if (Makefolder("output/moments/VxVy") != 0)
                cout<<"Warning: Folder 'output/moments/VxVy' exists" << endl;
        }
        if (Input::List().o_VyVy) {
            if (Makefolder("output/moments/VyVy") != 0)
                cout<<"Warning: Folder 'output/moments/VyVy' exists" << endl;
        }
        if (Input::List().o_VxVz) {
            if (Makefolder("output/moments/VxVz") != 0)
                cout<<"Warning: Folder 'output/moments/VxVz' exists" << endl;
        }
        if (Input::List().o_VyVz) {
            if (Makefolder("output/moments/VyVz") != 0)
                cout<<"Warning: Folder 'output/moments/VyVz' exists" << endl;
        }
        if (Input::List().o_VzVz) {
            if (Makefolder("output/moments/VzVz") != 0)
                cout<<"Warning: Folder 'output/moments/VzVz' exists" << endl;
        }
        if (Input::List().o_Vsq) {
            if (Makefolder("output/moments/Vsq") != 0)
                cout<<"Warning: Folder 'output/moments/Vsq' exists" << endl;
        }
        if (Input::List().o_Temperature) {
//                 if (Makefolder("output/moments/T_eV") != 0) 
//                     cout<<"Warning: Folder 'output/moments/T_eV' exists" << endl;
            if (Makefolder("output/moments/T") != 0)
                cout<<"Warning: Folder 'output/moments/T' exists" << endl;
        }
        if (Input::List().o_Pressure) {
            if (Makefolder("output/moments/P_Mbar") != 0)
                cout<<"Warning: Folder 'output/moments/P_Mbar' exists" << endl;
        }
        if (Input::List().o_ND) {
            if (Makefolder("output/moments/ND") != 0)
                cout<<"Warning: Folder 'output/moments/ND' exists" << endl;
        }
        if (Input::List().o_Nu) {
            if (Makefolder("output/moments/Nu") != 0)
                cout<<"Warning: Folder 'output/moments/Nu' exists" << endl;
        }
        if (Input::List().o_Qx) {
            if (Makefolder("output/moments/Qx") != 0)
                cout<<"Warning: Folder 'output/moments/Qx' exists" << endl;
        }
        if (Input::List().o_Qy) {
            if (Makefolder("output/moments/Qy") != 0)
                cout<<"Warning: Folder 'output/moments/Qy' exists" << endl;
        }
        if (Input::List().o_Qz) {
            if (Makefolder("output/moments/Qz") != 0)
                cout<<"Warning: Folder 'output/moments/Qz' exists" << endl;
        }
        if (Input::List().o_vNx) {
            if (Makefolder("output/moments/vNx") != 0)
                cout<<"Warning: Folder 'output/moments/vNx' exists" << endl;
        }
        if (Input::List().o_vNy) {
            if (Makefolder("output/moments/vNy") != 0)
                cout<<"Warning: Folder 'output/moments/vNy' exists" << endl;
        }
        if (Input::List().o_vNz) {
            if (Makefolder("output/moments/vNz") != 0)
                cout<<"Warning: Folder 'output/moments/vNz' exists" << endl;
        }
        if (Input::List().o_Ux) {
            if (Makefolder("output/moments/Ux") != 0)
                cout<<"Warning: Folder 'output/moments/Ux' exists" << endl;
        }
        if (Input::List().o_Uy) {
            if (Makefolder("output/moments/Uy") != 0)
                cout<<"Warning: Folder 'output/moments/Uy' exists" << endl;
        }
        if (Input::List().o_Uz) {
            if (Makefolder("output/moments/Uz") != 0)
                cout<<"Warning: Folder 'output/moments/Ux' exists" << endl;
        }
        if (Input::List().o_Z) {
            if (Makefolder("output/moments/Z") != 0)
                cout<<"Warning: Folder 'output/moments/Z' exists" << endl;
        }
        if (Input::List().o_ni) {
            if (Makefolder("output/moments/ni") != 0)
                cout<<"Warning: Folder 'output/moments/ni' exists" << endl;
        }
        if (Input::List().o_Ti) {
            if (Makefolder("output/moments/Ti") != 0)
                cout<<"Warning: Folder 'output/moments/Ti' exists" << endl;
        }

        if (Input::List().particlepusher) {
            if (Makefolder("output/particles") != 0)
                cout<<"Warning: Folder 'output/particles' exists" << endl;
            if (Makefolder("output/particles/prtx") != 0)
                cout<<"Warning: Folder 'output/particles/prtx' exists" << endl;
            if (Makefolder("output/particles/prtpx") != 0)
                cout<<"Warning: Folder 'output/particles/prtpx' exists" << endl;
            if (Makefolder("output/particles/prtpy") != 0)
                cout<<"Warning: Folder 'output/particles/prtpy' exists" << endl;
            if (Makefolder("output/particles/prtpz") != 0)
                cout<<"Warning: Folder 'output/particles/prtpz' exists" << endl;
        }
    }

    if (  Input::List().o_p1x1 || Input::List().o_p2x1 || Input::List().o_p3x1 ||
      Input::List().o_p1p3x1 || Input::List().o_p1p2x1 || Input::List().o_p2p3x1 ||
      Input::List().o_f0x1 ||  Input::List().o_f10x1 ||  Input::List().o_f11x1 
      ||  Input::List().o_f20x1 || Input::List().o_fl0x1) {
        if (Makefolder("output/distributions") != 0)
            cout<<"Warning: Folder 'output/distributions' exists" << endl;

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

// //  Moments
    part.push_back("prtx");
    part.push_back("prtpx");
    part.push_back("prtpy");
    part.push_back("prtpz");

//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px");
        pvsx.push_back( "py");
        pvsx.push_back( "pz");
    }

//  f-x
    for (size_t s(0); s < species; ++s) {
        fvsx.push_back( "f0");//+stringify(s) );
        fvsx.push_back( "f10");
        fvsx.push_back( "f11");
        fvsx.push_back( "f20");
        fvsx.push_back( "fl0");
    }

//  p-p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "pxpy");
        pvsx.push_back( "pypz");
        pvsx.push_back( "pxpz");
    }

}



//--------------------------------------------------------------


// Definition of the output axis
// Constructor 
// Export_Files::oAxis::oAxis() : label(""), units(""), min(0.0), max(1.0), sz(3) {}
// Export_Files::oAxis::oAxis( const float _m, const float _M,
//     const size_t _sz)
// : label(""), units(""), min(_m), max(_M), sz(_sz) {}
// Export_Files::oAxis::oAxis(const string _l, const string _u, const float _m, const float _M,
//    const size_t _sz)
// : label(_l), units(_u), min(_m), max(_M), sz(_sz) {}
// // Copy constructor 
// Export_Files::oAxis::oAxis(const oAxis& other) {
//     label = other.label;
//     units = other.units;
//     min   = other.min;
//     max   = other.max;
//     sz    = other.sz;
// }
//--------------------------------------------------------------

//--------------------------------------------------------------
Export_Files::Header::Header(string _axis_units, string _quantity_units, string _outputDir)
: axis_units(_axis_units), quantity_units(_quantity_units), outputDir(_outputDir){}
// 1D header constructor
// Export_Files::Header::Header(oAxis _x,
//  string _Ql, float _Qc,
//  string _tl, string _tu, float _tc,
//  string _oD)
// : title(_Ql), titleC(_Qc),
// time(_tl),  timeU(_tu),  timeC(_tc),
// oDir(_oD) {
//     xyz_axis.push_back(_x);
// }

// // 2D header constructor
// Export_Files::Header::Header(oAxis _x, oAxis _y,
//  string _Ql, float _Qc,
//  string _tl, string _tu, float _tc,
//  string _oD)
// : title(_Ql),  time(_tl), timeU(_tu),
// titleC(_Qc), timeC(_tc),
// oDir(_oD) {
//     xyz_axis.push_back(_x);
//     xyz_axis.push_back(_y);
// }

// // 3D header constructor
// Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z,
//  string _Ql, float _Qc,
//  string _tl, string _tu, float _tc,
//  string _oD)
// : title(_Ql), time(_tl), timeU(_tu),
// titleC(_Qc), timeC(_tc),
// oDir(_oD) {
//     xyz_axis.push_back(_x);
//     xyz_axis.push_back(_y);
//     xyz_axis.push_back(_z);
// }

// // 3D header constructor
// Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z, oAxis _imre,
//  string _Ql, float _Qc,
//  string _tl, string _tu, float _tc,
//  string _oD)
// : title(_Ql), time(_tl), timeU(_tu),
// titleC(_Qc), timeC(_tc),
// oDir(_oD) {
//     xyz_axis.push_back(_x);
//     xyz_axis.push_back(_y);
//     xyz_axis.push_back(_z);
//     xyz_axis.push_back(_imre);
// }

// // xD header constructor
// Export_Files::Header::Header(vector< oAxis > _xyz,
//  string _Ql, float _Qc,
//  string _tl, string _tu, float _tc,
//  string _oD)
// : xyz_axis(_xyz), title(_Ql), time(_tl), timeU(_tu),
// titleC(_Qc), timeC(_tc),
// oDir(_oD) {}

// number of header dimensions
// size_t Export_Files::Header::dim() {
//     return xyz_axis.size();
// }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// valarray<double> Export_Files::Header::axis(const size_t i) {
// //    return Algorithms::MakeAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz);
//     return Algorithms::MakeCAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz);
// }
// string Export_Files::Header::label(const size_t i) {
//     return xyz_axis[i].label;
// }
// string Export_Files::Header::units(const size_t i) {
//     return xyz_axis[i].units;
// }

// string Export_Files::Header::Title_label() { return title; }
// float Export_Files::Header::Title_conv()  { return titleC; }
// string Export_Files::Header::Time_label()  { return time; }
// float Export_Files::Header::Time_conv()   { return timeC; }
// string Export_Files::Header::Directory()   { return oDir; }
//--------------------------------------------------------------


//--------------------------------------------------------------
// Constructor of the export facility for data structures
Export_Files::Xport::Xport(const Algorithms::AxisBundle<double>& _axis,
   const vector< string > oTags,
   string homedir)

{
    size_t species(_axis.pdim());
    DefaultTags dTags(species);

// //  xyz Axis
    size_t xloc(0);                  // Find the location of the right tag
    while ( ( xloc < dTags.space.size()-1 ) &&
        ( find(oTags.begin(),oTags.end(), dTags.space[xloc]) == oTags.end() ) ) {
        ++xloc;
    }
    // xyz[0].label = "x["+ formulary().Label(dTags.space[xloc]) +"]";
    string axis_units = dTags.space[xloc];

//  Tags for Fields -->
    for (size_t i(0); i < dTags.fld.size(); ++i) 
    {
        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fld[i]) != oTags.end() ) 
        {
            string nounits = dTags.fld[i].substr(0, dTags.fld[i].find("_"));
            string folder = homedir + "output/fields/" + nounits + "/";
                
            Hdr[nounits] = Header(axis_units, dTags.fld[i], folder);
        }
    } // <--

//  Tags for Moments -->
    for (size_t i(0); i < dTags.mom.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.mom[i]) != oTags.end() ) {

            string nounits = dTags.mom[i].substr(0, dTags.mom[i].find("_"));
            string folder = homedir + "output/moments/" + nounits + "/";
            //         Generate a header file for this tag
            Hdr[nounits] = Header(axis_units, dTags.mom[i], folder);
        }
    } // <--

//  Tags for Particles -->
    for (size_t i(0); i < dTags.part.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.part[i]) != oTags.end() ) {

            string nounits = dTags.part[i].substr(0, dTags.part[i].find("_"));
            string folder = homedir + "output/particles/" + nounits + "/";
            //         Generate a header file for this tag
        }
    } // <--


//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx.size(); ++i) {

        // for (size_t k(0); k < oTags.size(); ++k)
        //         std::cout << " \n \n k = " << oTags[k];

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.pvsx[i]) != oTags.end() ) {

            string folder = homedir + "output/distributions/" + dTags.pvsx[i] + "/";
            Makefolder(folder);


//          Generate a header file for this tag
//          For each 9 you have a different species 
            // Hdr[dTags.pvsx[i]] = Header( pr[i/3], xyz[0],
            //  "f"+stringify(i/3), 1.0, tlabel, tunits, tconv, folder);
            // std::cout << " \n k = " << dTags.pvsx[i] << "\n";
            // std::cout << " \n \n 11 ";
            Hdr[dTags.pvsx[i]] = Header(axis_units, dTags.pvsx[i], folder);
        }
    } //<--

//  Tags for f-x -->
    for (size_t i(0); i < dTags.fvsx.size(); ++i) {

        //     If this tag is an output tag
        if ( find(oTags.begin(),oTags.end(), dTags.fvsx[i]) != oTags.end() ) {

            string folder = homedir + "output/distributions/" + dTags.fvsx[i] + "/";
            Makefolder(folder);

//          Generate a header file for this tag
            // Hdr[dTags.fvsx[i]] = Header( pr[i/5], xyz[0], imre[0],
            //  "f", 1.0, tlabel, tunits, tconv, folder);
            Hdr[dTags.fvsx[i]] = Header(axis_units, dTags.fvsx[i], folder);

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

    if (!rank) Makefolder(hdir+"restart/");

//        if (Makefolder(hdir+"restart/") != 0) cout<<"Warning: Folder "<< hdir+"restart/"<<" exists\n";
}

//--------------------------------------------------------------
//  Read restart file
void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, State1D& Y, double time_start) {

//      Generate filename 
    string   filename(hdir+"restart/re_1D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ifstream  fin(filename.c_str(), ios::binary);

    if (fin)
    {
        fin.read((char *) &time,sizeof(time_start));
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
void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, State1D& Y, double time_dump) {

//      Generate filename 
    string   filename(hdir+"restart/re_1D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ofstream  fout(filename.c_str(), ios::binary);

    fout.write((char *) &time_dump,sizeof(time_dump));

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
void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, State2D& Y, double time_start) {

//      Generate filename 
    string   filename(hdir+"restart/re_2D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ifstream  fin(filename.c_str(), ios::binary);

    if (fin)
    {
        fin.read((char *) &time,sizeof(time_start));
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
void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, State2D& Y, double time_dump) {

//      Generate filename 
    string   filename(hdir+"restart/re_2D_");
    filename.append(rFextension(rank,re_step));

//      Open file
    ofstream  fout(filename.c_str(), ios::binary);


    fout.write((char *) &time_dump,sizeof(time_dump));
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

    size_t sz(((Nm+1)*(2*Nl-Nm+2))/2);

//  Generate the structure to save the polynomials
    for (size_t i(0);i < sz; ++i)
    {
        plegendre.push_back(Array2D<double>(px.size(),py.size()));
    }
    
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
                for (size_t l(0); l < Nl + 1; ++l) 
                {
                    for (size_t m(0); m < ((Nm < l) ? Nm : l) + 1; ++m) 
                    {
                        (plegendre)[k](ipx, ipy) = vL(l, m);
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
    for (size_t i(0); i < other.dim(); ++i) {
        (plegendre).push_back( other(i) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre2D::~PLegendre2D(){
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
Array2D<double>  Output_Data::harmonicvsposition::operator()(DistFunc1D& df, size_t l, size_t m, size_t x0, size_t s) {
    
    Array2D<double> fx1(Np(s),2); 

    for (size_t ip(0); ip < Np(s); ++ip) {
        //  At each location in |p|
        fx1(ip,0) = static_cast<double>( (df(l,m))(ip,x0).real() );
        fx1(ip,1) = static_cast<double>( (df(l,m))(ip,x0).imag() );
    }

    return fx1;
}

Array2D<double>  Output_Data::harmonicvsposition::operator()(DistFunc2D& df, size_t l, size_t m, size_t x0, size_t y0, size_t s) {
    
    Array2D<double> fx1(Np(s),2); 

    for (size_t ip(0); ip < Np(s); ++ip) {
        //  At each location in |p|
        fx1(ip,0) = static_cast<double>( (df(l,m))(ip,x0,y0).real() );
        fx1(ip,1) = static_cast<double>( (df(l,m))(ip,x0,y0).imag() );
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
    
        pout1D_p1.push_back( valarray<double>( _G.axis.Npx(s))); 
        pout1D_p2.push_back( valarray<double>( _G.axis.Npy(s))); 
        pout1D_p3.push_back( valarray<double>( _G.axis.Npz(s)));         
        pout2D_p1p2.push_back( Array2D<double>( _G.axis.Npx(s),_G.axis.Npy(s))); 
        pout2D_p1p3.push_back( Array2D<double>( _G.axis.Npx(s),_G.axis.Npz(s))); 
        pout2D_p2p3.push_back( Array2D<double>( _G.axis.Npy(s),_G.axis.Npz(s))); 

        pout3D.push_back( Array3D<double>( _G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
    }
    
    

    double py_sq, pz_sqppy_sq;
    size_t pindlt;

    for (size_t s(0); s < _G.axis.pdim(); ++s) 
    {
        pradius.push_back(Array3D<double>(_G.axis.Npx(s),_G.axis.Npy(s),_G.axis.Npz(s)));
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
                        
                        phi[s](ipy,ipz) = (atan2((_G.axis.pz(s))[ipz],(_G.axis.py(s))[ipy]) + M_PI);
                    }
                    else
                    {
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
valarray<double>  Output_Data::fulldistvsposition::p1(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout1D_p1[s][ipx] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]*dpz[s][ipz]);
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

                            pout1D_p1[s][ipx] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re*dpy[s][ipy]*dpz[s][ipz]);
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
valarray<double>  Output_Data::fulldistvsposition::p1(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

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
                        pout1D_p1[s][ipx] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]*dpz[s][ipz]);
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

                            pout1D_p1[s][ipx] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re*dpy[s][ipy]*dpz[s][ipz]);
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
valarray<double>  Output_Data::fulldistvsposition::p2(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout1D_p2[s][ipy] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpz[s][ipz]);
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

                            pout1D_p2[s][ipy] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpz[s][ipz]);
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
valarray<double>  Output_Data::fulldistvsposition::p2(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

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
                        pout1D_p2[s][ipy] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpz[s][ipz]);
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

                            pout1D_p2[s][ipy] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpz[s][ipz]);
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
valarray<double>  Output_Data::fulldistvsposition::p3(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout1D_p3[s][ipz] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpy[s][ipy]);
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

                            pout1D_p3[s][ipz] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpy[s][ipy]);
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
valarray<double>  Output_Data::fulldistvsposition::p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

    pout1D_p3[s] = 0.;
    
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
                        pout1D_p3[s][ipz] += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]*dpy[s][ipy]);
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

                            pout1D_p3[s][ipz] += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re* dpx[s][ipx]*dpy[s][ipy]);
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
Array2D<double>  Output_Data::fulldistvsposition::p1p2(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout2D_p1p2[s](ipx,ipy) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpz[s][ipz]);
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

                            pout2D_p1p2[s](ipx,ipy) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpz[s][ipz]);
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
Array2D<double>  Output_Data::fulldistvsposition::p1p2(DistFunc2D& df, size_t x0,  size_t y0, size_t s) {

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
                        pout2D_p1p2[s](ipx,ipy) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpz[s][ipz]);
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

                            pout2D_p1p2[s](ipx,ipy) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpz[s][ipz]);
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
Array2D<double>  Output_Data::fulldistvsposition::p2p3(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout2D_p2p3[s](ipy,ipz) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]);
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

                            pout2D_p2p3[s](ipy,ipz) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpx[s][ipx]);
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
Array2D<double>  Output_Data::fulldistvsposition::p2p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

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
                        pout2D_p2p3[s](ipy,ipz) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpx[s][ipx]);
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

                            pout2D_p2p3[s](ipy,ipz) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpx[s][ipx]);
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
Array2D<double>  Output_Data::fulldistvsposition::p1p3(DistFunc1D& df, size_t x0, size_t s) {

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
                        pout2D_p1p3[s](ipx,ipz) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]);
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

                            pout2D_p1p3[s](ipx,ipz) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpy[s][ipy]);
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
Array2D<double>  Output_Data::fulldistvsposition::p1p3(DistFunc2D& df, size_t x0, size_t y0, size_t s) {

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
                        pout2D_p1p3[s](ipx,ipz) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy) * dpy[s][ipy]);
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

                            pout2D_p1p3[s](ipx,ipz) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re * dpy[s][ipy]);
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
Array3D<double>  Output_Data::fulldistvsposition::p1p2p3(DistFunc1D& df, size_t x0, size_t s){
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
                        pout3D[s](ipx,ipy,ipz) += static_cast<double>(YSH_re * (PL2D[s])(i_dist)(ipx,ipy));
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

                            pout3D[s](ipx,ipy,ipz) += static_cast<double>(2.0*(PL2D[s])(i_dist)(ipx,ipy)*YSH_re);//*dpy[s]*dpz[s];
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

void Output_Data::Output_Preprocessor::operator()(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    if (Input::List().o_Ex) {
        Ex( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ey) {
        Ey( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ez) {
        Ez( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Bx) {
        Bx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_By) {
        By( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Bz) {
        Bz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_x1x2) {
        n( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Temperature) {
        T( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jx) {
        Jx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jy) {
        Jy( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jz) {
        Jz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qx) {
        Qx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qy) {
        Qy( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qz) {
        Qz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNx) {
        vNx( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNy) {
        vNy( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNz) {
        vNz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ux) {
        Ux( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Uy) {
        Uy( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Uz) {
        Uz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Z) {
        Z( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_ni) {
        ni( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ux) {
        Ti( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().particlepusher) {
        particles_x( Y, grid, tout, time, dt, PE );
        particles_px( Y, grid, tout, time, dt, PE );
        particles_py( Y, grid, tout, time, dt, PE );
        particles_pz( Y, grid, tout, time, dt, PE );
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::distdump(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
   const Parallel_Environment_1D& PE) 
{

    if (Input::List().o_p1x1){
        px( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p2x1){
        py( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f0x1){
        f0( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f10x1){
        f10( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f11x1){
        f11( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f20x1)
    {
        f20( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_fl0x1){
        fl0( Y, grid, tout, time, dt, PE );
    }

}

void Output_Data::Output_Preprocessor::bigdistdump(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
   const Parallel_Environment_1D& PE) 
{

    if (Input::List().o_p1p2x1)
    {
        pxpy( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p1p3x1)
    {
        pxpz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p2p3x1)
    {
        pypz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p1p2p3x1)
    {
        // pxpypz( Y, grid, tout, time, dt, PE );
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::operator()(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_2D& PE) {

    if (Input::List().o_Ex) {
        Ex( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ey) {
        Ey( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Ez) {
        Ez( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Bx) {
        Bx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_By) {
        By( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Bz) {
        Bz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_x1x2) {
        n( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Temperature) {
        T( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jx) {
        Jx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jy) {
        Jy( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Jz) {
        Jz( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qx) {
        Qx( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qy) {
        Qy( Y, grid, tout, time, dt, PE );
    }

    if (Input::List().o_Qz) {
        Qz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNx) {
        vNx( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNy) {
        vNy( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_vNz) {
        vNz( Y, grid, tout, time, dt, PE );
    }

    // if (Input::List().o_Ux) {
    //     Ux( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().o_Uy) {
    //     Uy( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().o_Uz) {
    //     Uz( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().o_Z) {
    //     Z( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().o_ni) {
    //     ni( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().o_Ux) {
    //     Ti( Y, grid, tout, time, dt, PE );
    // }

    // if (Input::List().particlepusher) {
    //     particles_x( Y, grid, tout, time, dt, PE );
    //     particles_px( Y, grid, tout, time, dt, PE );
    //     particles_py( Y, grid, tout, time, dt, PE );
    //     particles_pz( Y, grid, tout, time, dt, PE );
    // }

}
//--------------------------------------------------------------
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::distdump(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
   const Parallel_Environment_2D& PE) 
{
    if (Input::List().o_p1x1){
        px( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p2x1){
        py( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f0x1){
        f0( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f10x1){
        f10( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f11x1){
        f11( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_f20x1)
    {
        f20( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_fl0x1){
        fl0( Y, grid, tout, time, dt, PE );
    }
}

void Output_Data::Output_Preprocessor::bigdistdump(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
   const Parallel_Environment_2D& PE) 
{
    if (Input::List().o_p1p2x1)
    {
        pxpy( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p1p3x1)
    {
        pxpz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p2p3x1)
    {
        pypz( Y, grid, tout, time, dt, PE );
    }
    if (Input::List().o_p1p2p3x1)
    {
        // pxpypz( Y, grid, tout, time, dt, PE );
    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ex(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); 
    
    double Exbuf[msg_sz];
    vector<double> ExGlobal(outNxGlobal,0.);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Exbuf[i] = static_cast<double>( Y.EMF().Ex()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Exbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ExGlobal[i] = Exbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Exbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Ex", xaxis, ExGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ey(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double Eybuf[msg_sz];
    
    vector<double> EyGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Eybuf[i] = static_cast<double>( Y.EMF().Ey()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Eybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EyGlobal[i] = Eybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Eybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Ey", xaxis, EyGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ez(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double Ezbuf[msg_sz];
    
    vector<double> EzGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Ezbuf[i] = static_cast<double>( Y.EMF().Ez()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Ezbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                EzGlobal[i] = Ezbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Ezbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Ez", xaxis, EzGlobal, tout, time, dt);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Bxbuf[msg_sz];
    vector<double> BxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Bxbuf[i] = static_cast<double>( Y.EMF().Bx()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BxGlobal[i] = Bxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Bxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Bx", xaxis, BxGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::By(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Bybuf[msg_sz];
    vector<double> ByGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Bybuf[i] = static_cast<double>( Y.EMF().By()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ByGlobal[i] = Bybuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Bybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("By", xaxis, ByGlobal, tout, time, dt);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Bzbuf[msg_sz];
    vector<double> BzGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Bzbuf[i] = static_cast<double>( Y.EMF().Bz()(Nbc+i).real() );
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bzbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                BzGlobal[i] = Bzbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Bzbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Bz", xaxis, BzGlobal, tout, time, dt);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ex(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal);
    double Exbuf[msg_sz];
    Array2D<double> ExGlobal(outNxGlobal,outNyGlobal);
    size_t i(0);
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Exbuf[i] = static_cast<double>( ( Y.EMF().Ex()(Nbc+ix,Nbc+iy).real() ));
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Exbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Exbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("Ex", xaxis, yaxis, ExGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ey(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); 
    double Eybuf[msg_sz];
    Array2D<double> EyGlobal(outNxGlobal,outNyGlobal);

    size_t i(0);
    int rankx,ranky; 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Eybuf[i] = static_cast<double>( Y.EMF().Ey()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Eybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Eybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("Ey", xaxis, yaxis, EyGlobal, tout, time, dt);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ez(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    double Ezbuf[msg_sz];
    Array2D<double> EzGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky; 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Ezbuf[i] = static_cast<double>( Y.EMF().Ez()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Ezbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Ezbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("Ez", xaxis, yaxis, EzGlobal, tout, time, dt);


}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    double Bxbuf[msg_sz];
    Array2D<double> BxGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky; 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bxbuf[i] = static_cast<double>( Y.EMF().Bx()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Bxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("Bx", xaxis, yaxis, BxGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::By(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    double Bybuf[msg_sz];
    Array2D<double> ByGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());
    size_t i(0);
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bybuf[i] = static_cast<double>( Y.EMF().By()(Nbc+ix,Nbc+iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Bybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("By", xaxis, yaxis, ByGlobal, tout, time, dt);

}
//--------------------------------------------------------------     
//--------------------------------------------------------------
////--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Bz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal);
    double Bzbuf[msg_sz];
    Array2D<double> BzGlobal(outNxGlobal,outNyGlobal);
    size_t i(0);
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(size_t ix(0); ix < outNxLocal; ++ix) {
        for(size_t iy(0); iy < outNyLocal; ++iy) {
            Bzbuf[i] = static_cast<double>( Y.EMF().Bz()(Nbc + ix, Nbc + iy).real() );
            ++i;
        }
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Bzbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                MPI_Recv(Bzbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

    if (PE.RANK() == 0) expo.Export_h5("Bz", xaxis, yaxis, BzGlobal, tout, time, dt);

}
//--------------------------------------------------------------   
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::px(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npx(grid.axis.Npx(s));
        int msg_sz(outNxLocal*Npx);
        
        Array2D<double> p1x1Global(Npx,outNxGlobal); 
        vector<double> p1axis(valtovec(grid.axis.px(s)));

        double pxbuf[Npx*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) 
        {
            valarray<double> data1D = p_x.p1( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npx; ++j) {
                pxbuf[j+i*Npx]=data1D[j];
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                    MPI_Recv(pxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("px", p1axis, xaxis, p1x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::py(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npy(grid.axis.Npy(s));
        int msg_sz(outNxLocal*Npy);
        Array2D<double> p2x1Global(Npy,outNxGlobal); 

        vector<double> p2axis(valtovec(grid.axis.py(s)));


        double pybuf[Npy*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) 
        {
            valarray<double> data1D = p_x.p2( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npy; ++j) {
                pybuf[j+i*Npy]=data1D[j];
            }       
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                    MPI_Recv(pybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < Npy; ++j) {
                            p2x1Global(j,i + outNxLocal*rr) = pybuf[j+i*Npy];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) 
            {
                for (size_t j(0); j < Npy; ++j) 
                {
                    p2x1Global(j,i) = pybuf[j+i*Npy];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("py", p2axis, xaxis, p2x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npz(grid.axis.Npz(s));
        int msg_sz(outNxLocal*Npz);
        
        Array2D<double> p3x1Global(Npz,outNxGlobal); 
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        double pzbuf[Npz*outNxLocal];

        for (size_t i(0); i < outNxLocal; ++i) {

            valarray<double> data1D = p_x.p3( Y.DF(s), i+Nbc, s);

            for (size_t j(0); j < Npz; ++j) {
                pzbuf[j+i*Npz]=data1D[j];
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pzbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                    MPI_Recv(pzbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("pz", p3axis, xaxis, p3x1Global, tout, time, dt, s);

    }

}
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::px(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npx(grid.axis.Npx(s));
        int msg_sz(outNxLocal*outNyLocal*Npx);
        Array3D<double> p1x1Global(Npx,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        double pxbuf[Npx*outNxLocal*outNyLocal];

        vector<double> p1axis(valtovec(grid.axis.px(s)));

        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {

                valarray<double> data1D = p_x.p1( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npx; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(pxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("px", p1axis,xaxis, yaxis, p1x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::py(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npy(grid.axis.Npy(s));
        int msg_sz(outNxLocal*outNyLocal*Npy);
        Array3D<double> p1x1Global(Npy,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        double pxbuf[Npy*outNxLocal*outNyLocal];
        vector<double> p2axis(valtovec(grid.axis.py(s)));
        

        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                valarray<double> data1D = p_x.p2( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npy; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) 
        {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(pxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("py", p2axis, xaxis, yaxis, p1x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//-----------------------------------------------------------------------------------
void Output_Data::Output_Preprocessor::pz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_2D& PE) {
    // std::cout << "0 \n";
    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(int s(0); s < Y.Species(); ++s) {
        size_t Npz(grid.axis.Npz(s));
        int msg_sz(outNxLocal*outNyLocal*Npz);
        Array3D<double> p1x1Global(Npz,outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        double pxbuf[Npz*outNxLocal*outNyLocal];
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        size_t counter(0);

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {

                valarray<double> data1D = p_x.p3( Y.DF(s), ix+Nbc, iy + Nbc, s);

                for (size_t j(0); j < Npz; ++j) {
                    pxbuf[counter]=data1D[j];
                    ++counter;
                }
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(pxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(pxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("pz", p3axis, xaxis, yaxis, p1x1Global, tout, time, dt, s);

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
void Output_Data::Output_Preprocessor::pxpy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) 
{

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    
    for(int s(0); s < Y.Species(); ++s) 
    {
        
        int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npy(s));        
        Array3D<double> pxpyGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal);
        
        double pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal];        
        vector<double> p1axis(valtovec(grid.axis.px(s)));
        vector<double> p2axis(valtovec(grid.axis.py(s)));

        ind = 0;    
        for (size_t i(0); i < outNxLocal; ++i) 
        {
            Array2D<double> data2D = p_x.p1p2( Y.DF(s), i+Nbc, s);
            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
            {

                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                {
                    pbuf[ind]=data2D(j,k);
                    ++ind;
                }

            }
            
        }

        if (PE.MPI_Processes() > 1) 
        {
           if (PE.RANK()!=0) {
               MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                            pxpyGlobal(j,k,i) = pbuf[ind];
                            ++ind;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    ind = 0;
                    for(size_t i(0); i < outNxLocal; i++) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                            {
                                pxpyGlobal(j,k,i + outNxLocal*rr) = pbuf[ind];
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
                       pxpyGlobal(j,k,i) = pbuf[ind];
                       ++ind;
                   }
               }
           }
       }

       if (PE.RANK() == 0) expo.Export_h5("pxpy", p1axis, p2axis, xaxis, pxpyGlobal, tout, time, dt, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    size_t counter;

    for(int s(0); s < Y.Species(); ++s) 
    {

        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<double> pxpyGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        double pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        vector<double> p1axis(valtovec(grid.axis.px(s)));
        vector<double> p2axis(valtovec(grid.axis.py(s)));

        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<double> data2D = p_x.p1p2( Y.DF(s), ix+Nbc, iy+Nbc, s);

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
                MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                                pxpyGlobal(j,k,ix,iy) = pbuf[counter];
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

                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    pxpyGlobal(j,k,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
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
                            pxpyGlobal(j,k,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pxpy", p1axis, p2axis, xaxis, yaxis, pxpyGlobal, tout, time, dt, s);

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
void Output_Data::Output_Preprocessor::pypz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(outNxLocal*grid.axis.Npy(s)*grid.axis.Npz(s));
        Array3D<double> pypzGlobal(grid.axis.Npy(s),grid.axis.Npz(s),outNxGlobal);
        double pbuf[grid.axis.Npy(s)*grid.axis.Npz(s)*outNxLocal];

        ind = 0;

        vector<double> p2axis(valtovec(grid.axis.py(s)));
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = p_x.p2p3( Y.DF(s), i+Nbc, s);

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
               MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                           pypzGlobal(j,k,i) = pbuf[ind];
                           ++ind;
                        }
                    }
                }
                // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    ind = 0;
                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < grid.axis.Npy(s); ++j) {
                            for (size_t k(0); k < grid.axis.Npz(s); ++k) {
                                pypzGlobal(j,k,i + outNxLocal*rr) = pbuf[ind];
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
                        pypzGlobal(j,k,i) = pbuf[ind];
                        ++ind;
                    }
                }
            }
        }

       if (PE.RANK() == 0) expo.Export_h5("pypz", p2axis, p3axis, xaxis, pypzGlobal, tout, time, dt, s);

   }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pypz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    size_t counter;


    for(int s(0); s < Y.Species(); ++s) 
    {
        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<double> dataGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); 
        double pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        vector<double> p2axis(valtovec(grid.axis.py(s)));
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<double> data2D = p_x.p2p3( Y.DF(s), ix+Nbc, iy+Nbc, s);

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
                MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                                dataGlobal(j,k,ix,iy) = pbuf[counter];
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

                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    dataGlobal(j,k,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
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
                            dataGlobal(j,k,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pypz", p2axis, p3axis, xaxis, yaxis, dataGlobal, tout, time, dt, s);

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
void Output_Data::Output_Preprocessor::pxpz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);
    vector<double> xaxis(valtovec(grid.axis.x(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npz(s));
        Array3D<double> pxpzGlobal(grid.axis.Npx(s),grid.axis.Npz(s),outNxGlobal); 
        double pbuf[grid.axis.Npx(s)*grid.axis.Npz(s)*outNxLocal];
        vector<double> p1axis(valtovec(grid.axis.px(s)));
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = p_x.p1p2( Y.DF(s), i+Nbc, s);
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
               MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                            pxpzGlobal(j,k,i) = pbuf[ind];
                            ++ind;
                        }
                    }
                }
               // Fill data for rank > 0
                for (int rr(1); rr < PE.MPI_Processes(); ++rr){
                    ind = 0;
                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    for(size_t i(0); i < outNxLocal; i++) 
                    {
                        for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                        {
                            for (size_t k(0); k < grid.axis.Npz(s); ++k) 
                            {
                                pxpzGlobal(j,k,i + outNxLocal*rr) = pbuf[ind];
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
                        pxpzGlobal(j,k,i) = pbuf[ind];
                        ++ind;
                    }
                }
            }
        }

       if (PE.RANK() == 0) expo.Export_h5("pxpz", p1axis, p3axis, xaxis, pxpzGlobal, tout, time, dt, s);

   }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::pxpz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_2D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    size_t counter;


    for(int s(0); s < Y.Species(); ++s) 
    {

        int msg_sz(outNxLocal*outNyLocal*grid.axis.Npx(s)*grid.axis.Npy(s));
        Array4D<double> dataGlobal(grid.axis.Npx(s),grid.axis.Npy(s),outNxGlobal,outNyGlobal); //, yglob_axis.dim());
        double pbuf[grid.axis.Npx(s)*grid.axis.Npy(s)*outNxLocal*outNyLocal];
        vector<double> p1axis(valtovec(grid.axis.px(s)));
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        counter = 0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                Array2D<double> data2D = p_x.p1p3( Y.DF(s), ix+Nbc, iy+Nbc, s);

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
                MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                                dataGlobal(j,k,ix,iy) = pbuf[counter];
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

                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                    
                    for(size_t ix(0); ix < outNxLocal; ++ix) 
                    {
                        for(size_t iy(0); iy < outNyLocal; ++iy)
                        {
                            for (size_t j(0); j < grid.axis.Npx(s); ++j) 
                            {
                                for (size_t k(0); k < grid.axis.Npy(s); ++k) 
                                {
                                    dataGlobal(j,k,ix + outNxLocal*rankx,iy + outNyLocal*ranky) = pbuf[counter];
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
                            dataGlobal(j,k,ix,iy) = pbuf[counter];
                            ++counter;
                        }
                    }
                }
            }
        }
       

       if (PE.RANK() == 0) expo.Export_h5("pxpz", p1axis, p3axis, xaxis, yaxis, dataGlobal, tout, time, dt, s);

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
void Output_Data::Output_Preprocessor::pxpypz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t ind(0);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(outNxLocal*grid.axis.Npx(s)*grid.axis.Npy(s)*grid.axis.Npz(s));
        Array4D<double> pxpypzGlobal(grid.axis.Npx(s),grid.axis.Npy(s),grid.axis.Npz(s),outNxGlobal); 
        double pbuf[msg_sz];
        vector<double> p1axis(valtovec(grid.axis.px(s)));
        vector<double> p2axis(valtovec(grid.axis.py(s)));
        vector<double> p3axis(valtovec(grid.axis.pz(s)));

        ind = 0;

        for (size_t i(0); i < outNxLocal; ++i) {

            Array3D<double> data3D = p_x.p1p2p3( Y.DF(s), i+Nbc, s);

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
                MPI_Send(pbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                    MPI_Recv(pbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("pxpypz", p1axis, p2axis, p3axis, xaxis, pxpypzGlobal, tout, time, dt, s);
    }

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f0(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {
        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<double> f0x1Global(outNxGlobal,f_x.Np(s),2); 

        vector<double> paxis(valtovec(grid.axis.p(s)));

        double f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = f_x( Y.DF(s),0,0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(i + outNxLocal*rr,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(i + outNxLocal*rr,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else 
        {
            for(size_t i(0); i < outNxGlobal; i++) 
            {
                for (size_t j(0); j < f_x.Np(s); ++j) 
                {
                    f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f0", xaxis, paxis, re_im_axis, f0x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f10(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    
    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);
    
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<double> f0x1Global(outNxGlobal,f_x.Np(s),2); //, yglob_axis.dim());
        vector<double> paxis(valtovec(grid.axis.p(s)));

        double f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = f_x( Y.DF(s), 1, 0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) 
            {
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(i + outNxLocal*rr,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(i + outNxLocal*rr,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) 
            {
                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f10", xaxis, paxis, re_im_axis, f0x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f11(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<double> f0x1Global(outNxGlobal,f_x.Np(s),2); 
        vector<double> paxis(valtovec(grid.axis.p(s)));

        double f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = f_x( Y.DF(s),1,1, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) {
                // std::cout << "\n f0(" << i << "," << j << ") = (" << data2D(j,1) << "," << data2D(j,2) <<")";
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);

                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(i + outNxLocal*rr,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(i + outNxLocal*rr,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                
                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("f11", xaxis, paxis, re_im_axis, f0x1Global, tout, time, dt, s);
    }


}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f20(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {
        if (Y.DF(s).l0() > 1)
        {
            int msg_sz(2*outNxLocal*f_x.Np(s));
            Array3D<double> f0x1Global(outNxGlobal,f_x.Np(s),2); //, yglob_axis.dim());
            vector<double> paxis(valtovec(grid.axis.p(s)));
            double f0xbuf[msg_sz];

            for (size_t i(0); i < outNxLocal; ++i) {

                Array2D<double> data2D = f_x( Y.DF(s), 2, 0, i+Nbc, s);

                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);
                    f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
                }

            }

            if (PE.MPI_Processes() > 1) {
                if (PE.RANK()!=0) {
                    MPI_Send(f0xbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    for(size_t i(0); i < outNxLocal; i++) {

                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(f0xbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                        for(size_t i(0); i < outNxLocal; i++) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                f0x1Global(i + outNxLocal*rr,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                                f0x1Global(i + outNxLocal*rr,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                            }
                        }
                    }
                }
            }
            else {
                for(size_t i(0); i < outNxGlobal; i++) {
                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f20", xaxis, paxis, re_im_axis, f0x1Global, tout, time, dt, s);

        }
    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::fl0(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {

    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        int msg_sz(2*outNxLocal*f_x.Np(s));
        Array3D<double> f0x1Global(outNxGlobal,f_x.Np(s),2);
        vector<double> paxis(valtovec(grid.axis.p(s)));

        double f0xbuf[msg_sz];

        for (size_t i(0); i < outNxLocal; ++i) {

            Array2D<double> data2D = f_x( Y.DF(s), Y.DF(s).l0(),0, i+Nbc, s);

            for (size_t j(0); j < f_x.Np(s); ++j) 
            {
                f0xbuf[2*j+   2*i*f_x.Np(s)]=data2D(j,0);
                f0xbuf[2*j+1+ 2*i*f_x.Np(s)]=data2D(j,1);
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(f0xbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {

                    for (size_t j(0); j < f_x.Np(s); ++j) {
                        f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                        f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                    }
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(f0xbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

                    for(size_t i(0); i < outNxLocal; i++) {
                        for (size_t j(0); j < f_x.Np(s); ++j) {
                            f0x1Global(i + outNxLocal*rr,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                            f0x1Global(i + outNxLocal*rr,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                        }
                    }
                }
            }
        }
        else {
            for(size_t i(0); i < outNxGlobal; i++) {
                for (size_t j(0); j < f_x.Np(s); ++j) {
                    f0x1Global(i,j,0) = f0xbuf[2*j+   2*i*f_x.Np(s)];
                    f0x1Global(i,j,1) = f0xbuf[2*j+1+ 2*i*f_x.Np(s)];
                }
            }
        }

        if (PE.RANK() == 0) expo.Export_h5("fl0", xaxis, paxis, re_im_axis, f0x1Global, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f0(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<double> global(outNxGlobal,outNyGlobal,f_x.Np(s),2); //, yglob_axis.dim());
        double buf[msg_sz];
        vector<double> paxis(valtovec(grid.axis.p(s)));

        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<double> data2D = f_x(Y.DF(s), 0, 0, ix + Nbc, iy + Nbc, s);

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
                    MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(ix, iy, j, 0) = buf[i];
                                global(ix, iy, j, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 0) = buf[i];
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 1) = buf[i+1];
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
                            global(ix, iy, j, 0) = buf[i];
                            global(ix, iy, j, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f0", xaxis, yaxis, paxis, re_im_axis, global, tout, time, dt, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
    void Output_Data::Output_Preprocessor::f10(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                                 const Parallel_Environment_2D& PE) {
        size_t Nbc(Input::List().BoundaryCells);
        MPI_Status status;
        size_t i(0);
        size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
        size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
        size_t outNxGlobal(grid.axis.Nxg(0));
        size_t outNyGlobal(grid.axis.Nxg(1));
        int rankx,ranky;
        vector<double> xaxis(valtovec(grid.axis.xg(0)));
        vector<double> yaxis(valtovec(grid.axis.xg(1)));

        vector<double> re_im_axis;
        re_im_axis.push_back(0.);
        re_im_axis.push_back(1.);


        for(int s(0); s < Y.Species(); ++s) {
            int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
            Array4D<double> global(outNxGlobal,outNyGlobal,f_x.Np(s),2);
            vector<double> paxis(valtovec(grid.axis.p(s)));

            double buf[msg_sz];
            i=0;
            for (size_t ix(0); ix < outNxLocal; ++ix) {
                for (size_t iy(0); iy < outNyLocal; ++iy) {
                    Array2D<double> data2D = f_x(Y.DF(s), 1, 0, ix + Nbc, iy + Nbc, s);
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
                        MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                    }
                    else {
                        // Fill data for rank = 0
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(ix, iy, j, 0) = buf[i];
                                    global(ix, iy, j, 1) = buf[i+1];
                                    i+=2;
                                }
                            }
                        }
                        // Fill data for rank > 0
                        for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                            MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                            rankx = rr % PE.MPI_X();
                            ranky = rr / PE.MPI_X();
                            i=0;
                            for (size_t ix(0); ix < outNxLocal; ++ix) {
                                for (size_t iy(0); iy < outNyLocal; ++iy) {
                                    for (size_t j(0); j < f_x.Np(s); ++j) {
                                        global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 0) = buf[i];
                                        global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 1) = buf[i+1];
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
                                global(ix, iy, j, 0) = buf[i];
                                global(ix, iy, j, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                }

                if (PE.RANK() == 0) expo.Export_h5("f10", xaxis, yaxis, paxis, re_im_axis, global, tout, time, dt, s);

            }

        }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f11(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<double> global(outNxGlobal,outNyGlobal,f_x.Np(s),2);
        double buf[msg_sz];
        vector<double> paxis(valtovec(grid.axis.p(s)));

        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<double> data2D = f_x(Y.DF(s), 1, 1, ix + Nbc, iy + Nbc, s);

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
                    MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(ix, iy, j, 0) = buf[i];
                                global(ix, iy, j, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 0) = buf[i];
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 1) = buf[i+1];
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
                            global(ix, iy, j, 0) = buf[i];
                            global(ix, iy, j, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f11", xaxis, yaxis, paxis, re_im_axis, global, tout, time, dt, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::f20(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

//    double buf[msg_sz];
//    Array2D<double> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<double> global(outNxGlobal,outNyGlobal,f_x.Np(s),2); //, yglob_axis.dim());
        double buf[msg_sz];
        vector<double> paxis(valtovec(grid.axis.p(s)));


        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<double> data2D = f_x(Y.DF(s), 2, 0, ix + Nbc, iy + Nbc, s);

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
                    MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(ix, iy, j, 0) = buf[i];
                                global(ix, iy, j, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 0) = buf[i];
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 1) = buf[i+1];
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
                            global(ix, iy, j, 0) = buf[i];
                            global(ix, iy, j, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("f20", xaxis, yaxis, paxis, re_im_axis, global, tout, time, dt, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::fl0(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    vector<double> re_im_axis;
    re_im_axis.push_back(0.);
    re_im_axis.push_back(1.);

//    double buf[msg_sz];
//    Array2D<double> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) {
        int msg_sz(2*outNxLocal*outNyLocal*f_x.Np(s));
        Array4D<double> global(outNxGlobal,outNyGlobal,f_x.Np(s),2); //, yglob_axis.dim());
        double buf[msg_sz];
        vector<double> paxis(valtovec(grid.axis.p(s)));

        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {

                Array2D<double> data2D = f_x(Y.DF(s), Y.DF(s).l0(), 0, ix + Nbc, iy + Nbc, s);

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
                    MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
                }
                else {
                    // Fill data for rank = 0
                    i=0;
                    for (size_t ix(0); ix < outNxLocal; ++ix) {
                        for (size_t iy(0); iy < outNyLocal; ++iy) {
                            for (size_t j(0); j < f_x.Np(s); ++j) {
                                global(ix, iy, j, 0) = buf[i];
                                global(ix, iy, j, 1) = buf[i+1];
                                i+=2;
                            }
                        }
                    }
                    // Fill data for rank > 0
                    for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                        MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
                        rankx = rr % PE.MPI_X();
                        ranky = rr / PE.MPI_X();
                        i=0;
                        for (size_t ix(0); ix < outNxLocal; ++ix) {
                            for (size_t iy(0); iy < outNyLocal; ++iy) {
                                for (size_t j(0); j < f_x.Np(s); ++j) {
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 0) = buf[i];
                                    global(ix + outNxLocal * rankx, iy + outNyLocal * ranky, j, 1) = buf[i+1];
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
                            global(ix, iy, j, 0) = buf[i];
                            global(ix, iy, j, 1) = buf[i+1];
                            i+=2;
                        }
                    }
                }
            }

            if (PE.RANK() == 0) expo.Export_h5("fl0", xaxis, yaxis, paxis, re_im_axis, global, tout, time, dt, s);

        }

    }
//-----------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::n(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
     
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double nbuf[msg_sz];
    vector<double> nGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
        valarray<double> pra( (grid.axis.p(s)));
        for(size_t i(0); i < msg_sz; ++i) {
            nbuf[i] = 4.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(nbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    nGlobal[i] = nbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(nbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("n", xaxis, nGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_x(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    double buf[msg_sz];
    vector<double> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().x(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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
    

    // if (PE.RANK() == 0) expo.Export_h5("prtx", pGlobal, tout, time, dt);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_px(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    double buf[msg_sz];
    vector<double> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().px(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    // if (PE.RANK() == 0) expo.Export_h5("prtpx", pGlobal, tout, time, dt);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_py(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    double buf[msg_sz];
    vector<double> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().py(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    // if (PE.RANK() == 0) expo.Export_h5("prtpy", pGlobal, tout, time, dt);

}
//--------------------------------------------------------------
//--------------------------------------------------------------    
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::particles_pz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    int msg_sz(Y.particles().numpar()) ; 
    double buf[msg_sz];
    vector<double> pGlobal(Y.particles().numpar()); 

    for (int ip(0); ip < Y.particles().numpar(); ++ip) {
        buf[ip] = Y.particles().pz(ip)* (double (Y.particles().ishere(ip)));
    }


    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < msg_sz; i++) {
                pGlobal[i] += buf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    // if (PE.RANK() == 0) expo.Export_h5("prtpz", pGlobal, tout, time, dt);

}
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::T(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    
    double tbuf[msg_sz];
    vector<double> tGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    double convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);

    for(int s(0); s < Y.Species(); ++s) {
        valarray<double> pra( (grid.axis.p(s)) );
        
        for(size_t i(0); i < msg_sz; ++i) {
            tbuf[i] = 4.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 4);
            tbuf[i] /= 3.0*4.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 2);

            tbuf[i] *= 1.0/Y.DF(s).mass();
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(tbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    tGlobal[i] = tbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(tbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("T", xaxis, tGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double Jxbuf[msg_sz];
    vector<double> JxGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {
        valarray<double> pra( (grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jxbuf[i] = Y.DF(s).q()*4.0/3.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JxGlobal[i] = Jxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Jxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Jx", xaxis, JxGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double Jybuf[msg_sz];
    vector<double> JyGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {
        valarray<double> pra( (grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            Jybuf[i] = Y.DF(s).q()*8.0/3.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jybuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JyGlobal[i] = Jybuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Jybuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Jy", xaxis, JyGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double Jzbuf[msg_sz];
    vector<double> JzGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {
        valarray<double> pra( (grid.axis.p(s)) );
        for(size_t i(0); i < msg_sz; ++i) {
            Jzbuf[i] = Y.DF(s).q()*-8.0/3.0*M_PI*Algorithms::moment(  vdouble_imag( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 3);
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Jzbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    JzGlobal[i] = Jzbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Jzbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Jz", xaxis, JzGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Qxbuf[msg_sz];

    vector<double> QxGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        for(size_t i(0); i < msg_sz; ++i) {

            Qxbuf[i] = 4.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_real( (Y.SH(s,1,0)).xVec(i+Nbc) ), pra, 5);
            Qxbuf[i] *= 0.5;
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Qx", xaxis, QxGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Qxbuf[msg_sz];
    vector<double> QxGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
        valarray<double> pra( (grid.axis.p(s)) );
        for(size_t i(0); i < msg_sz; ++i) {

            Qxbuf[i] = 8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_real( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
            Qxbuf[i] *= 0.5;

        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Qy", xaxis, QxGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;

    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);

    double Qxbuf[msg_sz];
    vector<double> QxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) 
    {

        valarray<double> pra( (grid.axis.p(s)) );
        
        for(size_t i(0); i < msg_sz; ++i) 
        {
            Qxbuf[i] = -8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_imag( (Y.SH(s,1,1)).xVec(i+Nbc) ), pra, 5);
            Qxbuf[i] *= 0.5;
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(Qxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    QxGlobal[i] = Qxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(Qxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("Qz", xaxis, QxGlobal, tout, time, dt, s);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNx(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double vNxbuf[msg_sz];
    vector<double> vNxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {

//        valarray<double> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<double> pra( (grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<double>( (1.0 / 6.0 * (Algorithms::moment(vdouble_real((Y.SH(s, 1, 0)).xVec(i + Nbc) ), pra, 6)
              / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))) );
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("vNx", xaxis, vNxGlobal, tout, time, dt, 0);

    }

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double vNxbuf[msg_sz];
    vector<double> vNxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<double> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );
        valarray<double> pra( (grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<double>( (2.0 / 6.0 * (Algorithms::moment(vdouble_real((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
              / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))) );
        }
        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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
        if (PE.RANK() == 0) expo.Export_h5("vNy", xaxis, vNxGlobal, tout, time, dt, 0);
    }
}
// --------------------------------------------------------------
// --------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
  const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal);
    double vNxbuf[msg_sz];
    vector<double> vNxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(int s(0); s < Y.Species(); ++s) {
//        valarray<double> pra( Algorithms::MakeAxis( f_x.Pmin(s), f_x.Pmax(s), f_x.Np(s) ) );

        valarray<double> pra( (grid.axis.p(s)) );

        for(size_t i(0); i < msg_sz; ++i) {
            vNxbuf[i] = static_cast<double>( (-2.0 / 6.0 * (Algorithms::moment(vdouble_imag((Y.SH(s, 1, 1)).xVec(i + Nbc) ), pra, 6)
               / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(i+Nbc) ), pra, 5))));
        }
        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(vNxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
            }
            else {
                // Fill data for rank = 0
                for(size_t i(0); i < outNxLocal; i++) {
                    vNxGlobal[i] = vNxbuf[i];
                }
                // Fill data for rank > 0
                for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                    MPI_Recv(vNxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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
        if (PE.RANK() == 0) expo.Export_h5("vNz", xaxis, vNxGlobal, tout, time, dt, 0);
    }
    
}
// --------------------------------------------------------------
// --------------------------------------------------------------
void Output_Data::Output_Preprocessor::n(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                            const Parallel_Environment_2D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc), outNyLocal(grid.axis.Nx(1) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0)), outNyGlobal(grid.axis.Nxg(1));

    int msg_sz(outNxLocal*outNyLocal); //*szy);
    double nbuf[msg_sz];
    Array2D<double> nGlobal(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));

    for(int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );

        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) {
            for(size_t iy(0); iy < outNyLocal; ++iy) {
                // std::cout << "f00n[" << i << "]=" << (Y.SH(s,0,0)).xVec(i+Nbc)[0] << "\n";
                nbuf[i] = 4.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 2);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(nbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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
                    MPI_Recv(nbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

        if (PE.RANK() == 0) expo.Export_h5("n", xaxis, yaxis, nGlobal, tout, time, dt, s);

    }

    

}
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor::T(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                            const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    int outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    int outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal*outNyLocal);

    double tbuf[msg_sz];
    Array2D<double> tGlobal(outNxGlobal,outNyGlobal);

    double convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);

    for(size_t s(0); s < Y.Species(); ++s) 
    {
        valarray<double> pra( (grid.axis.p(s)) );
        i=0;

        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                tbuf[i] = Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 4);
                tbuf[i] /= 3.0*Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 2);

                tbuf[i] *= 1.0/Y.DF(s).mass();

                // std::cout << "T1[" << ix << "," << iy << "] = " << tbuf[i] << "\n";
                
                ++i;
     
     
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(tbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(tbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("T", xaxis, yaxis, tGlobal, tout, time, dt, s);

    }

    

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal*outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal,outNyGlobal); //, yglob_axis.dim());

    for(int s(0); s < Y.Species(); ++s) 
    {
        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<double>(Y.DF(s).q()*4.0/3.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,1,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Jx", xaxis, yaxis, global, tout, time, dt, s);

    }

    

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal*outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal,outNyGlobal);

    for(int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for(size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for(size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<double>(Y.DF(s).q()*8.0/3.0*M_PI*Algorithms::moment(  vdouble_real( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK()!=0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Jy", xaxis, yaxis, global, tout, time, dt, s);

    }

    

}
// //--------------------------------------------------------------
// //--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Jz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {


    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<double>(Y.DF(s).q()*-8.0/3.0*M_PI*Algorithms::moment(  vdouble_imag( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 3));
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Jz", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*4.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_real( (Y.SH(s,1,0)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Qx", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_real( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Qy", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Qz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                             const Parallel_Environment_2D& PE) {
    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal);

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = 0.5*-8.0*M_PI/3.0*Y.DF(s).mass()*Algorithms::moment(  vdouble_imag( (Y.SH(s,1,1)).xVec(ix+Nbc,iy+Nbc) ), pra, 5);
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("Qz", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNx(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<double>( (1.0 / 6.0 * (Algorithms::moment(vdouble_real((Y.SH(s, 1, 0)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
              / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("vNx", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNy(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) {
            for (size_t iy(0); iy < outNyLocal; ++iy) {
                buf[i] = static_cast<double>( (2.0 / 6.0 * (Algorithms::moment(vdouble_real((Y.SH(s, 1, 1)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
              / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("vNy", xaxis, yaxis, global, tout, time, dt, s);

    }

    
}
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::vNz(const State2D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
                                              const Parallel_Environment_2D& PE) {

    size_t Nbc(Input::List().BoundaryCells);
    MPI_Status status;
    size_t i(0);
    size_t outNxLocal(grid.axis.Nx(0) - 2 * Nbc);
    size_t outNyLocal(grid.axis.Nx(1) - 2 * Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));
    size_t outNyGlobal(grid.axis.Nxg(1));
    int rankx,ranky;
    vector<double> xaxis(valtovec(grid.axis.xg(0)));
    vector<double> yaxis(valtovec(grid.axis.xg(1)));
    int msg_sz(outNxLocal * outNyLocal);

    double buf[msg_sz];
    Array2D<double> global(outNxGlobal, outNyGlobal); //, yglob_axis.dim());

    for (int s(0); s < Y.Species(); ++s) {

        valarray<double> pra( (grid.axis.p(s)) );
        i=0;
        for (size_t ix(0); ix < outNxLocal; ++ix) 
        {
            for (size_t iy(0); iy < outNyLocal; ++iy) 
            {
                buf[i] = static_cast<double>( (-2.0 / 6.0 * (Algorithms::moment(vdouble_imag((Y.SH(s, 1, 1)).xVec(ix + Nbc, iy + Nbc) ), pra, 6)
                        / Algorithms::moment(  vdouble_real( (Y.SH(s,0,0)).xVec(ix+Nbc, iy+Nbc) ), pra, 5))) );
                ++i;
            }
        }

        if (PE.MPI_Processes() > 1) {
            if (PE.RANK() != 0) {
                MPI_Send(buf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
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

                    MPI_Recv(buf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);

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

        if (PE.RANK() == 0) expo.Export_h5("vNz", xaxis, yaxis, global, tout, time, dt, s);

    }
}
//--------------------------------------------------------------
//--------------------------------------------------------------

//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ux(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); 
    double Uxbuf[msg_sz];
    vector<double> UxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) 
    {
        Uxbuf[i] = static_cast<double>(Y.HYDRO().vx(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Ux", xaxis, UxGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Uy(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); 
    double Uxbuf[msg_sz];
    vector<double> UxGlobal(outNxGlobal); 
    vector<double> xaxis(valtovec(grid.axis.xg(0)));  

    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<double>(Y.HYDRO().vy(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Uy", xaxis, UxGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Uz(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Uxbuf[msg_sz];
    vector<double> UxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<double>(Y.HYDRO().vz(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Uz", xaxis, UxGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Z(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
    const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Uxbuf[msg_sz];
    vector<double> UxGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));


    for(size_t i(0); i < msg_sz; ++i) {
        Uxbuf[i] = static_cast<double>(Y.HYDRO().Z(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Uxbuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                UxGlobal[i] = Uxbuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Uxbuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Z", xaxis, UxGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::ni(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double nibuf[msg_sz];
    vector<double> niGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        nibuf[i] = static_cast<double>(Y.HYDRO().density(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(nibuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                niGlobal[i] = nibuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(nibuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("ni",  xaxis, niGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
//
//--------------------------------------------------------------
//--------------------------------------------------------------
void Output_Data::Output_Preprocessor::Ti(const State1D& Y, const Grid_Info& grid, const size_t tout, const double time, const double dt,
 const Parallel_Environment_1D& PE) {


    size_t Nbc = Input::List().BoundaryCells;
    MPI_Status status;
    //  
    size_t outNxLocal(grid.axis.Nx(0) - 2*Nbc);
    size_t outNxGlobal(grid.axis.Nxg(0));

    int msg_sz(outNxLocal); //*szy);
    double Thydrobuf[msg_sz];
    vector<double> ThydroGlobal(outNxGlobal);
    vector<double> xaxis(valtovec(grid.axis.xg(0)));

    for(size_t i(0); i < msg_sz; ++i) {
        Thydrobuf[i] = static_cast<double>(511000.0/3.0*Y.HYDRO().temperature(i+Nbc));
    }

    if (PE.MPI_Processes() > 1) {
        if (PE.RANK()!=0) {
            MPI_Send(Thydrobuf, msg_sz, MPI_DOUBLE, 0, PE.RANK(), MPI_COMM_WORLD);
        }
        else {
            // Fill data for rank = 0
            for(size_t i(0); i < outNxLocal; i++) {
                ThydroGlobal[i] = Thydrobuf[i];
            }
            // Fill data for rank > 0
            for (int rr = 1; rr < PE.MPI_Processes(); ++rr){
                MPI_Recv(Thydrobuf, msg_sz, MPI_DOUBLE, rr, rr, MPI_COMM_WORLD, &status);
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

    if (PE.RANK() == 0) expo.Export_h5("Ti", xaxis, ThydroGlobal, tout, time, dt, 0);

}
//--------------------------------------------------------------
//--------------------------------------------------------------

//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
 std::vector<double> &axis1, std::vector<double> &data,
 const size_t  step, const double  time, const double  dt,
 const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

    //  Check Header file correctness
    // if (Hdr[tag].dim() != 1) {
    //     cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 1D structure\n";
    //     exit(1);
    // }

    //  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    // we create a new hdf5 file
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    // lets create a dataset of native double with the size of the vector
    // 'data'
    HighFive::DataSet dataset =
        file.createDataSet<double>(tag, HighFive::DataSpace::From(data));

    // lets write our vector of double to the HDF5 dataset
    dataset.write(data);

    add_attributes(dataset,tag,time,dt);

    HighFive::Group Axes = file.createGroup("Axes");

    HighFive::DataSet dataset_axis1 =
        Axes.createDataSet<double>("Axis1", HighFive::DataSpace::From(axis1));

    // lets write our vector of double to the HDF5 dataset
    dataset_axis1.write(axis1);

}
//--------------------------------------------------------------
// --------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
 std::vector<double> &axis1, std::vector<double> &axis2, Array2D<double> &dataA,
 const size_t  step, const double  time, const double  dt,
 const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

    //  Check Header file correctness
    // if (Hdr[tag].dim() != 2) {
    //     cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 2D structure\n";
    //     exit(1);
    // }

    //  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    // we create a new hdf5 file
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    // lets create a dataset of native double with the size of the Array2D
    // 'data'
    vector<double> dummyvec(dataA.dim2());
    vector<vector<double> > data;


    // exit(1);
    
    for (size_t i(0); i < dataA.dim1(); ++i)
    {
        for (size_t j(0); j < dataA.dim2(); ++j)
        {
            // std::cout << "(i,j) = " << i << "," << j << "\n";
            dummyvec[j] = dataA(i,j);
        }

        data.push_back(dummyvec);
    }

    

    std::vector<size_t> dims(2);
    dims[0] = dataA.dim1();
    dims[1] = dataA.dim2();
    
    
    HighFive::DataSet dataset =
        file.createDataSet<double>(tag, HighFive::DataSpace(dims));

    // lets write our vector of double to the HDF5 dataset
    dataset.write(data);

    add_attributes(dataset,tag,time,dt);

    HighFive::Group Axes = file.createGroup("Axes");

    HighFive::DataSet dataset_axis1 =
        Axes.createDataSet<double>("Axis1", HighFive::DataSpace::From(axis1));
    dataset_axis1.write(axis1);

    HighFive::DataSet dataset_axis2 =
    Axes.createDataSet<double>("Axis2", HighFive::DataSpace::From(axis2));
    dataset_axis2.write(axis2);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
 std::vector<double> &axis1, std::vector<double> &axis2, std::vector<double> &axis3,
 Array3D<double> &dataA,
 const size_t  step, const double  time, const double  dt,
 const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

    //  Check Header file correctness
    // if (Hdr[tag].dim() != 3) {
    //     cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n";
    //     exit(1);
    // }

    //  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    // we create a new hdf5 file
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    // lets create a dataset of native double with the size of the vector
    // 'data'
    
    vector<double> dummyvec(dataA.dim3());
    vector<vector<double> > dummyvec2(dataA.dim2());
    vector<vector<vector<double> > > data;

    for (size_t i(0); i < dataA.dim1(); ++i)
    {
        for (size_t j(0); j < dataA.dim2(); ++j)
        {
            for (size_t k(0); k < dataA.dim3(); ++k)
            {
                dummyvec[k] = dataA(i,j,k);
            }
            dummyvec2[j] = dummyvec;
        }
        data.push_back(dummyvec2);
    }

    std::vector<size_t> dims(3);
    dims[0] = dataA.dim1();
    dims[1] = dataA.dim2();
    dims[2] = dataA.dim3();
    
    HighFive::DataSet dataset =
        file.createDataSet<double>(tag, HighFive::DataSpace(dims));

    // lets write our vector of double to the HDF5 dataset
    dataset.write(data);

    /// Attributes
    /// Only a few defined for now
    /// dt, time
    add_attributes(dataset,tag,time,dt);

    HighFive::Group Axes = file.createGroup("Axes");

    HighFive::DataSet dataset_axis1 =
        Axes.createDataSet<double>("Axis1", HighFive::DataSpace::From(axis1));
    dataset_axis1.write(axis1);

    HighFive::DataSet dataset_axis2 =
        Axes.createDataSet<double>("Axis2", HighFive::DataSpace::From(axis2));
    dataset_axis2.write(axis2);

    HighFive::DataSet dataset_axis3 =
        Axes.createDataSet<double>("Axis3", HighFive::DataSpace::From(axis3));
    dataset_axis3.write(axis3);

}
//--------------------------------------------------------------
void Export_Files::Xport:: Export_h5(const std::string tag,
 std::vector<double> &axis1, std::vector<double> &axis2, std::vector<double> &axis3, std::vector<double> &axis4,
 Array4D<double> &dataA,
 const size_t  step, const double  time, const double  dt,
 const int spec){
//--------------------------------------------------------------
//  Export data to H5 file
//--------------------------------------------------------------

    string      filename(Hdr[tag].Directory());

    //  Check Header file correctness
    // if (Hdr[tag].dim() != 3) {
    //     cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n";
    //     exit(1);
    // }

    //  Open File
    filename.append(tag).append(oH5Fextension(step,spec));
    // we create a new hdf5 file
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    // lets create a dataset of native double with the size of the vector
    // 'data'
    
    vector<double> dummyvec(dataA.dim4());
    vector<vector<double> > dummyvec2(dataA.dim3());
    vector<vector<vector<double> > > dummyvec3(dataA.dim2());
    vector<vector<vector<vector<double> > > > data;

    for (size_t i(0); i < dataA.dim1(); ++i)
    {
        for (size_t j(0); j < dataA.dim2(); ++j)
        {
            for (size_t k(0); k < dataA.dim3(); ++k)
            {
                for (size_t l(0); k < dataA.dim4(); ++l)
                {
                    dummyvec[l] = dataA(i,j,k,l);
                }
                dummyvec2[k] = dummyvec;
            }
            dummyvec3[j] = dummyvec2;
        }
        data.push_back(dummyvec3);
    }

    std::vector<size_t> dims(4);
    dims[0] = dataA.dim1();
    dims[1] = dataA.dim2();
    dims[2] = dataA.dim3();
    dims[3] = dataA.dim4();
    
    HighFive::DataSet dataset =
        file.createDataSet<double>(tag, HighFive::DataSpace(dims));

    // lets write our vector of double to the HDF5 dataset
    dataset.write(data);

    /// Attributes
    /// Only a few defined for now
    /// dt, time
    add_attributes(dataset,tag,time,dt);

    HighFive::Group Axes = file.createGroup("Axes");

    HighFive::DataSet dataset_axis1 =
        Axes.createDataSet<double>("Axis1", HighFive::DataSpace::From(axis1));
    dataset_axis1.write(axis1);

    HighFive::DataSet dataset_axis2 =
        Axes.createDataSet<double>("Axis2", HighFive::DataSpace::From(axis2));
    dataset_axis2.write(axis2);

    HighFive::DataSet dataset_axis3 =
        Axes.createDataSet<double>("Axis3", HighFive::DataSpace::From(axis3));
    dataset_axis3.write(axis3);

    HighFive::DataSet dataset_axis4 =
        Axes.createDataSet<double>("Axis4", HighFive::DataSpace::From(axis4));
    dataset_axis4.write(axis4);

}
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Export_Files::Xport::add_attributes(HighFive::DataSet &dataset, const std::string tag, 
        const double  time, const double  dt) {
//--------------------------------------------------------------
//    Add initial attributes:
//    time step, iteration number, name of diagnostic,
//    physical time, time units, type, and max and min of axes ranges
//--------------------------------------------------------------

        // Now let's add a attribute on this dataset
        HighFive::Attribute adt = dataset.createAttribute<double>("dt", HighFive::DataSpace::From(dt));
        adt.write(dt);

        HighFive::Attribute at = dataset.createAttribute<double>("Time (c/\\omega_p)", HighFive::DataSpace::From(time));
        at.write(time);

        double timeps = time*formulary().Uconv("Time_ps");
        HighFive::Attribute at_ps = dataset.createAttribute<double>("Time (ps)", HighFive::DataSpace::From(timeps));
        at_ps.write(timeps);


        // HighFive::Attribute aname = dataset.createAttribute<std::string>("Name", HighFive::DataSpace::From(tag));
        // aname.write(tag);

        // std::string timeunits = "1/\\omega_0";
        // HighFive::Attribute atu = dataset.createAttribute<std::string>("Time Units", HighFive::DataSpace::From(timeunits));
        // atu.write(timeunits);
    }
//--------------------------------------------------------------


