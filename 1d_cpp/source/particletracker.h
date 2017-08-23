/*!\brief  Particle Tracker - Declarations
 * \author PICKSC
 * \file   particletracker.h
 *
 * Includes declarations for spatial advection, electric field advection, bfield and current
 *
 */

#ifndef DECL_PARTICLETRACKER_H
#define DECL_PARTICLETRACKER_H

/** \addtogroup vfp1d
 *  @{
 */
//--------------------------------------------------------------
//--------------------------------------------------------------
//  Spatial advection
class Particle_Pusher {
//--------------------------------------------------------------
public:
//      Constructors/Destructors
    Particle_Pusher(std::vector<double> xpos, 
                    std::vector<double> px,
                    std::vector<double> py,
                    std::vector<double> pz,
                    double xminLocal, double xmaxLocal, int NxLocal,
                    Particle1D& Pin);
//          Advance
    void push(State1D& Y, double dt);
    
private:

    void findnearestneighbortotheright(double position);
    
    size_t numpar;

    int nearestneighbortotheright;
    double normalizeddistancetoleft;

    double xmin, xmax, dx;

    double chargetomass;
};
//--------------------------------------------------------------
#endif