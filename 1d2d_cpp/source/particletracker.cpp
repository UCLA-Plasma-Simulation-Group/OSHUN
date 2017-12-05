/*!\brief  Particle Pusher Routines - Definitions
* \author  PICKSC
 * \date   March, 2017
 * \file   particletracker.cpp
 *
 * In here are the routines for the particle tracker object.
 *
 */

//  Standard libraries
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <cfloat>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

// Declarations
#include "input.h"
#include "state.h"
#include "particletracker.h"


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
//
Particle_Pusher::Particle_Pusher(std::vector<double> xpos, 
	std::vector<double> px, std::vector<double> py, std::vector<double> pz,
	double xminLocal, double xmaxLocal, int NxLocal,
	Particle1D& Pin): numpar(py.size()), 
	chargetomass(Pin.charge()/Pin.mass()),
	xmin(xminLocal),xmax(xmaxLocal),dx((xmaxLocal-xminLocal)/NxLocal)
{
		// std::cout << "\n xmax = " << xmax << "\n";
		// std::cout << "\n xmin = " << xmin << "\n";
	for (int ip(0); ip < numpar; ++ip){
		
		Pin.x(ip) = xpos[ip];
		if (Pin.x(ip) < xmax && Pin.x(ip) > xmin)
		{
			// std::cout << "\nParticle " << ip << " is in cell with xmin = " << xmin;
			Pin.ishere(ip) = 1;
		}
		else Pin.ishere(ip) = 0;

		Pin.px(ip) = px[ip];
		Pin.py(ip) = py[ip];
		Pin.pz(ip) = pz[ip];

		// std::cout << "min = " << xmin << ", and Pin.ishere(ip)" << Pin.ishere(ip) << "\n";
	}
}
/**
 * Pusher
 */
void Particle_Pusher::push(State1D& Y, double dt){

	for (int ip(0); ip < numpar; ++ip)
	{
		// std::cout << "min = " << xmin << ", position = " << Y.particles().x(ip) << ", Pin.ishere = " << Y.particles().ishere(ip) << "\n";
		if (Y.particles().ishere(ip))
		{	

			Y.particles().x(ip) = Y.particles().x(ip) + 0.5*dt*Y.particles().px(ip);
			// if (ip == 2)
			// std::cout << "\n ---- Particle Pusher for particle " << ip  << " \n";
			findnearestneighbortotheright(Y.particles().x(ip));

			// if (ip == 2)
			// std::cout << "nnr = " << nearestneighbortotheright << "\n";
			// if (ip == 2) 
			// std::cout << "xmax = " << xmin << "\n";


			Y.particles().px(ip) = Y.particles().px(ip) + dt*chargetomass*
			
			((normalizeddistancetoleft*(Y.EMF().Ex()(nearestneighbortotheright)).real() + 
			(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Ex()(nearestneighbortotheright-1)).real())) + 
			(
				Y.particles().py(ip)*((normalizeddistancetoleft*(Y.EMF().Bz()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Bz()(nearestneighbortotheright-1)).real())) 
			- 
				Y.particles().pz(ip)*((normalizeddistancetoleft*(Y.EMF().By()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().By()(nearestneighbortotheright-1)).real())) 
			);

			Y.particles().py(ip) = Y.particles().py(ip) + dt*chargetomass*
			
			((normalizeddistancetoleft*(Y.EMF().Ey()(nearestneighbortotheright)).real() + 
			(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Ey()(nearestneighbortotheright-1)).real())) + 
			(
				Y.particles().pz(ip)*((normalizeddistancetoleft*(Y.EMF().Bx()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Bx()(nearestneighbortotheright-1)).real())) 
			- 
				Y.particles().px(ip)*((normalizeddistancetoleft*(Y.EMF().Bz()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Bz()(nearestneighbortotheright-1)).real())) 
			);

			Y.particles().pz(ip) = Y.particles().pz(ip) + dt*chargetomass*
			
			((normalizeddistancetoleft*(Y.EMF().Ez()(nearestneighbortotheright)).real() + 
			(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Ez()(nearestneighbortotheright-1)).real())) + 
			(
				Y.particles().px(ip)*((normalizeddistancetoleft*(Y.EMF().By()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().By()(nearestneighbortotheright-1)).real())) 
			- 
				Y.particles().py(ip)*((normalizeddistancetoleft*(Y.EMF().Bx()(nearestneighbortotheright)).real() + 
				(1.0 - 	normalizeddistancetoleft)*(Y.EMF().Bx()(nearestneighbortotheright-1)).real())) 
			);

			// if (ip == 2) 
			// std::cout << "pxnew = " << Y.particles().px(ip) << "\n";


			Y.particles().x(ip) = Y.particles().x(ip) + 0.5*dt*Y.particles().px(ip);
			
			
			// exit(1);
			// dt*Y.particles().px(ip)
			if (Y.particles().x(ip) >= xmax )
			{	
				// std::cout << "\n Leaving right! \n";
				Y.particles().ishere(ip) = 0;
				Y.particles().goingright(ip) = 1;
			}
			else if (Y.particles().x(ip) < xmin)
			{
				// std::cout << "\n Leaving left! \n";
				Y.particles().ishere(ip) = 0;
				Y.particles().goingright(ip) = -1;
			}
			// if (ip == 2) 
			// std::cout << "xnew = " << Y.particles().x(ip) << "\n";
		}
	}
}
//
//
void Particle_Pusher::findnearestneighbortotheright(double position){
	nearestneighbortotheright = ceil((position - xmin) / dx);
	normalizeddistancetoleft = (position - xmin) / dx - (nearestneighbortotheright-1.0);
}