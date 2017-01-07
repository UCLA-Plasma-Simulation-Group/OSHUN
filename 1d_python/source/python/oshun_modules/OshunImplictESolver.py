# Cordon off all the implict stuff into this file for now then
# TODO: refactor timplict code

import numpy as np
import sys
import math


class RK2ForImplicit:
	def __init__(self, species, state, CFG, step_size):
		from Oshun1d import ElectricFieldOnDistribution1d, CurrentSolver, SpatialAdvection1d
		# maybe move these
		self.e_on_dist 					=		ElectricFieldOnDistribution1d	(species, state, CFG)
		#self.current_solver 				= 		CurrentSolver						(species, state, CFG)
		self.spatial_advection_solver = 		SpatialAdvection1d				(species, state, CFG)
		self.implicitESolver = ImplictECalcaTronicWonderMachine( species, state, CFG )

		self.h = step_size

		Yh = species.F
		self.trim_list = []
		self.do_filter = True
		cutoff = 20
		num_l = Yh.num_l
		self.time_step = 0
		self.sub_time_step = 0
		
		# this is the actuall results of micheal's filter...
		for l in xrange(0, num_l):
			if l > cutoff:
				threshold = cutoff
			elif l > 1:
				threshold = l-1
			else:
				threshold = 0
			#PRINT("%d,= %d" % (l, threshold))
			if(threshold > 0):
				self.trim_list.append( (l, 0, Yh.linear_index_for_harmonic( l), threshold) )

		pass
	
	# F: this is the function to evaluate (i.e. simulate)
	# results are of the form of adjustedments that are added directly to the ditribution-function-of-record
	#	for each species
	def rk2Core(self, F, spec, state, CFG):
		# initialize stuffs
		Y0 = spec.Y0; Yh = spec.Yh;
		# the 'real, final' state.
		Y  = spec.F				
		
		Y0.copy(Y);	

		# RK2 step 1
		# Yh = F(Y0)
		F(Y0, Yh, spec, state, CFG, Ein=state.E_field)

		# Y0 = Y0 + h*Yh
		Yh.scalar_mult(self.h);
		Y0.add(Yh);
		#Y = Y + (h/2)*Yh
		Yh.scalar_mult(0.5);
		Y.add(Yh);

		F(Y0, Yh, spec, state, CFG, Ein=state.E_field)
		# Y = y + (h/2)*Yh
		Yh.scalar_mult(0.5*self.h);
		Y.add(Yh);
	# the electric field is not used.
	def advance_p1(self, spec, state, CFG):
		self.rk2Core( self.Adv1, spec, state, CFG)
		pass
	# Ein is used (only as read-only)
	def advance_p2(self, spec, state, CFG):
		self.rk2Core( self.Adv2, spec, state, CFG)
		pass
	def advance_IEx(self, spec, state, CFG, Ein=None):
		self.rk2Core( self.ImpEx, spec, state, CFG)
		pass
	def advance_IEy(self, spec, state, CFG):
		# not neded for 1d
		#return Act.ImpEy(Yslope);
		self.rk2Core( self.ImpEy, spec, state, CFG)
		pass
	def advance_IEz(self, spec, state, CFG):
		# not neded for 1d
		# return Act.ImpEz(Yslope);
		self.rk2Core( self.ImpEz, F, spec, state, CFG)
		pass

	def Adv1(self, Y0, Yh, species, state, CFG, Ein=None):
		Yh.scalar_mult(0.0)
		self.spatial_advection_solver.calc(Y0, Yh, species, state, CFG)
	def Adv2(self, Y0, Yh, species, state, CFG, Ein=None):
		Yh.scalar_mult(0.0)
		self.e_on_dist.calc(Y0, Yh, species, state, CFG, Ein = Ein)
		for item in self.trim_list:
			idx = item[2]
			level_to_trim = item[3]
			for ix in xrange(0, Yh.num_total_cells[0]):
				Yh.cells[ix].harmonics[idx].momentum_projection[:level_to_trim] = 0.0
			#PRINT("zeros %d cells for %d,%d" % (  level_to_trim, item[0], item[1] ))
		pass
	def ImpEx(self, Y0, Yh, species, state, CFG, Ein=None):
		Yh.scalar_mult(0.0)
		self.e_on_dist.Implicit_Ex1d_calc(Y0, Yh, species, state, CFG, Ein = Ein)

	# called by main loop
	def advanceE(self, species, state, CFG, fokkerPlankImplictOperator, n, time):
		Yh = species.Yh;
		# the 'real, final' state.
		Y  = species.F				
		self.implicitESolver.advance( Y, Yh, species, state, CFG, fokkerPlankImplictOperator, self, n, time, self.h, Eout = None)

	#def ImpEy(self, Y0, Yh, species, state, CFG, Ein=None):
	#	Yh.scalar_mult(0.0)
	#	self.e_on_dist.Implicit_Ey1d_calc(Y0, Yh, species, state, CFG, Ein = Ein)
	#def ImpEz(self, Y0, Yh, species, state, CFG, Ein=None):
	#	Yh.scalar_mult(0.0)
	#	self.e_on_dist.Implicit_Ez1d_calc(Y0, Yh, species, state, CFG, Ein = Ein)



class Current_xyz:
	def __init__(self, species, state, CFG):
		from Oshun1d import Axis

		self.jay = state.E_field.clone()
		self.f11 = species.F.harmonic_centric_view()
		self.f10 = species.F.harmonic_centric_view()
		self.integration_temp = np.zeros ( species.F.num_total_cells[0] )

		# constructor copted from 'CurrentSolver' in main code..
		# TODO: unify 'Current_xyz' and 'CurrentSolver'

		# we'll do it micheal's way for now..
		self.p30g = Axis(axis=species.momentum_axis)
		self.delta_p = self.p30g.dx()
		self.small = species.momentum_axis.values[0]

		# calc p^3/gamma for all the gridded momentum values
		for p in xrange(0, self.p30g.num_points):
			p_val = self.p30g.values[p]
			self.p30g.values[p] = (p_val*p_val*p_val)/(math.sqrt(1.0+p_val*p_val))

		self.small *= self.small; self.small *= self.small
		self.small *= species.momentum_axis.values[0]
		self.small *= 0.2; self.small *= 1.0/( species.momentum_axis.values[1] ) 

		# not sure where current will end up living...
		# but in 1d it does not survive past this routine.
		#self.jayx = Field1d(STATE.spatial_axes[0].num_points, 1, state.E_field.num_boundry_cells)
		self.nump =  species.momentum_axis.num_points
		self.species = species
		self.species_charge = self.species.species_config.charge*(-1.0)
		self.final_jayx_factor = (4.0 * np.pi / 3.0) * self.species_charge

	# The input Distribution Function is Yin and it should not be written to.
	# The output is added to Yh. It Yh is initalized it it's proper state before entering.
	def calculateJ(self, Yin, Yh, species, state, CFG, Eout = None):
		Yin.copy_into_harmonic_centric_view(1,self.f10)
		self.f10.mpaxis(self.p30g.values)
		
		jayx = self.jay.components[0]
		jayx[:] = 0.0
		
		# integrate.
		jayx[:] =  np.sum(self.f10.data, axis=1)
		self.integration_temp[:] = self.f10.data[:, 0]
		self.integration_temp += self.f10.data[:, -1]
		self.integration_temp *= 0.5
		jayx[:] -= self.integration_temp[:]
		jayx *= self.delta_p
		
		# Add tiny contribution from 0-p0
		for x in xrange(0, Yin.num_total_cells[0]):
			jayx[x] += self.small * self.f10(1,x)
		
		#self.final_jayx_factor = (4.0 * np.pi / 3.0)*(charge of species * -1.0)
		# Difference form 'CurrentSolver'... Mickhail puts an explixit -1 here.. 
		# TODO: sort out the signs.
		jayx *= self.final_jayx_factor 
		
		# Difference from main 'CurrentSolver'.. this line is commented out.
		#Eout.components[0] += jayx

	# in 1D, this routine is the same as calculateJ
	def Calculate_Jx(self, Yin, Yh, species, state, CFG, Eout = None):
		self.calculateJ( Yin, Yh, species, state, CFG, Eout = Eout)

	def Calculate_Jy(self, Yin, Yh, species, state, CFG, Eout = None):
		# in 2d/3d, this routine is the 'y' portion of calculateJ
		pass

	def Calculate_Jz(self, Yin, Yh, species, state, CFG, Eout = None):
		# in 2d/3d, this routine is the 'z' portion of calculateJ
		pass
	def calculateJComponent(self,Yin, Yh, species, state, CFG, Eout = None, component=1):
		if( component == 1):
			self.Calculate_Jx(Yin, Yh, species, state, CFG, Eout = Eout)
		if( component == 2):
			self.Calculate_Jy(Yin, Yh, species, state, CFG, Eout = Eout)
		if( component == 3):
			self.Calculate_Jz(Yin, Yh, species, state, CFG, Eout = Eout)
	pass

class ImplictECalcaTronicWonderMachine:
	def __init__(self, species, state, CFG):
		# This is really memoery wastefull
		# TODO: refactor this to reuse/loft out the constant resources.. this is going to add-up 

		# The new (i.e. updated) current
		self.jayNew = Current_xyz(species, state, CFG)
		# Current as result of assumeing E=0
		self.jay0 = Current_xyz(species, state, CFG)
		# Current resulting from Ex
		self.jayEx = Current_xyz(species, state, CFG)

		self.Enew = state.E_field.clone()
		self.E0 = state.E_field.clone()
		self.DE = state.E_field.clone()

		self.f00 = species.F.harmonic_centric_view()
		self.f10 = species.F.harmonic_centric_view()
		self.f11 = species.F.harmonic_centric_view()
		self.f20 = species.F.harmonic_centric_view()
		self.f21 = species.F.harmonic_centric_view()
		self.f22 = species.F.harmonic_centric_view()

		# TODO: get this from the axis object.. I trust that more
		cell_size_x =  CFG.cell_size[0]; #cell_size_y =  CFG.cell_size[1]; cell_size_z =  CFG.cell_size[2];	
		self.idx = -1.0/(2.0 * cell_size_x); # self.idy = -1.0/(2.0 * cell_size_y); self.idz = -1.0/(2.0 * cell_size_z);

		self.executionAttempts = 0
		self.Eps = 16.0*sys.float_info.epsilon
		self.LargeEps = np.sqrt(self.Eps)

	# Ein: const
	def findDE(self, Ein):
		self.DE.copyFrom( Ein, component=0); 	#self.DE.copyFrom( Yin.E_field)
		
		#print "DE input "
		#print self.DE.components[0][:]

		self.DE.components[0][:] *= Ein.components[0][:]
		# get DE^2
 		# self.DE.comp(0)[:] += self.DE.comp(1)[:]; self.DE.comp(0)[:] += self.DE.comp(2)[:]

		np.sqrt( self.DE.comp(0), out=self.DE.comp(0))
		#print "DE squared "
		#print self.DE.components[0][:]

		self.DE.comp(0)[:] *= self.LargeEps
		self.DE.comp(0)[:] += self.Eps
		# self.DE.comp(1)[:] = self.DE.comp(0)[:]; self.DE.comp(2)[:] = self.DE.comp(0)[:];
		#print "DE output "
		#print self.DE.components[0][:]

	# Yin: the offical F of the species.. changes made to this are immediatly reflected into the simulation
	#  Yh:	the temperary 'slope' state storage variabee held for each species.
	def advance(self, Yin, Yh, species, state, CFG, fokkerPlankImplictOperator, implicitRK, n, time, dt, Eout = None):
		self.executionAttempts = 0
		Yin.copy_into_harmonic_centric_view(0,self.f00);
		Yin.copy_into_harmonic_centric_view(1,self.f10);#Yin.copy_into_harmonic_centric_view(1,self.f11);
		#Yin.copy_into_harmonic_centric_view(1,self.f20);Yin.copy_into_harmonic_centric_view(1,self.f21);Yin.copy_into_harmonic_centric_view(1,self.f22);
		self.findDE(state.E_field)


		zerosInDeterminate = 1
		while( zerosInDeterminate>0 and self.executionAttempts < 4):
			zerosInDeterminate = 0
			self.executionAttempts += 1

			# Assume the E=0. If this is true then find the effects on distribution function from Collision
			fokkerPlankImplictOperator.advance_1(Yin, species, state, state.config, n, time, dt)
			self.jay0.calculateJ(Yin, Yin, species, state, CFG, Eout = None)
			Yin.copy_from_harmonic_centric_view( 0,self.f00 ); Yin.copy_from_harmonic_centric_view( 1,self.f10 ); #Yin.copy_from_harmonic_centric_view( ,self.f11 );

			# Now set Ex to be the x component of the DE we found earlier.
			state.E_field.comp(0)[:] = self.DE.comp(0)[:]
			implicitRK.advance_IEx(species, state, CFG, Ein=None)
			fokkerPlankImplictOperator.advance_1(Yin, species, state, state.config, n, time, dt)
			self.jayEx.calculateJ(Yin, Yin, species, state, CFG, Eout = None)
			Yin.copy_from_harmonic_centric_view( 0, self.f00 ); Yin.copy_from_harmonic_centric_view(1, self.f10 ); Yin.copy_from_harmonic_centric_view( 2, self.f20 );
			#Yin.copy_from_harmonic_centric_view( self.f11 ); Yin.copy_from_harmonic_centric_view( self.f21 );
			#Yin.copy_from_harmonic_centric_view( self.f22 );
			
			# find the E!
			self.Enew.comp(0)[:] = self.DE.comp(0)[:]; self.Enew.comp(0)[:] *= self.jay0.jay.comp(0);
			self.jay0.jay.comp(0)[:] -= self.jayEx.jay.comp(0)[:]
			self.Enew.comp(0)[:] /= self.jay0.jay.comp(0)[:]
			#self.Enew.comp(1)[:] = 0.0; self.Enew.comp(2)[:] = 0.0; 

		# write out new 'E' to offical state.
		state.E_field.comp(0)[:] = self.Enew.comp(0)[:]
		#state.E_field.comp(1)[:] = self.Enew.comp(1)[:]
		#state.E_field.comp(2)[:] = self.Enew.comp(2)[:]

