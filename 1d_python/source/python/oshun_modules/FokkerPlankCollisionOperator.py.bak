import OshunGlobalState as OshunGlobalState
from OshunGlobalState import PRINT
USE_C_VERSION_GLOBAL = OshunGlobalState.USE_C_VERSION

import numpy as np
import math
import scipy.linalg
from OshunMpi import *
import traceback

if USE_C_VERSION_GLOBAL:
	import numpy.ctypeslib
	from OshunCInterface import *

# Fokker plank has zetho harmonic step that will ensure that 		
class FokkerPlankExplicit:
	
	def __init__(self, species, state, CFG, use_c_version = False):
		
		"""
			tried python cpp, python
			
		"""
		
		# configure the FP code (pure python or CPP or GPU etc..)
		self.explicit__execution_config__python 					= False
		self.explicit__execution_config__cpp    					= True
		
		self.implicit_calc_flm__run_all_in_cpp						= True
		self.explicit_calc_f00__run_all_in_cpp						= True
		
		self.implicit_advance__execution_config__python 		= False
		self.implicit_advance__execution_config__cpp    		= True
		self.test_cpp_vs_python_marix = False
		
		self.implicit_reset_coeff___execution_config__python	= False
		self.implicit_reset_coeff__execution_config__cpp    	= True
		
		if use_c_version:
			self.explicit_calc_f00__run_all_in_cpp = True
			self.implicit_calc_flm__run_all_in_cpp = True
			self.implicit_advance__execution_config__cpp = True
			self.implicit_reset_coeff__execution_config__cpp = True
			self.explicit__execution_config__cpp = True
		else:
			self.explicit_calc_f00__run_all_in_cpp = False
			self.implicit_calc_flm__run_all_in_cpp = False
			
			self.explicit__execution_config__python = True
			self.explicit__execution_config__cpp = False
			
			self.implicit_advance__execution_config__python = True
			self.implicit_advance__execution_config__cpp = False
			self.test_cpp_vs_python_marix = False
			
			self.implicit_reset_coeff___execution_config__python = True
			self.implicit_reset_coeff__execution_config__cpp = False
		
		if mpi_info.rank == 0:
			PRINT("")
			PRINT("Fokker Planck Execution ")
			PRINT("Implict Solve Mode Uses Tridiagonal Approximation: %s" %( str( CFG.if_tridiagonal ) ))
			PRINT("----------------------------------------------------------")
			
		if self.explicit__execution_config__python:
			self.set_run_type__fp_explicit_python()
			PRINT("\t'Explicit' Portion of Collsions: PYTHON")
		else:
			self.set_run_type__fp_explicit_cpp()
			PRINT("\t'Explicit' Portion of Collsions:  CPP")
		
		self.calc_f00 = self.calc_f00__python
		if self.explicit_calc_f00__run_all_in_cpp:
			self.calc_f00 = fokker_plank__calc_f00__cpp
			PRINT("\t'Explicit' Portion of Collsions: TAKEN ALL OVER BY CPP. All other flags invalidated.")
		
		self.calc_flm = self.calc_flm__python; self.advance_1 = self.euler_backward_solver_for_fp_advance_1;self.advance_flm = self.euler_backward_solver_for_fp_advance_lm
		if self.implicit_calc_flm__run_all_in_cpp:
			self.calc_flm = fokker_plank__calc_flm__cpp
			self.advance_1 = fokker_plank__advance_1__cpp
			self.advance_flm = fokker_plank__advance_flm__cpp
			PRINT("\t'Implicit' Portion of Collsions: TAKEN ALL OVER BY CPP. All other flags invalidated.")
		
		if self.implicit_reset_coeff___execution_config__python:
			self.set_run_type__fp_implicit_reset_coeff__python()
			PRINT("\t'Implict Coeff Reset' Portion of Collsions: PYTHON")
		else:
			self.set_run_type__fp_implicit_reset_coeff__cpp()
			PRINT("\t'Implict Coeff Reset' Portion of Collsions: CPP")
		
		if self.implicit_advance__execution_config__python:
			PRINT("\t'Implict Advance' Portion of Collsions: PYTHON")
		else:
			PRINT("\t'Implict Advance' Portion of Collsions': CPP")
		
		if self.implicit_calc_flm__run_all_in_cpp:
			self.printStatus = self.printStatusCpp  
		else:
			self.printStatus = self.printStatusPython
		PRINT("----------------------------------------------------------")
		PRINT("")
		PRINT("")
		
		self.ttt = 0
		# variables for the RK4 part.. 
		self.F0 = np.zeros ( species.momentum_axis.num_points )
		self.Fh = np.zeros ( species.momentum_axis.num_points )
		self.F1 = np.zeros ( species.momentum_axis.num_points )
		#TODO: get rid of this.. it just for debugging
		self.Ftemp = np.zeros ( species.momentum_axis.num_points )
		# these is a 'pointer' (view) to the |p| of the l=0,m=0 mode for the current spatial cell we are looking at
		# 	THIS is only for debugging....
		#self.data_for_current_spatial_cell = species.F.distribution_function_data
		#self.debug_f00_data_for_current_spatial_cell = self.species.F.distribution_function_data[
		
		
		#TODO: talk to Micheal and make this more robust.. this should be sperated out into
		#			A global, specialized module.. I hate how it is spread out in things like Osiris.
		#			Since this uses RK4, then all the code modules must use atleast 4 boundry cells.
		# TODO: read these from the CFG object.
		#self.explicit_solver_order = 3
		self.NB = int( species.F.num_boundry_cells[0][0] )
		self.low_P_transition_cell = 4
		
		self.num_plain_cells = np.copy( species.F.num_plain_cells )
		self.num_total_cells = np.copy( species.F.num_total_cells )
		self.num_bounrdry_cells = np.copy( species.F.num_total_cells )
		self.density_np = species.species_config.density_np
		self.zeta = species.species_config.zeta
		self.is_tridiagonal = CFG.if_tridiagonal
		self.is_implicit1D = CFG.if_implicit1D
		self.num_l = species.num_l
		self.num_m = species.num_m
		
		
		#TODO: talk to Micheal.. centralize the time handling.. porlly not right
		# to be simply reading from input object.
		self.num_subcycling_steps = int(CFG.dt / CFG.small_dt) + 1
		self.subcycle_dt = CFG.dt / float(self.num_subcycling_steps)
		
		# this is the timestep..
		self.h = self.subcycle_dt
		self.species = species
		
		
		# TODO: check boundry cells or else..
		
		# setup inital stuffs for the FP explicit calcuation...
		self.vr		= np.zeros( species.momentum_axis.num_points )
		self.U4 		= np.zeros( species.momentum_axis.num_points )
		self.U4m1	= np.zeros( species.momentum_axis.num_points )
		self.U2		= np.zeros( species.momentum_axis.num_points )
		self.U2m1	= np.zeros( species.momentum_axis.num_points )
		self.U1		= np.zeros( species.momentum_axis.num_points )
		self.U1m1	= np.zeros( species.momentum_axis.num_points )
		
		self.J1		= np.zeros( species.momentum_axis.num_points )
		#self.I2		= np.zeros( species.momentum_axis.num_points )
		self.I4		= np.zeros( species.momentum_axis.num_points )
		
		self.U3		= np.zeros( species.momentum_axis.num_points )
		self.Qn		= np.zeros( species.momentum_axis.num_points )
		self.Pn		= np.zeros( species.momentum_axis.num_points )
		
		self.re = 2.8179402894e-13
		self.kp = np.sqrt( 4.0*np.pi*self.density_np*self.re )
		self.c_kpre = 4.0*np.pi/3.0*self.re*self.kp
		
		#PRINT("FOKKERPLANK: kp=%e, c_kpre=%e" % (self.kp, self.c_kpre))
		
		# shorten some nmaes teroparlily to makes things prettier
		pr_axis =  species.momentum_axis
		vr = self.vr
		self.num_pr_cells = pr_axis.num_points
		
		
		#PRINT(species.momentum_axis.values)
		# initalize all the needed constants..
		for i in xrange(0, pr_axis.num_points):
			self.vr[i] = pr_axis.values[i]
			self.vr[i] = self.vr[i] /( np.sqrt(1 + self.vr[i]*self.vr[i]) ) 
			
		
		for i in xrange(1, pr_axis.num_points):
			self.U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1])     
			self.U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1])         
			self.U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1])         
			self.U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1])         
			self.U1[i]   = 0.5 * vr[i]            * (vr[i]-vr[i-1])         
			self.U1m1[i] = 0.5 * vr[i-1]          * (vr[i]-vr[i-1])
			
		for i in xrange(0, pr_axis.num_points):
			self.U3[i] = pow(vr[i],3)
			
		self.Qn[0] = 1.0 / ((vr[0]*vr[0]*vr[1])/2.0)
		for i in xrange(1, pr_axis.num_points-1):
			self.Qn[i] = 1.0 / (vr[i]*vr[i]*(vr[i+1]-vr[i-1])/2.0)

		for i in xrange(0, pr_axis.num_points-1):
			self.Pn[i] = 1.0 / ((vr[i+1]-vr[i])/2.0*(vr[i+1]+vr[i]))
		
		
		self.G_constant_1 = (vr[0]*vr[0])/(vr[1]*vr[1])
		self.G_constant_vr3 		= np.power( vr, 3) 
		self.G_constant_vr5 		= np.power( vr, 5) 
		self.G_constant_vr5_2 	= np.power( vr, 5) 
		self.G_constant_vr7	 	= np.power( vr, 7)
		
		self.G_constant_vr5 *= 0.2
		self.G_constant_vr5_2 *= (0.2/(vr[1]*vr[1]))
		self.G_constant_vr7	 *= (1.0/(vr[1]*vr[1]*7.0))
		
		
		#------ Stuffs for the implict, higher harmonic calcuations...
		self.kpre = self.re*self.kp
		self.J1m		= np.zeros( species.momentum_axis.num_points )
		self.I0		= np.zeros( species.momentum_axis.num_points )
		self.I2		= np.zeros( species.momentum_axis.num_points )
		self.Scattering_Term = np.zeros( species.momentum_axis.num_points )
		self.df0		= np.zeros( species.momentum_axis.num_points )
		self.ddf0	= np.zeros( species.momentum_axis.num_points )
		
		self.TriI1	= np.zeros( species.momentum_axis.num_points )
		self.TriI2	= np.zeros( species.momentum_axis.num_points )
		self.vr3 = np.zeros( species.momentum_axis.num_points )
		
		for i in xrange(0, pr_axis.num_points):
			self.vr3[i] = (self.vr[i]*self.vr[i]*self.vr[i])
		
		# TODO: can get rid of this..
		#self.fout = np.zeros ( species.momentum_axis.num_points )
		temp = np.zeros( (3, pr_axis.num_points))
		self.tridiagonal_solver_temp = np.matrix( temp )
		
		# matricies!
		# TODO: cehck for any hidden reason.. and make this pcart more tranparent.
		# why did I do the 'temp1=np.zeros(...,..) stuff?
		temp1 = np.zeros( (pr_axis.num_points, pr_axis.num_points) )
		self.Alpha_Tri 		= np.matrix( temp1 )
		temp2 = np.zeros( (pr_axis.num_points, pr_axis.num_points) )
		self.Alpha 				= np.matrix( temp2 )
		
		# these are purly denugging variables....
		# used to redirect the solving routines for comparison to python
		# can be taken out!
		if self.test_cpp_vs_python_marix:
			temp3 = np.zeros( (pr_axis.num_points, pr_axis.num_points) )
			self.Alpha_debug 				= np.matrix( temp3 )
				
		# some misc precomputations.....
		self.IvDnDm1 = np.zeros(  self.num_pr_cells )
		self.IvDnDp1 = np.zeros(  self.num_pr_cells )
		self.Ivsq2Dn = np.zeros(  self.num_pr_cells )
		self.f00_factor = ( (vr[0]*vr[0])/(vr[1]*vr[1]) )
		
		# These are calculated in the implict part of the FP code.. (reset_coeff) and are used
		# across function calls.. so if we using c++, we need to set these.
		self._LOGee = 0.0
		self._ZLOGei = 0.0

		# slush array for solver to put various stats into.
		self.solverStats = np.zeros( 10 )
		
		
		for i in xrange(1,  self.num_pr_cells-1):
			# some common 'constants' used in tthe implicy marriage cerimony
			self.IvDnDm1[i] = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]  -vr[i-1]))       	# ( v*D_n*D_{n-1/2} )^(-1)
			self.IvDnDp1[i] = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1]-vr[i]  ))       	#  ( v*D_n*D_{n+1/2} )^(-1)
			self.Ivsq2Dn[i] = 1.0 / (vr[i] * vr[i]                 * (vr[i+1]-vr[i-1]))       	#  ( v^2 * 2*D_n )^(-1)
		
		if use_c_version:
			# TODO move this.. i want to keep all ctypes things seperate...
			# create  helper objects that allow C/C++ to easily call pyhton functions (they are called 'callbacks')
			#		These will be used for now as a (lazy) compremise where I am using
			#		numpy matrix solvers rather then linking to pur c/c++ ones.. I'll
			#		get to thi later.. it effect perfromace at all.
			#
			self.tridiag_solver_callback = matrix_solver_callback_template( self.solve_tridiagonal_matrix_c_callback )
			self.full_matrix_solver_callback = matrix_solver_callback_template( self.solve_full_matrix_c_callback )
			self.debug_matrix_callback = debug_callback_template( self.debug_matrix_solver_callback )

			init__fokker_plank_implicit__cpp(	self.vr, self.vr3, 
														self.df0, self.ddf0, 
														self.Scattering_Term,
														self.Alpha, self.Alpha_Tri, self.J1m,
														self.TriI1, self.TriI2,
														self.IvDnDm1, self.IvDnDp1, self.Ivsq2Dn,
														self.I0,
														self.kpre, self.zeta, self.f00_factor,
														self.is_tridiagonal, self.is_implicit1D, self.species.index,
														self._LOGee, self._ZLOGei, 
														self.tridiag_solver_callback,
														self.full_matrix_solver_callback,
														self.debug_matrix_callback, self.solverStats)
			
			init__fokker_plank_explicit__cpp( 	self.vr, 
														self.U4, self.U4m1, 
														self.U2, self.U2m1, 
														self.U1, self.U1m1, 
														self.J1, self.I4, self.U3, 
														self.Qn, self.Pn, self.I2, 
														self.G_constant_1, self.G_constant_vr3, self.G_constant_vr5,
														self.G_constant_vr5_2, self.G_constant_vr7,
														self.density_np, self.NB, self.c_kpre, self.num_subcycling_steps, 
														self.species.index)
			# TAKE OUT NOW: This was one leve back
			# hack to make timer work well...
			fokker_plank__implicit__rest_coeff__cpp.implcit_FP_object = self
	
	def set_run_type__fp_explicit_python(self):
		self.lowest_harmonic_fp_slope = self.lowest_harmonic_fp_slope_python
	def set_run_type__fp_explicit_cpp(self):
		self.lowest_harmonic_fp_slope = fokker_plank__explicit__advance__cpp
	def set_run_type__fp_implicit_reset_coeff__python(self):
		self.reset_coeff = self.reset_coeff__python
	def set_run_type__fp_implicit_reset_coeff__cpp(self):
		self.reset_coeff = fokker_plank__implicit__advance__cpp   #self.reset_coeff__python
		
	
	def lowest_harmonic_fp_slope_G(self, n, fin):
		vr = self.vr	
		f00_ = (fin[0] - fin[1]*(vr[0]*vr[0])/(vr[1]*vr[1]))/ (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]))
		i2s = f00_*np.power(vr[n],3)/3.0 + (fin[1]-f00_)*np.power(vr[n],5)/(vr[1]*vr[1])*0.2
		i4s = f00_*np.power(vr[n],5)*0.2 + (fin[1]-f00_)*np.power(vr[n],7)/(vr[1]*vr[1]*7.0)
		
		temp = fin[n]*i4s + (np.power(vr[n],3)*fin[n]-3.0*i2s) * self.J1[n]
		return temp
	
	def lowest_harmonic_fp_slope_G_adam(self, n, fin):
		f00 = (fin[0] - fin[1]*self.G_constant_1) / (1.0 - self.G_constant_1)
       
		#i2s = f00*pow(vr[n],3)/3.0 + (fin[1]-f00)*pow(vr[n],5)/(vr[1]*vr[1])*0.2;
		#i4s = f00*pow(vr[n],5)*0.2 + (fin[1]-f00)*pow(vr[n],7)/(vr[1]*vr[1]*7.0);
		i2s = f00*self.G_constant_vr3[n]/3.0 + (fin[1]-f00)*self.G_constant_vr5_2[n]
		i4s = f00*self.G_constant_vr5[n] + (fin[1]-f00)*self.G_constant_vr7[n]
		
		#return fin[n]*i4s + (pow(vr[n],3)*fin[n]-3.0*i2s) * J1[n];		
		return fin[n]*i4s + (self.G_constant_vr3[n]*fin[n]-3.0*i2s) * self.J1[n]
	
	def LOGee(self, ne, Te):
		if ne > 0.000000001:
			Te /= (3.0*ne)
			Te *= 511000
			ne *= self.density_np 
			Te = np.log(Te)
			ne = np.log(ne)
			lnee = 23.5 - 0.5*ne + 1.25*Te - np.sqrt(0.00001+0.0625*(Te-2.0)*(Te-2.0))
			if lnee > 2.0:
				return lnee
		return 2.0
	def ZLOGei(self, ne, Te):
		if ne > 0.000000001:
			Te /= (3.0*ne)
			Te *= 511000
			ne *= self.density_np 
			Te = np.log(Te)
			ne = np.log(ne)
			lnei = 24.0 - 0.5*ne + Te
			if lnei > 2.0:
				return lnei * self.zeta
		return 2.0 * self.zeta

#     <><><><><><><><><><><><><><><><><><>
#     Tony's Energy Conserving Algorithm 
#     <><><><><><><><><><><><><><><><><><>
	def lowest_harmonic_fp_slope_python(self, fin, fh):
		# just for ease for now...
		vr = self.vr; U4=self.U4; U4m1=self.U4m1
		U2=self.U2; U2m1=self.U2m1; U1=self.U1; U1m1=self.U1m1
		J1=self.J1; I2=self.I2; I4=self.I4
		U3=self.U3; Qn=self.Qn; Pn=self.Pn
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Evaluate the integrals in a standard manner
		I4[0] = 0;
		for n in xrange(1, self.num_pr_cells):
			I4[n]  = U4[n]*fin[n]+U4m1[n]*fin[n-1]
			I4[n] += I4[n-1]
		
		
		I2[0] = 0;
		for n in xrange(1, self.num_pr_cells):
			I2[n]  = U2[n]*fin[n]+U2m1[n]*fin[n-1]; 
			I2[n] += I2[n-1];
		
		
		J1[ self.num_pr_cells-1] = 0;
		for n in xrange(self.num_pr_cells-2, -1, -1):
			J1[n]  = U1[n+1]*fin[n+1]+U1m1[n+1]*fin[n]; 
			J1[n] += J1[n+1];
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Use I2 and I4 to calculate the Coulomb logarithm now, the
		# arrays I2 I4 will be modified later and it will be too late.
		
		Ln_ee = self.LOGee(4.0*np.pi*I2[self.num_pr_cells-1],4.0*np.pi*I4[self.num_pr_cells-1]);
		#PRINT("lnee: %e" % Ln_ee)
		# <><><><><><><><><><><><><><><><><><>
		# Tony's Energy Conserving Algorithm 
		# <><><><><><><><><><><><><><><><><><>
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Evaluate G using the standard integrals
		# J1 needs to be used again later so it is not modified!
		for n in xrange(0, self.num_pr_cells):
			I2[n] *= J1[n];          # J1(n) * I2(n)
			I2[n] *= -3.0;           # -3 * J1(n) * I2(n)
			I4[n] += U3[n] * J1[n];  # I4(n) + u_n^3 * J1(n) 
			I4[n] *= fin[n];         # fn * I4(n) + u_n^3 * fn * J1(n)
			I4[n] += I2[n];          # Gn = fn * I4(n) + u_n^3 * fn * J1(n) - 3 * J1(n) * I2(n)
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Evaluate G assuming a parabolic f(v << vt)
		for n in xrange(0, self.low_P_transition_cell):
			self.I4[n] = self.lowest_harmonic_fp_slope_G(n, fin)
			#PRINT("%e " % self.I4[n],)
		#PRINT("")
		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Find -DG
		fh[0]  = (-1.0)*I4[0];
		for n in xrange(0, self.num_pr_cells-1):
			I4[n] -= I4[n+1];
		
		# Find -DG/(vDv)
		fh[0] *= 2.0/ (vr[0]*vr[0]);
		for n in xrange(0, self.num_pr_cells-1):
			I4[n] *= Pn[n];
		
		# Find DDG/(vDv)
		fh[0] -= I4[0]
		for n in xrange(0, self.num_pr_cells-1):
			I4[n] -= I4[n+1];

		# Find DDG/(v^3*DDv)
		fh[0] *= Qn[0];
		for n in xrange(0, self.num_pr_cells-1):
			fh[n+1] = Qn[n+1]*I4[n];

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Normalize
		fh *=  self.c_kpre *  Ln_ee		
	
	def quick_compare_dist(self, F1, F2):
		for ip in xrange(0,  self.num_pr_cells):
			if F1[ip] != F2[ip]:
				PRINT("They are unequal!")
				exit(-1)
		return
	
	def compare_against_save(self, F2):
		return self.quick_compare_dist(self.Ftemp, F2)
	
	def quick_save_a_dist(self , F1):
		self.Ftemp[:] = F1[:]
		return
	
	# TODO: make role of h,time step more explicit instead of burying it in 
	#		this object's state along with 20 other vars..
	def RK4_for_fokker_plank(self, eval_func, F):
		# F is the offical state used in the rest of the code.
		# e.g. updates to 'F' will be passed on to the rest of the code.
		self.F0[:] = F
		self.F1[:] = F
		
		# results of the evaluation are store inside of Fh, the input
		#		must remain unchanaged.
		
		#self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		#self.compare_against_save( self.F1 )
		
		# f1 = f1 + (h/2)*fh
		self.Fh *= (0.5*self.h)
		self.F1 += (self.Fh)
		# F = F + (h/6)*fh
		self.Fh *= (1.0 / 3.0)
		F += ( self.Fh )
		
		
		#------- Step 2
		#self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		#self.compare_against_save( self.F1 )
		
		self.F1[:] = self.F0[:]
		# f1 = f0 + (h/2)*fh
		self.Fh *= (0.5*self.h)
		self.F1 += (self.Fh)
		# F = F + (h/3)*Fh --> F + h/6*k1+h/2*k2
		self.Fh *= (2.0/3.0)
		F += ( self.Fh )
		#PRINT("result 2")
		#PRINT(F)
		
		#------- Step 3
		#self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		#self.compare_against_save( self.F1 )
		
		#f1 = f0 + h*Fh
		self.Fh *= (self.h)
		self.F0 += (self.Fh)
		# F = F + (h/3)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3
		self.Fh *= (1.0/3.0)
		F += (self.Fh)
		#PRINT("result 3")
		#PRINT(F)
		
		#------- Step 4
		#self.quick_save_a_dist( self.F0)
		eval_func(self.F0, self.Fh)
		#self.( self.F0 )
		
		# F = F + (h/6)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3*+1/6*k4
		self.Fh *= (self.h/6.0)
		F += (self.Fh)
		

#
#---------Higher Order harmonics calculations
#
	# this leaves 'fin' unaltered.. and sets up internal variables.
	#@timeit
	def reset_coeff__python(self, fin, Delta_t):
		# just for ease for now...
		vr = self.vr; U4=self.U4; U4m1=self.U4m1
		U2=self.U2; U2m1=self.U2m1; U1=self.U1; U1m1=self.U1m1
		J1m=self.J1m; I0=self.I0; I2=self.I2; Scattering_Term=self.Scattering_Term; 
		df0=self.df0; ddf0=self.ddf0; TriI1=self.TriI1; TriI2=self.TriI2;
		Alpha_Tri = self.Alpha_Tri; _LOGee = self._LOGee; _ZLOGei = self._ZLOGei
		
		Dt = Delta_t
		self.Dt = Dt
		
		# INTEGRALS
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# Temperature integral I2 = 4*pi / (v^2) * int_0^v f(u)*u^4du
		I2[0] = 0
		for k in xrange(1,  self.num_pr_cells):
			I2[k]  = U4[k]*fin[k]+U4m1[k]*fin[k-1]; 
			I2[k] += I2[k-1];
		I2 *= 4.0 * np.pi 
		I2_temperature = I2[self.num_pr_cells-1]
		for k in xrange(0,  self.num_pr_cells):
			I2[k] /=  vr[k]*vr[k];
		#PRINT("I2")
		#PRINT(I2)
		#PRINT("I2_temperature")
		#PRINT(I2_temperature)
		
		# Density integral I0 = 4*pi*int_0^v f(u)*u^2du 
		I0[0] = 0
		for k in xrange(1,  self.num_pr_cells):
		  I0[k]  = U2[k]*fin[k]+U2m1[k]*fin[k-1]
		  I0[k] += I0[k-1];
		
		I0 *= 4.0 * np.pi
		I0_density = I0[self.num_pr_cells-1];
		#PRINT("I0")
		#PRINT(I0)
		#PRINT("I0_density")
		#PRINT(I0_density)
		
		# Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
		J1m[self.num_pr_cells-1] = 0
		for k in xrange(self.num_pr_cells-2, -1, -1):
		  J1m[k]  = U1[k+1]*fin[k+1]+U1m1[k+1]*fin[k]; 
		  J1m[k] += J1m[k+1];
		
		for k in xrange(0,  self.num_pr_cells):
		  J1m[k] *= 4.0 * np.pi *vr[k];
		#PRINT("J1m")
		#PRINT(J1m)
		
		# COULOMB LOGARITHMS
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		_LOGee  = self.LOGee (I0_density,I2_temperature)
		_ZLOGei = self.ZLOGei(I0_density,I2_temperature)
		
		# BASIC INTEGRALS FOR THE TRIDIAGONAL PART
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
		# Temporary arrays
		for i in xrange(0,  self.num_pr_cells):
			TriI1[i] = I0[i] + (2.0*J1m[i] - I2[i]) / 3.0    	#(-I2 + 2*J_{-1} + 3*I0) / 3
			TriI2[i] = ( I2[i] + J1m[i] ) / 3.0             	#( I2 + J_{-1} ) / 3
			
		# SCATTERING TERM
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		Scattering_Term[:] = TriI1[:]
		Scattering_Term *= _LOGee
		Scattering_Term[0] = 0.0
		Scattering_Term += _ZLOGei * I0_density
		Scattering_Term /= self.vr3
		#np.power(vr, 3.0, Scattering_Term)
		Scattering_Term *=  (self.kpre * Dt)
		#PRINT(_LOGee)
		#PRINT(_ZLOGei)
		#PRINT("Scattering_Term")
		#PRINT(Scattering_Term)
		
		# MAKE TRIDIAGONAL ARRAY
		# TODO: is zerosing necessary?
		Alpha_Tri[:,:] = 0.0
		Alpha_Tri[0,0] = 8.0 * np.pi * fin[0]; 
		
		for i in xrange(1,  self.num_pr_cells-1):
			# July lame opt..
			#IvDnDm1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]  -vr[i-1]))       	# ( v*D_n*D_{n-1/2} )^(-1)
			#IvDnDp1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1]-vr[i]  ))       	#  ( v*D_n*D_{n+1/2} )^(-1)
			#Ivsq2Dn = 1.0 / (vr[i] * vr[i]                 * (vr[i+1]-vr[i-1]))       	#  ( v^2 * 2*D_n )^(-1)
			Alpha_Tri[i, i  ] = 8.0 * np.pi * fin[i] - TriI2[i] * (self.IvDnDm1[i] + self.IvDnDp1[i])
			Alpha_Tri[i, i-1] = TriI2[i] * self.IvDnDm1[i] - TriI1[i] * self.Ivsq2Dn[i]
			Alpha_Tri[i, i+1] = TriI2[i] * self.IvDnDp1[i] + TriI1[i] * self.Ivsq2Dn[i]
		
		factor = (-1.0) * _LOGee * self.kpre * self.Dt
		Alpha_Tri *=  factor        																	# (-1) because the matrix moves to the LHS in the equation
		#PRINT("TriI1")
		#PRINT(TriI1)
		#PRINT("TriI2")
		#PRINT(TriI2)
		
		#PRINT(self.Alpha_Tri)
		
		# CALCULATE DERIVATIVES
		# Evaluate the derivative
		for n in xrange(1,  self.num_pr_cells-1):
			df0[n]  = fin[n+1]-fin[n-1]
			df0[n] /= vr[n+1]-vr[n-1]
		# Evaluate the second derivative
		# -df/dv_{n-1/2}
		for n in xrange(1,  self.num_pr_cells):
			ddf0[n] = (fin[n-1]-fin[n])
			ddf0[n]  /= (vr[n] - vr[n-1])
		# D(df/dv)/Dv
		for n in xrange(1,  self.num_pr_cells-1):
			ddf0[n] -= ddf0[n+1]
			ddf0[n]  /= 0.5*(vr[n+1]-vr[n-1])
		# Calculate zeroth cell
		# july dumb opt
		#f00 = ( fin[0] - ( (vr[0]*vr[0])/(vr[1]*vr[1]) ) *fin[1] ) / (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]));
		f00 =  ( fin[0] - self.f00_factor*fin[1] ) / ( 1.0 -  self.f00_factor )
		ddf0[0] = 2.0 * (fin[1] - f00) / (vr[1]*vr[1]);
		df0[0] = ddf0[0] * vr[0];		
		# Calculate 1/(2v)*(d^2f)/(dv^2),  1/v^2*df/dv
		for n in xrange(0,  self.num_pr_cells-1):
			df0[n]  /= vr[n]*vr[n]
			ddf0[n] /= 2.0  *vr[n]
		
		#PRINT("df0")
		#PRINT(df0)
		#PRINT("ddf0")
		#PRINT(ddf0)
	
	# this update directly into the |p| array 'fin'
	def implicit_advance(self, fin, dt, current_L, current_x_index):
		Alpha = self.Alpha; vr = self.vr; df0 = self.df0; kpre = self.kpre; Scattering_Term = self.Scattering_Term
		ddf0 = self.ddf0; df0 = self.df0; _LOGee = self._LOGee;
		#dt = self.dt
		
		
		if self.test_cpp_vs_python_marix == False and self.implicit_advance__execution_config__cpp == True:
			#Alpha[:] = self.Alpha_Tri
			#if current_L > 1:
			#	Alpfa[0] = 0.0
			
			#if self.is_tridiagonal:
			#	fokker_plank__implicit__advance(fin, current_L, self._LOGee, dt)
			#else:
			#	fokker_plank__implicit__advance(fin, current_L, self._LOGee, dt)
			fokker_plank__implicit__advance__cpp(fin, current_L, _LOGee, dt, current_x_index)
		else:
			if self.test_cpp_vs_python_marix:
				# in this case, write into the Alpha_debug matrix... we write
				# into this to see if we get the same matrix from different codebases...
				Alpha = self.Alpha_debug
			
			Alpha[:] = self.Alpha_Tri
			if current_L > 1:
				Alpha[0,0] = 0.0
			
			if not self.is_tridiagonal:				
				#  CONSTRUCT COEFFICIENTS and then full array
				LL = current_L
				A1 = (LL+1.0)*(LL+2.0) / ((2.0*LL+1.0)*(2.0*LL+3.0))         
				A2 = (-1.0) *(LL-1.0)* LL      / ((2.0*LL+1.0)*(2.0*LL-1.0))
				B1 = (-1.0) *( 0.5 *LL*(LL+1.0) +(LL+1.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0))
				B2 = (       (-0.5)*LL*(LL+1.0) +(LL+2.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0))
				B3 = ( 0.5 *LL*(LL+1.0) +(LL-1.0) ) / ((2.0*LL+1.0)*(2.0*LL-1.0))
				B4 = ( 0.5 *LL*(LL+1.0) - LL      ) / ((2.0*LL+1.0)*(2.0*LL-1.0))				
				
				#TODO: check that ' self.num_pr_cells-1' referes to correct dim of
				#	Alpha.. since its almost always squar i havent seen problems
				# 	but a rectangualr case needs to run. (this is loop over dim1)
				for i in xrange(0,  self.num_pr_cells-1):
					t1 = A1*ddf0[i] + B1*self.df0[i]
					t1 *= (-1.0) * _LOGee * kpre * dt
					t2 = A1*ddf0[i] + B2*df0[i]
					t2 *= (-1.0) * _LOGee * kpre * dt
					t3 = A2*ddf0[i] + B3*df0[i]
					t3 *= (-1.0) * _LOGee * kpre * dt
					t4 = A2*ddf0[i] + B4*df0[i]
					t4 *= (-1.0) * _LOGee * kpre * dt
					#if( i < 5):
					#	PRINT("t's-------")
					#	PRINT(t1)
					#	PRINT(t2)
					#	PRINT(t3)
					#	PRINT(t4)
					
					Alpha[i,0] += t1 * ( 2.0*np.pi * np.power(vr[0]/vr[i],current_L+2)*vr[0]*vr[0]*(vr[1]-vr[0]) )
					Alpha[i,0] += t3 * ( 2.0*np.pi * np.power(vr[0]/vr[i],current_L)  *vr[0]* vr[0]*(vr[1]-vr[0]) )
					
					#if(i==0):
					#	PRINT("alook1")
					#	PRINT(Alpha[i,0])
					
					for j in xrange(1, i):
						Alpha[i,j] += t1 * ( 2.0*np.pi*np.power(vr[j]/vr[i],current_L+2)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
						Alpha[i,j] += t3 * ( 2.0*np.pi*np.power(vr[j]/vr[i],current_L)  *vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
					#if(i==0):
					#	PRINT("alook1.5")
					#	PRINT(Alpha[i,i])
					#	PRINT("t1: %e" % t1)
					#	PRINT("t3: %e" % t3)
					#	PRINT("vr[i]: %e" % vr[i])
					#	PRINT("f1: %e" % ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) ))
					#	PRINT("diff11: %e" % (vr[i]-vr[i-1]))
					
					Alpha[i,i] += t1 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) )
					Alpha[i,i] += t3 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) )
					#if(i==0):
					#	PRINT("alook2")
					#	PRINT(Alpha[i,i])
					
					Alpha[i,i] += t2 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i+1]-vr[i]) )
					Alpha[i,i] += t4 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i+1]-vr[i]) )
					#if(i==0):
					#	PRINT("alook3")
					#	PRINT(Alpha[i,i])
					
					#TODO: check that ' self.num_pr_cells-1' referes to correct dim of
					#	Alpha.. since its almost always squar i havent seen problems
					# 	but a rectangualr case needs to run. (this is loop over  dim2)
					for j in xrange(i+1, self.num_pr_cells-1):
						Alpha[i,j] += t2 * ( 2.0*np.pi*np.power(vr[j]/vr[i],-current_L-1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
						Alpha[i,j] += t4 * ( 2.0*np.pi*np.power(vr[j]/vr[i],-current_L+1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
				pass
		
		
		#PRINT("Alpha 0,0 %e" % (Alpha[0,0], ))
		#PRINT("Alpha")
		#PRINT(Alpha)
		# INCLUDE SCATTERING TERM
		ll1 = float(current_L)
		ll1 *= (-0.5)*(ll1 + 1.0);
		#PRINT(current_L)
		#PRINT(ll1)
		
		# this num_pr_cells == Alpha.dim1
		for i in xrange(0,  self.num_pr_cells):
			temp_calc = 1.0 - ll1 * Scattering_Term[i]
			#PRINT("\t PYTHON diagonal cleanup [%d,%d]=%e  and scattering[i]=%e and additive=%e" % (i, i, Alpha[i,i],  Scattering_Term[i], temp_calc))
			Alpha[i,i] += 1.0 - ll1 * Scattering_Term[i]
		
		"""
		if current_x_index==0 and current_L == 1:
			PRINT("Data f1 x=0, l=1\n------------")
			PRINT("Scatterign term")
			PRINT(Scattering_Term)
			PRINT("Alpha")
			PRINT(Alpha)
		"""
		
		#PRINT("Alpha2")
		#alpha_out = []
		#for ip in xrange(0,  self.num_pr_cells-1):
		#	alpha_out.append(Alpha[ip,ip-1])
		#	alpha_out.append(Alpha[ip,ip])
		#	alpha_out.append(Alpha[ip,ip+1])
		#	#PRINT("%e,%e,%e" % ( Alpha[ip,ip-1], Alpha[ip,ip],Alpha[ip,ip+1] ))
		#PRINT(Alpha[self.num_pr_cells-1, self.num_pr_cells-2])
		#PRINT(Alpha[self.num_pr_cells-1,self.num_pr_cells-1])
		#PRINT("%e,%e,%e" % ( Alpha[self.num_pr_cells-1, self.num_pr_cells-2], Alpha[self.num_pr_cells-1,self.num_pr_cells-1],0.0 ))
		#alpha_out.append(Alpha[self.num_pr_cells-1, self.num_pr_cells-2])
		#alpha_out.append(Alpha[self.num_pr_cells-1, self.num_pr_cells-1])
		#alpha_out.append(0.0);
		#np.save('alpha1.npy', Alpha)
		
		if ( self.is_tridiagonal ):
			#is_tri_diagonal_ok = Thomas_Tridiagonal(Alpha, fin, fout)
			#if is_tri_diagonal_ok == False:
			#	PRINT( "WARNING: Matrix is not diagonally dominant!!!")
			#PRINT("TRI")
			temp = self.solve_tridiagonal_matrix( Alpha, fin )
			#np.save('fout1.npy', temp)
			#np.save('fin1.npy', fin)
			fin[:] = temp[:]
		else:
			# fully solve A
			#self.invert_matrix(Alpha, fin, fout)
			#fin[:] = fout[:]
			#PRINT("FULL!")
			temp = self.solve_matrix( Alpha, fin )
			fin[:] = temp[:]
		#pass
		#PRINT("Fout!!!")
		#PRINT(temp)
		pass
	
	
	# TODO: make this nicer/faster....
	# This is a 'c' callback... which a Python term for a Python function
	#		that is packaged up to be easily called from C.
	# It solves Ax=b.. the vector b is has the values for the |p| projection for a specific harmonic 
	#		at a specific spatial cell. It is copied directly from the current distribution function. 
	# The solution is the 'x' vector. This vector is written into the 'b' vector. Thus we update the
	#		distribution function in-place avoid uncessary copies etc...
	def solve_tridiagonal_matrix_c_callback(self, right_hand_side, current_x_index, current_L, dt):
		# right_hand_side  comes in as a normal C double pointer.. this needed to
		#		converted (not so much converted.. just 'dressed up') as a Numpy array..
		#		luckily this is easy and cheap (althought there are faster ways.. but
		#		worth dealing with since this not a code bottleneck)
		
		# the following line caused a memoery leak upon each call and was replaced.. I need to satify myself why.
		#b =  numpy.ctypeslib.as_array( right_hand_side, shape =  (self.num_pr_cells, ) )
		import ctypes
		b = np.ctypeslib.as_array((ctypes.c_double * self.num_pr_cells).from_address(ctypes.addressof(right_hand_side.contents)))
		try:
			temp = self.solve_tridiagonal_matrix( self.Alpha, b)
			b[:] = temp[:]
		
		except:
			PRINT("Error in Solve Tridiagonal Matrix Callback (spatial cell=%d, L=%d, dt=%0.3e) \n"  % ( current_x_index,current_L,  dt))
			PRINT("\t The Exception thrown was: \n --------------------------------- \n")
			PRINT(traceback.format_exc())
			PRINT("")
			exit(-1)
			return -1

		return 0
		
		"""
		if self.test_cpp_vs_python_marix:
			# redo the whole matrix construction and then solve.. just for comparison..
			#PRINT("\t\t redoing Matrix calc in python....self.is_tridiagonal=%d" % (self.is_tridiagonal , ))
			python_answer = np.copy( b )
			self.implicit_advance(python_answer, dt, current_L, current_x_index)
			
			for x in xrange(0,self.num_pr_cells):
				for y in xrange(0,self.num_pr_cells):
					delta = np.fabs(self.Alpha_debug[x,y] - self.Alpha[x,y])
					if delta > 1e-5:
						PRINT("alpha matrix differs at %d,%d by %e" %( self.Alpha_debug[x,y], self.Alpha[x,y], delta))
						exit(-1)
			
			PRINT("\t\t solving the c++ version")
		
		
		temp = self.solve_tridiagonal_matrix( self.Alpha, b)
		b[:] = temp[:]
		
		
		if self.test_cpp_vs_python_marix:
			if not np.allclose( python_answer, temp ):
				PRINT("Error!!!! x=%d, L=%d dt=%f" % (current_x_index, current_L, dt))
				PRINT("--------------------------------------")
				PRINT("Python gets -- c++ gets")
				for i in xrange(0, temp.shape[0]):
					PRINT("%e     %e" % ( temp[i], python_answer[i] ))
			#else:
			#	PRINT("\t\t compare ok!!!")
		
		
		return 0
		"""
	
	def solve_full_matrix_c_callback(self, right_hand_side, current_x_index, current_L, dt):
		#PRINT("in solve_full_matrix_c_callback  FULL!!!!!!!"\)
		import ctypes
		try:
			# the following line caused a memoery leak upon each call and was replaced.. I need to satify myself why.
			#b =  numpy.ctypeslib.as_array( right_hand_side, shape =  (self.num_pr_cells, ) )
			b = np.ctypeslib.as_array((ctypes.c_double * self.num_pr_cells).from_address(ctypes.addressof(right_hand_side.contents)))
			temp = self.solve_matrix( self.Alpha, b)
			b[:] = temp
		except:
			PRINT("Error in Solve Full Solver Matrix Callback (spatial cell=%d, L=%d, dt=%0.3e) \n"  % ( current_x_index,current_L,  dt))
			PRINT("\t The Exception thrown was: \n --------------------------------- \n")
			PRINT(traceback.format_exc())
			PRINT("")
			exit(-1)
			return -1
		return 0
		
	
	def solve_tridiagonal_matrix(self, matrix, right_hand_side):
		#ud = np.insert(np.diag(A,1), 0, 0) 					# upper diagonal
		#d  = np.diag(A) 											# main diagonal
		#ld = np.insert(np.diag(A,-1), len(d)-1, 0)		# lower diagonal
		
		size = self.tridiagonal_solver_temp.shape[1]
		self.matrix_debug = matrix
		self.tridiagonal_solver_temp[0,1:] = np.diag(matrix, 1)
		self.tridiagonal_solver_temp[1,:]  = np.diag(matrix, 0)
		self.tridiagonal_solver_temp[2,0:(size-1)]  = np.diag(matrix, -1)
		self.tridiagonal_solver_temp[0,0] = 0.0
		self.tridiagonal_solver_temp[2,size-1] = 0.0

		return scipy.linalg.solve_banded( (1,1), self.tridiagonal_solver_temp, right_hand_side)
		
	
	def solve_matrix( self, matrix, right_hand_side):
		self.matrix_debug = matrix
		self.right_hand_side_debug = right_hand_side
		return scipy.linalg.solve( matrix, right_hand_side )
	
	# This is a debugging function called from C when errors are detected when comparing the C based matrix solver
	#	to the python one.
	def debug_matrix_solver_callback(self, arg):
		if self.is_tridiagonal:
			print "-------- upper diagonal (first elment ignored)"
			print self.tridiagonal_solver_temp[0,:]
			print "--------       diagonal"
			print self.tridiagonal_solver_temp[1,:]
			print "-------- lower diagonal (last element ignored)"
			print self.tridiagonal_solver_temp[2,:]
			print "-------------------"
			print "The first,2nd matrix rows"
			print self.matrix_debug[0,:]
			print self.matrix_debug[1,:]
			print "The second to last,last matrix rows"
			print self.matrix_debug[-2,:]
			print self.matrix_debug[-1,:]

		else:
			print "The first,2nd matrix rows"
			print self.matrix_debug[0,:]
			print self.matrix_debug[1,:]
			print "The second to last,last matrix rows"
			print self.matrix_debug[-2,:]
			print self.matrix_debug[-1,:]
			print "right hand side"
			print self.right_hand_side_debug
			pass
		return 0

	# result put inplace into Yin
	def f1_loop(self, Yin, dt, max_m):
		# loop over spaitial cells....
		#f00 = self.F0;  fc = self.F1;
		for x in xrange(0, self.num_plain_cells[0]):
			# reset the integrals an ceofficients..
			# reset_coeff does not alter the input.
			#PRINT("Starting cell x=%d (fl_loop) dt=%e" % (x, dt))
			temp_1 = Yin.cells[x + self.NB].harmonics[0].momentum_projection[1]
			self.reset_coeff( Yin.cells[x + self.NB].harmonics[0].momentum_projection, dt)
			temp_12 = Yin.cells[x + self.NB].harmonics[0].momentum_projection[1]
			
			temp_2 = Yin.get_harmonic_at_location(x + self.NB, 1).momentum_projection[1]
			
			# loop over harmonics
			#for m in xrange(0, max_m):
			#	self.implicit_advance(Yin.get_harmonic_at_location(x + self.NB, 1, m).momentum_projection, 1)
			self.implicit_advance(Yin.get_harmonic_at_location(x + self.NB, 1).momentum_projection, dt, 1, x)
			#PRINT("\t\tl=1 (fl_loop)")
			temp_3 = Yin.get_harmonic_at_location(x + self.NB, 1).momentum_projection[1]
			
			#PRINT("For PYTHON annoying 'x=%d+%d,L=1' %e, %e, %e, %e" %(x, self.NB, temp_1, temp_12, temp_2, temp_3))
	
	def flm_loop(self, Yin, dt, max_m):
		
		for x in xrange(0, self.num_plain_cells[0]):
			#PRINT("Starting cell x=%d (flm_loop) dt=%e" % (x, dt))
			# reset_coeff does not alter the input.
			# TODO: Assume you did this when you called f1_loop ???? Yes or no?
			#F_start_ix = (x + self.NB)*Yin.num_l*self.num_pr_cells
			#F_end_ix = (x + self.NB)*Yin.num_l*self.num_pr_cells + self.num_pr_cells
			#t2 = Yin.distribution_function_data[ F_start_ix:F_end_ix ]
			
			"""
			# debugging assertion test that the addressing of the dist function is correct..
			t1 = Yin.cells[x + self.NB].harmonics[0].momentum_projection;
			t2 = Yin.distribution_function_data[ x + self.NB, 0, :] 
			for tme in xrange(0, self.num_pr_cells ):
				diff = np.fabs( t1[tme] - t2[tme] )
				if diff > 1e-8:
					PRINT("unequal!!!!!")
					PRINT("buffer 1")
					PRINT(t1)
					PRINT("buffer 1")
					PRINT(t2)
			"""
			
			self.reset_coeff( Yin.cells[x + self.NB].harmonics[0].momentum_projection, dt)
			
			# loop over harmonics for this spatial point.
			for l in xrange(2, self.num_l):
				temp_1 = Yin.get_harmonic_at_location(x + self.NB, l).momentum_projection[0]
				
				self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l).momentum_projection, dt, l, x)
				temp_2 = Yin.get_harmonic_at_location(x + self.NB, l).momentum_projection[0]
				
				#PRINT("In python x=%d, l=%d   %e,%e" % (x, l, temp_1, temp_2);)
				
				# we don't m>0 yet, so comment the following block out...
				"""
				if max_m < l:
					for m in xrange(0, max_m):
						self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l, m).momentum_projection, l)					
				else:
					for m in xrange(0, l):
						self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l, m).momentum_projection, l)
				"""
		pass
	
	
	# this does the integration for f10, f11
	# num_m is not necessary for now as it's always 1 in 1D
	# TODO: compact and refactor... dont need 3 levels of functions here...
	def euler_backward_solver_for_fp_advance_1(self, Yin, dt, num_m):
			if self.is_implicit1D:
				# if 1D, do calcualtion just for f10 (f11 not needed)
				self.f1_loop(Yin, dt, 1)
			else:
				# Do calcualtion just for both f10 and f11
				self.f1_loop(Yin, dt, 2)
	
	# num_m is not necessary for now as it's always 1 in 1D
	def euler_backward_solver_for_fp_advance_lm(self, Yin, dt, num_m):
			if self.is_implicit1D:
				self.flm_loop(Yin, dt, 1)
			else:
				self.flm_loop(Yin, dt, num_m)
	
	# main entry point for higher harmonic collsions.. results will be written
	#	directly into the input distribution function.
	#@timeit
	def calc_flm__python(self, Yin, species, state, CFG, n_step, time, dt):
		
		# this does the integration for f10, f11
		self.euler_backward_solver_for_fp_advance_1(Yin, dt, self.num_m)
		
		# this does the remaing higher harmonics..
		self.euler_backward_solver_for_fp_advance_lm(Yin, dt, self.num_m)
		#exit(-1)
	
	
	# main entry point for f00 collsions.. results will be written
	#	directly into the input distribution function.
	def calc_f00__python(self, Yin, species, state, CFG, n_step, time, dt):
		eval_func = self.lowest_harmonic_fp_slope
		#PRINT("self.NB %d with %d subcycle steps (each with dt=%e)" %  (self.NB, self.num_subcycling_steps, self.h))
		
		# loop over the useful (non-boundry) cells.. (we don;t need to calculate th FP operator
		#	on the boundry cells.
		for x in xrange(0, self.num_plain_cells[0]):
			# grab |p| projections from the  L=0, M=0 harmonic at this spatial location
			momentum_projection = Yin.cells[x+self.NB].harmonics[0].momentum_projection
			# update these values using the FP-operator. This is repeated many times
			#	as a "subcycling" operation.. that is, many small FP update steps are performed
			#	for 1 step of main simulation.
			for n in xrange(0, self.num_subcycling_steps):
				# results of Rk4 will be written directly into Yin's L=0, M=0 momentum projection.
				self.RK4_for_fokker_plank( eval_func, momentum_projection )
		pass
	def printStatusPython(self):
		pass
	def printStatusCpp(self):
		if self.is_tridiagonal == False:
			if mpi_info.is_I_root():
				print("\t Max itterations for implicit collsion calculation this timestep: %d (max ever: %d) " % (int(self.solverStats[0]),int(self.solverStats[0])))

	