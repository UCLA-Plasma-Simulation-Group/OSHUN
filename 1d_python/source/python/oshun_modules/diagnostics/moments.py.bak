#TODO:		THERE ARE TWO VAIRABLES f00_ f00__.. I did this because they are both commonly used in a dumb initialization routine with the moments
#				But I don't know what mood I was in (probably angry) `cause these are terrible variable names

# TODO: 	'vz' is not used in the Q_x calculation

# TODO: 	Vy relies on the L=1, M=1 harmonic... Vz relies on L=1, M=-1
# 				check these.. might be something hiding... as it stands it is unchecked

# TODO:		Do the relativistic momentum moments

# TODO: get to the bottom of comments in the various integraters that look like:
			#I am redoing the calculation here because the pr(0)=0.. but it might not numerically..
			#	This is because this portion of the calcuation is operating on a higher (l=2) harmonic.
			#		To make the intergraton of F consistent, we assume that there are no higher harmonics assocaited
			#	with the lowest |p| cell.
			#   so this is just ot be safe.
			#  TODO: include case for this so do not have to redo.
			#
			#	I have an idea.. its that F20, F21 anf F22, becasue of the harmonic filterting (in fact all harmonics with l>1),
			#		have zero as the weight for the p(pmin) cell.
			#		And look again becasue I though L=1 would also have same thing.. but no
			#		And look at exact defintion of P axis.. where does it start.. where are values defined.. i forgot!
			#

#TODO:	put in the Q_y diagnpstic

#pointer, read_only_flag = self.axis_temp.__array_interface__['data']; print "before pointer"; print pointer
import numpy as np
#from pylab import *
#print dir()
#TODO: don't do calculations of guard cells,,, it's a waste!

class MomentQuantity:
	def __init__(self, name):
		self.name = name
		self.data = None
		self.last_update_timestep = -1
		
class MomentCalculationsState:
	def __init__(self, species, state):
		self.MomentCalculations = MomentCalculations(species, state)
		self.moments = {}
		
		# 0th smoment		
		self.moments['density'] 			= MomentQuantity('density')
		# 1st velocity moments
		self.moments['vx'] 				= MomentQuantity('vx')
		self.moments['vy'] 				= MomentQuantity('vy')
		self.moments['vz'] 				= MomentQuantity('vz')
		# 2nd moments diagonal (isotropic)
		self.moments['vxvx'] 				= MomentQuantity('vxvx')
		self.moments['vyvy'] 				= MomentQuantity('vyvy')
		self.moments['vzvz'] 				= MomentQuantity('vzvz')
		self.moments['moment1_v_v_'] 	= MomentQuantity('moment1_v_v_')
		self.moments['moment2_v_v_'] 	= MomentQuantity('moment2_v_v_')
		# 2nd moments tensor
		self.moments['vxvy'] 				= MomentQuantity('vxvy')
		self.moments['vxvz'] 				= MomentQuantity('vxvz')
		self.moments['vyvz'] 				= MomentQuantity('vyvz')
		# tempaerature-ish and collsions
		self.moments['v2'] 				= MomentQuantity('v2')
		self.moments['temp_ev'] 			= MomentQuantity('temp_ev')
		self.moments['pressure'] 		= MomentQuantity('pressure')
		self.moments['ND'] 				= MomentQuantity('ND')
		self.moments['nu'] 				= MomentQuantity('nu')
		# heat conduction and flux
		self.moments['qx'] 				= MomentQuantity('qx')
		self.moments['kx'] 				= MomentQuantity('kx')
		
		self.current_time = 0.0
		self.current_step = 0
		self.species = species
		self.state = state
	
		
	def is_dirty_var(self, name):
		try:
			var = self.moments[name]
		except:
			raise Exception("The moment '%s' is not valid." % (str(name),) )
		if var.data == None:
			#print "\t Moment '%s' is NULL dirty" % (name,)
			return True
		if var.last_update_timestep == self.current_step:
			#print "\t Moment '%s' is clean" % (name,)
			return False
		#print "\t Moment '%s' is dirty" % (name,)
		return True
	
	def data(self, name ):
		var = self.checkout_var(name)
		return var.data
	
	def checkout_var(self, name ):
		try:
			var = self.moments[name]
		except:
			raise Exception("The moment '%s' is not valid." % (str(name),) )
		
		if var.data == None:
			var.data = self.create_var(var)
			var.last_update_timestep = self.current_step
		return var
	
	def mark_var_updated(self, name):
		try:
			var = self.moments[name]
		except:
			raise Exception("The moment '%s' is not valid." % (str(name),) )
		var.last_update_timestep = self.current_step
	
	def create_var(self, var):
		return np.zeros ( (self.MomentCalculations.num_data_cells_x,) )
		 
	def notify_simulation_time(self, n, time):
		self.current_time = time
		self.current_step = n
	
	def get_density(self, force_update = False):
		if self.is_dirty_var('density') or force_update:
			self.MomentCalculations.calculate_density1(self.species, self.state, density=self.data('density') )
			self.mark_var_updated('density')
		return self.data('density')

	def get_vx(self, force_update = False):
		if self.is_dirty_var('vx'):
			self.MomentCalculations.calculate_vx1(self.species, self.state, density = self.get_density(),
																vx=self.data('vx') )
			self.mark_var_updated('vx')
		return self.data('vx')
	def get_vy(self, force_update = False):
		if self.is_dirty_var('vy'):
			self.MomentCalculations.calculate_vy1(self.species, self.state, density = self.get_density(),
																vy=self.data('vy') )
			self.mark_var_updated('vy')
		return self.data('vy')	
	def get_vz(self, force_update = False):
		if self.is_dirty_var('vz'):
			self.MomentCalculations.calculate_vz1(self.species, self.state, density = self.get_density(),
																vz=self.data('vz') )
			self.mark_var_updated('vz')
		return self.data('vz')
	
	def get_vxvx(self, force_update = False):
		if self.is_dirty_var('vxvx'):
			self.MomentCalculations.calculate_vxvx(self.species, self.state, density = self.get_density(),
																vxvx=self.data('vxvx'),  moment1=self.data('moment1_v_v_'))
			self.mark_var_updated('vxvx')
			self.mark_var_updated('moment1_v_v_')
		return self.data('vxvx')
	
	def get_moment1_v_v_(self, force_update = False):
		if self.is_dirty_var('moment1_v_v_'):
			self.MomentCalculations.calculate_vxvx(self.species, self.state, density = self.get_density(),
																vxvx=self.data('vxvx'),  moment1=self.data('moment1_v_v_'))
			self.mark_var_updated('vxvx')
			self.mark_var_updated('moment1_v_v_')
		return self.data('moment1_v_v_')
	
	def get_vyvy(self, force_update = False):
		if self.is_dirty_var('vyvy'):
			self.MomentCalculations.calculate_vyvy(self.species, self.state, density = self.get_density(),
																vyvy=self.data('vyvy'),  moment2=self.data('moment2_v_v_'), moment1=self.get_moment1_v_v_() )
			self.mark_var_updated('vyvy')
			self.mark_var_updated('moment2_v_v_')
		return self.data('vyvy')
	
	def get_moment2_v_v_(self, force_update = False):
		if self.is_dirty_var('moment2_v_v_'):
			self.MomentCalculations.calculate_vyvy(self.species, self.state, density = self.get_density(),
																vyvy=self.data('vyvy'),  moment2=self.data('moment2_v_v_'), moment1=self.get_moment1_v_v_() )
			self.mark_var_updated('vyvy')
			self.mark_var_updated('moment2_v_v_')
		return self.data('moment2_v_v_')
	
	def get_vzvz(self, force_update = False):
		if self.is_dirty_var('vzvz'):
			self.MomentCalculations.calculate_vzvz(self.species, self.state, density = self.get_density(),
																vzvz=self.data('vzvz'),  moment2=self.get_moment2_v_v_(), moment1=self.get_moment1_v_v_() )
			self.mark_var_updated('vzvz')
		return self.data('vzvz')
	
	def get_vxvy(self, force_update = False):
		if self.is_dirty_var('vxvy'):
			self.MomentCalculations.calculate_vxvy(self.species, self.state, density = self.get_density(),
																vxvy=self.data('vxvy'))
			self.mark_var_updated('vxvy')
		return self.data('vxvy')
	def get_vxvz(self, force_update = False):
		if self.is_dirty_var('vxvz'):
			self.MomentCalculations.calculate_vxvz(self.species, self.state, density = self.get_density(),
																vxvz=self.data('vxvz'))
			self.mark_var_updated('vxvz')
		return self.data('vxvz')
	def get_vxvz(self, force_update = False):
		if self.is_dirty_var('vyvz'):
			self.MomentCalculations.calculate_vyvz(self.species, self.state, density = self.get_density(),
																vyvz=self.data('vyvz'))
			self.mark_var_updated('vyvz')
		return self.data('vyvz')
	
	def get_v2(self, force_update = False):
		if self.is_dirty_var('v2'):
			self.MomentCalculations.calculate_v2(self.species, self.state, density = self.get_density(),
																v2=self.data('v2'), vx=self.get_vx(), vy=self.get_vy(), vz=self.get_vz(),
																vxvx=self.get_vxvx(), vyvy=self.get_vyvy(), vzvz=self.get_vzvz())
			self.mark_var_updated('v2')
		return self.data('v2')
	
	def get_temp_ev(self, force_update = False):
		if self.is_dirty_var('temp_ev'):
			self.MomentCalculations.calculate_temp_ev(self.species, self.state, density=self.get_density(),
																	temp_ev=self.data('temp_ev'), v2=self.get_v2() )
			
			self.mark_var_updated('temp_ev')
		return self.data('temp_ev')
	
	def get_pressure(self, force_update = False):
		if self.is_dirty_var('pressure'):
			self.MomentCalculations.calculate_presure(self.species, self.state, density=self.get_density(),
																	pressure=self.data('pressure'), temp_ev=self.get_temp_ev() )
			
			self.mark_var_updated('pressure')
		return self.data('pressure')
	
	def get_ND(self, force_update = False):
		if self.is_dirty_var('ND'):
			self.MomentCalculations.calculate_ND(self.species, self.state, density=self.get_density(),
																	ND=self.data('ND'), temp_ev=self.get_temp_ev() )
			
			self.mark_var_updated('ND')
		return self.data('ND')
	
	def get_qx(self, force_update = False):
		if self.is_dirty_var('qx'):
			self.MomentCalculations.calculate_qx(self.species, self.state, 
																self.get_vx(), self.get_vy(),  self.get_vz(), self.get_vxvy(), 
																self.get_vxvx(), vyvy=self.get_vyvy(), vzvz=self.get_vzvz(),
																qx=self.data('qx'), density=self.get_density() )
			self.mark_var_updated('qx')
		return self.data('qx')
	
	def get_nu(self, force_update = False):
		if self.is_dirty_var('nu'):
			self.MomentCalculations.calculate_collision_frequency(self.species, self.state, 
																temp_ev=self.get_temp_ev(), v2=self.get_v2(),
																nu=self.data('nu'), ND_=self.get_ND(), density=self.get_density() )
			self.mark_var_updated('nu')
		return self.data('nu')
	
	def get_heat_conductivity_x(self, force_update = False):
		if self.is_dirty_var('kx'):
			self.MomentCalculations.heat_conductivity_x(self.species, self.state, self.get_v2(), self.get_nu(),
																		self.get_qx(), density=self.get_density(), kx=self.data('kx'))
			
			self.mark_var_updated('kx')
		return self.data('kx')
		
		
class MomentCalculations:
	def __init__(self, species, state):
		from ..Oshun1d import Axis,DistributionFunction, SHarmonic, DistributionFunctionSpatialCell, Species, HarmonicCentricView, Config, SpeciesConfig, DaState
		
		F = species.F
		self.momentum_axis = F.momentum_axis
		self.pr = F.momentum_axis.values
		self.pr_num_points = F.momentum_axis.num_points
		self.pr__dp = F.momentum_axis.dp()
		self.state = state
		
		self.p0p1_sq = self.pr[0]*self.pr[0] / (self.pr[1]*self.pr[1])
		self.inv_mp0p1_sq = 1.0/(1.0-self.p0p1_sq)
		self.f00 = 0.0
		self.f00_ = np.zeros( F.num_total_cells[0] )
		self.f00__ = np.zeros( F.num_total_cells[0] )
		self.moment = np.zeros( F.num_total_cells[0] )
		self.moment_temp = np.zeros( F.num_total_cells[0] )
	
		self.gamma_axis = Axis(axis=F.momentum_axis)
		self.inverse_gamma_axis = Axis(axis=F.momentum_axis)
		self.axis_temp = np.array ( np.zeros ( self.gamma_axis.num_points ) )
		for p in xrange(0, self.gamma_axis.num_points):
			self.gamma_axis.values[p] = np.sqrt( (1.0 + self.pr[p]*self.pr[p] ) )
			self.inverse_gamma_axis.values[p] = 1.0/self.gamma_axis.values[p]
		
		self.num_data_cells_x = F.num_plain_cells[0] # - F.num_boundry_cells[0][0] - F.num_boundry_cells[0][1]
		#self.num_data_cells_y = F.num_cells[1] - F.num_boundry_cells[1][0] - F.num_boundry_cells[1][1]
		#self.num_data_cells_z = F.num_cells[2] - F.num_boundry_cells[2][0] - F.num_boundry_cells[2][1]
		
		self.num_data_cells_x 		= F.num_total_cells[0]
		self.num_plain_cells 		= F.num_plain_cells
		self.num_boundry_cells 		= F.num_boundry_cells
		
		self.moment_temp__no_guard_cells = np.zeros ( (self.num_data_cells_x) )
		self.F0  = F.harmonic_centric_view()
		self.F10 = F.harmonic_centric_view()
		self.F11 = F.harmonic_centric_view()
		self.F20 = F.harmonic_centric_view()
		self.F21 = F.harmonic_centric_view()
		self.F22 = F.harmonic_centric_view()
		
		#readability
		self.num_bdry_cells_left__x = F.num_boundry_cells[0][0]
		self.num_bdry_cells_left__y = F.num_boundry_cells[0][1]
		
		self.FLT_MIN = 1.175494351e-38
		self.FLT_MIN4 = 4.0*1.175494351e-38
		self.DBL_MIN = 2.2250738585072014e-308
		# random factor for various momentus that are low hanging fruit.

		# stuff for MPI
		pass
	
	# This calculates the density at each spatial cell (It does data cells as well as all the guard cells..
	#	(It does guard cells and data cells.. a small optimazation is to not calculate ensity for the guard cells)
	def calculate_density1(self, spec, F, density=None):
		pass
		from ..Oshun1d import DistributionFunction, SHarmonic, DistributionFunctionSpatialCell, Species, HarmonicCentricView, Config, SpeciesConfig, DaState
		F = spec.F
		if density == None:
			density = np.zeros( (self.num_data_cells_x) )
		

		F.copy_into_harmonic_centric_view(0,self.F0)
		
		# generate momentum^2 for all momentum cells..
		#	and apply this (i.e. multiply in 1d strips) to the 0th Harmonic's data 
		#	(via .mpaxis(xxx) which is equivalent to: for x in num_spatial_points: harmonic_centric_data[x, :] *= self.axis_temp
		self.axis_temp[:] = self.pr[:]
		self.axis_temp *= self.pr[:]
		self.F0.mpaxis( self.axis_temp )
		
		# Do average of first/last |p| cells (At all points in space simultaneously)
		#	(we will call first cell |p|(0) and last cell |p|(-1) )
		self.moment[:] = self.F0.data[:, 0]
		self.moment += self.F0.data[:, -1]
		self.moment_temp[:] = self.moment[:]
		self.moment *= 0.5
		
		# sum over |p| to needed value at each spatial point... but the following sum will also include
		#	the 0th and last cell of |p|.. which we just already accounted for by adding in 1/2*(|p|(0)+|p|(-1))
		#	so after the sum over all points, the first/last points will be double counted
		#	To be exact, the contribution of the first/last cells will be 3/2*(|p|(0)+|p|(-1)) after the sum so we will
		#	subtract out (|p|(0)+|p|(-1)) so that we are left with the desired 1/2*(|p|(0)+|p|(-1)) contribution for the first/last cells
		#	This is done this way becasue its much faster in numpy to use unrestructed sums.
		self.moment += np.sum( self.F0.data, axis = 1)
		self.moment -= self.moment_temp
		
		# Take care of the edge interpolation ( we decided this is overkill and will prolly be removed )
		#	See the latexx mathematical documentation for deatials about this.. It is usually a tiny effect.
		self.f00_ [:] = self.F0.data[:,  1]
		self.f00_ *= self.p0p1_sq
		self.f00__[:] = self.F0.data[:,  0]
		self.f00__ -= self.f00_ 
		factor = ( self.inv_mp0p1_sq*(1.0/3.0-1.0/5.0*self.p0p1_sq) )
		self.f00__ *= factor
		
		self.f00_ /= 5.0
		self.f00__ += self.f00_
		
		factor = ( self.pr[0]*(self.pr[0]/self.momentum_axis.dx())*self.pr[0] )
		self.f00__ *= factor
		
		self.moment += self.f00__

		# now compute the normalzation factor for the intergral...
		density_norm_factor = 4.0 * np.pi * self.momentum_axis.dx()
		self.moment *= density_norm_factor

		density[:] = self.moment[0:self.num_data_cells_x]
		
		#TODO: Mikhail has a strange lower bound check here that I didn't put in. See if it matters

		return density
		
	def do_moment(self, spec, F, functs, target_harmonics, result=None):
		from ..Oshun1d import DistributionFunction, SHarmonic, DistributionFunctionSpatialCell, Species, HarmonicCentricView, Config, SpeciesConfig, DaState
		if result == None:
			result = np.zeros( (self.num_data_cells_x) )
		
		# Don't do moments if 'target_harmonics' is invalid or if 'target_harmonics' specifies that the moment
		#	being calculated depends on a higher order harmonics than exists in the system.
		#	In these cases, return zeros as the result (since all higher then specfied haronics have magnitude zero.. ie.e we have truncated the sphereical harmonic series.
		if( target_harmonics[0][0] >= spec.F.num_l ):
			result = 0.0
			return result
		if( target_harmonics[0][1] >= spec.F.num_m ): # this was > before, but was changes becsue those terms were giving zero anyway.
			result = 0.0
			return result
		
		spec.F.copy_into_harmonic_centric_view(target_harmonics[0][0], target_harmonics[0][2])
		target_harmonic_data = target_harmonics[0][2]

		functs['p_axis_prep']()
		target_harmonic_data.mpaxis( self.axis_temp )
		
		
		self.moment[:] = target_harmonic_data.data[:, 0]
		self.moment += target_harmonic_data.data[:, -1]
		self.moment_temp[:] = self.moment[:]
		
		if 'prep_first_last_p_cell' in functs.keys():
			functs['prep_first_last_p_cell']()
		# sum over |p| to needed value at spatial point... but we have also included 
		#	the 0,1 cells of |p|.. they are double counted.. so we subtract them out.
		self.moment += np.sum( target_harmonic_data.data, axis = 1)
		self.moment -= self.moment_temp
		
		if 'f00_factor' in functs.keys():
			factor = functs['f00_factor']()
			self.f00_[:] = target_harmonic_data.data[:,  1]
			self.f00_ *= factor
			self.moment += self.f00_
		if 'f00_need_pr1_and_pr0' in functs.keys():
			self.f00_[:] = target_harmonic_data.data[:,  1]
			self.f00__[:] = target_harmonic_data.data[:, 0]
			functs['f00_need_pr1_and_pr0'](self.f00_, self.f00__)
			self.moment += self.f00_
			
		result[:] = self.moment[0:self.num_data_cells_x]	
		if 'overall_factor' in functs.keys():
			factor = functs['overall_factor']()
			# NEW I THINK FIXES BUG
			result *= factor

		if 'final_norm' in functs.keys():
			functs['final_norm'](result)

		return result
		
		#if 'functs['compenstate_first_last_p_cell' in functs.keys():
		#	functs['compenstate_first_last_p_cell']()
		
	def calculate_vx1(self, spec, F, vx=None, density=None, calc_density=False):
		if vx == None:
			vx = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		functs = {}
		
		target_harmonics = [(1,0, self.F10)]
		def p_axis_prep():
			np.power(self.pr, 3, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		
		def prep_first_last_p_cell():
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		
		#def compenstate_first_last_p_cell(self):
		#functs['compenstate_first_last_p_cell'] = prep_first_last_p_cell
		def f00_factor():
			factor = self.pr[0]*self.pr[0]
			factor *= factor
			factor *= self.pr[0]
			factor /= ( self.momentum_axis.dx()*self.pr[1]*5.0 )
			return factor
		functs['f00_factor'] = f00_factor
		
		def overall_factor():
			return (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		def final_norm(result):
			result /= density
		functs['final_norm'] = final_norm
		
		vx = self.do_moment(spec, F, functs, target_harmonics, result=vx)
		return vx
	
	
	def calculate_vy1(self, spec, F, vy=None, density=None, calc_density=False):
		if vy == None:
			vy = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		functs = {}
		
		# Vy relies on the L=1, M=1 harmonic...
		# TODO: check this.. might be somehting hdding... as it stands it is unchecked
		target_harmonics = [(1,1, self.F10)]
		def p_axis_prep():
			np.power(self.pr, 3, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		
		def prep_first_last_p_cell():
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		
		
		def f00_factor():
			factor = self.pr[0]*self.pr[0]
			factor *= factor
			factor *= self.pr[0]
			factor /= ( self.momentum_axis.dx()*self.pr[1]*5.0 )
			return factor
		functs['f00_factor'] = f00_factor
		
		def overall_factor():
			return (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		def final_norm(result):
			result /= density
		functs['final_norm'] = final_norm
		
		vy = self.do_moment(spec, F, functs, target_harmonics, result=vy)
		return vy

	def calculate_vz1(self, spec, F, vz=None, density=None, calc_density=False):
		if vz == None:
			vz = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		functs = {}
		
		# Vz relies on the L=1, M=-1 harmonic...
		# TODO: check this.. might be somehting hdding... as it stands it is unchecked
		target_harmonics = [(1,1, self.F10)]
		def p_axis_prep():
			np.power(self.pr, 3, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		
		def prep_first_last_p_cell():
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		
		
		def f00_factor():
			factor = self.pr[0]*self.pr[0]
			factor *= factor
			factor *= self.pr[0]
			factor /= ( self.momentum_axis.dx()*self.pr[1]*5.0 )
			return factor
		functs['f00_factor'] = f00_factor
		
		def overall_factor():
			return (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		def final_norm(result):
			result /= density
		functs['final_norm'] = final_norm
		
		vz = self.do_moment(spec, F, functs, target_harmonics, result=vz)
		return vz
		
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
	
	def ghetto_conductivity(self, spec, temp_ev, density, q_x, kx=None):
		if kx == None:
			kx = np.empty( (self.num_data_cells_x) )
		
		gamma0 = 3.0
		ion_z = spec.species_config.zeta
		denisty_np = spec.species_config.density_np
		kb__J_K = 1.380648e-23
		kb__ev_K =  8.617332478e-5 
		vt__cm_sec = 4.19e7 * np.sqrt( t_ev )
		
		
		for ix in xrange(0, self.num_plain_cells):
			
			density 			= density[num_bdry_cells_left__x+ ix]
			temp_eq__ev 	= temp_ev[num_bdry_cells_left__x + ix]
			Q_x 				= q_x[num_bdry_cells_left__x + ix]
			dx 				= self.state.spatial_axes[0].dx()
			
			density__1_cm3 				= density*denisty_np
			wp__1_sec 						= u.wpe(density__1_cm3)
			skindepth__cm 					= 3e10 / wp__1_sec
			
			x__cm 							= x*skindepth__cm
			dx__cm							= dx*skindepth__cm
			ee_collision_freq 			= u.electron_collision_freq(temp_eq__ev, density__1_cm3)
			electron_mfp_cm 				= u.electron_mfp__cm(temp_eq__ev, density__1_cm3, ion_z)
			ee_collsion_time 				= u.ee_collsion_time(temp_eq__ev, density__1_cm3)
			factor_convert_wp_to_ee 	= ee_collision_freq / wp__1_sec
			
			# delta T (central difference)
			deltaT__ev_cm 					= (temp_eq__ev[ix+1] - temp_eq__ev[ix-1])/(2.0*dx__cm)
			k_measured__1_cm				= np.abs( Q_x/deltaT__ev_cm)
			
			#spitzer calcs....
			v_thermal__cm_s 				= 4.19e7*math.sqrt(temp_eq__ev)
			q_free 							= density__1_cm3*temp_eq__ev*v_thermal__cm_s
			spitzer_heat_flow 			= electron_mfp_cm / (temp_eq__ev/abs(slope)) * q_free 
			spitzer_k 						=  spitzer_heat_flow / abs(slope)
			
			other_spitzer_k1 				= density__1_cm3*8.617e-5*ee_collsion_time*v_thermal__cm_s*v_thermal__cm_s
			other_spitzer_k2				= density__1_cm3*ee_collsion_time*v_thermal__cm_s*v_thermal__cm_s
		
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
	def heat_conductivity_x(self, spec, F, v2, nu, qx, density=None, kx = None):
		if kx == None:
			kx = np.empty( (self.num_data_cells_x) )
		
		kx[:] = v2[:]
		kx *= v2
		kx = self.derivitive_1d_dx(kx)
		kx *= (-0.5 / (self.state.spatial_axes[0].dx()) )
		# put some value for boundry
		kx[ 0] = kx[ 1]
		kx[-1] = kx[-2]
		kx *= density
		kx /= nu
		
		
		"""
		clf();plot(kx); title("debug conducitvity "); show();		
		temp = ones( self.num_data_cells_x)
		np.divide( temp, kx, temp)
		clf(); plot(temp); title("debug conducitvity "); show();
		"""


		kx *= (-1.0)/18.0
		# TODO: filter small values
		#print qx
		for ix in xrange(0, self.num_data_cells_x):
			if (kx[ix] < self.FLT_MIN4 and kx[ix] > -1.0*self.FLT_MIN4):
				kx[ix] = 1.0
			kx[ix] = qx[ix]/ kx[ix]
		
		"""
		clf(); plot(kx); title("debug conducitvity "); show();		
		#np.divide(qx, kx, kx)
		"""
		
		return kx
	
	
	def derivitive_1d_dx(self, vals):
		for i in range(0, len(vals)-2):
			vals[i] -= vals[i+2]
		for i in range(len(vals)-3, -1, -1):
			vals[i+1] = vals[i]
		return vals
	
	# TODO: 'vz' is not used in the Q_x calculation
	def calculate_qx(self, spec, F, vx, vy, vz, vxvy, vxvx, vyvy, vzvz, qx=None, density=None, calc_density=False):
		if qx == None:
			qx = np.empty( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)		
		
		qx[:] = 0.0
		
		# does the -<v^2>*vx term
		temp = np.add(vxvx, vyvy)
		np.add(temp, vzvz, temp)
		temp *= vx
		qx -= temp
		
		# does the -2*(VxVx*Vx+Vxy*Vy) term
		np.multiply(vxvx, vx, temp)
		temp += (vxvy*vy)
		temp *= (2.0)
		qx -= temp
		
		"""
		# this is the 1D way.. it is equailiivent..
		crap = np.zeros( (self.num_data_cells_x) )
		# other way...
		temp[:] = 0.0
		temp = np.add(vxvx, vyvy)
		np.add(temp, vzvz, temp)
		temp *= vx		
		crap -= temp
		
		temp = np.multiply(vx,vx,temp)
		temp -= vxvx
		temp *= vx
		temp *= 2.0
		crap += temp
		print "---------------------1D way-------------"
		print crap
		"""
		
		
		# does the +2*(Vx*Vx + Vy*Vy)*Vx term
		np.multiply(vx, vx, temp)
		temp += (vy*vy)
		temp *= vx
		temp *= 2.0
		qx += temp
		
		qx *= density
		
		#--------------------
		# get <v^2*vx>
		temp = self.calculate_qx_part_2(spec, F, temp, density)
		qx += temp
		qx *= 0.5
		
		return qx
		
		
	def calculate_qx_part_2(self, spec, F, qx2, density):
		functs = {}
		target_harmonics = [(1,0, self.F10)]
		
		def p_axis_prep():
			np.power(self.pr, 5, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		
		def prep_first_last_p_cell():
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		
		#def compenstate_first_last_p_cell(self):
		#functs['compenstate_first_last_p_cell'] = prep_first_last_p_cell
		def f00_factor():
			factor = np.power( self.pr[0], 7)/ (  self.momentum_axis.dx()* self.pr[1] * 7.0 )
			return factor
		functs['f00_factor'] = f00_factor
		
		def overall_factor():
			return (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		
		qx2 = self.do_moment(spec, F, functs, target_harmonics, result=qx2)
		return qx2
	
		
	"""
		Calculates the noralized electron-ion collision frequency (i.e. colision frequnecy relative to the plasma frequnecy)
			Major usage if in the heat conduction coefficicent calcualtion.

			TODO: update the log lamda to handle all cases and use centralized function to calculate

			Following Tzoufras 2011 (and standard lititure), equation 44 says:
				lambda_MFP=v_thermal*t_ei = 9* ND / (sqrt(2/pi) * log_lambda ) * lambda_debye
			now divide by lambda_debye on both sides (and use v_thermal/lambda_debye = w_p)
				v_thermal*t_ei/lambda_debye = w_p*t_ie = 9* ND / (sqrt(2/pi) * log_lambda )
			Now ( w_p*t_ie ) is the relaxtion time, in normalized simulation units. So the normalized frequency is one over this:
				1/( w_p*t_ie ) = frequency = (sqrt(2/pi) * log_lambda ) /  9* ND

			First, Log lambda is calucalted using log_lambda_ei = 24 - ln( n_e^1/2*Te^1/2) (Te is elctron temperature in ev, n_e is electron nuber density in non-normalizey)
			The Zeta (charge of ion species) factor is multiped in with log lambda as it enhances the colsionality linearly in Z)
				Multiplying log lambda is the primary place that the ion charge state (Zeta) comes into play
	"""
	def calculate_collision_frequency(self, spec, F, nu=None, temp_ev=None, v2=None, ND_=None, density=None, calc_density=False):
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		if ND_ == None:
			ND_ = np.zeros( (self.num_data_cells_x) )
			self.calculate_ND(self, spec, F, temp_ev=temp_ev, v2=v2, ND=ND_, density=density, calc_density=False)

		if temp_ev != None:
			nu[:] = temp_ev[:]
		elif v2 != None:
			# in this case, quickly calcuate the temperature (in ev) from the v2 calcualtion
			#	TODO: this can be removed.
			nu[:] = v2[:]
			nu *= 170880.776
		else:
			# in this case, quickly calcuate the temperature (in ev) from the v2 calcualtion
			#	and also calcualte v2 from scratch
			#	TODO: this can be removed.
			v2 = np.zeros( (self.num_data_cells_x) )
			self.calculate_v2(spec, F, v2 = v2, density = density)
			nu[:] = v2[:]
			nu *= 170880.776
		np.log( nu, nu)
		temp = np.multiply(density, spec.species_config.density_np)
		np.log(temp, temp); temp *= 0.5; nu -= temp; nu += 24.0
		nu *= spec.species_config.zeta
		factor = np.sqrt( 2.0/np.pi) / 9.0; nu *= factor
		nu /= ND_
		
	

	"""
		Calculates the number of particles in a debye sphere at each spatial point.
			This is a dimensionless parameter.
			The calculation method is pulled straight from the plasma formulary
				calculate_ND 	=4/3*PI*number_density*debye_length^3 = 1.72e9*temperature_in_ev^(3/2)/sqrt(number_density)
								=1.72e9*sqrt(temperature_in_ev^3/number_density)

	"""
	def calculate_ND(self, spec, F, ND=None, temp_ev=None, v2=None, density=None, calc_density=False):
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		if temp_ev != None:
			ND[:] = temp_ev[:]
		elif v2 != None:
			ND[:] = v2[:]
			ND *= 170880.776
		else:
			v2 = np.zeros( (self.num_data_cells_x) )
			self.calculate_v2(spec, F, v2 = v2, density = density)
			ND[:] = v2[:]
			ND *= 170880.776
			
		np.power( ND , 3, ND)
		ND /= spec.species_config.density_np
		ND /= density
		
		np.sqrt(ND, ND)
		ND *=  1.72e+9
		return ND
		
		
	"""
		Calculates the pressure in Mbar
			The basic equation is:	p = nKT
			Convert 1/cm^3 to 1/m^3 and eV to joules to get Pascals:
				p 	= (n_in_1_cm^3*10^6)*(T_in_ev*1.6022e-19J) = n_in_1_cm^3*T_in_ev*1.6022e-13 Pascals
			There are 10^5 pascals in a bar so
				p 	= n_in_1_cm^3*T_in_ev*1.6022e-18 bar = n_in_1_cm^3*T_in_ev*1.6022e-24 Megabar
					= 2.73785179e-19 * n_in_1_cm^3 * <v^2> (v in rest frame)
	"""
	def calculate_presure(self, spec, F, pressure=None, temp_ev=None, v2=None, density=None, calc_density=False):
		if pressure == None:
			pressure = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		if v2 == None:
			if temp_ev != None:
				v2 = np.zeros( (self.num_data_cells_x) )
				v2[:] = temp_ev[:]
				v2 /= 170880.776
			else:
				v2 = np.zeros( (self.num_data_cells_x) )
				self.calculate_v2(spec, F, v2 = v2, density = density)
		
		pressure[:] = v2[:]
		temp = 2.73785179e-19 * spec.species_config.density_np
		pressure *= temp
		pressure *= density
		return pressure
	
	# TODO: temperature calcualtion is not done yet.. need to subtract out drift velopcites
	#		<v^2 > = <Vx^2> + <Vy^2> + <Vz^2> - <Vx>^2 - <Vy>^2 - <Vz>^2
	#		And I am missing the last 3 terms (which have always been zeros so far)
	def calculate_temp_ev(self, spec, F, temp_ev=None, v2=None, density=None, calc_density=False):
		
		
		if temp_ev == None:
			temp_ev = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		if v2 == None:
			v2 = np.zeros( (self.num_data_cells_x) )
			self.calculate_v2(spec, F, v2 = v2, density = density)
			
		
		
		
		# setup the momentum integrating factor (momentum^4 * Delta_p)
		axis_temp = np.multiply( self.pr, self.pr)
		axis_temp *= axis_temp
		axis_temp *= self.pr__dp 		# multipy by the velocity space cell size.
		
		
		# apply the temperature momentum integrating factors the the L=0, M=0 harmonic weights simultaniously.
		#	(	.mpaxis( temp_ev ) = self.F0[x, :] *= temp_ev[:]	where x in range(0, self.num_data_cells_x)
		spec.F.copy_into_harmonic_centric_view(0,self.F0)
		self.F0.mpaxis( axis_temp )
		
		# do the integral....
		temp_ev[:] = np.sum( self.F0.data, axis = 1)
		
		# normalize by the density cell by cell
		temp_ev /= density
		
		# take care of other normalizations
		temp_ev *=  4.0*np.pi/3.0
		
		# conversion factor is c^2*me/qe
		ev_convert_factor = (2.99792458e8)*(2.99792458e8)*(9.1093829e-31)/(1.602176565e-19);
		
		# calculate, into a temp variable:  1+p^2  (i.e. 1/gamma^2 )
		f00_ = np.ones(  (self.num_data_cells_x) )
		f00_ += temp_ev
		
		# temp_ev will now have the (non-relativisitic) velocity squared which we will use to get the (non relativisitic) temperature.
		temp_ev /= f00_
		temp_ev *= ev_convert_factor
		
		return temp_ev
		
		
	"""
		# old version of temperature
		# Assumes ths < v^2 > is the non-relatisitic total velocity squared.
		temp_ev[:] = v2[:]
		
		temp_ev *= 170880.776
		return temp_ev
		
	"""
	
	def temp_reflexsive( self, spec, F, density, num_data_cells, num_p_cell, pr):
		temp_ev = np.zeros( (self.num_data_cells_x) )
		axis_temp = np.multiply( pr, pr)
	
			
	def calculate_v2(self, spec, F, vx, vy, vz, vxvx, vyvy, vzvz, v2=None, density=None, calc_density=False):
		if v2 == None:
			v2 = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		# v2 = <vx2>+<vy2>+<vz2> - <vx>^2 - <vy>^2 - <vz>^2
		v2[:] = vxvx[:]
		v2 += vyvy
		v2 += vzvz
		temp = np.multiply(vx, vx)
		v2 -= temp
		np.multiply(vy, vy, temp)
		v2 -= temp
		np.multiply(vz, vz, temp)
		v2 -= temp
		
		
		#print "vxvx"; print vxvx[0:10]; print ""; print ""
		#print "vyvy"; print vyvy[0:10]; print ""; print ""
		#print "vzvz"; print vzvz[0:10]; print ""; print ""
		#print "vx"; print vx[0:10]; print ""; print ""
		#print "vy"; print vy[0:10]; print ""; print ""
		
		return v2

	# non relativisitic velocity offdiagonal tensor terms......
	def calculate_vyvz(self, spec, F, vyvz = None, density=None, calc_density=False):
		if vyvz == None:
			vyvz = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		target_harmonics = [(2,2, self.F22)]
		normalization = (-1.0) * 16.0 * np.pi / 5.0 * self.momentum_axis.dx()
		vyvz = self.calculate_V_V_(spec, F, target_harmonics, normalization, vyvz, density)
		return vyvz
	
	def calculate_vxvz(self, spec, F, vxvz = None, density=None, calc_density=False):
		if vxvz == None:
			vxvz = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		target_harmonics = [(2,1, self.F21)]
		normalization = (-1.0) * 8.0 * np.pi / 5.0 * self.momentum_axis.dx()
		vxvz = self.calculate_V_V_(spec, F, target_harmonics, normalization, vxvz, density)
		return vxvz
	
	def calculate_vxvy(self, spec, F, vxvy = None, density=None, calc_density=False):
		if vxvy == None:
			vxvy = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		target_harmonics = [(2,1, self.F21)]
		normalization = 8.0 * np.pi / 5.0 * self.momentum_axis.dx()
		vxvy = self.calculate_V_V_(spec, F, target_harmonics, normalization, vxvy, density)
		return vxvy
		
	def calculate_V_V_(self, spec, F, target_harmonics, overall_factor_, result, density):
	
		functs = {}		
		#target_harmonics = [(2,1, self.F21)]
		
		def p_axis_prep():
			np.power(self.pr, 4, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		def prep_first_last_p_cell():
			target_harmonic_data = target_harmonics[0][2]
			# For higher harmonics, OSHUN has a filter that zeros the |p|(0) coeffecient. This is a quick check to warn the user if this is not zero.
			#	It should not ever occur unless this filter is disabled.
			if( (target_harmonic_data.data[:, 0] != 0.0).any() ):
				print "WARNING!!! p(p_min)!=0.0 in a moment calculation on a higher order harmonin (L=%d,m=%d)" % ( target_harmonics[0][0], target_harmonics[0][1] )
			# Just in case the |p|(0) mode is not zero, still do not include it in the overall intergral.
			# 	so that the overall contribution form the lowest and highest |p| terms just has highest |p|: 1/2*(|p|(-1))
			#	So we set self.moment to  1/2*(|p|(-1))
			#	In the do_moment routine (where the calcuation is done), there is a variable moment_temp that has value  |p|: (|p|(0)+|p|(-1))
			#	Also in do_moment, a sum over all modes is done and added to self.moment.. then moment_temp is subtracted to undo oversounting of these modes
			#	The upshot of all this is is any non-zzero |p|(0) valuw will be subtracted from the overall sum.
			self.moment[:] = target_harmonic_data.data[:, -1]
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		def overall_factor():
			return overall_factor_
			#return (8.0 * np.pi / 5.0) * self.momentum_axis.dx()
		functs['overall_factor'] = overall_factor			

		result = self.do_moment(spec, F, functs, target_harmonics, result=result)
		result /= density
		return result
		
		
	def calculate_vzvz(self, spec, F, vzvz=None, moment1=None, moment2=None, density=None, calc_density=False):
		
		if vzvz == None:
			vzvz = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		vzvz[:] = moment1
		vzvz -= moment2
		vzvz /= density
		
		
		return vzvz
		
	def calculate_vyvy(self, spec, F, vyvy=None, moment1=None, moment2=None, density=None, calc_density=False):
	
		if vyvy == None:
			vyvy = np.zeros( (self.num_data_cells_x) )
		if moment2 == None:
			moment2 = np.zeros( (self.num_data_cells_x) )
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		# moment1 is needed for now to get VyVy
		if moment1 == None:
			moment1 = np.zeros( (self.num_data_cells_x) )
			calculate_vxvx(self, spec, F, moment1=moment1, density=density)
		
		vyvy = self.calculate_vyvy_1(spec, F, vyvy, moment1, density)
		moment2[:] = vyvy		# To be used for VzVz

		
		#finsih up vyvy
		vyvy += moment1		#  This is the integral 16pi/5f22 + 4pi/3f00 - 4pi/15*f20 
		vyvy /= density
		
		
		return vyvy
	
	def calculate_vyvy_1(self, spec, F, vyvy, moment1, density):
		functs = {}		
		target_harmonics = [(2,0, self.F20)]
		
		def p_axis_prep():
			np.power(self.pr, 4, self.axis_temp)
			#FINK: when fixing temperature, I commented out following 2 lines.. I think they should be there.. keep this comment in until sure + tested.
			# devide by 1/gamma^2 to convert to non-relativistic velocity
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		def prep_first_last_p_cell():
			target_harmonic_data = target_harmonics[0][2]
			# For higher harmonics, OSHUN has a filter that zeros the |p|(0) coeffecient. This is a quick check to warn the user if this is not zero.
			#	It should not ever occur unless this filter is disabled.
			if( (target_harmonic_data.data[:, 0] != 0.0).any() ):
				print "WARNING!!! p(p_min)!=0.0 in a moment calculation on a higher order harmonin (L=%d,m=%d)" % ( target_harmonics[0][0], target_harmonics[0][1] )
			# Just in case the |p|(0) mode is not zero, still do not include it in the overall intergral.
			# 	so that the overall contribution form the lowest and highest |p| terms just has highest |p|: 1/2*(|p|(-1))
			#	So we set self.moment to  1/2*(|p|(-1))
			#	In the do_moment routine (where the calcuation is done), there is a variable moment_temp that has value  |p|: (|p|(0)+|p|(-1))
			#	Also in do_moment, a sum over all modes is done and added to self.moment.. then moment_temp is subtracted to undo oversounting of these modes
			#	The upshot of all this is is any non-zzero |p|(0) valuw will be subtracted from the overall sum.
			self.moment[:] = target_harmonic_data.data[:, -1]
			self.moment *= 0.5

		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		def overall_factor():
			return (16.0 * np.pi / 5.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor			
		
		vyvy = self.do_moment(spec, F, functs, target_harmonics, result=vyvy)
		return vyvy		
	
	def calculate_vxvx(self, spec, F, vxvx=None, moment1=None, density=None, calc_density=False):
	
		if vxvx == None:
			vxvx = np.zeros( (self.num_data_cells_x) )
		if moment1 == None:
			moment1 = np.zeros( (self.num_data_cells_x) )		
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		functs = {}		
		target_harmonics = [(2,0, self.F20)]
		
		
		moment1 = self.calculate_vxvx_part1(spec, F, vxvx_1=moment1)
		self.moment_temp__no_guard_cells  = self.calculate_vxvx_part2(spec, F, vxvx_2=self.moment_temp__no_guard_cells )
		
		#print "vxvx part1"; print moment1[0:10]; print ""; print ""
		#print "vxvx part2"; print self.moment_temp__no_guard_cells[0:10]; print ""; print ""
		
		vxvx[:] = moment1
		vxvx += self.moment_temp__no_guard_cells 
		vxvx /= density
		
		# move on....
		self.moment_temp__no_guard_cells  *= 0.5;
		moment1 -= self.moment_temp__no_guard_cells ;  # This is the integral 4pi/3f00 - 4pi/15*f20
				
		return vxvx
		
	def calculate_vxvx_part2(self, spec, F, vxvx_2=None, calc_density=False):
		if vxvx_2 == None:
			vxvx_2 = np.zeros( (self.num_data_cells_x) )
		functs = {}		
		target_harmonics = [(2,0, self.F20)]
		def p_axis_prep():
			np.power(self.pr, 4, self.axis_temp)
			#FINK: when fixing temperature, I commented out following 2 lines.. I think they should be there.. keep this comment in until sure + tested.
			# devide by 1/gamma^2 to convert to non-relativistic velocity
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep
		def prep_first_last_p_cell():
			target_harmonic_data = target_harmonics[0][2]
			# For higher harmonics, OSHUN has a filter that zeros the |p|(0) coeffecient. This is a quick check to warn the user if this is not zero.
			#	It should not ever occur unless this filter is disabled.
			if( (target_harmonic_data.data[:, 0] != 0.0).any() ):
				print "WARNING!!! p(p_min)!=0.0 in a moment calculation on a higher order harmonin (L=%d,m=%d)" % ( target_harmonics[0][0], target_harmonics[0][1] )
			# Just in case the |p|(0) mode is not zero, still do not include it in the overall intergral.
			# 	so that the overall contribution form the lowest and highest |p| terms just has highest |p|: 1/2*(|p|(-1))
			#	So we set self.moment to  1/2*(|p|(-1))
			#	In the do_moment routine (where the calcuation is done), there is a variable moment_temp that has value  |p|: (|p|(0)+|p|(-1))
			#	Also in do_moment, a sum over all modes is done and added to self.moment.. then moment_temp is subtracted to undo oversounting of these modes
			#	The upshot of all this is is any non-zzero |p|(0) valuw will be subtracted from the overall sum.
			self.moment[:] = target_harmonic_data.data[:, -1]
			self.moment *= 0.5

		functs['prep_first_last_p_cell'] = prep_first_last_p_cell		
		def overall_factor():
			return (8.0 * np.pi / 15.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		vxvx_2 = self.do_moment(spec, F, functs, target_harmonics, result=vxvx_2)
		return vxvx_2
		
	def calculate_vxvx_part1(self, spec, F, vxvx_1=None, calc_density=False):
		
		if vxvx_1 == None:
			vxvx_1 = np.zeros( (self.num_data_cells_x) )
		functs = {}
		
		target_harmonics = [(0,0, self.F0)]
		def p_axis_prep():
			# create sampled values of p_x^4/gamma^2
			np.power(self.pr, 4, self.axis_temp)
			#FINK: when fixing temperature, I commented out following 2 lines.. I think they should be there.. keep this comment in until sure + tested.
			# devide by 1/gamma^2 to convert to non-relativistic velocity
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
			np.multiply(self.axis_temp, self.inverse_gamma_axis.values, self.axis_temp)
		functs['p_axis_prep'] = p_axis_prep		
		def prep_first_last_p_cell():
			self.moment *= 0.5
		functs['prep_first_last_p_cell'] = prep_first_last_p_cell
		def f00_need_pr1_and_pr0(F_pr1, F_pr0):
			#temp0 = np.power(self.pr[0], 4 )*self.inverse_gamma_axis.values[0]*self.inverse_gamma_axis.values[0]
			#temp1 = np.power(self.pr[1], 4 )*self.inverse_gamma_axis.values[1]*self.inverse_gamma_axis.values[1]
			F_pr1 /= self.axis_temp[1]
			F_pr0 /= self.axis_temp[0]
			
			# undo the axis scalling..
			F_pr1 *= self.p0p1_sq; F_pr0 -= F_pr1; F_pr0 *= self.inv_mp0p1_sq;
			factor = ( (1.0/5.0) - ( (1.0/7.0)*self.p0p1_sq) ); F_pr0 *= factor;
			
			F_pr1 /= 7.0; F_pr1 += F_pr0
			factor = np.power( self.pr[0], 5) / self.momentum_axis.dx()
			F_pr1 *= factor
			
			#print ("low bound %e .. main %e" % ( F_pr1[2], self.moment[2]) )
		functs['f00_need_pr1_and_pr0'] = f00_need_pr1_and_pr0
		def overall_factor():
			return (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		functs['overall_factor'] = overall_factor
		
		vxvx_1 = self.do_moment(spec, F, functs, target_harmonics, result=vxvx_1)
		return vxvx_1		

		"""
	# unused.. here for legacy.. erase it.
	def calculate_vx(self, spec, F, vx=None, density=None, calc_density=False):
		
		spec.F.copy_into_harmonic_centric_view(0,0,self.F10)
		if vx == None:
			vx = np.zeros( (self.num_data_cells_x) )
		
		#If density provided, use it.. unless calc_density=True
		#If density not prvided, so ahead and calculate it.
		if density == None:
			density = self.calculate_density1(spec, F)
		elif calc_density:
			density = self.calculate_density1(spec, F, density=density)
		
		self.axis_temp[:] = self.pr[:] 
		self.axis_temp *= self.pr[:]
		self.axis_temp *= self.pr[:]
		self.axis_temp *= self.inverse_gamma_axis.values[:]
		#TODO: mpaxis without boundry cells..
		self.F10.mpaxis( self.axis_temp )
		
		
		#Do first and last |p| cell (At all points in space)
		self.moment = self.F10.data[:, 0]
		self.moment += self.F10.data[:, -1]
		self.moment_temp[:] = self.moment[:]
		self.moment *= 0.5

		# sum over |p| to needed value at spatial point... but we have also included 
		#	the 0,1 cells of |p|.. s they are double counted.. so we subtract them out.
		self.moment += np.sum( self.F10.data, axis = 1)
		self.moment -= self.moment_temp
		
		factor = self.pr[0]*self.pr[0]
		factor *= factor
		factor *= self.pr[0]
		factor /= ( self.momentum_axis.dx()*self.pr[1]*5.0 )
		
		self.f00_[:] = self.F10.data[:,  1]
		self.f00_ *= factor
		self.moment += self.f00_
		
		vx[:] = self.moment[ self.num_bdry_cells_left__x:(self.num_data_cells_x+self.num_bdry_cells_left__x) ]
		factor = (4.0 * np.pi / 3.0) * self.momentum_axis.dx();
		vx *= factor
		
		#divide by density
		vx /= density
		
		
		#print "---------------- old-----------"
		#print vx
		#print "---------------- new-----------"
		self.calculate_vx1(spec, spec.F)
		#print "---------------- --------------"
		return vx
	


	MomentCalculations: old an anot used
	def calculate_density(self, spec, F, density=None):
		# yes.. this is bad Python form.. but it's here to avoid cicural dependancies
		# during development.. later I will offically pslit the file in a nice way.
		from ..Oshun1d import DistributionFunction, SHarmonic, DistributionFunctionSpatialCell, Species, HarmonicCentricView, Config, SpeciesConfig, DaState
		F = spec.F
		if density == None:
			density = np.zeros( (self.num_data_cells_x) )
		F.copy_into_harmonic_centric_view(0,self.F0)

		for ix in xrange(0, self.num_data_cells_x):

			denisity_norm_factor = 4.0 * np.pi * self.momentum_axis.dx()
			# strange correction...
			self.f00  =  self.F0.data[ix+self.num_bdry_cells_left__x, 1]/5.0*self.p0p1_sq
			self.f00 += ( (self.F0.data[ix+self.num_bdry_cells_left__x, 0] - self.F0.data[ix+self.num_bdry_cells_left__x, 1]*self.p0p1_sq)
							*self.inv_mp0p1_sq*(1.0/3.0-1.0/5.0*self.p0p1_sq) )
			self.f00 *= denisity_norm_factor
			
			
			#normal-ish
			self.f00 *= ( self.pr[0]*(self.pr[0]/self.momentum_axis.dx())*self.pr[0] )
			density[ix] =  self.F0.data[ix+self.num_bdry_cells_left__x,  0] * self.pr[ 0] * self.pr[ 0]
			density[ix] += self.F0.data[ix+self.num_bdry_cells_left__x, -1] * self.pr[-1] * self.pr[-1]
			density[ix] *= 0.5
			
			#normal plus add in strange correction
			for p in xrange(1, self.momentum_axis.num_points - 1):
				density[ix] += self.pr[p] * self.pr[p] * self.F0.data[ix+self.num_bdry_cells_left__x, p]*denisity_norm_factor
				density[ix] += self.f00
		#print np.sum(density)
		return density
	"""