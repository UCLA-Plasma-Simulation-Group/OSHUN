"""
Module that calcuates the scutal numerical distribution function by evaluating the 
	spherical harmonic expansion (in velocity space) on a cartesian velocity grid.
	(Space is already on a cartestian grid.. so in essense the velocity projection is simply repeated
		at each spatial cell to get a the full cartesian ditribution function)
"""

import numpy as np
import numpy.polynomial.legendre as legendre
import scipy
from ..OshunMpi import *

class ActualDistrubtionFunctionPScafolding:
	def __init__(self, abs_p_value, abs_p_index):
		self.abs_p_value = abs_p_value
		self.abs_p_index = abs_p_index
	
	
class ActualDistrubtionFunction:
	def __init__(self, F):
		
		self.num_harmonics_per_cell = F.get_num_harmonics_per_cell()
		self.final_data = np.zeros( ( (F.num_plain_cells[0] ), 2*F.momentum_axis.num_points ) )
		self.abs_p_temp = np.zeros( F.num_l )
		self.abs_p_temp_with_m = np.zeros( self.num_harmonics_per_cell )
		
		
		self.da_real_distribution_funcion_last_update = -1
		self.nump = F.momentum_axis.num_points
		
		self.full_distbution_function = None
		
		self.momenum_axis_sq = np.power( F.momentum_axis.values[:], 2)
		self.momenum_axis_inverse = np.zeros( self.nump )
		self.momenutm_axis_deltas = np.zeros( self.nump )
		
		for i in xrange(0, self.nump):
			self.momenum_axis_inverse[i] = 1.0/F.momentum_axis.values[i]
		
		#momenum_axis_inverse= 1.0 / F.momentum_axis.values[:]
		self.momenutm_axis_deltas[:] = F.momentum_axis.values[:]
		self.momenutm_axis_deltas[1:] -= F.momentum_axis.values[0:self.nump-1]
		
		
		# init these once
		#self.vals_plm_at_thetas 							= zeros( F.momentum_axis.num_points )
		#self.vals_plm_at_thetas_minus 					= zeros( F.momentum_axis.num_points )
		self.rsintheta_delta 								= np.zeros ( F.momentum_axis.num_points )
		self.rsintheta_cell 									= np.zeros ( F.momentum_axis.num_points )
		self.values_at_a_px_temp 							= np.zeros( (F.momentum_axis.num_points) )
		self.values_at_a_px_temp_minus 					= np.zeros( (F.momentum_axis.num_points) )
		self.values_at_a_px_accumlator 					= np.zeros( (F.momentum_axis.num_points) )
		self.values_at_a_px_accumlator_minus 			= np.zeros( (F.momentum_axis.num_points) )
		self.values_at_a_px_accumlator_temp 					= np.zeros( (F.momentum_axis.num_points) )
		self.valid_absp_cells								= np.zeros( (F.momentum_axis.num_points) )
		self.valid_absp_cells_temp								= np.zeros( (F.momentum_axis.num_points) )
		
		self.F_temp_buffer = np.zeros( (F.num_plain_cells[0], F.num_l, F.momentum_axis.num_points) )  #DistributionFunction(F.num_plain_cells[0], F.num_l, F.momentum_axis, F.num_boundry_cells)
		self.F_minus_temp_buffer  = np.zeros( (F.num_plain_cells[0], F.num_l, F.momentum_axis.num_points) )  #DistributionFunction(F.num_plain_cells[0], F.num_l, F.momentum_axis, F.num_boundry_cells)
		self.F_temp_buffer_save = np.zeros( (F.num_plain_cells[0], F.num_l, F.momentum_axis.num_points) )
		
		self.setup_calc_structure( F )
		
	def get_full_F( self, F, n, t, get_mb=False):
		if self.da_real_distribution_funcion_last_update == n:
			return self.full_distbution_function
		num_bndry_cells_lower = F.num_boundry_cells[0][0]
		num_cells = F.num_plain_cells[0]
		
		for x in range( 0, num_cells ):
			f_at_cell_x = np.array(F.momentum_axis.num_points)
			for p in xrange( 0, F.momentum_axis.num_points ):
				# for each |p|, grab  all the harmonics l at this cell!
				for l in xrange(0, F.num_l ):
					index = F.linear_index_for_harmonic(l,0)
					self.abs_p_temp[l] = F.cells[num_bndry_cells_lower+x].harmonics[index].momentum_projection[p]
				vals = legendre.legval([-1.0, 1.0], self.abs_p_temp)
				self.final_data[x, self.nump-p-1] = vals[0]
				self.final_data[x, self.nump+p] = vals[1]
		if mpi_gather.gather_2d( "out_f_slice", 	STATE.spatial_axes__global[0].num_points, 
											2*F.momentum_axis.num_points, self.final_data):
			# only rank 0 will get here.
			self.full_distbution_function = mpi_gather.get_2d_data("out_f_slice")
		else:
			self.full_distbution_function = None
		
		return self.full_distbution_function
	
	def setup_calc_structure(self, F):
		num_bndry_cells_lower = F.num_boundry_cells[0][0]
		num_cells = F.num_plain_cells[0]
		self.abs_p_calc_scafolding = [None]*F.momentum_axis.num_points
		m = 0
		self.twopi = 2.0*np.pi 
		
		self.all_l_at_once_buffer = np.zeros( (F.num_l, F.momentum_axis.num_points) )
		self.all_l_at_once_buffer_minus = np.zeros( (F.num_l, F.momentum_axis.num_points) )
		
		for p in xrange( 0, F.momentum_axis.num_points ):
			scafolding_for_p = ActualDistrubtionFunctionPScafolding(  F.momentum_axis.values[p], p)
			self.abs_p_calc_scafolding[p] = scafolding_for_p
			
			current_value_px = F.momentum_axis.values[p]
			current_value_px2 = current_value_px*current_value_px
			
			scafolding_for_p.vals_cartestian_p_cos_thetas 				= np.zeros( F.momentum_axis.num_points )	
			scafolding_for_p.vals_cartestian_p_transverse_step_size	= np.zeros( F.momentum_axis.num_points )	
			scafolding_for_p.vals_cartestian_p_thetas_minus				= np.zeros( F.momentum_axis.num_points )	
			
			scafolding_for_p.vals_cartestian_p_cos_thetas[p] 								= 1.0
			scafolding_for_p.vals_cartestian_p_thetas_minus[p] 							= -1.0
			scafolding_for_p.vals_cartestian_p_transverse_step_size[p]					= 0.0

			scafolding_for_p.vals_cartestian_p_cos_thetas[p+1:] 						= 	current_value_px*self.momenum_axis_inverse[p+1:]
			scafolding_for_p.vals_cartestian_p_thetas_minus[p+1:] 					=	scafolding_for_p.vals_cartestian_p_cos_thetas[p+1:]
			scafolding_for_p.vals_cartestian_p_thetas_minus[p+1:] 					*=	-1.0
			
			# no arccos because we are oging to immediatly calulate cos(theta) anyeway (for the Ylm)
			scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:] 			= self.momenum_axis_sq[p+1:]
			scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:] 			-= current_value_px2
			np.sqrt( scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:], scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:])
			
			# and also precompute all the legenre coeffecient we will need....
			scafolding_for_p.vals_plm_at_theta_for_lm 													= np.zeros( (F.num_l,F.momentum_axis.num_points) ) #[None] * F.num_l
			scafolding_for_p.vals_plm_at_thetas_minus_for_lm 											= np.zeros( (F.num_l,F.momentum_axis.num_points) ) #[None] * F.num_l
			for l in xrange(0,F.num_l):
				#scafolding_for_p.vals_plm_at_theta_for_lm[l]							=	zeros( F.momentum_axis.num_points )
				#scafolding_for_p.vals_plm_at_thetas_minus_for_lm[l]				=	zeros( F.momentum_axis.num_points )
				vals_plm_at_thetas 			= scafolding_for_p.vals_plm_at_theta_for_lm[l]
				vals_plm_at_thetas_minus	= scafolding_for_p.vals_plm_at_thetas_minus_for_lm[l]
				
				vals_plm_at_thetas[p:] 											= scipy.special.lpmv(m,l, scafolding_for_p.vals_cartestian_p_cos_thetas[p:])
				vals_plm_at_thetas_minus[p:] 									= scipy.special.lpmv(m,l, scafolding_for_p.vals_cartestian_p_thetas_minus[p:])
				
	#@profile_func(filename="testme.grind")
	def get_full_lm( self, F, n, t, STATE, get_mb=False):
		if self.da_real_distribution_funcion_last_update == n:
			return self.full_distbution_function
		
		num_bndry_cells_lower = F.num_boundry_cells[0][0]
		num_bndry_cells_upper = F.num_boundry_cells[0][1]
		
		num_cells = F.num_plain_cells[0]
		
		# in the scaffolding...
		#vals_cartestian_p_cos_thetas = zeros( F.momentum_axis.num_points )
		#vals_cartestian_p_transverse_step_size = zeros( F.momentum_axis.num_points )
		#vals_cartestian_p_thetas_minus = zeros( F.momentum_axis.num_points )
		
		# init these once
		#vals_plm_at_thetas = zeros( F.momentum_axis.num_points )
		#vals_plm_at_thetas_minus = zeros( F.momentum_axis.num_points )
		#rsintheta_delta = np.zeros ( F.momentum_axis.num_points )
		#rsintheta_cell = np.zeros ( F.momentum_axis.num_points )
		#values_at_a_px_temp = zeros( (F.momentum_axis.num_points) )
		#values_at_a_px_temp_minus = zeros( (F.momentum_axis.num_points) )
		#values_at_a_px_accumlator = zeros( (F.momentum_axis.num_points) )
		#values_at_a_px_accumlator_minus = zeros( (F.momentum_axis.num_points) )
		#thetas_temp = np.zeros(  F.momentum_axis.num_points )
		#thetas_delta  = np.zeros(  F.momentum_axis.num_points )
		
		
		last_p_seen = -1
		
		#print "starting dist calc."
		current_x = num_bndry_cells_lower
		
		self.F_temp_buffer_save[:,:,:] =  F.distribution_function_data[num_bndry_cells_lower:-num_bndry_cells_upper,:,:]
		
		for p in xrange( 0, F.momentum_axis.num_points ):
			nump_valid_values = self.nump - p
			self.F_temp_buffer[:,:,:nump_valid_values] 			= self.F_temp_buffer_save[:,:,p:]
			self.F_minus_temp_buffer[:,:,:nump_valid_values] 	= self.F_temp_buffer_save[:,:,p:]
			
			
			scafolding_for_p 														= 	self.abs_p_calc_scafolding[p] 
			vals_cartestian_p_cos_thetas 										=	scafolding_for_p.vals_cartestian_p_cos_thetas
			vals_cartestian_p_transverse_step_size 						=	scafolding_for_p.vals_cartestian_p_transverse_step_size
			vals_cartestian_p_thetas_minus 									=	scafolding_for_p.vals_cartestian_p_thetas_minus

			#print "p: %d" % p
			#print	" \t\t Big calc1"
			# this will do all x and l at same time in numpy...
			self.F_temp_buffer[:,:,:nump_valid_values] *= scafolding_for_p.vals_plm_at_theta_for_lm[:,p:]
			temp = np.sum( self.F_temp_buffer, axis=1)
			#self.values_at_a_px_accumlator[0:nump_valid_values] = temp[p:]
			
			#print	" \t\t Big calc2"
			self.F_minus_temp_buffer[:,:,:nump_valid_values] *= scafolding_for_p.vals_plm_at_thetas_minus_for_lm[:,p:]
			#self.F_minus_temp_buffer[:] *= scafolding_for_p.vals_plm_at_thetas_minus_for_lm
			temp_minus = np.sum( self.F_minus_temp_buffer, axis=1)
			#self.values_at_a_px_accumlator_minus[0:nump_valid_values] = temp[p:]
			
			#print	" \t\t Litterstuff."
			
			# we are integrating in shpereical coords... so we need a volume element p^2sin(theta) d(theta) d(phi) (dp)
			#		the "vals_cartestian_p_transverse_step_size" is [p*sin(theta)] dp is the itervals defined betwen out theta points.
			#		d_phi the the constant angle step size. So we need 
			
			offset = 0
			
			# calculate the delta(P*sin(theta)) between each evaluation point.
			self.rsintheta_delta[0:nump_valid_values-1] 						= 		scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:]
			self.rsintheta_delta[0:nump_valid_values] 						-=  	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
			
			
			# calculate the cells sizes that emcompass each evaluation point...
			self.rsintheta_cell[0] 													=  	self.rsintheta_delta[0]
			self.rsintheta_cell[nump_valid_values-1] 							=  	self.rsintheta_delta[nump_valid_values-2]
			if nump_valid_values > 1:
				self.rsintheta_cell[1:nump_valid_values-1] 					= 		self.rsintheta_delta[0:nump_valid_values-2]
				self.rsintheta_cell[1:nump_valid_values-1] 					+= 	self.rsintheta_delta[1:nump_valid_values-1]
			self.rsintheta_cell /= 2.0
			
			# multiply that interesting SH's value at the interction point by this cell size
			temp[:,0:nump_valid_values] 																	*=  	self.rsintheta_cell[0:nump_valid_values]
			temp_minus[:,0:nump_valid_values] 														*=		self.rsintheta_cell[0:nump_valid_values]
			
			temp[:,0:nump_valid_values] 												*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
			temp_minus[:,0:nump_valid_values] 									*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
			
			#temp[:,p:] 																	*=  	self.rsintheta_cell[0:nump_valid_values]
			#temp_minus[:,p:] 														*=		self.rsintheta_cell[0:nump_valid_values]
			#
			#temp[:,p:] 												*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
			#temp_minus[:,p:] 									*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
			
			#print "temp...."
			#print temp.shape
			# Do the integration sum (integration is just sum at this point since all the wighting factors were already multiplied in)
			value 																		= 		np.sum( temp[:, 0:nump_valid_values] , axis = 1)
			value_minus																	= 		np.sum( temp_minus[:, 0:nump_valid_values], axis = 1 )
			
			
			# save the final value... and don't forget to multiply by the 2*PI that comes from integration over the lonely azimuthal angle.
			value *= self.twopi; value_minus *= self.twopi;
			self.final_data[:, self.nump-p-1] = value[:]
			self.final_data[:, self.nump+p] = value_minus[:]
			
			#self.final_data[:, self.nump-p-1] 									= value_minus * 2.0*np.pi
			#self.final_data[:, self.nump+p]										= value * 2.0*np.pi
			
			"""
			for x in range( 0, num_cells ):
				# at this |p|, find angles for this |p| and lower....
				# do this |p|
				scafolding_for_p 														= 	self.abs_p_calc_scafolding[p] 
				vals_cartestian_p_cos_thetas 										=	scafolding_for_p.vals_cartestian_p_cos_thetas
				vals_cartestian_p_transverse_step_size 						=	scafolding_for_p.vals_cartestian_p_transverse_step_size
				vals_cartestian_p_thetas_minus 									=	scafolding_for_p.vals_cartestian_p_thetas_minus
				
				#self.values_at_a_px_accumlator[:] 								= 0.0
				#self.values_at_a_px_accumlator_minus[:] 						= 0.0
				
				# so now thetas_temp has all the theta angles of one 'ray' of the points...
				#	and theta_deltas has the d(theta) for each step (ignoing first which ois on the Px axis so need sperate treatment.
				nump_valid_values = self.nump - p
								
				# now apply this grid to each spherical harmonic...
				sum_over_all_harmonics = 0
				sum_over_all_harmonics_minus = 0
				
				lin_harmonic_idx = 0
				np.multiply(scafolding_for_p.vals_plm_at_theta_for_lm, F.distribution_function_data[current_x], out = self.all_l_at_once_buffer)
				np.multiply(scafolding_for_p.vals_plm_at_thetas_minus_for_lm, F.distribution_function_data[current_x], out = self.all_l_at_once_buffer_minus)
				np.sum(  self.all_l_at_once_buffer, axis= 0, out= self.values_at_a_px_accumlator_temp )
				self.values_at_a_px_accumlator[:nump_valid_values] = self.values_at_a_px_accumlator_temp[p:]
				
				np.sum(  self.all_l_at_once_buffer_minus, axis= 0, out= self.values_at_a_px_accumlator_temp )
				self.values_at_a_px_accumlator_minus[:nump_valid_values] = self.values_at_a_px_accumlator_temp[p:]
				
				
				for l in xrange(0,F.num_l):
					np.mult(scafolding_for_p.vals_plm_at_theta_for_lm, F.distribution_function_data[current_x], out = self.all_l_at_once_buffer)
					np.mult(scafolding_for_p.vals_plm_at_thetas_minus_for_lm, F.distribution_function_data[current_x], out = self.all_l_at_once_buffer_minus)
					
					
					#self.all_l_at_once_buffer_minus = np.zeros
					#scafolding_for_p.vals_plm_at_theta_for_lm[l][p:]
					
					
					m = 0
					index = F.linear_index_for_harmonic(l)
					# 
					# numpy's Associated legendre function (lpmv) uses the Condon-Shortly phase (-1)^m (checked)
					#
					# notice destintion index shift to 0->(nump_valid_values-1)
					#
					# can prolly speed up by bunchinng all 
					#self.vals_plm_at_thetas[p:] 											= scipy.special.lpmv(m,l, vals_cartestian_p_cos_thetas[p:])
					self.values_at_a_px_temp[0:nump_valid_values] 					= scafolding_for_p.vals_plm_at_theta_for_lm[l][p:] #self.vals_plm_at_thetas[p:]
					
					#self.vals_plm_at_thetas_minus[p:] 									= scipy.special.lpmv(m,l, vals_cartestian_p_thetas_minus[p:])
					self.values_at_a_px_temp_minus[0:nump_valid_values] 			= scafolding_for_p.vals_plm_at_thetas_minus_for_lm[l][p:]	#self.vals_plm_at_thetas_minus[p:]
					
					#self.valid_absp_cells[:nump_valid_values] 						= F.distribution_function_data[current_x, index, p:] #F.cells[current_x].harmonics[index].momentum_projection[p:]
					#self.values_at_a_px_temp[0:nump_valid_values] 					*= valid_absp_cells[p:] #valid_absp_cells
					#self.values_at_a_px_temp_minus[0:nump_valid_values] 			*= valid_absp_cells[p:]
					self.valid_absp_cells_temp[:] 										= F.distribution_function_data[current_x, index, :] #F.cells[current_x].harmonics[index].momentum_projection[p:]
					self.valid_absp_cells[:nump_valid_values] 						= self.valid_absp_cells_temp[p:]
					self.values_at_a_px_temp												*= self.valid_absp_cells
					self.values_at_a_px_temp_minus										*= self.valid_absp_cells
					
					# values_at_a_px_temp/minus now contains the value of the current l,m harmonic at the grid points that are valid for the current p (abs_p)
					# now copy values in with the other harmonics at this 
					self.values_at_a_px_accumlator[0:nump_valid_values] 			+= self.values_at_a_px_temp[0:nump_valid_values]
					self.values_at_a_px_accumlator_minus[0:nump_valid_values] 	+= self.values_at_a_px_temp_minus[0:nump_valid_values]
				
					
				# we are integrating in shpereical coords... so we need a volume element p^2sin(theta) d(theta) d(phi) (dp)
				#		the "vals_cartestian_p_transverse_step_size" is [p*sin(theta)] dp is the itervals defined betwen out theta points.
				#		d_phi the the constant angle step size. So we need 
				
				offset = 0
				
				# calculate the delta(P*sin(theta)) between each evaluation point.
				self.rsintheta_delta[0:nump_valid_values-1] 						= 		scafolding_for_p.vals_cartestian_p_transverse_step_size[p+1:]
				self.rsintheta_delta[0:nump_valid_values] 						-=  	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
				
				
				# calculate the cells sizes that emcompass each evaluation point...
				self.rsintheta_cell[0] 													=  	self.rsintheta_delta[0]
				self.rsintheta_cell[nump_valid_values-1] 							=  	self.rsintheta_delta[nump_valid_values-2]
				if nump_valid_values > 1:
					self.rsintheta_cell[1:nump_valid_values-1] 					= 		self.rsintheta_delta[0:nump_valid_values-2]
					self.rsintheta_cell[1:nump_valid_values-1] 					+= 	self.rsintheta_delta[1:nump_valid_values-1]
				self.rsintheta_cell /= 2.0
				
				# multiply that interesting SH's value at the interction point by this cell size
				self.values_at_a_px_accumlator[0:nump_valid_values] 			*=  	self.rsintheta_cell[0:nump_valid_values]
				self.values_at_a_px_accumlator_minus[0:nump_valid_values] 	*=		self.rsintheta_cell[0:nump_valid_values]
				
				self.values_at_a_px_accumlator[0:nump_valid_values] 			*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
				self.values_at_a_px_accumlator_minus[0:nump_valid_values] 	*= 	scafolding_for_p.vals_cartestian_p_transverse_step_size[p:]
				
				# Do the integration sum (integration is just sum at this point since all the wighting factors were already multiplied in)
				value 																		= 		np.sum( self.values_at_a_px_accumlator[0:nump_valid_values] )
				value_minus																	= 		np.sum( self.values_at_a_px_accumlator_minus[0:nump_valid_values] )
				
				# save the final value... and don't forget to multiply by the 2*PI that comes from integration over the lonely azimuthal angle.
				self.final_data[x, self.nump-p-1] 									= value_minus * 2.0*np.pi
				self.final_data[x, self.nump+p]										= value * 2.0*np.pi
			current_x += 1
			"""
		if mpi_gather.gather_2d( "out_f", 	STATE.spatial_axes__global[0].num_points, 
											2*F.momentum_axis.num_points, self.final_data):
			# only rank 0 will get here.
			#print "RANK 0 got here!!!"
			self.full_distbution_function = mpi_gather.get_2d_data("out_f")
		else:
			self.full_distbution_function = None
		#print ""
		
		return self.full_distbution_function
