import numpy as np
#from ctypes import *
#import ctypes as np__c
#import  numpy.ctypeslib
#import numpy.polynomial.legendre as legendre
#import math
from .diagnostics.moments import *
from .diagnostics.DistributionFunction import *
from OshunMpi import *
#from oshun_c import *
#import OshunGlobalState as OrshunGlobalState
from OshunGlobalState import NP_SAVE,SAVEFIG,PRINT
from matplotlib.pyplot import *
import matplotlib.pyplot as pyplot
import numpy.fft as fft
import time as time_utils


"""
		Was in the 1 time initalization in main_loop
"""
class OutputDiagnostics:

	def __init__(self, STATE):

		# Enumerate momentum space axes ( used to label axes when making distribution functions plots/views/slices)
		self.f_slice_axis = np.zeros( 2* STATE.species[0].F.momentum_axis.num_points)
		for ip in xrange(0, STATE.species[0].F.momentum_axis.num_points):
			pval=  STATE.species[0].F.momentum_axis.values[ip]
			self.f_slice_axis[STATE.species[0].F.momentum_axis.num_points-1-ip] = -1.0* pval
			self.f_slice_axis[STATE.species[0].F.momentum_axis.num_points+ip ] = pval 
		self.actual_dist_functions = []
		self.initial_dist = None
		
		# take care of output intervals..
		try:
			self.collsions_enabled = STATE.config.collsions_enabled
		except:
			self.collsions_enabled = False
		

		try:
			self.F_output_interval = STATE.config.F_output_interval
			if self.F_output_interval > 0:
				for i in range(0,  STATE.config.num_species):
					self.actual_dist_functions.append( ActualDistrubtionFunction(STATE.species[i].F) )
		except:
			STATE.config.F_output_interval = -1
		pass
		
		try:
			self.density_interval = STATE.config.density_interval
		except:
			self.density_interval = -1
		
		try:
			self.E_interval = STATE.config.E_interval
		except:
			self.E_interval = -1
		
		try:
			self.F_output_interval = STATE.config.F_output_interval
		except:
			self.F_output_interval = -1
		
		try:
			self.SH_diag_interval = STATE.config.SH_diag_interval
		except:
			self.SH_diag_interval = -1
		
		try:
			self.Q_interval = STATE.config.Q_interval
		except:
			self.Q_interval = -1
		
		try:
			self.temperature_interval = STATE.config.temperature_interval
		except:
			self.temperature_interval = -1
		
		try:
			self.heat_conduction_interval = STATE.config.heat_conduction_interval
		except:
			self.heat_conduction_interval = -1
		
		try:
			self.v2_interval = STATE.config.v2_interval
		except:
			self.v2_interval = -1
		
		try:
			self.F_output_velocity_interval = STATE.config.F_output_range
		except:
			self.F_output_velocity_interval = None
		try:
			self.F_output_range = STATE.config.F_output_range
		except:
			self.F_output_range = None
		try:
			self.F_output_log_scale = STATE.config.F_output_log_scale
		except:
			self.F_output_log_scale = False
		try:
			self.F_output_lineouts = STATE.config.F_output_lineouts
		except:
			self.F_output_lineouts = None

class SimpleTimer:
	def __init__(self):
		self.timers = {}
	
	def start_timer(self, time_label):
		try:
			if time_label in self.timers:
				timer_val = self.timers[time_label]
			else:
				timer_val = {}
				self.timers[time_label] = timer_val
				timer_val['end_val'] = []
				timer_val['interval'] = []
				timer_val['accumlated_time'] = 0.0
				timer_val['val'] = []
			timer_val['val'].append(time_utils.time())
		except:
			pass
	
	def end_timer(self, time_label):
		#try:
		end_time = time_utils.time()
		timer_val = self.timers[time_label]	
		interval = end_time - timer_val['val'][-1]
		timer_val['end_val'].append( end_time )		
		timer_val['interval'].append( interval )
		timer_val['accumlated_time'] += interval
		#except:
		#	pass
		return
	def getTimer(self, time_label):
		try:
			return self.timers[time_label]
		except:
			return None


class GeneralPerfMonitor:
	def __init__(self, STATE, simple_timer):
		"""
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS = []
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__EXCLUDE_IO = []
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__TIME_IN_IO = []

		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL_START_TIME = 0
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___TIMER_FOR_IO_ACCUMULATOR = 0.0
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL = 100
		TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL_TO_PLOT_MULTIPLIER = 100*2
		"""
		try:
			self.perf_measure_interval = STATE.config.perf_measure_interval
		except:
			self.perf_measure_interval = -1
		
		try:
			self.perf_measure_num_measures_before_output = STATE.config.perf_measure_num_measures_before_output
		except:
			self.perf_measure_num_measures_before_output = -1
		
		try:
			self.simple_timer = None
			self.simple_timer = simple_timer
		except:
			pass
		
	def take_measurement( self, STATE, n):
		if STATE.simple_perf_timer.perf_measure_interval > 0 and n%STATE.simple_perf_timer.perf_measure_interval == 0:
			self.simple_timer.end_timer("main_loop_timing")
			self.simple_timer.start_timer("main_loop_timing")
		
		if self.perf_measure_num_measures_before_output > 0 and n%self.perf_measure_num_measures_before_output == 0 and n > 0:
			pass
			"""
			temp = time_utils.time() - TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL_START_TIME
			
			TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS.append(  temp )
			TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__EXCLUDE_IO.append( temp - TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___TIMER_FOR_IO_ACCUMULATOR)
			TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__TIME_IN_IO.append( TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___TIMER_FOR_IO_ACCUMULATOR )
			
			TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL_START_TIME = time_utils.time() 
			TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___TIMER_FOR_IO_ACCUMULATOR = 0.0

			if mpi_info.rank == 0 and n%TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL_TO_PLOT_MULTIPLIER == 0:
				clf(); timer_events = np.arange( 0, len(TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS) );
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS, label="Time");
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__EXCLUDE_IO, label="Time w/o IO");
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__TIME_IN_IO, label="Time in IO");						
				legend(prop={'size':6})
				title('Time elapsed for %d loops as simulation progresses'% TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL); xlabel('Ordinal of a Chunk of %d loops'% TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL); ylabel('Seconds required to execute %d loops' % TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL )
				SAVEFIG("timing__timesteps")
				clf(); output_time1 = dt*float(TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL); timer_events *= output_time1;
				print dt
				print output_time1
				print timer_events[0]; print timer_events[1];
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS, label="Time");
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__EXCLUDE_IO, label="Time w/o IO");
				plot( timer_events, TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS__TIME_IN_IO, label="Time in IO");						
				legend(prop={'size':6})
				title('Time elapsed for %d loops versus simulation time.'%TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL); xlabel('Simulation Time'); ylabel('Seconds required to execute %d loops' % TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___INTERVAL )
				SAVEFIG("timing__simulation_time")
			"""
		pass
	pass
	

class DiagDoer:
	def __init__(self):
		self.simple_E_buffer = None
		self.moments = []
		
		pass
	
	def do_it_for_this_timestep(self, STATE, n, time, moments):
		#if((n%output_interval) == 0):
		
		# the following line magically fixes errors with pyploy state being garbage collected.
		# TODO: understand why		
		#pyplot.clf()

		if ( True ):
			for i in range(0,  STATE.config.num_species):
				species_name = STATE.species[i].species_config.name
				species_name_prefix = ""
				if  STATE.config.num_species > 1:
					species_name_prefix = species_name + "_"
				
				if True and STATE.out_diags.F_output_interval> 0 and n % STATE.out_diags.F_output_interval == 0:
					distribution_function = STATE.out_diags.actual_dist_functions[i].get_full_lm( STATE.species[i].F, n, time, STATE)
					f_slice_axis = STATE.out_diags.f_slice_axis
					if(distribution_function != None):
						if STATE.out_diags.initial_dist == None:
							STATE.out_diags.initial_dist = np.copy ( distribution_function )
							if( STATE.out_diags.F_output_log_scale):
								STATE.out_diags.initial_dist = np.log(STATE.out_diags.initial_dist)
							pass
						else:
							"""
							# for debugging, dump stupid slice of sittrubtion....
							clf()
							ax = subplot(111)
							slut=np.zeros( distribution_function.shape[1])
							slut1=np.zeros( distribution_function.shape[1])
							slut[:] = distribution_function[20,:]
							slut1[:] = distribution_function[40,:]
							slut -= initial_dist[20,:]
							slut1 -= initial_dist[40,:]
							
							
							clf()
							title("slices of F at 20 cheap at time=%f, step=%d" % (time, n) )
							plot(f_slice_axis, np.log( np.fabs(slut)), color='blue' )
							#plot(f_slice_axis, np.log(np.fabs(initial_dist[20,:])), color='black', alpha=.2)
							SAVEFIG("./%sf_cheap_1_cross_section_%06d.png" % (species_name_prefix, n))
							
							clf()
							title("slices of F at 40 cheap at time=%f, step=%d" % (time, n) )
							plot(f_slice_axis, np.log(np.fabs(slut1)), color='green' )
							#plot(f_slice_axis, np.log(np.fabs(initial_dist[40,:])), color='black', alpha=.2)
							SAVEFIG("./%sf_cheap_2_cross_section_%06d.png" % (species_name_prefix, n))
							
							distribution_function -= initial_dist
							
							clf()
							spot = distribution_function.shape[1] / 5
							plot( distribution_function[:, spot-5] )
							plot( distribution_function[:, spot+5] )
							title("slices of F along x at cells +-5 at time=%f, step=%d" % (time, n) )
							SAVEFIG("./%sf_long_cheap_cross_section_%06d.png" % (species_name_prefix, n))
							"""
							pass
						"""
						clf()
						ax = subplot(111)
						logme_dist_slice = False
						if logme_dist_slice:
							plot(f_slice_axis, np.log(np.fabs( distribution_function[20,:] ) ))
							plot(f_slice_axis, np.log(np.fabs( distribution_function[40,:] ) ))
							vth_intersection = np.max( np.log(np.fabs( distribution_function[20,:] ) ) ) - 1.0 
							ylim([-20,10])
						else:
							#max =  np.max( distribution_function[20,:] )
							#max1 =  np.max( distribution_function[40,:] )
							#if max1 > max:
							#	max = max1
							
							# normalize
							#plot(f_slice_axis, distribution_function[20,:]/max )
							#plot(f_slice_axis, distribution_function[40,:]/max1 )
							
							# dont normalized...
							plot(f_slice_axis, distribution_function[20,:])
							plot(f_slice_axis, distribution_function[40,:])

							#vth_intersection = 1.0 - np.exp(-1.0)
							vth_intersection = .5
							
							#ylim([0.0,1.1])
						try:
							vth = STATE.species[i].species_config.temp_profile.val
							
						except:
							vth = 0
						
						#axvline(x=vth)
						#if STATE.config.E_field_profile != None:
						vp = None; vth = None; vth1 = None
						if DaWaveDrivah != None and False:
							
							vp = DaWaveDrivah.v_phase
							vth = DaWaveDrivah.vth
							vth1 = DaWaveDrivah.vth * np.sqrt(1.5) #1.224744871	#TODO: sprt out this sqrt(1.5) and 1p vs 3p ness.
						elif STATE.config.E_field_profile != None:
							try:
								vth1 = STATE.species[i].species_config.temp_profile.val  * np.sqrt(1.5)
								wavenumber = STATE.config.E_field_profile.wavenumber
								wp_e = 1.0
								complex_freq = landau_damping_frequency(wp_e, vth1, wavenumber, maxwellian_convention_factor = 1.0 )
								vp = np.real( complex_freq ) / wavenumber
							
							except:
								pass
						
						if vp != None and vth1 != None:
							try:
								if mpi_info.rank == 0:
									print "V phase = %f with vth=%f anf k=%f" % ( vp, vth1, wavenumber )
								
								axvline(x=vp, color='black')
								
								axvline(x=vth1, color='red')
								#axhline(y=vth_intersection, color='purple')
								#axvline(x=vth1, color='purple')
								loc = ax.yaxis.get_ticklocs()[5]
								text(vp,loc,'$V_\\phi$',rotation=90)
								
								
								#xlim( (-0.2, 0.2) )
								#xlim( (f_slice_axis[0], f_slice_axis[-1]) )
								loc = ax.yaxis.get_ticklocs()[5]
								text(vth1,loc,'$V_{th}$',rotation=90, color='red',  weight='bold')
								
								#val = STATE.species[i].species_config.temp_profile.val * 1.224744871
								#axvline(x=vp)
								#loc = ax.yaxis.get_ticklocs()[2]
								#text(vp,loc,'$V_\\phi$',rotation=90)
									
								#except:
								#	print "Daig messup bail"
								#	exit(-1)
								#	pass
								#xlim([-1.0*STATE.species[i].F.momentum_axis.values[-1],  STATE.species[i].F.momentum_axis.values[-1]])
								
								xlim([-5.0*vth1, 5.0*vth1] )
							except:
								pass
						try:
							title("Dist function versus absolute value of Px\nstep %f, time %e" % (n, time))
							xlabel("Momentum along X")
							SAVEFIG("./%sf_cross_section_%06d.png" % (species_name_prefix, n) )
						except:
							pass
						"""
						
						# save the numpy file for entire distuibution function.					
						NP_SAVE("px_vs_x_%06d.npy" % n, distribution_function)

						# now just do any requested PNGs
						# do log scale, if requested.
						if( STATE.out_diags.F_output_log_scale):
							distribution_function_to_plot = np.log(distribution_function)
						else:
							distribution_function_to_plot = distribution_function

						# output full distribution function
						clf()
						imshow( ( distribution_function_to_plot) ,  aspect='auto', cmap='PRGn' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						colorbar()
						title("X versus Px\nStep %f, Time %e" % (n, time))
						ylabel("Space (X)")
						xlabel("Momentum (Px)")
						SAVEFIG("dist/%spx_vs_x_%06d.png" % (species_name_prefix, n) )

						# output any slices requested.
						if( STATE.out_diags.F_output_lineouts ):
							#for slice_cell in STATE.out_diags.F_output_lineouts:
							#	title("Distribution Function at x=%.3e (cell %d)" % (slice_cell*STATE.spatial_axes__global[0].dx(),slice_cell) )
							clf()
							title("Distribution Function at fixed x values")
							xlabel("Momentum (Px)")
							ylabel("Distribution Function Value")
							for slice_cell in STATE.out_diags.F_output_lineouts:
								plot( f_slice_axis, distribution_function_to_plot[ slice_cell] )
							if( STATE.out_diags.F_output_velocity_interval ):
								xlim(STATE.out_diags.F_output_velocity_interval)
							if( STATE.out_diags.F_output_range ):
								ylim( STATE.out_diags.F_output_range )
							SAVEFIG("dist/%spx_slices_%06d.png" % (species_name_prefix, n) )

							# do maxwell Botlzmann
							clf()
							title("Distribution Function of V2 (at fixed x values")
							xlabel("Momentum (Px)")
							ylabel("Distribution Function Value")
							for slice_cell in STATE.out_diags.F_output_lineouts:
								slice_values = distribution_function_to_plot[ slice_cell]
								abs_p = slice_values[0:slice_values.shape[0]/2][::-1] + slice_values[slice_values.shape[0]/2:]
								abs_p *= STATE.species[i].F.momentum_axis.values[:]
								abs_p *= STATE.species[i].F.momentum_axis.values[:]
								#values_to_plot = scipy.integrate.simps( abs_p, dx = STATE.species[i].F.momentum_axis.dp() , axis=0)
								plot( f_slice_axis[f_slice_axis.shape[0]/2:], abs_p  )
							if( STATE.out_diags.F_output_velocity_interval ):
								xlim(STATE.out_diags.F_output_velocity_interval)
							if( STATE.out_diags.F_output_range ):
								ylim( STATE.out_diags.F_output_range )
							SAVEFIG("dist/%spx_slices_maxwell_boltzmann_%06d.png" % (species_name_prefix, n) )
							
							#clf()
							#stem(  STATE.species[0].F.cells[10].harmonics[0].momentum_projection[:] )
							#print STATE.species[0].F.cells[10].harmonics[0].momentum_projection[:]
							#SAVEFIG("dist/%szzzz_%06d.png" % (species_name_prefix, n) )


						# output full delta distribution function
						#distribution_function_to_plot -= STATE.out_diags.initial_dist
						#clf()
						#imshow( ( distribution_function_to_plot) ,  aspect='auto', cmap='PRGn' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						#colorbar()
						#title("Change from Inital: X versus Px\nStep %f, Time %e" % (n, time))
						#SAVEFIG("dist/%spx_vs_x_delta%06d.png" % (species_name_prefix, n) )

						"""
						clf()
						title("Distribution Function Slice (integrated along x direction")
						ylabel("F");xlabel("Momentum allong X");
						integrated_slice = scipy.integrate.simps( distribution_function, dx = STATE.spatial_axes__global[0].dx() , axis=0)
						plot( f_slice_axis, integrated_slice)
						if( STATE.out_diags.F_output_velocity_interval ):
							xlim(STATE.out_diags.F_output_velocity_interval)
						if( STATE.out_diags.F_output_range ):
							ylim( STATE.out_diags.F_output_range )
						SAVEFIG("dist/%sp1x1_slice_int%06d.png" % (species_name_prefix, n) )

						clf()
						title("Distribution Function Slice")
						ylabel("F");xlabel("Momentum allong X");	
						plot( f_slice_axis, distribution_function[ 20] )
						plot( f_slice_axis, distribution_function[ 40] )
						plot( f_slice_axis, distribution_function[ STATE.spatial_axes__global[0].num_points/2 ] )
						if( STATE.out_diags.F_output_velocity_interval ):
							xlim(STATE.out_diags.F_output_velocity_interval)
						if( STATE.out_diags.F_output_range ):
							ylim( STATE.out_diags.F_output_range )
						SAVEFIG("dist/%sp1x1_slice_%06d.png" % (species_name_prefix, n) )
						
						clf()
						# cmap='PRGn'
						#imshow(np.log(np.fabs(distribution_function)), aspect='auto', cmap="copper", extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						#imshow(np.log(np.fabs(distribution_function)), cmap='Paired', aspect='auto', extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						imshow( np.log( distribution_function) ,  aspect='auto', cmap='Paired' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						colorbar()
						#clim( (-4.5e-4, 4.5e4) )
						#xlim( (50,135) )
					
						#X, Y = np.meshgrid(f_slice_axis, np.arange(0,distribution_function.shape[0]) )
						#contourf( X,Y,distribution_function, 40, cmap=cm.Paired)
						
						# xlim( (-0.075, 0.075) )
						# xlim( (-0.5, 0.5) )
						# xlim( (f_slice_axis[0], f_slice_axis[-1]) )
						
						title("X versus Px\nStep %f, Time %e" % (n, time))
						xlim(STATE.out_diags.F_output_range)
						SAVEFIG("dist/%slog_p1x1_%06d.png" % (species_name_prefix, n) )
						xlim([f_slice_axis[0],f_slice_axis[-1]] )
						SAVEFIG("dist/%slogall_p1x1_%06d.png" % (species_name_prefix, n) )
						NP_SAVE("p1x1_%06d.npy" % n, distribution_function)
						clf()
						imshow( ( distribution_function) ,  aspect='auto', cmap='Paired' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						colorbar()
						title("X versus Px\nStep %f, Time %e" % (n, time))
						xlim( STATE.out_diags.F_output_range )
						SAVEFIG("dist/%sp1x1_%06d.png" % (species_name_prefix, n) )
						
						clf()
						distribution_function -= STATE.out_diags.initial_dist
						imshow( distribution_function,  aspect='auto', cmap='PRGn' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						colorbar()
						title("F1 (change from maxwellian) X versus Px\nStep %f, Time %e" % (n, time))
						SAVEFIG("dist/%sp1x1_f1_%06d.png" % (species_name_prefix, n) )
						NP_SAVE("p1x1_f1_%06d.npy" % n, distribution_function)
						
						#xlim( (-2.0, 2.0) )
						#SAVEFIG("./%sp1x1_zoom_%06d.png" % (species_name_prefix, n) )
						"""

						"""
						clf()
						subplot(211)
						labels = []
						stop = STATE.species[i].F.num_l
						if stop > 6:
							stop = 6
						for l in xrange(0, stop):
							m_stop = l+1
							if m_stop > STATE.species[i].F.num_m:
								m_stop = STATE.species[i].F.num_m
							for m in xrange(0,m_stop):
								labels.append( "(%d,%d)"% (l, m) )
						
						xlocations = np.array( range(len(labels))) + 0.1; width = 0.1;
						bar( xlocations, np.log(actual_dist_functions[i].diag__harmonic_content[:len(xlocations)]), width=width )
						xticks( xlocations + width/2.0, labels, size="xx-small")
						xlim(0, xlocations[-1]+width*2.0)
						ylim([-10.0, 10.0])
						gca().get_yaxis().tick_left()
						#SAVEFIG("./harmonics_%06d.png" % n)
						subplot(212)
						bar( xlocations, np.log(actual_dist_functions[i].diag__harmonic_content_minus[:len(xlocations)]), width=width )
						xticks( xlocations + width/2.0, labels,  size="xx-small")
						xlim(0, xlocations[-1]+width*2.0)
						ylim([-10.0, 10.0])
						gca().get_yaxis().tick_left()
						SAVEFIG("./harmonics_%06d.png" % n)
						"""
				
				if STATE.out_diags.SH_diag_interval > 0 and  n%STATE.out_diags.SH_diag_interval == 0:
					# package the data on all nodes...
					if tt_crap_vuf == None:
						tt_crap_vuf = np.zeros( STATE.species[i].F.num_total_cells )
					hhh_numl =   STATE.species[i].F.num_l - 1
					for gggg in xrange(  STATE.species[i].F.num_total_cells):
						hhh_temp_last =  STATE.species[i].F.distribution_function_data[gggg,hhh_numl,:]
						hhh_temp_all =  np.sum( STATE.species[i].F.distribution_function_data[gggg,:,:] )
						tt_crap_vuf[ gggg] = np.sum( hhh_temp_last / hhh_temp_all )
					
					# only the root node will execute this code block.
					if mpi_gather.gather( tt_crap_vuf, remove_guard_cells=True):
						#if simple_hhh_buffer == None:
						#	simple_hhh_buffer = np.zeros( [100,mpi_gather.stiching_buffer.shape[0]] )
						#	simple_hhh_buffer__index = 0
						#simple_hhh_buffer[ simple_hhh_buffer__index, : ] = mpi_gather.stiching_buffer[:]
						#simple_hhh_buffer__index += 1
						#if simple_hhh_buffer__index == simple_hhh_buffer.shape[0]:
						#	NP_SAVE("hhh_%06d.npy" % (n, ) , simple_hhh_buffer)
						#	simple_hhh_buffer__index = 0
						NP_SAVE("hhh_%06d.npy" % (n, ) , tt_crap_vuf)
						if True: 
							clf()
							#subplot(211)
							plot(mpi_gather.stiching_buffer, marker="*", ms=.8)
							axhline( y = 0)
							title("Last L field at step %f, time %e" % (n, time));
							SAVEFIG("./hhh_%06d.png" % (n, ) )
						pass
				
				if i == 0:
					# This is the 'raw' E-field data from the simulation... so we need to remove the guard cells
					
					#if mpi_gather.gather( STATE.E_field__external.components[0], remove_guard_cells=True):
					if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
						if self.simple_E_buffer == None:
							simple_E_buffer = np.zeros( [100,mpi_gather.stiching_buffer.shape[0]] )
							simple_E_buffer__index = 0
						simple_E_buffer[simple_E_buffer__index, : ] = mpi_gather.stiching_buffer[:]
						simple_E_buffer__index  += 1
						if simple_E_buffer__index == simple_E_buffer.shape[0]:
							NP_SAVE("E_%06d.npy" % (n, ) , simple_E_buffer)
							simple_E_buffer__index = 0
							#print "SAVED 500 E-feild entries!!!!"
						
						if n > 260 and n < 310:
							ass = True
						else:
							ass = False
						ass = False
						
				if STATE.out_diags.E_interval > 0 and n%STATE.out_diags.E_interval == 0: # and time > 240.0 and time < 300.0:
					if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
						clf()
						subplot(211)
						plot(mpi_gather.stiching_buffer, marker="*", ms=.8)
						#ylim( [-0.0015, 0.0015] )
						#ylim( [-5e-5, 5e-5] )
						axhline( y = 0)
						title("E vs X (step %f, time %e)" % (n, time));
						ylabel("Electric Field")
						subplot(212)
						size = np.arange(0, 2.0, 2.0/mpi_gather.stiching_buffer.shape[0])
						size -= 1.0
						#fft_data = np.absolute(fft.fftshift(fft.fft(mpi_gather.stiching_buffer)))
						fft_data = np.absolute((fft.fftshift( fft.fft(mpi_gather.stiching_buffer ))))   		#, n=256)))
						
						plot(np.log10(fft_data))
						title("Spectral Density (FFT) of E field")
						ylabel("Log spctral Density")

						#ylim( [0.0, .04] )
						#ylim( [0.0, .005] )
						SAVEFIG("E_%06d.png" % (n, ) )
						NP_SAVE("E_%06d.npy" % (n, ) , mpi_gather.stiching_buffer)

						#print "Signal to noise: %e" % scipy.stats.signaltonoise( mpi_gather.stiching_buffer)
						#print "length of fft buffer %d" % fft_data.shape[0]
				
				
				
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if True and STATE.out_diags.density_interval>0 and n%STATE.out_diags.density_interval == 0 and mpi_gather.gather(  moments[i].get_density(), remove_guard_cells=True ):
					try:
						"""
						print "-------------"
						print (moments[i].get_density()).shape
						print (moments[i].get_density())
						print "---------------"
						"""
						
						if initial_density == None:
							initial_density = np.zeros ( mpi_gather.stiching_buffer.shape[0] )
							initial_density[:] = mpi_gather.stiching_buffer[:]
							initial_density[:] = F_from_External_E_Solver.init_f0_density[:] #mpi_gather.stiching_buffer[:]
						else:
							mpi_gather.stiching_buffer -= initial_density
							pass
					except:
						
						initial_density = None
						pass
					
					clf(); 
					#subplot(211)
					plot(mpi_gather.stiching_buffer); title("Density %.02f, time %.02e" % (n, time));
					xlabel('Space (x)')
					ylabel('Normalized Density')
					#draw(); 
					#subplot(212)
					#size = np.arange(0, 2.0, 2.0/mpi_gather.stiching_buffer.shape[0]); size -= 1.0
					#fft_data = np.absolute(fft.fftshift(fft.fft(mpi_gather.stiching_buffer)))
					#fft_data += 1e-30
					#plot(size, np.log10(fft_data))
					SAVEFIG("%sdensity_%06d.png" % (species_name_prefix, n) )
					NP_SAVE("%sdensity_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
					
					#print fft_data
					print "\t Total charge %e" % (np.sum(moments[i].get_density()), )
				
				
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if STATE.out_diags.temperature_interval > 0 and n%STATE.out_diags.temperature_interval==0 and mpi_gather.gather(  moments[i].get_temp_ev(), remove_guard_cells = True):
					tempp = moments[i].get_temp_ev()					
					clf(); plot( mpi_gather.stiching_buffer); 
					title("Temperture (ev) at step %f, time %e" % (n, time)); draw();
					SAVEFIG("%stemp_ev_%06d.png" % (species_name_prefix, n) )
					NP_SAVE("%stemp_ev_%06d.npy" % (species_name_prefix,n ), mpi_gather.stiching_buffer)
				
				
				if STATE.out_diags.v2_interval > 0 and n%STATE.out_diags.v2_interval == 0 and mpi_gather.gather(  moments[i].get_v2(), remove_guard_cells=True ):
				# the daignostics emmited from Moments already have the gaurd cells trimmed
					tempp = moments[i].get_v2()
					clf(); plot(mpi_gather.stiching_buffer); title("$v^2$ verus $x$ at step %f, time %e" % (n, time));
					#ylim( (-1e-8, 1e-8) ); 
					axhline(y=0, color='black'); #lines_temp_qx[0].set_ydata( moments[i].get_qx())
					SAVEFIG("%sv2_%06d.png" % (species_name_prefix, n) )
					NP_SAVE("%sv2_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
				
				if STATE.out_diags.Q_interval > 0 and n%STATE.out_diags.Q_interval == 0 and mpi_gather.gather(  moments[i].get_qx(), remove_guard_cells=True ):
				# the daignostics emmited from Moments already have the gaurd cells trimmed
					tempp = moments[i].get_qx()
					clf(); plot(mpi_gather.stiching_buffer); title("Qx at step %f, time %e" % (n, time));
					ylim( (-1e-8, 1e-8) ); axhline(y=0, color='black'); #lines_temp_qx[0].set_ydata( moments[i].get_qx())
					SAVEFIG("%sqx_%06d.png" % (species_name_prefix, n) )
					NP_SAVE("%sqx_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
				
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if STATE.out_diags.heat_conduction_interval > 0 and n%STATE.out_diags.heat_conduction_interval ==0 and mpi_gather.gather( moments[i].get_heat_conductivity_x(), remove_guard_cells=True ):
					clf(); plot(mpi_gather.stiching_buffer ); k_ave =  mpi_gather.stiching_buffer.mean()
					title("Heat Conduction at step %d, time %0.2e (ave: %0.3e)" % (n, time, k_ave));
					axhline(y=0, color='black'); #lines_temp_kx[0].set_ydata(moments[i].get_heat_conductivity_x())
					ylim( (-4, +4) ); draw();
					SAVEFIG("%skx_%06d.png" % (species_name_prefix, n) )
					NP_SAVE("%skx_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
				
				
				#output_f( STATE.species[i].F, n , time ,plot=plotme1.plot_ps)
				#if mpi_info.rank == 0:
				#PRINT("Output written for n=%d" % n)
				#exit(-1)
			pass
		pass