import numpy as np
import numpy.polynomial.legendre as legendre
import math
import matplotlib
matplotlib.use('Agg')
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.linalg
import os
import scipy.integrate

#import timeit

# for mayavi....
#from mayavi import mlab


import time
from scipy.stats import norm


import threading
#import msvcrt
from multiprocessing import Pipe
import PyQt4.QtCore as QtCore
from PyQt4 import QtGui

from oshun_moments import *
from oshun_mpi import *
from oshun_c import *
import scipy.fftpack as fft
import scipy.special
import scipy.integrate

import h5py
import inspect
import cPickle as pickle

import scipy.special
import scipy.optimize
import scipy.stats

def plasma_dispersion( value):
	a = scipy.special.wofz(value)
	a *= np.sqrt(np.pi)*1j
	return a
def plasma_dispersion_prime( value):
	z = plasma_dispersion( value )
	return -2.0*(1.0+value*z)

def landau_damping_frequency(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0, root_guess = 1.01 ):
	#
	# we are gonna solve for the roots of the electrostatic warm plasma dispersion relation
	#	e.g. 	solve:  	eplison(w,k) = 0
	#			where:		eplison(w,k) = 1 + sum_of_terms_for_all_species_of_form( chi_s(w,k) ) 
	#			here:				( for example, chi_e = electron susceptibility or chi_i = ion susceptibility etc....)
	#			now:						chi_s(x,s) = - ( wp_s/(2*vth_s*k_wave) )^2 * plasma_dispersion_prime( w_wave / (sqrt(2)*vth_s*k_wave) )
	#												for any species s.
	#														wp_s = plasma frew for species s
	#														vth_s = thermal speed f species s
	#														k_wave = the wavenumber of the wave traveling in the plasma.
	#														w_wave = the frequency of the wave traveling in the plasma.
	#
	#		We will use this function to, given a wave's wavenumber as input, solve the dipersion realtion to find the (complex) frequency (i.e. pole
	#			of the dispersion equation)
	#
	#		The default guess_root is 1.01 just because we expect that the roots we want have a real frequency that is always > 1 (you can see this from
	#			the Bohm-Gross relation, for example...) We'll see how good it works to see if we need to be any smarter...
	#
	# example:
	#					wp_e = 1.0						
	#					vth_e = sqrt( .1 )
	#					wavenumber = 1.0				# setup so that vth_e*wavenumber = sqrt(.1)  
	#														#		any other combinations will work.. only the product vth_e*wavenumber matters
	#														#		Note: vth_e*wavenumber is often refered to as k*lambda_debye... this is an apporximation
	#														#			that assumes that the wave frequency is very close in value to wp_e (w_wave ~= wp_e)
	#														#			A rule of thumb is that any value of k*lambda_debye greater then .33 or so is pushing it.. 
	#														#			so dont't go much higher.. [also notice that's not an appoximation when in the coefficent chi..
	#														#			that reduces to 1/ k*lambda_debye... but it is an appoximation inside the eponential. ]
	#
	#														# for vth_e*wavenumber = sqrt(.1), wp_e=1.0, the verifed answer is 1.179161-0.0184479j
	#
	#					print landau_damping_frequency( wp_e, vth_e, wavenumber)
	#														
	chi_e = np.power( ( wp_e/(vth__e*wavenumber) ), 2.0) / maxwellian_convention_factor
	def plasma_epsilon1( x ):
		val = 1.0 - chi_e*plasma_dispersion_prime( x )
		return val
	epsilon_root = scipy.optimize.newton( plasma_epsilon1, root_guess )
	return epsilon_root * wavenumber * vth__e * np.sqrt(maxwellian_convention_factor)

#from profiling import *

#matplotlib.use('GTKAgg')
pipe1main = None
pipe1upd = None
mypass = None
t0 = None
runthread = None

class SimMetadata:
	def __init(self):
		self.header = "Oshun1d metadata"
		self.dirty = True
	def set_dirty(self):
		self.dirty = True
	def is_dirty(self):
		return self.dirty
	def save(self):
		if mpi_info.rank == 0:
			pickle.dump( self, open( "oshun_metadata.p", "wb" ) )
	
def DEBUG_E(E, bail=False, indent=False):
	if mpi_info.rank == 0:
		limit = 50
		if E.num_plain_cells[0] < limit:
			limit = E.num_plain_cells[0]
		num_boundry_cells = E.num_boundry_cells
		
		for ix in xrange(20, 30):
			#print "%d: %e" %( ix, E.components[0][num_boundry_cells[0][0] + ix])
			if indent:
				print "\t\t",
			print "%d: %.8e" %( ix,  E.components[0][ix])
		if bail:
			exit(-1)

def DEBUG_BAIL():
	if mpi_info.rank == 0:
		exit(-1)
	pass
def DEBUG(YYY=None, L=0, M=0, bail=False, indent=False):
	if mpi_info.rank == 0:
		temp_species = STATE.species[0]
		if YYY == None:
			print "Setting DEBUG Y to defualt for species"
			YYY = temp_species.F
		
		index_ml = YYY.linear_index_for_harmonic(L,M)
		num_boundry_cells = YYY.num_boundry_cells
		#print "INDEX OF HARM:" + str(index_ml)
		x_points = [0,1,2,3,4,5,6]
		#print id(YYY)
		if indent:
			print "\t\t",
		
		for ix in x_points:
			print "\t%.9e" % STATE.spatial_axes[0].values[ix],
		print ""
		if indent:
			print "\t\t",	
		print "--------------------------------------------------------"
		
		#print "beam temp %f" % beam_temp
		#print "EXP arg factor is: %f" % exp_arg_coeff
		#print "FROT factor is: %f" % factor_coeff
		#for ip in xrange(0, temp_species.momentum_axis.num_points):
		
		for ip in xrange(20, 30):
			if indent:
				print "\t\t",
			
			print "%d: "%( ip,) ,
			for ix in x_points:
				#print "\t%e"%( YYY.cells[num_boundry_cells[0][0]+ix].harmonics[index_ml].momentum_projection[ip] ),
				if( ix == num_boundry_cells[0][0]):
					print "",
				print "\t%.8e"%( YYY.cells[ix].harmonics[index_ml].momentum_projection[ip] ),
			print ""
		if bail:
			exit(-1)

def DEBUG_HCV(HCV, L=0, M=0, bail=False, indent=False):
	if mpi_info.rank == 0:
		x_points = [0,1,2,3,4,5,6]
		if indent:
			print "\t\t",
		for ix in x_points:
			print "\t%.8e" % STATE.spatial_axes[0].values[ix],
		print ""
		temp_species = STATE.species[0]
		nump = temp_species.momentum_axis.num_points
		num_boundry_cells = temp_species.F.num_boundry_cells
		#for ip in xrange(0, nump):
		for ip in xrange(20, 30):
			if indent:
				print "\t\t",
			
			print "%d: "%( ip,) ,
			for ix in x_points:
				if( ix == num_boundry_cells[0][0]):
					print "",
				#print "\t%e" % ( HCV.data[num_boundry_cells[0][0]+ix,ip],) ,
				print "\t%.8e" % ( HCV.data[ix,ip],) ,
			print ""
		if bail:
			exit(-1)

def dump_F_at_xlm(F, x_points, l, m, title=None):
	try:
		temp = len(x_points)
	except:
		x_points = [x_points]
	
	num_boundry_cells = F.num_boundry_cells
	sh_index = F.linear_index_for_harmonic(l,m)
	
	print ""
	if title != None:
		print title
		print ""
	for ip in xrange(0, F.momentum_axis.num_points):
		print "dump_F_at_xlm: %d: "%( ip,) ,
		for ix in x_points:
			print "\t%f"%( F.cells[num_boundry_cells[0][0]+ix].harmonics[sh_index].momentum_projection[ip] ),
		print ""
	pass	
	
def heaviside_step(x):
	if x< 0:
		return 0.0
	return 1.0;

def timeit_basic(f):

    def timed(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()

        print 'func:%r args:[%r, %r] took: %2.4f sec' % \
          (f.__name__, args, kw, te-ts)
        return result

    return timed
	
class timings_info:
	
	def __init__(self, fragment_name):
		self.name = fragment_name
		self._num_times_run = 0
		self._total_time = 0.0
		self._min_time = 1e12
	def report_timing(self, time_taken):
		self._num_times_run += 1
		self._total_time += time_taken
		if time_taken < self._min_time:
			self._min_time = time_taken
	def min_time(self):
		return self._min_time
	def average_time(self):
		return self._total_time / float(self._num_times_run)
	def total_time(self):
		return self._total_time
	
	
class timings:
	def __init__(self):
		self.fragments = {}
		self.order_list = []
	def add_timing(self, func, time_taken):
		class_of_func =  self.get_class_that_defined_method(func)
		if class_of_func == None:
			name = func.__name__
		else:
			name = "%s->%s" % (func.__name__, class_of_func.__name__)
		
		if name in self.fragments:
			timing_info = self.fragments[name]
		else:
			timing_info = timings_info(name)
			self.order_list.append( name)
			self.fragments[name] = timing_info
		timing_info.report_timing(time_taken)
	def add_timing_by_name( self, name, time_taken):
		if name in self.fragments:
			timing_info = self.fragments[name]
		else:
			timing_info = timings_info(name)
			self.order_list.append( name)
			self.fragments[name] = timing_info
		timing_info.report_timing(time_taken)
		
	
	def get_class_that_defined_method1(self, func):
		# see if it is a method call of a class..
		if func.__self__:
			classes = [func.__self__.__class__]
		else:
			#in this case, it is an unbound function...
			classes = [method.im_class]
		
		while classes:
			c = classes.pop()
			if func.__name__ in c.__dict__:
				return c
			else:
				classes = list(c.__bases__) + classes
		return None	
	def get_class_that_defined_method(self, meth):
		for cls in inspect.getmro(meth.im_class):
			if meth.__name__ in cls.__dict__: return cls
		return None
	
	def print_results(self):
		for name in self.order_list:
			entry = self.fragments[name]
			print "------------------\n\t\t%s\n\t\t\t: min: %e " % (name, entry.min_time() )
		print ""
	

TIMINGS = timings()

# function decorator version of the timing...
def timeit(f, id=None):
	
	def time_clousure(*args, **kw):
		_id = id
		if _id == None:
			stack_frames = inspect.getouterframes(inspect.currentframe())
			try:
				calling_line_index = stack_frames[1][5]
				calling_line_source_code = stack_frames[1][4][calling_line_index]
				_id = calling_line_source_code.strip()
			finally:
				del stack_frames
		
		start_time = time.time()
		result = f(*args, **kw)
		end_time   = time.time()
		
		#TIMINGS.add_timing(f, end_time-start_time)
		TIMINGS.add_timing_by_name(_id, end_time-start_time)
		
		return result
		
	return time_clousure

def launch_plots():
	global pipe1main, pipe1upd	
	pipe1main, pipe1upd = Pipe()
	
def pause_sim(interval = .1):
	while( KEY.pause_sim ):
		show(block=False)
		draw()
		#time.sleep(interval)
	return

class GuiThread(threading.Thread):
	def __init__(self):
		super(GuiThread, self).__init__()
	def run(self):
		print("Starting GUI Loop")
		show()
		print("ending GUI Loop")
"""	
class KeyboardInterrupt(threading.Thread):
	def __init__(self):
		super(KeyboardInterrupt, self).__init__()
		self._stop = False
		self.quit_sim = False
		self.pause_sim = False
	def run(self):
		#set initial 'ch', alternatively check for ch != None in the while loop
		ch = '1'	
 
		#loop until 
		while (not self._stop):
			#msvcrt.getch is blocking so wrap it in a hit test loop
			while (msvcrt.kbhit()):
				ch = msvcrt.getch()
				if (ch in ['q']):
					print "Keyboard interrupt triggered, exiting program"
					self.quit_sim = True
				if(ch in [' ']):
					if self.pause_sim == False:
						self.pause_sim = True
					else:
						self.pause_sim = False
			
			time.sleep(0.030)
 
	def stop (self):
		self._stop = True
"""
#matplotlib.interactive(True)
# standard central difference....
def central_difference_derivitive(vals, dx):
	result = np.ones( len(vals) )
	for i in range(1,len(vals)-1):
		result[i] = vals[i+1] - vals[i-1]
		result[i] /= 2.0*dx
	return result

# Micheal's central difference from the C++ code:
# 		is -1 times the normal one
#		does not do the division by "2h"
#		defines values that the normal central differnece does not:
#			element at 0: set to be the same as element 1
#			element at last position: unchanged.
def test_central_difference_derivitive_micheal(vals, dx):
	for i in range(0, len(vals)-2):
		vals[i] -= vals[i+2]
	for i in range(len(vals)-3, -1, -1):
		vals[i+1] = vals[i]
	print vals

# convert standrd central differnce results to the ones from Micheal just
#		as a test that I is doing it right.
def test_central_difference_boring_against_micheal(vals, dx):
	results = central_difference_derivitive(vals, dx)
	results *= (2*dx)
	results *= -1.0
	results[0] = results[1]
	results[len(vals)-1] = vals[len(vals)-1]
	return results

class WaveDrivah:
	def __init__(self, dt, species, num_wavelengths, num_cells_wave_travels_in_dt = 0):
		pass
		# ok see dmaping, but not lots...
		#self.ending_amp = .001
		#self.amp_rate = 2.0
		#self.starting_amp = .0001
		"""
		self.ending_amp = 1e-5
		self.amp_rate = 2.0
		self.starting_amp =  1e-6 # 1e-4   #.0001
		self.rise_time = 2.0
		self.sustain= 0.0
		"""
		
		# 1e-16 to 1e-10 over 2.0   num=25 was not spiked.
		# ok some non-lin 1e-14 to 1e-9 over 10.0 num 28
		#nice 1e-24 to 1e-20 over 10 num 28
		
		self.echo_2nd_pulse_time = 200.0 #200.0
		self.driver_on = False 
		
		""" 
		# GOOD WORKING 192 x cells, 150 p
		self.ultra_pure_harmonic_mode = False
		self.ending_amp = 1e-30
		self.final_amp = 1e-30
		self.amp_rate = 1.0
		self.starting_amp =  1e-36 # 1e-4   #.0001
		self.rise_time = 10.0
		self.sustain= 0.0
		self.fall_time = 0.0
		
		self.wp = 1.0   # WATCH OUT FOR THIS   				WARNING WARNING!!!!!!!!!!!!!!!
		"""
		
		self.ultra_pure_harmonic_mode = False
		self.ending_amp = 1e-19
		self.final_amp = 1e-19
		self.amp_rate = 1.0
		self.starting_amp =  1e-21 #1e-40 # 1e-4   #.0001
		self.rise_time = 2.5
		self.sustain= 0.0
		self.fall_time = 0
		self.wp = 1.0   # WATCH OUT FOR THIS   				WARNING WARNING!!!!!!!!!!!!!!!		
		
		self.current_amp = self.starting_amp
		self.num_wavelengths = 21.0 #16 # now try 32 /4 #30 # num_wavelengths
		#18 is .4999
		#16 is .44
		#10 is .2777
		# 9 is .2499
		spatial_axis_global = STATE.spatial_axes__global[0]
		
		dx = spatial_axis_global.dx()
		
		box_size = (spatial_axis_global.values[-1] - spatial_axis_global.values[0] )
		self.num_cells_in_box = spatial_axis_global.values.shape[0]
		self.wavelength = box_size/ (float( self.num_wavelengths + dx))
		self.wavelength = box_size/ (float( self.num_wavelengths))
		#self.wavenumber = 2.0*np.pi/(self.wavelength + dx) +
		self.wavenumber = 2.0*np.pi/( self.wavelength )
		self.discrete_factor = np.pi*2.0*float( self.num_wavelengths)/spatial_axis_global.values.shape[0]
		print self.discrete_factor
		
		
		self.species = species
		self.vth = species.species_config.temp_profile.val
		self.buffah = STATE.E_field.clone()
		self.buffah_accumulate = STATE.E_field.clone()
		
		self.e_profiles = []
		
		for labmdaa in xrange(int(self.num_wavelengths),int(self.num_wavelengths)+1): #(25, 26):
			self.e_profiles.append( EFieldSinusodialProfile(self.starting_amp, num_wavelengths=labmdaa) )
			#self.e_profiles[-1].apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah )
		
		
		#self.e_profile = EFieldSinusodialProfile(self.starting_amp, num_wavelengths=self.num_wavelengths)
		#self.e_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah )
		
		
		#self.wavenumber = self.e_profile.wavenumber
		self.w_wave__complex = landau_damping_frequnecy__complex = landau_damping_frequency(self.wp, self.vth, self.wavenumber)
		self.w_wave = np.real( self.w_wave__complex )
		#self.w_wave = np.sqrt( 1 + 3*self.wavenumber*self.wavenumber*self.vth*self.vth)
		self.v_phase =  (self.w_wave/ self.wavenumber)
		# interperate the rise_time as the number of driver wave cycles....
		self.rise_time *= (np.pi * 2.0 / self.w_wave)
		

		self.last_seen_multiple = 0
		self.update_interval = 2*np.pi/self.w_wave
		self.update_interval_steps = int(2*np.pi/self.w_wave/dt)
		self.driver_off = False
		vth_in_ev = self.normalized_vth_to_ev( self.vth )
		
		
		# if num_cells_wave_travels_in_dt > 0, then we are attempting to enforce a time step
		#		that will ensure pure harmonic content in the driver.
		self.phase_step_per_cell = 2.0*np.pi/ float(spatial_axis_global.values.shape[0])
		self.phase_step_per_dt = 0.0
		self.pure_harmonic_dt = -1.0
		self.num_cells_wave_travels_in_dt = num_cells_wave_travels_in_dt
		
		# HACK TAKE OUT....
		num_cells_wave_travels_in_dt = 0
		if self.w_wave > 0.0:
			self.pure_harmonic_dt =  (2.0*np.pi/ 50.0) * self.w_wave
			dt = self.pure_harmonic_dt
		else:
			self.pure_harmonic_dt =  (2.0*np.pi/ 50.0) * self.w_wave
			dt = self.pure_harmonic_dt
			
		#static wave
		self.w_wave = 0
		
		if num_cells_wave_travels_in_dt > 0:
			self.phase_step_per_dt = self.phase_step_per_cell * float(self.num_cells_wave_travels_in_dt)
			self.pure_harmonic_dt = self.phase_step_per_dt / self.w_wave
			print  spatial_axis_global.values.shape[0]
			print self.pure_harmonic_dt
			print self.w_wave
			print num_cells_wave_travels_in_dt
			
		
		if num_cells_wave_travels_in_dt > 0:
			self.amp_rate = (self.ending_amp - self.starting_amp) / (self.rise_time / self.pure_harmonic_dt )
			self.amp_fall_rate = (self.final_amp - self.ending_amp) / (self.fall_time / self.pure_harmonic_dt )
		else:
			self.amp_rate = (self.ending_amp - self.starting_amp) / (self.rise_time / dt )
			try:
				self.amp_fall_rate = 0.0
				#self.amp_fall_rate = (self.final_amp - self.ending_amp) / (self.fall_time / dt )			
			except:
				self.amp_fall_rate = 0.0
		
		f = open('!sim_info.txt','a')
		f.write("\n\n")
		f.write("Thermal velocity (%s) (ev): %f \r\n" % (self.species.species_config.name, vth_in_ev) )
		f.write( "expect wavelength (imposed by driver such that %d spatial periods fit in the box of size %e):  %f \r\n" % (self.num_wavelengths, box_size, self.wavelength) )
		f.write( "assigned wavenumber (2*pi/( wavelength)):  %f \r\n" % (self.wavenumber,) )
		f.write( "freq (from Bohm-gross dispersion, t=%f (normalized) wavenumber %f, wp~=1 (normalized) ): %e \r\n" % (self.vth, self.wavenumber, self.w_wave) )
		f.write( "vphase ((bohm gross w)/(wavenumber)): %e  \r\n" % self.v_phase)
		f.write( "l_debye ( (thermal velocity (normalized)) / (wp (normalized)) ): %f c/wp \r\n" % (self.vth / self.wp) )
		f.write( "k*vth: %f \r\n" % ( (self.wavenumber*self.vth) ) )
		#f.write( "T_bounce: %f \r\n" % ( (  ) ) )
		f.write("\n\n")
		f.close()
		
		self.start_of_2nd_pulse = True
		STATE.metadata.driver = {}
		STATE.metadata.driver['driver_stop_time'] = self.sustain + self.rise_time
			
		STATE.metadata.driver['starting_amp'] =  self.starting_amp
		STATE.metadata.driver['ending_amp_after_rise'] =  self.ending_amp
		STATE.metadata.driver['amp_rate'] = self.amp_rate
		STATE.metadata.driver['final_amp_after_fall'] =  self.final_amp
		STATE.metadata.driver['rise_time'] = self.rise_time
		STATE.metadata.driver['sustain']  = self.sustain
		STATE.metadata.driver['fall_time'] = self.fall_time
		
		
		
		STATE.metadata.driver['vth__ev'] = vth_in_ev
		STATE.metadata.driver['num_wavelengths'] = self.num_wavelengths
		STATE.metadata.driver['wavelength'] = self.wavelength
		STATE.metadata.driver['wavenumber'] = self.wavenumber
		STATE.metadata.driver['w_wave'] = self.w_wave
		STATE.metadata.driver['v_phase'] = self.v_phase
		STATE.metadata.driver['lambda_debye'] = (self.vth / self.wp)
		STATE.metadata.driver['k__vth'] = (self.wavenumber*self.vth)
		
		STATE.metadata.driver['num_cells_wave_travels_in_dt'] = self.num_cells_wave_travels_in_dt
		STATE.metadata.driver['phase_step_per_dt'] = self.phase_step_per_dt
		STATE.metadata.driver['pure_harmonic_dt'] = self.pure_harmonic_dt
		STATE.metadata.set_dirty()
		
		
		STATE.use_external_E_field = True
		STATE.E_field__external = STATE.E_field.clone()
		STATE.E_field__external.scalar_mult(0.0)
		
		
	def normalized_vth_to_ev(self, normalized_vth):
		vth = normalized_vth 
		ev1 = vth/np.sqrt(vth*vth +1)*3e10/4.19e7
		return ev1*ev1*self.species.species_config.mass
	
	def update_boundry(self, time, n, dt):
		
		return
		
		if mpi_info.rank > 0:
			return
		
		if self.current_amp > ( self.ending_amp + self.sustain):
			return
		
		
		direct_driving_region_start = 0
		direct_driving_region_end = 1
		
		#self.current_amp = .001
		# in the full box drivah, we drive waves using a known K and let w take care of itself...
		#		here, we are frorced in theother direction.. We will tap out an omega at edge of box
		#		and let K take length as dictacted by the dispersion (Bohm-Gross) relation.
		wt = time*self.w_wave 
		#diver_new_value = self.current_amp * np.sin(discrete_factor*global_cell_index + phase_cell_offset + phase_me )
		diver_new_value = self.current_amp * np.sin( wt ) #self.current_amp * np.sin( wt )
		STATE.E_field.components[0][STATE.E_field.num_boundry_cells[0][0]] = diver_new_value
		if mpi_info.rank == 0:
			print "Done with drivah (amp=%e)" % self.current_amp
		if time < self.ending_amp:
			self.current_amp += self.amp_rate
		
	# working (but Bulgarian) full box traveling wave...
	def update(self, time, n, dt):
		self.driver_on = False
		STATE.use_external_E_field = False
		#if self.current_amp > ( self.ending_amp + self.sustain):
		go_away = False
		if time > ( self.rise_time + self.sustain + self.fall_time):
			go_away = True
		if time >= self.echo_2nd_pulse_time and time < (self.echo_2nd_pulse_time+self.rise_time + self.sustain+ self.fall_time):
			go_away = False
			time = time - self.echo_2nd_pulse_time
			if self.start_of_2nd_pulse:
				self.current_amp = self.starting_amp
				k2 = 21.0 # now try 48 /4 #45 #self.num_wavelengths+25
				self.e_profiles[0] = EFieldSinusodialProfile(self.starting_amp, num_wavelengths=k2)	
			self.start_of_2nd_pulse = False
			print "ECHO PUSLEEEE!!!!!!!!!!!!!!"
		if go_away:
			
			return
		
		self.driver_on = True
		STATE.use_external_E_field = True
		
		#self.e_profile.amp = self.current_amp
		
		#calculate the time phase factor (w*time):
		if self.ultra_pure_harmonic_mode == True:
			# generate the k dependence.....
			self.num_wavelengths = 2
			#a = zeros( self.num_cells_in_box+1, dtype='complex128')
			a = zeros( self.num_cells_in_box )
			phaseshift1 = np.exp(1j*2.00*np.pi / float(self.num_cells_in_box) * float(self.num_cells_wave_travels_in_dt) * 0.0)
			phaseshift2 = np.conjugate( phaseshift1) # np.exp(1j*2.00*np.pi / self.num_cells_in_box * self.num_cells_wave_travels_in_dt + np.pi)
			#phaseshift1 = 1.0
			#phaseshift2 = 1.0
			clf()
			a[1] = 1.0
			ddd = fft.irfft( a ) 
			plot(ddd )
			print fft.rfft( ddd)
			savefig('nose.png')
			exit(-1)			
			
			print phaseshift1
			print np.absolute( phaseshift1 )
			print phaseshift2
			print np.absolute( phaseshift2 )
			a[ self.num_wavelengths ] = phaseshift1# *self.num_cells_in_box/2.0*np.cos(phaseshift1)/2.0/np.pi
			a[ self.num_cells_in_box - self.num_wavelengths ] = phaseshift2 #*self.num_cells_in_box/2.0*np.cos(phaseshift1)/2.0/np.pi
			a[0] = 0.0
			#print a
			waveform = ( fft.ifft ( a[0: self.num_cells_in_box]  ) )
			print waveform
			#print a
			clf()
			#plot(np.log1p(abs(fftshift(fft.fft( waveform ) ) )))
			plot ( waveform )
			print ((fftshift(fft.rfft( np.real(waveform) ) ) ))
			savefig('nose.png')
			exit(-1)
			wt = -1.0*time*self.w_wave
		else:
			wt = -1.0*time*self.w_wave
		
		#wt = 0.0
		#self.e_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah, phase_me=wt )
		#STATE.E_field.add( self.buffah)
		
		self.buffah_accumulate.scalar_mult(0.0)
		
		for i in xrange( 0,len(self.e_profiles) ):
			self.e_profiles[i].amp = self.current_amp
			self.e_profiles[i].apply( STATE.E_field__external, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah, phase_me=wt )
			self.buffah_accumulate.add( self.buffah)
		
		#print "wt %f %f" % (wt, self.w_wave) 
		#np.save( "field_%06d.npy" % n, self.buffah_accumulate.components[0] )
		STATE.E_field__external.copy( self.buffah_accumulate )
		
		#STATE.E_field.copy( self.buffah_accumulate )
		
		if time < self.rise_time:
			self.current_amp += self.amp_rate
			pass
		
		if time > ( self.rise_time+ self.sustain):
			#self.current_amp += self.amp_fall_rate
			if mpi_info.rank == 0:
				print "lowered driver to drivah (amp=%e) (rate=%e) " % (self.current_amp,  self.amp_fall_rate)
		else:
			if mpi_info.rank == 0:
				print "Done with drivah (amp=%e) (rate=%e) " % (self.current_amp,  self.amp_rate)
	
	
	def update1(self, time, n, dt):
		# attempt at simple static wave amping...
		if n % self.update_interval_steps != 0:
			if mpi_info.rank == 0:
				if self.driver_off:
					print "E field driver is off."
				else:
					print "%d steps until E field push. (cuurent amp=%e) " % ( (self.update_interval_steps - (n % self.update_interval_steps)), self.current_amp)
			return
		
		"""
		multiple = time / self.update_interval 
		if multiple <= self.last_seen_multiple:
			return
		self.last_seen_multiple = multiple
		"""
		
		self.current_amp *= self.amp_rate
		if self.current_amp >= self.ending_amp:
			if mpi_info.rank == 0:
				print "wave driver upped wave to %f already and is OFF." % self.current_amp
				self.driver_off = True
			return
		
		self.e_profile.amp = self.current_amp
		self.e_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah )
		STATE.E_field.add( self.buffah)
		
		if mpi_info.rank == 0:
			print "wave driver upped wave to %f" % self.current_amp
		
		#vth = CFG.species_config[0].temp_profile.val
		#w_wave = np.sqrt( 1 + 3*CFG.E_field_profile.wavenumber*CFG.E_field_profile.wavenumber*vth*vth)
		#print "freq (BG): %e " % w_wave
		#print "vphase: %e" % (w_wave/ CFG.E_field_profile.wavenumber)
		#print "K l_debye: %f" % (CFG.E_field_profile.wavenumber*vth)
		#CFG.E_field_profile.v_phase = (w_wave/ CFG.E_field_profile.wavenumber)


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
		self.momenum_axis_inverse = zeros( self.nump )
		self.momenutm_axis_deltas = zeros( self.nump )
		
		for i in xrange(0, self.nump):
			self.momenum_axis_inverse[i] = 1.0/F.momentum_axis.values[i]
		
		#momenum_axis_inverse= 1.0 / F.momentum_axis.values[:]
		self.momenutm_axis_deltas[:] = F.momentum_axis.values[:]
		self.momenutm_axis_deltas[1:] -= F.momentum_axis.values[0:self.nump-1]
		
		
		
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
	
	#@profile_func(filename="testme.grind")
	def get_full_lm( self, F, n, t, get_mb=False):
		if self.da_real_distribution_funcion_last_update == n:
			return self.full_distbution_function
		
		num_bndry_cells_lower = F.num_boundry_cells[0][0]
		num_cells = F.num_plain_cells[0]
		
		vals_cartestian_p_cos_thetas = zeros( F.momentum_axis.num_points )
		vals_cartestian_p_transverse_step_size = zeros( F.momentum_axis.num_points )
		vals_plm_at_thetas = zeros( F.momentum_axis.num_points )
		vals_plm_at_thetas_minus = zeros( F.momentum_axis.num_points )
		vals_cartestian_p_thetas_minus = zeros( F.momentum_axis.num_points )
		rsintheta_delta = np.zeros ( F.momentum_axis.num_points )
		rsintheta_cell = np.zeros ( F.momentum_axis.num_points )
		
		values_at_a_px_temp = zeros( (F.momentum_axis.num_points) )
		values_at_a_px_temp_minus = zeros( (F.momentum_axis.num_points) )
		values_at_a_px_accumlator = zeros( (F.momentum_axis.num_points) )
		values_at_a_px_accumlator_minus = zeros( (F.momentum_axis.num_points) )
		thetas_temp = np.zeros(  F.momentum_axis.num_points )
		thetas_delta  = np.zeros(  F.momentum_axis.num_points )
		
		
		last_p_seen = -1
		
		#print "starting dist calc."
		for x in range( 0, num_cells ):
			for p in xrange( 0, F.momentum_axis.num_points ):
				# at this |p|, find angles for this |p| and lower....
				# do this |p|
				vals_cartestian_p_cos_thetas[p] = 1.0
				vals_cartestian_p_thetas_minus[p] = -1.0
				vals_cartestian_p_transverse_step_size[p] = 0.0
				values_at_a_px_accumlator[:] = 0.0
				values_at_a_px_accumlator_minus[:] = 0.0
				
				current_value_px = F.momentum_axis.values[p]
				current_value_px2 = current_value_px*current_value_px
				
				vals_cartestian_p_cos_thetas[p+1:] = current_value_px*self.momenum_axis_inverse[p+1:]
				
				vals_cartestian_p_thetas_minus[p+1:] =vals_cartestian_p_cos_thetas[p+1:]
				vals_cartestian_p_thetas_minus[p+1:] *= -1.0
				# nor arccos because we are oging to immediatly calulate cos(theta) anyeway (for the Ylm)
				#vals_cartestian_p_thetas[p+1:] = np.arccos(vals_cartestian_p_thetas[p+1:])
				vals_cartestian_p_transverse_step_size[p+1:] = self.momenum_axis_sq[p+1:]
				vals_cartestian_p_transverse_step_size[p+1:] -= current_value_px2
				np.sqrt( vals_cartestian_p_transverse_step_size[p+1:], vals_cartestian_p_transverse_step_size[p+1:])
				thetas_temp[p:] = np.arccos(  vals_cartestian_p_cos_thetas[p:] )
				thetas_temp[p] = 0.0
				# thetha_delta will have [p] = 0.0... and the last value withh be theta(N-(N-1))
				thetas_delta[p:] 	= thetas_temp[p:]
				thetas_delta[p+1:] 	-= thetas_temp[p:self.nump-1]
				
				# so now thetas_temp has all the theta angles of one 'ray' of the points...
				#	and theta_deltas has the d(theta) for each step (ignoing first which ois on the Px axis so need sperate treatment.
				
				nump_valid_values = self.nump - p
				
				#if p < F.momentum_axis.num_points-1:
				#	vals_cartestian_p_transverse_step_size[p] = vals_cartestian_p_transverse_step_size[p+1]/2.0
				
				# now apply this grid to each spherical harmonic...
				sum_over_all_harmonics = 0
				sum_over_all_harmonics_minus = 0
				
				lin_harmonic_idx = 0
				for l in xrange(0,F.num_l):
					m = 0
					index = F.linear_index_for_harmonic(l)
					# 
					# numpy's Associated legendre function (lpmv) uses the Condon-Shortly phase (-1)^m (checked)
					#
					# notice destintion index shift to 0->(nump_valid_values-1)
					vals_plm_at_thetas[p:] = scipy.special.lpmv(m,l, vals_cartestian_p_cos_thetas[p:])
					values_at_a_px_temp[0:nump_valid_values] = vals_plm_at_thetas[p:]
					
					vals_plm_at_thetas_minus[p:] = scipy.special.lpmv(m,l, vals_cartestian_p_thetas_minus[p:])
					values_at_a_px_temp_minus[0:nump_valid_values] = vals_plm_at_thetas_minus[p:]
					
					valid_absp_cells = F.cells[num_bndry_cells_lower+x].harmonics[index].momentum_projection[p:]
					values_at_a_px_temp[0:nump_valid_values] *= valid_absp_cells
					values_at_a_px_temp_minus[0:nump_valid_values] *= valid_absp_cells
					# values_at_a_px_temp/minus now contains the value of the current l,m harmonic at the grid points that are valid for the current p (abs_p)
					# now copy values in with the other harmonics at this 
					values_at_a_px_accumlator[0:nump_valid_values] += values_at_a_px_temp[0:nump_valid_values]
					values_at_a_px_accumlator_minus[0:nump_valid_values] += values_at_a_px_temp_minus[0:nump_valid_values]
					
					
				# we are integrating in shpereical coords... so we need a volume element p^2sin(theta) d(theta) d(phi) (dp)
				#		the "vals_cartestian_p_transverse_step_size" is [p*sin(theta)] dp is the itervals defined betwen out theta points.
				#		d_phi the the constant angle step size. So we need 
				
				offset = 0
				
				# calculate the delta(P*sin(theta)) between each evaluation point.
				rsintheta_delta[0:nump_valid_values-1] = vals_cartestian_p_transverse_step_size[p+1:]
				rsintheta_delta[0:nump_valid_values] -=  vals_cartestian_p_transverse_step_size[p:]
				
				
				# calculate the cells sizes that emcompass each evaluation point...
				rsintheta_cell[0] =  rsintheta_delta[0]
				rsintheta_cell[nump_valid_values-1] =  rsintheta_delta[nump_valid_values-2]
				if nump_valid_values > 1:
					rsintheta_cell[1:nump_valid_values-1] = rsintheta_delta[0:nump_valid_values-2]
					rsintheta_cell[1:nump_valid_values-1] += rsintheta_delta[1:nump_valid_values-1]
				rsintheta_cell /= 2.0
				
				# multiply that interesting SH's value at the interction point by this cell size
				values_at_a_px_accumlator[0:nump_valid_values] *=  rsintheta_cell[0:nump_valid_values]
				values_at_a_px_accumlator_minus[0:nump_valid_values] *=  rsintheta_cell[0:nump_valid_values]
				
				values_at_a_px_accumlator[0:nump_valid_values] *= vals_cartestian_p_transverse_step_size[p:]
				values_at_a_px_accumlator_minus[0:nump_valid_values] *= vals_cartestian_p_transverse_step_size[p:]
				
				# Do the integration sum (integration is just sum at this point since all the wighting factors were already multiplied in)
				value = np.sum( values_at_a_px_accumlator[0:nump_valid_values] )
				value_minus = np.sum( values_at_a_px_accumlator_minus[0:nump_valid_values] )
				
				# save the final value... and don't forget to multiply by the 2*PI that comes from integration over the lonely azimuthal angle.
				self.final_data[x, self.nump-p-1] = value_minus * 2.0*np.pi
				self.final_data[x, self.nump+p] = value * 2.0*np.pi
		
		if mpi_gather.gather_2d( "out_f", 	STATE.spatial_axes__global[0].num_points, 
											2*F.momentum_axis.num_points, self.final_data):
			# only rank 0 will get here.
			print "RANK 0 got here!!!"
			self.full_distbution_function = mpi_gather.get_2d_data("out_f")
		else:
			self.full_distbution_function = None
		print ""
		
		return self.full_distbution_function

class Axis:
	def __init__(self,  min=None, step=None, num_points=None, max=None, axis=None):
		self.init(min, step, num_points, max, axis)
	
	def init(self, min, step=None, num_points=None, max=None, axis=None):
		if(axis != None):
			# copy the axis into this object...
			self.min = axis.min
			self.max = axis.max
			self.step = axis.step
			self.num_points = axis.num_points
			self.values = np.array(axis.values, copy=True)
			return
		if min != None and max != None and num_points != None:
			self.min = min
			self.max = max
			self.num_points = num_points
			self.values = np.zeros( num_points )
			for i in xrange(0, num_points):
				self.values[i] = i
			self.values *= (max - min)/(num_points - 1)
			self.values += min
			self.step = self.values[1] - self.values[0]
			#print "AXIS MADE!!!"
			#print min
			#print max
			#print num_points
			#print self.values
			#(self.values, self.step) = np.linspace(min, max, num=num_points, endpoint=False, retstep=True)
			return
		if min != None and step != None and num_points != None:
			self.min = min
			self.num_points = num_points
			self.max = num_points*step + min
			self.values = [(x*step+min) for x in range(0, num_points)]
			return
		raise Exception("Invalid Axis specification..")
		
		self.num_cells = slef.num_points
	# these are just syntax suger
	def dx(self):
		return self.step
	def dp(self):
		return self.step
	

# Fields will not be interleaved for now...
# Axes and distances in this class are all local to the node..
class Field1d:
	def __init__(self, num_cells, num_components, num_boundry_cells):
		
		self.num_boundry_cells = num_boundry_cells
		self.num_plain_cells = np.array( [num_cells,0,0] )
		self.num_total_cells = np.array( [0,0,0] )
		self.num_total_cells[0] = self.num_plain_cells[0] + self.num_boundry_cells[0][0] + self.num_boundry_cells[0][1]
		#self.num_total_cells[1] = self.num_plain_cells[1] + self.num_boundry_cells[1][0] + self.num_boundry_cells[1][1]
		#self.num_total_cells[2] = self.num_plain_cells[1] + self.num_boundry_cells[2][0] + self.num_boundry_cells[2][1]
		
		self.num_components = num_components
		self.components = []
		for i in xrange(0, self.num_components):
			self.components.append( np.zeros( self.num_total_cells[0] ) )
		
		# sometime an axis that covers the guard cells too is needed...
		new_min = STATE.spatial_axes[0].min-(self.num_boundry_cells[0][1]*STATE.spatial_axes[0].dx())
		
		self.spatial_axes_with_gc = [None, None, None]
		self.spatial_axes_with_gc[0] = Axis(min=new_min, num_points = self.num_total_cells[0], step=STATE.spatial_axes[0].dx() )
	
	def read(self, x, component):
		return (self.components[component])[x];
	def clone(self,copy_data = False):
		new_guy = Field1d( self.num_plain_cells[0],  self.num_components, self.num_boundry_cells)
		if copy_data:
			for i in xrange(0, self.num_components):
				new_guy.components[i][:] = self.components[i][:]
		return new_guy
	def scalar_mult(self, scalar):
		for i in xrange(0, self.num_components):
				self.components[i] *= scalar
		pass		
	def add(self, to_add):
		for i in xrange(0, self.num_components):
				self.components[i] += to_add.components[i]
		pass		
	def madd(self, scalar, to_add):
		for i in xrange(0, self.num_components):
				temp = np.copy( to_add.components[i] )
				temp *= scalar
				self.components[i] += temp
		pass		
	def copy(self, f_to_copy):
		for i in xrange(0, self.num_components):
			self.components[i][:] = f_to_copy.components[i][:]
	def Dx(self, component):
		data = self.components[component]
		for x in range(0, self.num_total_cells[0]-2):
			data[x] -= data[x+2]
		for x in range(self.num_total_cells[0]-3, -1, -1):
			data[x+1] = data[x]
		# normalize.
		data /= STATE.spatial_axes[component].dx()
		
class SHarmonic:
	def __init__(self, l, m, linear_index=-1, momentum_axis = None):
		self.my_l = l
		self.my_m = m
		self.my_linear_index = linear_index
		self.init(momentum_axis)
	def init(self, momentum_axis):
		self.momentum_axis = momentum_axis
		if momentum_axis != None:
			self.momentum_projection = np.zeros( (momentum_axis.num_points) )
		pass
	def copy(self, sh_to_copy):
		if self.my_l != sh_to_copy.my_l:
			raise Exception("Cannot copy sphereical harmonic because 'my_l' paramters differs.")
		if self.my_m != sh_to_copy.my_m:
			raise Exception("Cannot copy sphereical harmonic because 'my_m' paramters differs.")
		if self.my_linear_index != sh_to_copy.my_linear_index:
			raise Exception("Cannot copy sphereical harmonic because 'my_linear_index' paramters differs.")
		if self.momentum_axis != None:
			if sh_to_copy == None:
				raise Exception("Cannot copy sphereical harmonic because 'momentum_axis' is Null is source.")
			if self.momentum_axis.num_points != sh_to_copy.momentum_axis.num_points:
				raise Exception("Cannot copy sphereical harmonic because 'momentum_axis' is is differnt size.")
			self.momentum_axis = sh_to_copy.momentum_axis
		elif sh_to_copy != None:
				raise Exception("Cannot copy sphereical harmonic because 'moentum_axis' is Null is source.")
		
		self.momentum_axis = sh_to_copy.momentum_axis
		#np.copyto(self.momentum_projection, sh_to_copy.momentum_projection)
		#TODO: replace with faster
		self.momentum_projection[:] = sh_to_copy.momentum_projection[:]


class DistributionFunctionSpatialCell:
	def __init__(self, num_l, momentum_axis = None):
		if(num_l < 1):
			raise Exception("Invalid value for 'num_l'.. is needs to be >= 1")		
		num_harmonics = 0
		
		# 1d code restricted to num_m = 1
		num_m = 1
		
		# ultra lazy way to count the number of harmonics we will have.
		for l in range(0, num_l):
			if(l < num_m):
				num_harmonics += l + 1
			else:
				num_harmonics += num_m
		
		# so num_harmonics now has the proper number of harmonics.. let's just check..
		if(num_harmonics < 1):
			raise Exception("Invalid number of harmonics!")
		
		self.num_l = num_l
		self.num_m = 1
		
		self.num_harmonics = num_harmonics
		self.harmonics = []
		
		# make list of harmonics...
		for l in xrange(0, num_l):
			self.harmonics.append( SHarmonic(l,0,len(self.harmonics),momentum_axis) )
		pass
		
	def init(self):
		pass
	
	def copy(self, f_to_copy):
		if self.num_harmonics != f_to_copy.num_harmonics:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_harmonics' parameters")
		for i in xrange(0, self.num_harmonics):
			self.harmonics[i].copy( f_to_copy.harmonics[i] )
		pass


# cells and distances in this class are all local to the node...
class DistributionFunction:
	
	def __init__(self, num_cells, num_l, momentum_axis, num_boundry_cells):
		
		self.num_plain_cells 	= np.array( [num_cells,0,0] )
		self.num_boundry_cells 	= num_boundry_cells
		self.num_total_cells 	= np.array( (3,) )
		self.num_total_cells[0] = self.num_plain_cells[0] + self.num_boundry_cells[0][0] + self.num_boundry_cells[0][1]
		#self.num_total_cells[1] = self.num_plain_cells[1] + self.num_boundry_cells[1][0] + self.num_boundry_cells[1][1]
		#self.num_total_cells[2] = self.num_plain_cells[2] + self.num_boundry_cells[2][0] + self.num_boundry_cells[2][1]
		#self.num_total_cells[1] = 0
		#self.num_total_cells[2] = 0
		
		
		self.dim_x = self.num_plain_cells[0]; 
		self.dim_y = self.num_plain_cells[1]; 
		self.dim_z = self.num_plain_cells[2];
		self.momentum_axis = momentum_axis
		
		# create spatial cell entries...
		self.cells= [DistributionFunctionSpatialCell(num_l, momentum_axis) for x in range(0, self.num_total_cells[0])]
		
		self.num_l = self.cells[0].num_l
		self.num_m = 1
		
		# create a mapping to quickly link (l,m) coords to linear..
		self.harmonic_ordinal_map = [None] * (self.num_l)
		for i in range(0, self.cells[0].num_harmonics):
			l = self.cells[0].harmonics[i].my_l
			m = 0
			linear_index = self.cells[0].harmonics[i].my_linear_index
			self.harmonic_ordinal_map[l] = linear_index;
	
	def copy(self, f_to_copy):
		if f_to_copy.num_l != self.num_l:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_l' parameters")
		if f_to_copy.num_m != self.num_m:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_m' parameters")
		# SHIT LOOKIN THIS ADAM BUG
		for x in xrange(0, self.num_total_cells):
			self.cells[x].copy(f_to_copy.cells[x])
		pass
	
	def get_harmonic_at_location(self, x, l):
		idx = self.harmonic_ordinal_map[l]
		return self.cells[x].harmonics[idx]
	
	def linear_index_for_harmonic(self, l):
		return self.harmonic_ordinal_map[l]
	def get_harmonic(self, l, x):
		idx = self.harmonic_ordinal_map[l*num_m + m]
		return cells[x].harmonics[idx]
	
	def scalar_mult_spatial_cell( self, scalar, spatial_cell_index):
		for i in xrange(0, self.cells[spatial_cell_index].num_harmonics):
			self.cells[spatial_cell_index].harmonics[i].momentum_projection *= scalar
	
	def scalar_mult(self, scalar):
		for x in xrange(0, self.num_total_cells[0]):
			for i in xrange(0, self.cells[x].num_harmonics):
				self.cells[x].harmonics[i].momentum_projection *= scalar
		pass
	def max(self):
		# mainly for debugging
		temp = []
		for x in xrange(0, self.num_total_cells[0]):
			for i in xrange(0, self.cells[x].num_harmonics):
				temp.append(np.max(self.cells[x].harmonics[i].momentum_projection))
		tempp = np.array( temp)
		
		
	def add(self, to_add):
		for x in xrange(0, self.num_total_cells[0]):
			for i in xrange(0, self.cells[x].num_harmonics):
				self.cells[x].harmonics[i].momentum_projection += to_add.cells[x].harmonics[i].momentum_projection
		pass
	def madd(self, scalar, to_add):
		for x in xrange(0, self.num_total_cells[0]):
			for i in xrange(0, self.cells[x].num_harmonics):
				val = scalar*to_add.cells[x].harmonics[i].momentum_projection
				self.cells[x].harmonics[i].momentum_projection += val
		pass
	
	def harmonic_centric_view(self):
		return HarmonicCentricView( self.num_total_cells, self.momentum_axis)
	
	def copy_into_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			hcv.data[x,:] = self.cells[x].harmonics[idx].momentum_projection[:]
	
	def get_num_harmonics_per_cell(self):
		return self.cells[0].num_harmonics
	def sub_from_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			self.cells[x].harmonics[idx].momentum_projection[:] -= hcv.data[x,:]		
	def add_from_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			self.cells[x].harmonics[idx].momentum_projection[:] += hcv.data[x,:]
	def copy_from_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			cells[x].harmonics[idx].momentum_projection[:] = hcv.data[x,:]
	def mult_by_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			cells[x].harmonics[idx].momentum_projection[:] *= hcv.data[x,:]
	
	def compare(self, Fother):
		if Fother.num_total_cells[0] != self.num_total_cells[0]:
			return False
		for ix in xrange(0, self.num_total_cells[0]):
			my_cell = self.cells[ix]; other_cell = Fother.cells[ix];
			if my_cell.num_harmonics != other_cell.num_harmonics:
				return False
			for sh in xrange(0,  my_cell.num_harmonics):
				my_harmonic = my_cell.harmonics[sh]; other_harmonic = other_cell.harmonics[sh];
				for ip in xrange(0, my_harmonic.momentum_projection.shape[0]):
					if my_harmonic.momentum_projection[ip]!= other_harmonic.momentum_projection[ip]:
						return False
				pass
		pass
		return True
	

class HarmonicCentricView:
	
	def __init__(self, num_total_cells, momentum_axis):
		self.momentum_axis = momentum_axis
		self.num_total_cells = num_total_cells
		self.data = np.zeros( (self.num_total_cells[0], self.momentum_axis.num_points) )
	
	def scalar_mult(self, val):
		self.data *= val
	
	def Dp(self):
		#TODO: this can be done better using numpy magic...
		abs_p_num_points = self.momentum_axis.num_points
		
		# do a derivitive along |p|... set zeroth element to 0
		#		for last element, pretend f(nump) = f(nump-1) and calc the deritivitive. 
		self.temp = np.zeros(self.num_total_cells[0] )
		self.temp[:] = self.data[:,-2]
		self.temp[:] -= self.data[:,-1]
			
		for x in xrange(0, self.num_total_cells[0]):
			for i in range(0, abs_p_num_points-2):
				self.data[x,i] -= self.data[x,i+2]
			for i in range(abs_p_num_points-3, -1, -1):
				self.data[x,i+1] = self.data[x, i]
		
		self.data[:, 0] = 0.0
		self.data[:, -1] = self.temp[:]
	
	def __call__(self, abs_p, x):
		return self.data[x,abs_p]
	
	"""
		Repeat over all spatial cells:
			Element-wise multiply |p| with the array 'mul_array' passed in
				(so 'mul_array' must be the same size as momentum_axis.num_points)
	"""
	def mpaxis(self, mul_array):
		# TDO: make faster with broadcasting
		for x in xrange(0, self.num_total_cells[0]):
			self.data[x,:] *= mul_array
	def mpaxis__constant_at_each_spatial_cell(self, values):
		"""
		if mpi_info.rank == 0:
			print "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
			print ""
			print ""
		"""
		# values[x] evaluates to a constant.
		for x in xrange(0, self.num_total_cells[0]):
			"""
			if mpi_info.rank == 0:
				print "---------------ix=%d" % x
				print self.data[x,:]
				print values[x]
			"""
			self.data[x,:] *= values[x]
		
	#TODO: this can be done better using numpy magic...
	def Dx(self):
		abs_p_num_points = self.momentum_axis.num_points
		num_total_cells = self.num_total_cells[0]
		
		for x in range(0, num_total_cells-2):
			self.data[x,:] -= self.data[x+2,:]
		for x in range(num_total_cells-3, -1, -1):
			self.data[x+1,:] = self.data[x,:]
		pass
	
	def copy(self, hcv):
		self.data[:,:] = hcv.data[:,:]
	def copy_from(self, hcv):
		self.data[:,:] = hcv.data[:,:]
	def copy_into(self, hcv):
		hcv.data[:,:] = self.data[:,:]
	
class Species:
	
	def __init__(self, num_cells, num_l, momentum_axis = None, species_config = None, num_boundry_cells = None):
		
		self.species_config = species_config
		self.momentum_axis = momentum_axis
		
		self.F  = DistributionFunction(num_cells, num_l, momentum_axis, num_boundry_cells)
		self.Y0 = DistributionFunction(num_cells, num_l, momentum_axis, num_boundry_cells)
		self.Y1 = DistributionFunction(num_cells, num_l, momentum_axis, num_boundry_cells)
		self.Yh = DistributionFunction(num_cells, num_l, momentum_axis, num_boundry_cells)
		
		# this is just for  debugging.. take out
		self.Ytemp = DistributionFunction(num_cells, num_l, momentum_axis, num_boundry_cells)
		self.J = None

		#self.num_total_cells = self.F.num_total_cells
		self.num_l = self.F.num_l
		self.num_m = 1
		
		# TODO: MOVE THIS
		self.moving_frame_u_field = STATE.E_field.clone()   #Field1d np.zeros( self.F.num_total_cells[0] )
		
		# other species specific data.. maybe fields etc...
		pass

class ConstantProfile:
	def __init__(self, val, x_start=None,  x_end=None):
		self.val = val
		self.limits = False
		
		if x_start != None and x_end != None:
			self.limits = True
			self.x_start = x_start
			self.x_end = x_end
	
	def apply(self, temperature_grid, spatial_axis_global, spatial_axis__local_in_global, pt_amb, species):
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			if self.limits:
				if x_val >= self.x_start and x_val <= self.x_end:
					temperature_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
				else:
					temperature_grid[species.F.num_boundry_cells[0][0] + n] *= pt_amb
			else:
				temperature_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
		pass
	
	def apply_density( self, density_grid, spatial_axis_global, spatial_axis__local_in_global, density_amb, species):
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			if self.limits:
				if x_val >= self.x_start and x_val <= self.x_end:
					density_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
				else:
					density_grid[species.F.num_boundry_cells[0][0] + n] *= density_amb
			else:
				density_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
		#for i in xrange(0, species.F.num_boundry_cells[0][0]):
		#	density_grid[n] = density_amb
		#for i in xrange(0, species.F.num_boundry_cells[0][1]):
		#	density_grid[-i] = density_amb

		
class SinusodialProfile:
	def __init__(self, pt_sinmax, num_wavelengths=1.0, phase=0):
		self.pt_sinmax = pt_sinmax
		self.num_wavelengths = num_wavelengths
		self.phase = phase
	
	def apply(self, temperature_grid, spatial_axis_global, spatial_axis__local_in_global, pt_amb, species, num_wavelengths = 1.0):
		
		self.umax = pt_amb
		self.Aver = self.pt_sinmax
		self.Ampl = 0.5*(self.umax*self.umax - self.Aver*self.Aver)
		self.Aver = 0.5*(self.umax*self.umax + self.Aver*self.Aver)
		
		box_size = (spatial_axis_global.values[-1] - spatial_axis_global.values[0] )/ num_wavelengths
		dx = spatial_axis_global.dx()
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			temperature_grid[species.F.num_boundry_cells[0][0] + n] = self.Aver+self.Ampl*np.sin(2.0*np.pi*x_val/(box_size+dx) )
			temperature_grid[species.F.num_boundry_cells[0][0] + n] = np.sqrt( temperature_grid[species.F.num_boundry_cells[0][0] + n]) 
		pass
		
class EFieldSolve1d:
	def get_charge_density_for_species(self, species_index):
		Y00 = STATE.species[species_index].F.harmonic_centric_view()
		STATE.species[species_index].F.copy_into_harmonic_centric_view( 0, Y00)
		integration_temp = np.zeros ( STATE.species[species_index].F.num_total_cells[0] )
		charge_density = np.zeros ( STATE.species[species_index].F.num_total_cells[0] )
		
		# we need to mutiply all values by (4*pi)*p^2*dp
		# take care of multiplying by p^2 (must do this before summing over p since value differs point by point)
		Y00.mpaxis(STATE.species[species_index].F.momentum_axis.values)
		Y00.mpaxis(STATE.species[species_index].F.momentum_axis.values)
		# sum (e.g. integrate) and weight 1st and last cells by .5
		charge_density[:] =  np.sum(Y00.data, axis=1)
		integration_temp[:] = Y00.data[:, 0]; integration_temp[:] += Y00.data[:, -1]; integration_temp *= 0.5
		charge_density[:] -= integration_temp[:]
		# multiply by cell size (since cell size constant can multiply after summation)
		charge_density *= STATE.species[species_index].F.momentum_axis.dp()
		# multiply by 4*pi to take care of integration of solid angle (no p dependence so can multiply after summation over p)
		charge_density *= 4.0*np.pi
		
		return charge_density
	
	# given an E field, this creates the consistent charge density so that E and the charge density are properly matched.
	def solve(self, E_field):
		# take d(Ex)/dx
		E = E_field.clone( copy_data = True)
		extant_net_charge_density = np.zeros(  E_field.num_total_cells[0] )
		
		E.Dx(0); net_charge_density = E.components[0]
		
		extant_charge_density = []
		for i in xrange(0, STATE.config.num_species):
			extant_charge_density.append( self.get_charge_density_for_species(i) )
			#print "%d out of %d" % ( i, STATE.config.num_species)
			#print extant_charge_density[-1][:]
			extant_net_charge_density += ( extant_charge_density[-1][:]* STATE.species[i].species_config.charge)
		
		for ix in xrange(E_field.num_boundry_cells[0][0], E_field.num_boundry_cells[0][0] + E_field.num_plain_cells[0]):
			net_charge_density_needed = net_charge_density[ix] - extant_net_charge_density[ix]
			amount_needed_per_species = net_charge_density_needed / float(STATE.config.num_species)
			for i in xrange(0, STATE.config.num_species):
				factor =  extant_charge_density[i][ix] - amount_needed_per_species*STATE.species[i].species_config.charge
				factor /= extant_charge_density[i][ix]
				STATE.species[i].F.cells[ix].harmonics[0].momentum_projection *= factor
		
		pass
		#E_field.components[0] *= 0.0
		#exit(-1)
	
	# one species for now...
	#def sinusodial_density_perturbation( self, ampltitude, num_wavelenths, mode_add = False):
	#	charge_density = self.get_charge_density_for_species(0)
	#	charge_density *= STATE.species[0].species_config.charge
		
		
		

class EFieldSinusodialProfile:
	def __init__(self, amp, phase=0.0, num_wavelengths=None, wavelength=None):
		self.amp = amp
		self.phase = phase
		self.num_wavelengths = num_wavelengths
		self.wavelength = wavelength
		if num_wavelengths == None and wavelength==None:
			self.num_wavelengths = 1.0
		
		
	def apply_frequency(self, E_field, w, time, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset, buffah=None, phase_me = 0.0):
		box_size = spatial_axis_global.values[-1] - spatial_axis_global.values[0]
		
		"""
		data = E_field.components[0]
		if buffah != None:
			data = buffah.components[0]
			dx = spatial_axis_global.dx()
			self.wavelength = box_size/ (float( self.num_wavelengths))
			self.wavenumber = 2.0*np.pi/(self.wavelength + dx)
			discrete_factor = np.pi*2*float( self.num_wavelengths)/spatial_axis_global.values.shape[0]
			
			data = E_field.components[0]
			if buffah != None:
				data = buffah.components[0]
			
			#for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
		pass
		"""
	
	
	# should prolly move this...
	def apply_perturbation(self, F, E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset ):
		# tkae incoming F and make it into (1+amp*cos(k*x))*F
		#print "NO!!!"
		#exit(-1)
		# first make the cosine...
		box_size = spatial_axis_global.values[-1] - spatial_axis_global.values[0]
		if self.num_wavelengths != None:
			buffer = E_field.clone( )
			self.apply( E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset, buffah=buffer )
			
			"""
			if( mpi_gather.gather(  buffer.components[0], remove_guard_cells=True ) ):
				print  mpi_gather.stiching_buffer
				print np.mean ( mpi_gather.stiching_buffer )
				exit(-1)
			"""
			
			buffer.components[0] += 1.0
			for ix in xrange( 0, F.num_total_cells[0] ):
				perturbattion_value = buffer.components[0][ix]
				F.scalar_mult_spatial_cell( perturbattion_value, ix )
		
		
			
	def apply(self, E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset, buffah=None, phase_me = 0.0):
		
		box_size = spatial_axis_global.values[-1] - spatial_axis_global.values[0]
		if self.num_wavelengths != None:
			# we will setup the profile so that this many wavelength are in the box!
			
			box_size = spatial_axis_global.values[-1] - spatial_axis_global.values[0]
			dx = spatial_axis_global.dx()
			#self.wavelength = box_size/ (float( self.num_wavelengths))
			#self.wavenumber = 2.0*np.pi/(self.wavelength + dx)
			#self.wavenumber = 2.0*np.pi/(self.wavelength)
			
			
			#discrete_factor = np.pi*2*float( self.num_wavelengths)/spatial_axis_global.values.shape[0]
			#discrete_factor = np.pi*2*float( self.num_wavelengths)/( (float(spatial_axis_global.values.shape[0]+0.0)))  ##BIG VHNSGR
			
			
			discrete_factor = np.pi*2.0*float( self.num_wavelengths)/( (float(spatial_axis_global.values.shape[0]+0.0)))
			phase_cell_offset = int( self.phase / discrete_factor)
			phase_cell_offset =  discrete_factor*float(phase_cell_offset)
			
			#print "IN DRIVER..... %f" % box_size #spatial_axis_global.values.shape[0]
			#print spatial_axis_global.values.shape[0]
			#print dx
			
			
			if buffah != None:
				data = buffah.components[0]
			else:
				data = E_field.components[0]
			
			indic = []
			#for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			#	global_cell_index = float( global_cell_offset[0] + n )
			#	data[ E_field.num_boundry_cells[0][0] + n] = self.amp * np.sin(discrete_factor*global_cell_index + phase_cell_offset + phase_me )
			
			num = len(spatial_axis__local_in_global.values) + E_field.num_boundry_cells[0][0] + E_field.num_boundry_cells[0][1]
			for n in xrange(0, num):
				global_cell_index = float( global_cell_offset[0] + n - E_field.num_boundry_cells[0][0])
				data[n] = self.amp * np.cos(2.0*np.pi*float( self.num_wavelengths)/(box_size +dx )*global_cell_index*dx) # + phase_cell_offset + phase_me )
				#global_cell_index = float( global_cell_offset[0] + n - E_field.num_boundry_cells[0][0])
				##indic.append( global_cell_index )
				#data[n] = self.amp * np.cos(discrete_factor*global_cell_index) # + phase_cell_offset + phase_me )
			
				
			#self.v_phase = 1.0/self.wavenumber
			#print phase_cell_offset
			#print phase_me
			#print "FOR %d (%d): \n %s" % ( mpi_info.rank, len( indic), np.array( indic ) )
			
			self.wavelength = box_size / float( self.num_wavelengths)
			self.wavenumber = 2*np.pi / self.wavelength
			#self.v_phase = 1.0/self.wavenumber
			return
		else:
			# use a specifc wavelength...
			self.wavenumber = 2.0*np.pi/(self.wavelength)
			self.num_wavelengths = box_size/self.wavenumber 
		data = E_field.components[0]
		if buffah != None:
			data = buffah.components[0]
		
		print "CCCCCCUNT"
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			data[ E_field.num_boundry_cells[0][0] + n] = self.amp * np.sin(self.wavenumber*x_val + self.phase )
		self.v_phase = 1.0/self.wavenumber
		#print "if wp=1, v_phase = wp/k, vphase=%f" % ( self.v_phase )
	
class GaussianProfile:
	def __init__(self, ptx, stdd, center ):
		self.ptx = ptx
		self.stdd = stdd
		self.center = center
		
	#def apply(self, temperature_grid, spatial_axis, pt_amb, species):
	def apply(self, temperature_grid, spatial_axis_global, spatial_axis__local_in_global, pt_amb, species):
		# temperature_grid is initalized to al 1.0
		for (n,x) in enumerate(spatial_axis__local_in_global.values):
			arg = (x - self.center) / self.stdd
			val = self.ptx*self.ptx*math.exp(-1.0*arg*arg)
			#print temperature_grid[n]
			temperature_grid[n + species.F.num_boundry_cells[0][0] ] *= val
		pass
		
	
		
class SpeciesConfig:
	def __init__(self):
		
		self.name = "No name"
		self.charge = -1.0
		self.mass = 1.0
		self.temperature_profile = None
		

class Config:
	
	def __init__(self, num_cells__global, system_size_min__global, system_size_max__global, num_species):
				
		self.dim = 1
		self.num_spatial_dims = 1
		self.num_cells__global 			= np.array( num_cells__global , dtype='int32')
		self.system_size_min__global 	= np.array(system_size_min__global)
		self.system_size_max__global 	= np.array(system_size_max__global)
		self.system_size__global 		= self.system_size_max__global - self.system_size_min__global
		
		
		#if len(num_cells__global) == 1:
		self.total_num_cells__global 	= self.num_cells__global[0]
		self.cell_size 					= np.divide(self.system_size__global, self.num_cells__global[0])
		self.cell_volume 					= self.cell_size[0]
		
		self.num_l 					= -1
		self.num_m 					= 1
		self.num_species 			= num_species
		
		self.dt 						= 0.0
		self.t_end 					= 0.0
		
		self.E_num_guard_cells = np.array((3,2))
		self.B_num_guard_cells = np.array((3,2))
		self.J_num_guard_cells = np.array((3,2))
		
		self.species_config = [None] * num_species
		

# This holds various dimensions of the problem... make this object the object-of-merit in these matters...
#		If no pre/in-fix is given, assume quantity id local to node. Other quantites will be be label global etc...
class DaState:
	def __init__(self, CFG ):
		self.species = [None]*CFG.num_species
	
		self.E_field 	= None
		self.B_field 	= None
		self.config 	= CFG
		self.metadata = SimMetadata()
		self.num_timesteps_between_outputs = -1
		
		self.use_external_E_field = False
		self.E_field__external = None
		
		
		# These are the main dimension of the problem....
		system_size_global_min 			= CFG.system_size_min__global
		system_size_global_max 			= CFG.system_size_max__global
		num_cells__global 				= CFG.num_cells__global 
		system_size__global				= CFG.system_size__global

		
		# create the spatial_axes__global axis.
		axis_x = Axis(	system_size_global_min[0], max=system_size_global_max[0], 
									num_points=num_cells__global[0])
		
		axis_y = None; axis_z = None;
		
		self.spatial_axes__global = [axis_x, axis_y, axis_z]

		# this is local to the node.. it extends from cell=(0, num_cells on th node)  x=(0.0, max_x_on_node)
		#	most claculations should make use of this... any explicit knowledge of gloabl positions should
		#	be used sparingly.
		# TODO: think a bit and talk to Micheal.. I am split.. part of me doesnt want to put these local quantites
		#			in the state.. I want to have flexability to have different object have different node layouts etc..
		#			I think that access to those things should occur throught the objects in question.
		#	It thik I will leave spatital_axes.. but then not have any specic boundry cell stuff.
		self.spatial_axes 		= [None, None, None]
		#self.num_total_cells 	= None
		#self.num_plain_cells 	= None
		#self.num_boundry_cells 	= None		
		# self.spatial_axes_with_gc = None
	
		
	def make_array(self, val):
		try:
			ftry = len(val)
			return val
		except:
			pass
		return [val]
	
	def init(self, system_size_global_min, system_size_global_max, system_size__cells_global, number_species, num_l, config=None):
		
		self.system_size_global_min =  self.make_array(system_size_global_min)
		self.system_size_global_max = self.make_array(system_size_global_max)
		self.system_size__cells_global = self.make_array(system_size__cells_global)
		self.number_species = number_species
		self.num_l = num_l
		self.num_m = 1
		
		#TODO: make general..
		#axis_x_node = Axis(0.0, max=system_size_node[0], num_points=system_size__cells_node[0])
		if config == None:
			self.config = Config(self.system_size__cells_global,  self.system_size_global_min, self.system_size_global_max, self.number_species)
		self.species = [None] * number_species
		
	
#
#
#
#		Start of Operatations....
#
#
class ElectricFieldOnDistribution1d:
	
	def __init__(self, species, state, CFG):
		num_l = species.num_l
		num_m = 1

		self.momentum_axis = Axis(axis=species.momentum_axis)
		
		# 'declare' the member vars...
		self.A1 = np.zeros( (num_l) )
		self.A2 = np.zeros( (num_l) )
		self.Hp0 = np.ones( (num_l) )
		
		self.H =  species.F.harmonic_centric_view()
		self.G =  species.F.harmonic_centric_view()
		self.TMP = species.F.harmonic_centric_view()
		
		self.num_l = num_l
		self.num_m = 1
		
		#now fill 'em up...
		for l in xrange(1, num_l):
			self.A1[l] = (l + 1.0) * l / ( (2.0*l + 1.0)*(l+1.0) )
			self.A2[l] = -(l) / (2.0 * l + 1.0)
		self.A1[0] = -1.0
		self.A1 *= (0.5/self.momentum_axis.dx() )
		self.A2 *= (0.5/self.momentum_axis.dx() )
		self.A2[0] = 1.0
		
		self.Hp0 *= (1.0/self.momentum_axis.values[0])
		ratio = self.momentum_axis.values[0]/self.momentum_axis.values[1]
		for l in xrange(1, num_l):
			self.Hp0[l] = self.Hp0[l-1]*ratio * (2.0*l+1.0)/(2.0*l-1)
		self.Hp0 *= (-2.0*self.momentum_axis.dx())
		
		# inverted momentum axis...
		self.invpr = Axis(axis=self.momentum_axis)
		for i in xrange(0, self.momentum_axis.num_points):
			self.invpr.values[i] = ( 1.0/self.momentum_axis.values[i] )
	
		self.species = species
		self.species_charge = self.species.species_config.charge*(-1.0)
		
	def G00(self, G, Y):
		
		pr = self.momentum_axis.values
		dx = self.momentum_axis.dx()
		
		p0p1_sq = pr[0]*pr[0]/(pr[1]*pr[1])
		inv_mp0p1_sq = 1.0/(1.0-p0p1_sq)
		g_r = -4.0*dx*pr[0]/(pr[1]*pr[1])
		
		# copy the zeroth harmonic into a view.
		# 	and diff it along |p|
		#TODO: put this together.. it will be much faster.
		#Y.copy_into_harmonic_centric_view(0,0,G)
		self.TMP.copy_from( G )
		G.Dp()
		
		"""
		print "p0p1_sq: %f" % p0p1_sq
		print "inv_mp0p1_sq: %f" % inv_mp0p1_sq
		print "g_r: %f" % g_r
		"""
		for x in xrange(0, Y.num_total_cells[0]):
			# now go and fix up the first |p| cell in each spitial cell.
			p0 = self.TMP.data[x,0]
			p1 = self.TMP.data[x,1]
			f00 = ( p0-p1*p0p1_sq)*inv_mp0p1_sq			
			G.data[x, 0] = (p1-f00)*g_r
		# G is the result (but is in same scope has the consuming function)
		pass
	
	def MakeGH(self, l0, Y, tester=None):
		harmonic_linear_index = Y.linear_index_for_harmonic(l0)
		inxpax = Axis(axis = self.invpr)
		inxpax.values *= (-2.0*(l0+1.0)*self.momentum_axis.dx())
		
		# TODO: can speed up at least with array copy.
		Y.copy_into_harmonic_centric_view(l0,self.H)
		Y.copy_into_harmonic_centric_view(l0,self.G)
		#if tester != None:
		#	self.G.copy_from( tester)
		#	self.H.copy_from( tester)
		
		self.G.Dp()
		self.H.mpaxis( inxpax.values )
		self.H.data += self.G.data
		self.G.data *= (-2.0*l0 - 1.0) / l0
		self.G.data += self.H.data
		
		for x in xrange(0, Y.num_total_cells[0]):
			# TODO: put array interface on harmonic view
			self.G.data[x,0] = 0.0
			#if tester != None:
			#	self.H.data[x,0] = tester.data[x,1]*self.Hp0[l0]
			#else:
			self.H.data[x,0] = Y.cells[x].harmonics[harmonic_linear_index].momentum_projection[1]*self.Hp0[l0]
		pass
		
	# TODO: make faster with numpy indexing...
	def state_harmonic_plus_Ex_x_scalar_mxy_matrix(self, target_l, target_distribution_function, scalar, E_field, 
																	_operator):
		try:
			sh_linear_index = target_distribution_function.linear_index_for_harmonic(target_l)
		except:
			return
		
		num_abs_p_cells = self.momentum_axis.num_points
		operator_data = _operator.data
		
		for x in xrange(0, target_distribution_function.num_total_cells[0]):
			momentum_projection = target_distribution_function.cells[x].harmonics[sh_linear_index].momentum_projection
			value_Ex = E_field.read(x, 0) * self.species_charge
			multiplier = value_Ex * scalar
			operator_data[x,:] *= multiplier
			np.add( momentum_projection[:], operator_data[x,:], momentum_projection[:])
		pass
	
	def copy_to_TMP(self, what_to_copy):
		self.TMP.copy_from( what_to_copy )
	
	# this is the 1D calcuation of the efect of E on the distribution function...
	# the 'result' is stored in the distribution funciton variable Yh
	# 
	# The input is Yin and it should not be written to.
	# It is assumed that Yh is initalized properly before entering this fine function.		
	def calc(self, Yin, Yh, species, state, CFG, Ein=None):
		G = self.G
		H = self.H
		TMP = self.TMP
		
		Yin.copy_into_harmonic_centric_view(0,self.G)
		self.G00(self.G, Yin)
		
		# Ex *= A1(0,0); TMP=G; Yh.SH(1,0) += G.mxy_matrix(Ex);
		# so target (l,m) is 1,0; results stored in target Yh; Ex pre multiplied by scalar A1[0,0]; operator is G.
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	1, Yh, self.A1[0], Ein, self.G)
		
		# m = 0, l = 1
		self.MakeGH(1,Yin)
		# Ex *= A2(1,0); Yh.SH(0,0) += H.mxy_matrix(Ex);
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	0, Yh, self.A2[1], Ein, self.H)
		# Ex *= A1(1,0); Yh.SH(2,0) += G.mxy_matrix(Ex);
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	2, Yh, self.A1[1], Ein, self.G)
		
		# m = 0, 1 < l < l0
		for l in xrange(2, species.num_l-1):
			self.MakeGH(l,Yin)
			# Ex *= self.A2[l,0]; Yh.SH(l-1,0) += H.mxy_matrix(Ex)
			self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(l-1, Yh, self.A2[l], Ein, self.H)
			# Ex *= self.A1[l,0]; Yh.SH(l+1,0) += G.mxy_matrix(Ex)
			self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(l+1, Yh, self.A1[1], Ein, self.G)
		
		#m = 0,  l = l0-1
		self.MakeGH(species.num_l-1, Yin)
		# Ex *= self.A2[num_l-1,0]; Yh.SH(l-2,0) += H.mxy_matrix(Ex)
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(species.num_l-2, Yh, self.A2[species.num_l-1], Ein, self.H)
		
		return
	

class Cudddd:
	def __init__(self):
		pass
	
		
class CurrentSolver:
	def __init__(self, species, state, CFG):
		self.f11 = species.F.harmonic_centric_view()
		self.f10 = species.F.harmonic_centric_view()
		self.integration_temp = np.zeros ( species.F.num_total_cells[0] )
				
		# we'll do it micheal's way for now..
		self.p30g = Axis(axis=species.momentum_axis)
		self.delta_p = self.p30g.dx()
		self.small = species.momentum_axis.values[0]
		
		for p in xrange(0, self.p30g.num_points):
			p_val = self.p30g.values[p]
			self.p30g.values[p] = (p_val*p_val*p_val)/(math.sqrt(1.0+p_val*p_val))
		
		self.small *= self.small; self.small *= self.small
		self.small *= species.momentum_axis.values[0]
		self.small *= 0.2; self.small *= 1.0/( species.momentum_axis.values[1] ) 
	
		# not sure where current will end up living...
		# but in 1d it does not survive past this routine.
		self.jayx = Field1d(STATE.spatial_axes[0].num_points, 1, state.E_field.num_boundry_cells)
		self.nump =  species.momentum_axis.num_points
		self.species = species
		self.species_charge = self.species.species_config.charge*(-1.0)
		self.final_jayx_factor = (4.0 * np.pi / 3.0) * self.species_charge
		
		
	# output is Ex...
	def calc(self, Yin, Yh, species, state, CFG, Eout = None):
		
		Yin.copy_into_harmonic_centric_view(1,self.f10)
		self.f10.mpaxis(self.p30g.values)
		
		jayx = self.jayx.components[0]
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
		jayx *= self.final_jayx_factor 
		
		#print "-------"
		#print "----%d \n %s" %( mpi_info.rank, str(jayx) )
		
		Eout.components[0] += jayx
		
# FLUID IONS
class VlasovMovingFrameCalculation:

	def __init__(self, species, state, CFG):
		
		self.momentum_axis = species.momentum_axis
		 
		num_l = species.num_l
		# constant sphereical harmonic-based coefficients for each harmonic
		# These are the additional adiabtibatic compression term
		A1 = np.zeros( (num_l) )
		A2 = np.zeros( (num_l) )
		A3 = np.zeros( (num_l) )
		A4 = np.zeros( (num_l) )
		self.Hp0 = np.ones( (num_l) )
		self.A__G = np.zeros( (num_l) )
		self.A__H = np.zeros( (num_l) )
		
		for l in xrange(0, num_l):
			L = float(l)
			A1[l] = L*(L-1) /( (2*l-3)*(2*l-1) )
			A2[l] = (L+1)*(L+1) / ( (2*L+3)*(2*L+1) )
			A3[l] = (L*L)/( (2*L+1)*(2*L-1) )
			A4[l] = (L+1)*(L+2)/( (2*L+3)*(2*L+5) )
		pass
		
		# now combine these constants to relfect the number of harmnonics we have..
		#
		self.A__G[0] = -1.0*A2[0]; self.A__H[0] = 0.0
		if species.num_l > 2:
			self.A__G[0] += A1[2]
		
		for l in xrange(1, num_l):
			self.A__G[l] = A2[l]
			self.A__H[l] = A3[l]
			
			if l > 2 and (l+2) < num_l:
				self.A__G[l] += A1[l+2]
			
			if l > 1 and (l+2) < num_l:
				self.A__H[l] += A1[l-2]
			
			self.A__G[l] *= l
			self.A__G[l] /= (l+1)
			
		# now multiply with factor to clean up the deritivitvies taken when making G,H
		self.A__G *= .5/species.momentum_axis.dx()
		self.A__H *= (-0.5)/species.momentum_axis.dx()
		
		
		self.U_temp = species.F.harmonic_centric_view()
		
		self.H =  species.F.harmonic_centric_view()
		self.G =  species.F.harmonic_centric_view()
		self.TMP = species.F.harmonic_centric_view()
		# inverted momentum axis...
		self.invpr = Axis(axis=self.momentum_axis)
		for i in xrange(0, self.momentum_axis.num_points):
			self.invpr.values[i] = ( 1.0/self.momentum_axis.values[i] )
		self.Hp0 *= (1.0/self.momentum_axis.values[0])
		ratio = self.momentum_axis.values[0]/self.momentum_axis.values[1]
		for l in xrange(1, num_l):
			self.Hp0[l] = self.Hp0[l-1]*ratio * (2.0*l+1.0)/(2.0*l-1)
		self.Hp0 *= (-2.0*self.momentum_axis.dx())	
	#
	# TODO: These are copied from the E-field-acting-on-F function.. unite them..
	#	
	#		G_00 = - 1 / (2*dp) * result
	#
	def G00(self, G, Y):
		
		pr = self.momentum_axis.values
		dx = self.momentum_axis.dx()
		
		p0p1_sq = pr[0]*pr[0]/(pr[1]*pr[1])
		inv_mp0p1_sq = 1.0/(1.0-p0p1_sq)
		g_r = -4.0*dx*pr[0]/(pr[1]*pr[1])
		
		# copy the zeroth harmonic into a view.
		# 	and diff it along |p|
		#TODO: put this together.. it will be much faster.
		#Y.copy_into_harmonic_centric_view(0,0,G)
		self.TMP.copy_from( G )
		G.Dp()
		
		for x in xrange(0, Y.num_total_cells[0]):
			# now go and fix up the first |p| cell in each spitial cell.
			p0 = self.TMP.data[x,0]
			p1 = self.TMP.data[x,1]
			f00 = ( p0-p1*p0p1_sq)*inv_mp0p1_sq			
			G.data[x, 0] = (p1-f00)*g_r
		# G is the result (but is in same scope has the consuming function)
		pass
	
	def MakeGH(self, l0, Y, tester=None):
		harmonic_linear_index = Y.linear_index_for_harmonic(l0)
		inxpax = Axis(axis = self.invpr)
		inxpax.values *= (-2.0*(l0+1.0)*self.momentum_axis.dx())
		
		# TODO: can speed up at least with array copy.
		Y.copy_into_harmonic_centric_view(l0,self.H)
		Y.copy_into_harmonic_centric_view(l0,self.G)
		
		self.G.Dp()
		self.H.mpaxis( inxpax.values )
		self.H.data += self.G.data
		self.G.data *= (-2.0*l0 - 1.0) / l0
		self.G.data += self.H.data
		
		for x in xrange(0, Y.num_total_cells[0]):
			self.G.data[x,0] = 0.0
			self.H.data[x,0] = Y.cells[x].harmonics[harmonic_linear_index].momentum_projection[1]*self.Hp0[l0]
		pass
	
	def advection_corrections(self, u_field, Yin, Yh, species, state, CFG, Ein=None, debug=False):
		# u_field is a field1d object that hold spaitial velocity data...
		
		#if mpi_info.rank == 0:
		#	print u_field.components[0]
		#	print u_field.components[0].shape
			
		# compute -u_x * d(f_lm)/dx
		for l in xrange(0, species.num_l):
			Yin.copy_into_harmonic_centric_view(l,self.U_temp)
			self.U_temp.Dx()
			self.U_temp.mpaxis__constant_at_each_spatial_cell( u_field.components[0] )
			Yh.sub_from_harmonic_centric_view(l, self.U_temp)
		#if mpi_info.rank == 0:
		#	exit(-1)
	def adibatic_compressions_corrections(self, u_field, Yin, Yh, species, state, CFG, Ein=None, debug=False):
		
		# find du/dx
		dudx = u_field.clone( copy_data = True )
		dudx.Dx(0)
		
		#
		# The general calculation (for a given L) involves:
		#			G_(L-2), G_L
		#			H_L and H(L+2)
		# m = 0, l = 0
		#
		
		# -----G00 term
		Yin.copy_into_harmonic_centric_view(0,self.G)
		self.G00(self.G, Yin)
		# calculate: p*G_00
		self.G.mpaxis( self.momentum_axis.values )
		# calculate: p*[d(ux)/d(x)]*G_00
		self.G.mpaxis__constant_at_each_spatial_cell( dudx.components[0] )
		self.G.scalar_mult( self.A__G[0] )
		
		#print "adam1"
		#print "------------------------"
		#print np.max( self.G.data )
		#exit(-1)
		Yh.add_from_harmonic_centric_view(0, self.G )
		
		
		for l in xrange(1, species.num_l):
			self.MakeGH( l, Yin )
			#print np.max( self.G.data )
			#print np.max( self.H.data )
			
			self.G.mpaxis( self.momentum_axis.values )
			self.G.mpaxis__constant_at_each_spatial_cell( dudx.components[0] )
			self.G.scalar_mult( self.A__G[l] )
			
			self.H.mpaxis( self.momentum_axis.values )
			self.H.mpaxis__constant_at_each_spatial_cell( dudx.components[0] )
			self.H.scalar_mult( self.A__H[l] )
			
			Yh.add_from_harmonic_centric_view(l, self.G )
			Yh.add_from_harmonic_centric_view(l, self.H )
		
	
	def calc(self, Yin, Yh, species, state, CFG, Ein=None, debug=False):
		# TODO: decide final resting spot of the fluid velocity....
		#self.advection_corrections( species.moving_frame_u_field, Yin, Yh, species, state, CFG, Ein=Ein, debug=debug)
		self.adibatic_compressions_corrections( species.moving_frame_u_field, Yin, Yh, species, state, CFG, Ein=Ein, debug=debug)
	
		


class SpatialAdvection1d:
	
	def __init__(self, species, state, CFG):
		num_l = species.num_l
		num_m = 1
		
		self.num_l = species.num_l
		self.num_m = 1
		self.fd1 = species.F.harmonic_centric_view()
		self.fd2 = species.F.harmonic_centric_view()
		self.db = species.F.harmonic_centric_view()
		
		self.A1 = np.zeros( (num_l) )
		self.A2 = np.zeros( (num_l) )
		adv_cell_size =  CFG.cell_size[0]
		
		self.idx = -1.0/(2.0 * adv_cell_size)
		
		# setup a local |p| momentum grid.
		self.vr = Axis(axis=species.momentum_axis)
		for p in xrange(0,self.vr.num_points):
			self.vr.values[p] = self.vr.values[p] / ( math.sqrt(1.0 + self.vr.values[p] * self.vr.values[p]) )
		
		for l in xrange(0, self.num_l):
			self.A1[l] = self.idx *(-1.0) * (l+1.0) / (2.0*l+1.0)
			self.A2[l] = self.idx *(-1.0) * (l)     / (2.0*l+1.0)
		
	
	def calc(self, Yin, Yh, species, state, CFG, Ein=None, debug=False):
		
		if debug:
			if mpi_info.rank > 0:
				self.debug = False
				debug = False
			else:
				self.debug = True
				debug = True
		
				
		# Do the X advection!!!!!!
		self.vt = Axis(axis=self.vr)
		if debug:
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AVECT THIS MUTHA FUCKA~~~~~~~~~~~~~~~~~~~"
		
		Yin.copy_into_harmonic_centric_view(0,self.fd1)
		self.fd1.Dx()
		
		if debug:
			print "TOP m=0"
			self.dump_hcv( self.fd1, "Derivitive of Y(0,0)" )
					
		self.vt.values *= self.A1[0]
		self.fd1.mpaxis(self.vt.values)
		Yh.add_from_harmonic_centric_view(1, self.fd1)
		
		if(debug):
			self.dump_sh( Yh, 1, 0, "Victim (1,0)" )
		
		
		#TODO: suspect index
		for l in xrange(1, self.num_l-1):
			Yin.copy_into_harmonic_centric_view(l,self.fd1)
			if debug:
					print "In L loop L=%d before 1st derivitive" % ( l,)
					self.dump_hcv( self.fd1, "")
			
			self.fd1.Dx()
			if debug:
				print "In L loop L=%d" % ( l,)
				self.dump_hcv( self.fd1, "Derivitive of Y(%d,0)" % (l,) )
			
			self.vt.values *= ( self.A2[l]/self.A1[l-1] )
			self.fd2.copy( self.fd1 )
			self.fd1.mpaxis(self.vt.values)
			Yh.add_from_harmonic_centric_view(l-1, self.fd1)
			
			self.vt.values *= ( self.A1[l]/self.A2[l] )
			self.fd2.mpaxis(self.vt.values)
			Yh.add_from_harmonic_centric_view(l+1, self.fd2)
			
			if debug:
				self.dump_sh( Yh, l-1, 0, "Victim (%d,0)" % (l-1,) )
				self.dump_sh( Yh, l+1, 0, "Victim (%d,0)" % (l+1,) )
							
		# The highest l
		Yin.copy_into_harmonic_centric_view(self.num_l - 1, self.fd1)
		self.fd1.Dx()
		if debug:
			self.dump_hcv( self.fd1, "After L loop Deritiveitve of Y(%d,0)" % (self.num_l - 1, ) )
		
		self.vt.values *= ( self.A2[self.num_l - 1]/self.A1[self.num_l - 2] )
		self.fd1.mpaxis(self.vt.values)
		Yh.add_from_harmonic_centric_view(self.num_l - 2, self.fd1)

		if debug:
			self.dump_sh(Yh, self.num_l - 2, "After L loop Deritiveitve of Y(%d,0)" % (self.num_l - 2,) )
		self.debug = False
		
	
# Fokker plank has zetho harmonic step that will ensure that 		
class FokkerPlankExplicit:
	
	def __init__(self, species, state, CFG ):
		
		# variables for the RK4 part.. 
		self.F0 = np.zeros ( species.momentum_axis.num_points )
		self.Fh = np.zeros ( species.momentum_axis.num_points )
		self.F1 = np.zeros ( species.momentum_axis.num_points )
		#TODO: get rid of this.. it just for debugging
		self.Ftemp = np.zeros ( species.momentum_axis.num_points )
		
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
		self.if_tridiagonal = CFG.if_tridiagonal
		self.if_implicit1D = CFG.if_implicit1D
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
		
		#print "FOKKERPLANK: kp=%e, c_kpre=%e" % (self.kp, self.c_kpre)
		
		# shorten some nmaes teroparlily to makes things prettier
		pr_axis =  species.momentum_axis
		vr = self.vr
		self.num_pr_cells = pr_axis.num_points
		
		
		#print species.momentum_axis.values
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
		temp1 = np.zeros( (pr_axis.num_points, pr_axis.num_points) )
		self.Alpha_Tri 		= np.matrix( temp1 )
		temp2 = np.zeros( (pr_axis.num_points, pr_axis.num_points) )
		self.Alpha 				= np.matrix( temp2 )
		
		
		init__fokker_plank_implicit__c(	self.vr, self.df0, self.ddf0, self.Scattering_Term, self.Alpha, self.kpre)
		init__fokker_plank_explicit__c( 	self.vr, self.U4, self.U4m1, self.U2, self.U2m1, self.U1, self.U1m1, 
													self.J1, self.I4, self.U3, self.Qn, self.Pn, self.I2, 
													self.G_constant_1, self.G_constant_vr3, self.G_constant_vr5,
													self.G_constant_vr5_2, self.G_constant_vr7,
													self.density_np, self.NB, self.c_kpre )
		
		
		# setup code type paratmeters to python defaults...
		self.set_run_type__fp_explicit_python()
		#self.set_run_type__fp_explicit_cpp()
		
		pass
		
	def set_run_type__fp_explicit_python(self):
		self.lowest_harmonic_fp_slope = self.lowest_harmonic_fp_slope_python
	def set_run_type__fp_explicit_cpp(self):
		self.lowest_harmonic_fp_slope = fokker_plank__explicit__advance
	
	
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
		#print "lnee: %e" % Ln_ee
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
			#print "%e " % self.I4[n],
		#print ""
       
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
				print "They are unequal!"
				exit(-1)
		return
	
	def compare_against_save(self, F2):
		return self.quick_compare_dist(self.Ftemp, F2)
	
	def quick_save_a_dist(self , F1):
		self.Ftemp[:] = F1[:]
		return
	
	def RK4_for_fokker_plank(self, eval_func, F):
		# F is the offical state used in the rest of the code.
		# e.g. updates to 'F' will be passed on to the rest of the code.
		self.F0[:] = F
		self.F1[:] = F
		
		# results of the evaluation are store inside of Fh, the input
		#		must remain unchanaged.
		
		self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		self.compare_against_save( self.F1 )
		
		# f1 = f1 + (h/2)*fh
		self.Fh *= (0.5*self.h)
		self.F1 += (self.Fh)
		# F = F + (h/6)*fh
		self.Fh *= (1.0 / 3.0)
		F += ( self.Fh )
		
		
		#------- Step 2
		self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		self.compare_against_save( self.F1 )
		
		self.F1[:] = self.F0[:]
		# f1 = f0 + (h/2)*fh
		self.Fh *= (0.5*self.h)
		self.F1 += (self.Fh)
		# F = F + (h/3)*Fh --> F + h/6*k1+h/2*k2
		self.Fh *= (2.0/3.0)
		F += ( self.Fh )
		#print "result 2"
		#print F
		
		#------- Step 3
		self.quick_save_a_dist( self.F1)
		eval_func(self.F1, self.Fh)
		self.compare_against_save( self.F1 )
		
		#f1 = f0 + h*Fh
		self.Fh *= (self.h)
		self.F0 += (self.Fh)
		# F = F + (h/3)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3
		self.Fh *= (1.0/3.0)
		F += (self.Fh)
		#print "result 3"
		#print F
		
		#------- Step 4
		self.quick_save_a_dist( self.F0)
		eval_func(self.F0, self.Fh)
		self.compare_against_save( self.F0 )
		
		# F = F + (h/6)*Fh --> F+ h/6*k1+h/3*k2+h/3*k3*+1/6*k4
		self.Fh *= (self.h/6.0)
		F += (self.Fh)
		

#
#---------Higher Order harmonics calculations
#
	# this leaves 'fin' unaltered.. and sets up internal variables.
	#@timeit
	def reset_coeff(self, fin, Delta_t):
		# just for ease for now...
		vr = self.vr; U4=self.U4; U4m1=self.U4m1
		U2=self.U2; U2m1=self.U2m1; U1=self.U1; U1m1=self.U1m1
		J1m=self.J1m; I0=self.I0; I2=self.I2; Scattering_Term=self.Scattering_Term; 
		df0=self.df0; ddf0=self.ddf0; TriI1=self.TriI1; TriI2=self.TriI2;
		
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
		#print "I2"
		#print I2
		#print "I2_temperature"
		#print I2_temperature
		
		# Density integral I0 = 4*pi*int_0^v f(u)*u^2du 
		I0[0] = 0
		for k in xrange(1,  self.num_pr_cells):
		  I0[k]  = U2[k]*fin[k]+U2m1[k]*fin[k-1]
		  I0[k] += I0[k-1];
		
		I0 *= 4.0 * np.pi
		I0_density = I0[self.num_pr_cells-1];
		#print "I0"
		#print I0
		#print "I0_density"
		#print I0_density
		
		# Integral J_(-1) = 4*pi * v * int_0^v f(u)*u^4du
		J1m[self.num_pr_cells-1] = 0
		for k in xrange(self.num_pr_cells-2, -1, -1):
		  J1m[k]  = U1[k+1]*fin[k+1]+U1m1[k+1]*fin[k]; 
		  J1m[k] += J1m[k+1];
		
		for k in xrange(0,  self.num_pr_cells):
		  J1m[k] *= 4.0 * np.pi *vr[k];
		#print "J1m"
		#print J1m
		  
		# COULOMB LOGARITHMS
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		self._LOGee  = self.LOGee (I0_density,I2_temperature)
		self._ZLOGei = self.ZLOGei(I0_density,I2_temperature)
		
		# BASIC INTEGRALS FOR THE TRIDIAGONAL PART
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		
		# Temporary arrays
		for i in xrange(0,  self.num_pr_cells):
			TriI1[i] = I0[i] + (2.0*J1m[i] - I2[i]) / 3.0    	#(-I2 + 2*J_{-1} + 3*I0) / 3
			TriI2[i] = ( I2[i] + J1m[i] ) / 3.0             	#( I2 + J_{-1} ) / 3
			
		# SCATTERING TERM
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		Scattering_Term[:] = TriI1[:]
		Scattering_Term *= self._LOGee
		Scattering_Term[0] = 0.0
		Scattering_Term += self._ZLOGei * I0_density
		Scattering_Term /= self.vr3
		#np.power(vr, 3.0, Scattering_Term)
		Scattering_Term *=  (self.kpre * Dt)
		#print _LOGee
		#print _ZLOGei
		#print "Scattering_Term"
		#print Scattering_Term
		
		# MAKE TRIDIAGONAL ARRAY
		# TODO: is zerosing necessary?
		self.Alpha_Tri[:,:] = 0.0
		self.Alpha_Tri[0,0] = 8.0 * np.pi * fin[0]; 
		
		for i in xrange(1,  self.num_pr_cells-1):
			IvDnDm1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i]  -vr[i-1]))       	# ( v*D_n*D_{n-1/2} )^(-1)
			IvDnDp1 = 1.0 / (vr[i] * 0.5*(vr[i+1]-vr[i-1]) * (vr[i+1]-vr[i]  ))       	#  ( v*D_n*D_{n+1/2} )^(-1)
			Ivsq2Dn = 1.0 / (vr[i] * vr[i]                 * (vr[i+1]-vr[i-1]))       	#  ( v^2 * 2*D_n )^(-1)
			self.Alpha_Tri[i, i  ] = 8.0 * np.pi * fin[i] - TriI2[i] * (IvDnDm1 + IvDnDp1)
			self.Alpha_Tri[i, i-1] = TriI2[i] * IvDnDm1 - TriI1[i] * Ivsq2Dn
			self.Alpha_Tri[i, i+1] = TriI2[i] * IvDnDp1 + TriI1[i] * Ivsq2Dn
		
		factor = (-1.0) * self._LOGee * self.kpre * self.Dt
		self.Alpha_Tri *=  factor        																	# (-1) because the matrix moves to the LHS in the equation
		#print "TriI1"
		#print TriI1
		#print "TriI2"
		#print TriI2
		
		#print self.Alpha_Tri
		
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
		f00 = ( fin[0] - ( (vr[0]*vr[0])/(vr[1]*vr[1]) ) *fin[1] ) / (1.0 - (vr[0]*vr[0])/(vr[1]*vr[1]));
		ddf0[0] = 2.0 * (fin[1] - f00) / (vr[1]*vr[1]);
		df0[0] = ddf0[0] * vr[0];		
		# Calculate 1/(2v)*(d^2f)/(dv^2),  1/v^2*df/dv
		for n in xrange(0,  self.num_pr_cells-1):
			df0[n]  /= vr[n]*vr[n]
			ddf0[n] /= 2.0  *vr[n]
		
		#print "df0"
		#print df0
		#print "ddf0"
		#print ddf0
		
	# this update directly into the |p| array 'fin'
	#@timeit
	def implicit_advance(self, fin, el):
		Alpha = self.Alpha; vr = self.vr; df0 = self.df0; kpre = self.kpre; Scattering_Term = self.Scattering_Term
		ddf0 = self.ddf0
		df0 = self.df0
		_LOGee = self._LOGee
		Dt = self.Dt
		
		Alpha[:,:] = self.Alpha_Tri[:,:]
		
		#print "Alpha 0,0 %e" % (Alpha[0,0], )
		# ZEROTH CELL FOR TRIDIAGONAL ARRAY
		if el > 1:
			Alpha[0,0] = 0.0;
		
		if self.if_tridiagonal:
			fokker_plank__implicit__advance(fin, el, _LOGee, Dt, 1)
		else:
			fokker_plank__implicit__advance(fin, el, _LOGee, Dt, 0)
		"""
		print "oncomng!!!"
		for i in xrange(0,  self.num_pr_cells):
			print "%e " % Alpha[i,i] ,
		print "outgoing"
		if not self.if_tridiagonal:
			#  CONSTRUCT COEFFICIENTS and then full array
			LL = el
			A1 = (LL+1.0)*(LL+2.0) / ((2.0*LL+1.0)*(2.0*LL+3.0))         
			A2 = (-1.0) *(LL-1.0)* LL      / ((2.0*LL+1.0)*(2.0*LL-1.0))
			B1 = (-1.0) *( 0.5 *LL*(LL+1.0) +(LL+1.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0))
			B2 = (       (-0.5)*LL*(LL+1.0) +(LL+2.0) ) / ((2.0*LL+1.0)*(2.0*LL+3.0))
			B3 = ( 0.5 *LL*(LL+1.0) +(LL-1.0) ) / ((2.0*LL+1.0)*(2.0*LL-1.0))
			B4 = ( 0.5 *LL*(LL+1.0) - LL      ) / ((2.0*LL+1.0)*(2.0*LL-1.0))
			#print "A1:    %e" % A1
			#print "A2:    %e" % A2
			#print "B1:    %e" % B1
			#print "B2:    %e" % B2
			#print "B3:    %e" % B3
			#print "B4:    %e" % B4
			#print
			
			
			for i in xrange(0,  self.num_pr_cells-1):
				t1 = A1*ddf0[i] + B1*self.df0[i]
				t1 *= (-1.0) * _LOGee * kpre * Dt
				t2 = A1*ddf0[i] + B2*df0[i]
				t2 *= (-1.0) * _LOGee * kpre * Dt
				t3 = A2*ddf0[i] + B3*df0[i]
				t3 *= (-1.0) * _LOGee * kpre * Dt
				t4 = A2*ddf0[i] + B4*df0[i]
				t4 *= (-1.0) * _LOGee * kpre * Dt
				#if( i < 5):
				#	print "t's-------"
				#	print t1
				#	print t2
				#	print t3
				#	print t4
					
				Alpha[i,0] += t1 * ( 2.0*np.pi * np.power(vr[0]/vr[i],el+2)*vr[0]*vr[0]*(vr[1]-vr[0]) )
				Alpha[i,0] += t3 * ( 2.0*np.pi * np.power(vr[0]/vr[i],el)  *vr[0]* vr[0]*(vr[1]-vr[0]) )
				
				#if(i==0):
				#	print "alook1"
				#	print Alpha[i,0]
				
				for j in xrange(1, i):
					Alpha[i,j] += t1 * ( 2.0*np.pi*np.power(vr[j]/vr[i],el+2)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
					Alpha[i,j] += t3 * ( 2.0*np.pi*np.power(vr[j]/vr[i],el)  *vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
				#if(i==0):
				#	print "alook1.5"
				#	print Alpha[i,i]
				#	print "t1: %e" % t1
				#	print "t3: %e" % t3
				#	print "vr[i]: %e" % vr[i]
				#	print "f1: %e" % ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) )
				#	print "diff11: %e" % (vr[i]-vr[i-1])
					
				Alpha[i,i] += t1 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) )
				Alpha[i,i] += t3 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i]-vr[i-1]) )
				#if(i==0):
				#	print "alook2"
				#	print Alpha[i,i]
				
				Alpha[i,i] += t2 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i+1]-vr[i]) )
				Alpha[i,i] += t4 * ( 2.0*np.pi *vr[i]*vr[i]*(vr[i+1]-vr[i]) )
				#if(i==0):
				#	print "alook3"
				#	print Alpha[i,i]
				
				for j in xrange(i+1, self.num_pr_cells-1):
					Alpha[i,j] += t2 * ( 2.0*np.pi*np.power(vr[j]/vr[i],-el-1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
					Alpha[i,j] += t4 * ( 2.0*np.pi*np.power(vr[j]/vr[i],-el+1)*vr[j]*vr[j]*(vr[j+1]-vr[j-1]) )
			pass
		
		
		#print "Alpha 0,0 %e" % (Alpha[0,0], )
		#print "Alpha"
		#print Alpha
		# INCLUDE SCATTERING TERM
		ll1 = float(el)
		ll1 *= (-0.5)*(ll1 + 1.0);
		print el
		print ll1
		for i in xrange(0,  self.num_pr_cells):
			Alpha[i,i] += 1.0 - ll1 * Scattering_Term[i]
			print "%e " % Alpha[i,i] ,
		exit(-1)
		"""
		#print "Alpha2"
		#alpha_out = []
		#for ip in xrange(0,  self.num_pr_cells-1):
		#	alpha_out.append(Alpha[ip,ip-1])
		#	alpha_out.append(Alpha[ip,ip])
		#	alpha_out.append(Alpha[ip,ip+1])
		#	#print "%e,%e,%e" % ( Alpha[ip,ip-1], Alpha[ip,ip],Alpha[ip,ip+1] )
		#print Alpha[self.num_pr_cells-1, self.num_pr_cells-2]
		#print Alpha[self.num_pr_cells-1,self.num_pr_cells-1]
		#print "%e,%e,%e" % ( Alpha[self.num_pr_cells-1, self.num_pr_cells-2], Alpha[self.num_pr_cells-1,self.num_pr_cells-1],0.0 )
		#alpha_out.append(Alpha[self.num_pr_cells-1, self.num_pr_cells-2])
		#alpha_out.append(Alpha[self.num_pr_cells-1, self.num_pr_cells-1])
		#alpha_out.append(0.0);
		#np.save('alpha1.npy', Alpha)
		
		if ( self.if_tridiagonal ):
			#is_tri_diagonal_ok = Thomas_Tridiagonal(Alpha, fin, fout)
			#if is_tri_diagonal_ok == False:
			#	print  "WARNING: Matrix is not diagonally dominant!!!"
			#print "TRI"
			temp = self.solve_tridiagonal_matrix( Alpha, fin )
			#np.save('fout1.npy', temp)
			#np.save('fin1.npy', fin)
			fin[:] = temp[:]
		else:
			# fully solve A
			#self.invert_matrix(Alpha, fin, fout)
			#fin[:] = fout[:]
			#print "FULL!"
			temp = self.solve_matrix( Alpha, fin )
			fin[:] = temp[:]
		#pass
		#print "Fout!!!"
		#print temp
		pass
	
	def solve_tridiagonal_matrix(self, matrix, right_hand_side):
		#ud = np.insert(np.diag(A,1), 0, 0) 					# upper diagonal
		#d  = np.diag(A) 											# main diagonal
		#ld = np.insert(np.diag(A,-1), len(d)-1, 0)		# lower diagonal
		size = self.tridiagonal_solver_temp.shape[1]
		self.tridiagonal_solver_temp[0,1:] = np.diag(matrix, 1)
		self.tridiagonal_solver_temp[1,:]  = np.diag(matrix, 0)
		shit = np.diag(matrix, -1)
		self.tridiagonal_solver_temp[2,0:(size-1)]  = np.diag(matrix, -1)
		self.tridiagonal_solver_temp[0,0] = 0.0
		self.tridiagonal_solver_temp[0,size-1] = 0.0
		return scipy.linalg.solve_banded( (1,1), self.tridiagonal_solver_temp, right_hand_side)
	#@timeit
	def solve_matrix( self, matrix, right_hand_side):
		return scipy.linalg.solve( matrix, right_hand_side )
	
	# result put inplace into Yin
	def f1_loop(self, Yin, Dt, max_m):
		# loop over spaitial cells....
		f00 = self.F0;  fc = self.F1;
		for x in xrange(0, self.num_plain_cells[0]):
			# reset the integrals an ceofficients..
			# reset_coeff does not alter the input.
			self.reset_coeff( Yin.cells[x + self.NB].harmonics[0].momentum_projection, Dt)
			
			# loop over harmonics
			#for m in xrange(0, max_m):
			#	self.implicit_advance(Yin.get_harmonic_at_location(x + self.NB, 1, m).momentum_projection, 1)
			self.implicit_advance(Yin.get_harmonic_at_location(x + self.NB, 1).momentum_projection, 1)
		pass
	
	def flm_loop(self, Yin, Dt, max_m=None):
		if max_m == None:
			max_m = self.num_m
		
		for x in xrange(0, self.num_plain_cells[0]):
			# reset_coeff does not alter the input.
			# TODO: Assume you did this when you called f1_loop ???? Yes or no?
			self.reset_coeff( Yin.cells[x + self.NB].harmonics[0].momentum_projection, Dt)
			
			# loop over harmonics for this spatial point.
			for l in xrange(2, self.num_l):
				self.implicit_advance( Yin.get_harmonic_at_location(x + self.NB, l).momentum_projection, l)
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
	def euler_backward_solver_for_fp_advance_1(self, Yin, Dt):
			if self.if_implicit1D:
				# if 1D, do calcualtion just for f10 (f11 not needed)
				self.f1_loop(Yin, Dt, 1)
			else:
				# Do calcualtion just for both f10 and f11
				self.f1_loop(Yin, Dt, 2)
	
	def euler_backward_solver_for_fp_advance_lm(self, Yin, Dt):
			if self.if_implicit1D:
				self.flm_loop(Yin, Dt, 1)
			else:
				self.flm_loop(Yin, Dt)
	
	# main entry point for f00 collsions.. results will be written
	#	directly into the input distribution function.
	#@timeit
	def calc_f00(self, Yin, species, state, CFG, n_step, time):
		eval_func = self.lowest_harmonic_fp_slope
		#print "self.NB %d with %d subcycle steps (each with dt=%e)" %  (self.NB, self.num_subcycling_steps, self.h)
		
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
	
	# main entry point for higher harmonic collsions.. results will be written
	#	directly into the input distribution function.
	#@timeit
	def calc_flm(self, Yin, species, state, CFG):
		# this does the integration for f10, f11
		self.euler_backward_solver_for_fp_advance_1(Yin, CFG.dt)
		# this does the remaing higher harmonics..
		self.euler_backward_solver_for_fp_advance_lm(Yin, CFG.dt)
	
class UpdateBoundryCells:
	def __init__(self, state, CFG):
		pass
	
	def calc(self, state, CFG):
		
		# so let's update spitial cells in the momenum distribution...
		for i in xrange(0, len(state.species)):
			
			F = state.species[i].F
			num_boundry_cells_left = F.num_boundry_cells[0][0]
			num_boundry_cells_right = F.num_boundry_cells[0][1]
			
			# fix up the right side..
			base_index_right = F.num_total_cells[0]-num_boundry_cells_right
			base_index_left = num_boundry_cells_left
			num_harmonics = F.cells[0].num_harmonics
			
			for bc in xrange(0,num_boundry_cells_right):
				cell_src = F.cells[ base_index_left  + bc]
				cell_dst = F.cells[ base_index_right + bc]					
				for sh in xrange(0,num_harmonics):
					cell_dst.harmonics[sh].momentum_projection[:] = cell_src.harmonics[sh].momentum_projection[:]
			
			# fix up the left side..
			base_index_right = F.num_total_cells[0]-num_boundry_cells_right-num_boundry_cells_left 
			base_index_left  = 0
			for bc in xrange(0,num_boundry_cells_left):
				cell_src = F.cells[ base_index_right + bc]
				cell_dst = F.cells[ base_index_left + bc]	
				for sh in xrange(0,num_harmonics):
					cell_dst.harmonics[sh].momentum_projection[:] = cell_src.harmonics[sh].momentum_projection[:]
			pass
		
		# now let's update the fields...
		data = state.E_field.components[0]
		num_cells = state.E_field.num_total_cells[0]
		num_boundry_cells_left = F.num_boundry_cells[0][0]
		num_boundry_cells_right = F.num_boundry_cells[0][1]
		
		# fix up right...
		src_start_index = num_boundry_cells_left
		src_end_index   = num_boundry_cells_left+num_boundry_cells_right
		dst_start_index = num_cells - num_boundry_cells_right
		dst_end_index = num_cells
		data[ dst_start_index:dst_end_index ] = data[ src_start_index:src_end_index ]
		
		# fix up left
		src_start_index = state.E_field.num_total_cells[0] - num_boundry_cells_right - num_boundry_cells_left
		src_end_index   = state.E_field.num_total_cells[0] - num_boundry_cells_right
		dst_start_index = 0
		dst_end_index = num_boundry_cells_left
		data[ dst_start_index:dst_end_index ] = data[ src_start_index:src_end_index ]
	
		# fix external field too, if present...
		if STATE.use_external_E_field:
			self.reconcile_field( STATE.E_field__external.components[0], state, num_cells, num_boundry_cells_left, num_boundry_cells_right)
		
	def reconcile_field(self, data, state, num_cells, num_boundry_cells_left, num_boundry_cells_right):
		# fix up right...
		src_start_index = num_boundry_cells_left
		src_end_index   = num_boundry_cells_left+num_boundry_cells_right
		dst_start_index = num_cells - num_boundry_cells_right
		dst_end_index = num_cells
		data[ dst_start_index:dst_end_index ] = data[ src_start_index:src_end_index ]
		
		# fix up left
		src_start_index = state.E_field.num_total_cells[0] - num_boundry_cells_right - num_boundry_cells_left
		src_end_index   = state.E_field.num_total_cells[0] - num_boundry_cells_right
		dst_start_index = 0
		dst_end_index = num_boundry_cells_left
		data[ dst_start_index:dst_end_index ] = data[ src_start_index:src_end_index ]
		
		
		
	

class RK3:
	def __init__(self, step_size):
		self.h = step_size
		self.E0 = STATE.E_field.clone()
		self.Eh = STATE.E_field.clone()
		self.time_step = 0
		
	def debug_me(self):
		return False
	
	def step(self, F, spec, state, CFG):
		
		LL = 1
		MM = 0
		
		Y0 = spec.Y0
		Yh = spec.Yh
		Y  = spec.F				# the 'real, final' state.
		E = STATE.E_field		# the 'real, final' Electric field.
		
		Y0.copy(Y);															self.E0.copy (E)
		
		# account for an external E-Field is present...
		#print STATE.use_external_E_field
		if STATE.use_external_E_field:
			#print "------"
			#print STATE.E_field.components[0]
			self.E0.add( STATE.E_field__external)
			pass
			
		
		if self.debug_me():
			print "Step1 start: Y0"; DEBUG(YYY=Y0, L=LL, M=MM)
		
		# step 1
		# Yh = F(Y0)
		F.sub_time_step = 0
		F(Y0, Yh, spec, state, CFG, Ein=self.E0, Eout=self.Eh)
		if self.debug_me():
			print "Step1 end rigth after: Yh"; DEBUG(YYY=Yh, L=LL, M=MM)
			print "Step1 end rigth after: Y0"; DEBUG(YYY=Y0, L=LL, M=MM)
			print "Step1 end rigth after: Eh"; DEBUG_E(self.Eh)
			
		# Y0 = Y0 + h*Yh
		Yh.scalar_mult(self.h);											self.Eh.scalar_mult(self.h);	
		Y0.add(Yh);															self.E0.add(self.Eh);
		
		if self.debug_me():
			print "Step2 start: Yh"; DEBUG(YYY=Yh, L=LL, M=MM)
			print "Step2 start: Y0"; DEBUG(YYY=Y0, L=LL, M=MM)
			print "Step2 start: Eh"; DEBUG_E(self.Eh)
        
		#step 2
		# Yh = F(Y0)
		F.sub_time_step = 1
		F(Y0, Yh, spec, state, CFG, Ein=self.E0, Eout=self.Eh)
		# Y0 = 1/4 * ( 3*Y + (Y0 + h*Yh) )
		Yh.scalar_mult(self.h);											self.Eh.scalar_mult(self.h);	
		Y0.add(Yh);															self.E0.add(self.Eh);
		Y0.madd(3.0, Y);													self.E0.madd(3.0, E)
		Y0.scalar_mult(0.25);											self.E0.scalar_mult(0.25);
		
		if self.debug_me():
			print "Step3		start: Y0"; DEBUG(YYY=Y0, L=LL, M=MM)
			print "Step3		start: E0";  DEBUG_E(self.E0)
			
		# step 3
		# Yh = F(Y0)
		F.sub_time_step = 2
		F(Y0, Yh, spec, state, CFG, Ein=self.E0, Eout=self.Eh)
		# Y  = 1/3 * ( Y + 2 * (Y0+h*Yh) )
		
		Yh.scalar_mult(self.h);											self.Eh.scalar_mult(self.h);	
		Y0.add(Yh);															self.E0.add(self.Eh);
		Y0.scalar_mult(2.0);												self.E0.scalar_mult(2.0);
		Y.add(Y0);															E.add(self.E0);
		Y.scalar_mult(1.0/3.0);											E.scalar_mult(1.0/3.0);	
		
		
		if self.debug_me():
			print "AT END OF RK3, DIST is (1,0):"; DEBUG(YYY=Y, L=LL, M=MM);
			print "AT END OF RK3, Ex is:"; DEBUG_E(E);
			if F.time_step == 1:
				DEBUG_BAIL()
		
		F.time_step = F.time_step + 1
		self.time_step += 1

class F_explicit:
	def __init__(self, species, state, CFG):
		self.e_on_dist 					=		ElectricFieldOnDistribution1d	(species, state, CFG)
		self.current_solver 				= 		CurrentSolver						(species, state, CFG)
		self.spatial_advection_solver = 		SpatialAdvection1d				(species, state, CFG)
		#self.moving_frame					=		VlasovMovingFrameCalculation  (species, state, CFG)
		
		# be paranoid like Micheal!
		# if  (l+m-1 > #p) then #p = 0.0
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
			#print "%d,= %d" % (l, threshold)
			if(threshold > 0):
				self.trim_list.append( (l, 0, Yh.linear_index_for_harmonic( l), threshold) )
		
		if self.do_filter:
			if mpi_info.rank == 0:
				print "Trimming filter will occur!!!"
	
	def do_db(self):
		return False
		if self.time_step == 0 and self.sub_time_step == 0:
			return True
		return False
	
	def __call__(self, Y0, Yh, species, state, CFG, Ein=None, Eout=None):
		if self.do_db():
			print "\t\t****************START EVAL Timestep %d,%d" % (self.time_step, self.sub_time_step)
		
		Yh.scalar_mult(0.0)
		Eout.scalar_mult(0.0)
		LL = 1
		MM = 0
		
		if self.do_db():
			print "\t\tIN EVAL: Yin In at start!!"; DEBUG(YYY=Y0, L=LL, M=MM, indent=True);
			print "\t\tIN EVAL:  Ein IN at start!!"; DEBUG_E(Ein, indent=True);
		
		self.e_on_dist.calc								(Y0, Yh, species, state, CFG, Ein = Ein)
		
		if self.do_db():
			print "\t\tIN EVAL: Yh AFTER E effect!!"; DEBUG(YYY=Yh, L=LL, M=MM, indent=True);
		
		"""
		if DaWaveDrivah != None:
			if STATE.time >= (DaWaveDrivah.rise_time + DaWaveDrivah.sustain):
				self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
			else:
				if mpi_info.rank == 0:
					print "Field solver skipped."
				
		else:
			self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
		"""
		#if DaWaveDrivah.driver_on == False:
		self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
		
		debug_avection = False
		if self.do_db():
			print "\t\tIN EVAL: Yh EX AFTERJ!!"; DEBUG_E(Eout, indent=True);
			debug_avection = False
		
		
		self.spatial_advection_solver.calc			(Y0, Yh, species, state, CFG, debug = debug_avection)
		
		#self.moving_frame.calc			(Y0, Yh, species, state, CFG, debug = debug_avection)
		
		if self.do_db():
			print "\t\tIN EVAL: Yh AFTER Avectection!!"; DEBUG(YYY=Yh, L=LL, M=MM, indent=True);
		
		if self.do_filter:
			for item in self.trim_list:
				idx = item[2]
				level_to_trim = item[3]
				for ix in xrange(0, Yh.num_total_cells[0]):
					Yh.cells[ix].harmonics[idx].momentum_projection[:level_to_trim] = 0.0
				#print "zeros %d cells for %d,%d" % (  level_to_trim, item[0], item[1] )
			pass
		
		
		
		
		if self.do_db():
			print "\t\t****************END EVAL Timestep %d,%d" % (self.time_step, self.sub_time_step)
		
		#if self.time_step == 1 and self.sub_time_step == 2:
		#	DEBUG_BAIL()
		


plotme = None
plotme1 = None
moments = []
DaWaveDrivah = None


def main_loop():
	global DaWaveDrivah
	# setup the calculation..
	# TODO: move this, of course.
	time = 0.0
	n = 0
	STATE.n = n
	STATE.time = time
	
	callback_mode = False
	
	#setup the simulation notes output file....
	if mpi_info.rank == 0:
		f = open('!sim_info.txt','w')
		f.write("\n\n")
		f.close()
	
	DaWaveDrivah = None
	try:
		if STATE.config.use_wave_driver:
			DaWaveDrivah = WaveDrivah(STATE.config.dt, STATE.species[0], 1.0, num_cells_wave_travels_in_dt=1.0)
			exit(-2)
	except:
		pass
	
	# set the dt....
	#STATE.config.dt = 2.0*np.pi / DaWaveDrivah.w_wave / 64.0
	#STATE.config.dt = 2*np.pi*30.0/ DaWaveDrivah.w_wave / 192.0
	#STATE.config.dt = 2*np.pi / 192.0 /  DaWaveDrivah.w_wave * 4.0

	
	
	# reconcile the boundry cells for the T=0 evaluatation.
	update_boundries = UpdateBoundryCells(STATE, STATE.config)
	if mpi_info.num_nodes == 1:
		update_boundries.calc(STATE, STATE.config)
	else:
		mpi_update.update(STATE, STATE.config)
	moments[0].notify_simulation_time(n, time)
	
	"""
	if STATE.config.E_field_profile != None and False:
		
		solva = EFieldSolve1d()
		print "Gonna solve!"
		exit(-1)
		solva.solve( STATE.E_field )
		print "Done solvaing"
		# redo boundry conditions...
		if mpi_info.num_nodes == 1:
			update_boundries.calc(STATE, STATE.config)
		else:
			mpi_update.update(STATE, STATE.config)		
	"""
	
	dt = STATE.config.dt;
	if DaWaveDrivah != None:
		if DaWaveDrivah.pure_harmonic_dt > 0.0:
			dt = DaWaveDrivah.pure_harmonic_dt
	
	if mpi_info.rank == 0:
		print "Delta Time %f" % dt
	
	t_end = STATE.config.t_end;
	num_species = STATE.config.num_species;

	solva = RK3(dt)
	
	f_slice_axis = zeros( 2* STATE.species[0].F.momentum_axis.num_points)
	for ip in xrange(0, STATE.species[0].F.momentum_axis.num_points):
		pval=  STATE.species[0].F.momentum_axis.values[ip]
		f_slice_axis[STATE.species[0].F.momentum_axis.num_points-1-ip] = -1.0* pval
		f_slice_axis[STATE.species[0].F.momentum_axis.num_points+ip ] = pval 
	actual_dist_functions = []
	eval_functions = []
	fokker_planks = []
	
	
	for i in range(0, num_species):
		eval_functions.append( F_explicit(STATE.species[i], STATE, STATE.config) )
		
		# play with frame...
		#TODO: TAKE OUT
		#STATE.species[i].moving_frame_u_field[:] = DaWaveDrivah.v_phase
		#STATE.species[i].moving_frame_u_field[:] = DaWaveDrivah.v_phase
		#STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset)
		
		
		#fred = EFieldSinusodialProfile( 0.02, 	num_wavelengths = 3.0 )
		#fred.apply(STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset,  buffah= STATE.species[i].moving_frame_u_field)
		
		#mpi_info.mpi_axes[0].reconcile_boundries_field( STATE.species[i].moving_frame_u_field.components[0] )
		#print STATE.species[i].moving_frame_u_field.components[0]
		
		if DaWaveDrivah != None:
			STATE.species[i].moving_frame_u_field.components[0][:] = DaWaveDrivah.v_phase
		#if mpi_info.rank > 0:
		#	print "------------ Fludid velocity---------------------"
		#	print STATE.species[i].moving_frame_u_field.components[0][:]
		
		
		fokker_planks.append( FokkerPlankExplicit	(STATE.species[i], STATE, STATE.config) )
		actual_dist_functions.append( ActualDistrubtionFunction(STATE.species[i].F) )
	
	if test_mode:
		testing1_main(eval_functions)
	
	#------------------ inital F perturbation and making E consistent...
	init_f0_density = None
	init_f1_density = None
	Esolva = EFieldSolve1d()
	if STATE.config.E_field_profile != None and True:
	
		# use moments to get density..
		if( mpi_gather.gather(  moments[0].get_density(), remove_guard_cells=True ) ):
			init_f0_density = np.copy ( mpi_gather.stiching_buffer )
		# save inital F
		local_f0__initial = Esolva.get_charge_density_for_species(0)
		# do perturbation
		STATE.config.E_field_profile.apply_perturbation( STATE.species[0].F, STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset )
		# record new F and subtract form unperturbed F to get 'f1'
		local_f1__initial = Esolva.get_charge_density_for_species(0)
		local_f1__initial -= local_f0__initial
		
		# now try to update E to reflect this perturbation...
		# It already does.. just one change... the perturbation is electrons.. so high electron density is large, then E should be negitive..
		#	(ie if init_f1_density > 0 then E should be negitive and if init_f1_density<0 then E should be positive
		#STATE.E_field.components[0][:] *= -1.0
		
		# args are: Yh, species, state, CFG, Eout = None
		#STATE.E_field.components[0][:] = 0.0
		#eval_functions[0].current_solver.calc( STATE.species[0].F, None, None, None, None, Eout = STATE.E_field )
		
		number_field_data_cells = STATE.E_field.num_plain_cells[0]
		num_boundry_cells_left  = STATE.E_field.num_boundry_cells[0][0]
		num_boundry_cells_right = STATE.E_field.num_boundry_cells[0][1]
		num_cells = STATE.spatial_axes__global[0].num_points
		box_size = 	STATE.spatial_axes__global[0].max
		#print "box size: \t\t %f" % float(box_size)
		#print "num_cells: \t\t %f" % float(num_cells)
		print "delta_x: \t\t %f (offical from axis: %f)" % ( (float(box_size) / float(num_cells)), STATE.spatial_axes__global[0].dx()  )
		delta_x = float(box_size) / ( float(num_cells) + 1.0 )
		delta_x_2 = ( float(box_size) + STATE.spatial_axes__global[0].dx() ) / float ( num_cells )
		print "adjusted dx  \t\t %f" % delta_x
		print "adjusted dx2 \t\t %f" % delta_x_2
		recieve_buffer_for_E_field = np.zeros (  STATE.spatial_axes__global[0].values.shape )
		
		# use Moments to solve for E
		moments[0].current_step -= 1
		if(True and mpi_gather.gather(  moments[0].get_density(), remove_guard_cells=True ) ):
			init_f1_density = np.copy ( mpi_gather.stiching_buffer )
			init_f1_density -= init_f0_density
			init_f1_density *= -1.0
			
			init_f1_density__integrate_buffer = np.zeros(  (STATE.spatial_axes__global[0].num_points + 2 ) )
			#init_f1_density__integrate_buffer[1:] =  init_f1_density[:]
			#init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-1]
			init_f1_density__integrate_buffer[0:-2] =  init_f1_density[:]
			init_f1_density__integrate_buffer[-2] =  init_f1_density__integrate_buffer[0]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			
			#init_f1_density__integrate_buffer[1:-1] =  init_f1_density[:]
			init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-2]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			#print init_f1_density__integrate_buffer
			
			# we fed in a vector of length N+1 and will get out a vector of length N
			temp = scipy.integrate.cumtrapz( init_f1_density__integrate_buffer, dx=( (float(box_size) / float(num_cells) ) ) )
			
			recieve_buffer_for_E_field[0] = temp[-2]
			recieve_buffer_for_E_field[1:] = temp[0:-2]
			print "MEAN!!!!!!! %f" % ( np.mean( recieve_buffer_for_E_field ) )
			
			recieve_buffer_for_E_field *= (4.0*np.pi)
			clf(); plot( recieve_buffer_for_E_field ); savefig('ass_cunt_E__moments')
			
			"""
			# basiclly works.. just need to hunt down normalization factors!
			# and try an FFT!!!
			dx = STATE.spatial_axes__global[0].dx()
			k_axis = np.arange(0.0, float(num_cells) ) * dx
			k_axis *= np.pi
			k_axis *= ( 2.0*np.pi / ( float(box_size) +  dx ) )
			fft_n = fft.fft( init_f1_density )
			
			symetic_k_axis = np.arange( float(-1*num_cells/2), float(num_cells/2) )
			symetic_k_axis *= np.pi
			symetic_k_axis *= dx
			symetic_k_axis /= ( 2.0*np.pi / ( float(box_size) +  dx ) )
			
			
			fft_n = fftshift( fft_n )
			
			for ik in xrange( 1 , num_cells ):
				if ik == (num_cells/2):
					print " skipped %f" % symetic_k_axis[ ik ]
					continue
				
				fft_n[ik] /= np.complex( 0, -1.0* symetic_k_axis[ ik ] )
				#fft_n[ik] /= np.complex( 0, -1.0* k_axis[ ik ] )
			
			E = fft.ifft( ifftshift( fft_n ) )
			E1 = E / dx
			E2 = E / np.sqrt(dx)
			print E1
			print E2
			"""
		moments[0].current_step += 1
			
			
		
		# This gets E from local densoty...
		if( False and  mpi_gather.gather(  local_f1__initial, remove_guard_cells=True ) ):
			# only rank==0 can get in here.. and we have all 'local_f1__initial' buffers stiched togther!
			init_f1_density__integrate_buffer = np.zeros(  (STATE.spatial_axes__global[0].num_points + 2 ) )
			#init_f1_density__integrate_buffer[1:] =  mpi_gather.stiching_buffer[:]
			#init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-1]
			init_f1_density__integrate_buffer[0:-2] =  mpi_gather.stiching_buffer[:]
			init_f1_density__integrate_buffer[-2] =  init_f1_density__integrate_buffer[0]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			
			#init_f1_density__integrate_buffer[1:-1] =  mpi_gather.stiching_buffer[:]
			init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-2]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			#print init_f1_density__integrate_buffer
			
			# we fed in a vector of length N+1 and will get out a vector of length N
			temp = scipy.integrate.cumtrapz( init_f1_density__integrate_buffer, dx=( delta_x_2 ) )
			
			
			recieve_buffer_for_E_field[0] = temp[-2]
			recieve_buffer_for_E_field[1:] = temp[0:-2]
			#print "MEAN!!!!!!! %f" % ( np.mean( recieve_buffer_for_E_field ) )
			
			recieve_buffer_for_E_field *= (4.0*np.pi)
			clf(); plot( recieve_buffer_for_E_field ); savefig('ass_cunt_E')
			#print recieve_buffer_for_E_field
			#print "MEAN final!!!!!!! %f" % ( np.mean( recieve_buffer_for_E_field ) )
			#exit(-1)
		data_start_index = mpi_info.mpi_axes[0].global_cell_offset[0]
		number_valid_data_cells = STATE.E_field.num_plain_cells[0]
		data_end_index = data_start_index + number_valid_data_cells
		
		#print "broadcast the integrated E field!"
		recieve_buffer_for_E_field = mpi_info.comm.bcast( recieve_buffer_for_E_field, root=0)
		#print "USE data from others to get what E is in our Local neighborhood"		
		appropiate_data = recieve_buffer_for_E_field[data_start_index:data_end_index ]
		STATE.E_field.components[0][ num_boundry_cells_left:num_boundry_cells_left+number_valid_data_cells] = appropiate_data[:]

		#print "NOw updates guard cells."
		# update everyones's gaurd cells...
		if mpi_info.num_nodes == 1:
			update_boundries.calc(STATE, STATE.config)
		else:
			mpi_update.update(STATE, STATE.config)		
		#print "Everyone's guard cells are good now!"
		
		
		
		# HACK BUG BAD FIX ADAM
		
		#if (mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True) ):
		#	my_E =  np.copy ( mpi_gather.stiching_buffer )
		#	print "ROOT GOT the full E from everyone."
		#
		#mpi_gather.stiching_buffer[:] = 666.0
		"""
		moments[0].current_step -= 1
		if( mpi_gather.gather(  moments[0].get_density(), remove_guard_cells=True ) ):
			init_f1_density = np.copy ( mpi_gather.stiching_buffer )
			init_f1_density -= init_f0_density
			
			#clf()
			#plot( my_E )
			#savefig("crap_start__E.png")
			clf()
			plot( init_f1_density )
			savefig("crap_start__f1.png")
			
			init_f1_density__integrate_buffer = np.zeros(  mpi_gather.stiching_buffer.shape[0] + 1)
			init_f1_density__integrate_buffer[0:-1] =  mpi_gather.stiching_buffer[:]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[0]
			
			# we fed in a vector of length N+1 and will get out a vector of length N
			recieve_buffer_for_E_field[:] = scipy.integrate.cumtrapz( init_f1_density__integrate_buffer, dx=(STATE.spatial_axes__global[0].dx() ) )
			recieve_buffer_for_E_field *= (4.0*np.pi)
			
			recieve_buffer_for_E_field[:] = init_f1_density[:]
		"""
		"""
		# BAD HACK TAKE OUT
		moments[0].current_step += 1
		print "Past ugly Adam moment hack"
		
		# e_field_integration  now has the full empirically integrated number density for all real cells along X (no guard cells)
		#		TODO: make faster
		#		As first pass (and since this is used in only 2 time steps of a simuation), send all nodes a copy of this full empirical E-field.
		
		#if mpi_info.rank > 0:
		#	recieve_buffer_for_E_field = None
		recieve_buffer_for_E_field = mpi_info.comm.bcast( recieve_buffer_for_E_field, root=0) 
		
		#print "Rank %d: %s" % ( mpi_info.rank, str( recieve_buffer_for_E_field ) )
		#mpi_info.comm.Barrier()	
		#print "after scatter"
		
		# now all modes have this E-field buffer. Uses it to update your local E-field piece.
		num_boundry_cells_left  = STATE.E_field.num_boundry_cells[0][0]
		num_boundry_cells_right = STATE.E_field.num_boundry_cells[0][1]
		
		data_start_index = mpi_info.mpi_axes[0].global_cell_offset[0]
		number_valid_data_cells = STATE.E_field.num_plain_cells[0]
		data_end_index = data_start_index + number_valid_data_cells
		
		print "USE data from others to get what E is in our Local neighborhood"
		
		appropiate_data = recieve_buffer_for_E_field[data_start_index:data_end_index ]
		STATE.E_field.components[0][ num_boundry_cells_left:num_boundry_cells_left+number_valid_data_cells] = appropiate_data[:]
		"""
	"""
	print "RANK %d is wating at block" % mpi_info.rank
	mpi_info.comm.Barrier()	

	print "RANK %d over block" % mpi_info.rank
	print "NOw updates guard cells."
	# update everyones's gaurd cells...
	if mpi_info.num_nodes == 1:
		update_boundries.calc(STATE, STATE.config)
	else:
		mpi_update.update(STATE, STATE.config)		
	
	print "Ueveryone's guardc ells are good now!"
	
	print STATE.E_field.components[0]
	"""
	"""
	clf()
	plot( init_f0_density )
	savefig( "f0__initial", dpi=200)
	
	clf()
	#plot( ( init_f1_density) )
	plot( ( my_E - init_f1_density ) )
	savefig( "f1__initial", dpi=200)
	
	x_data = np.arange(0.0, 96.0); x_data *= (10.0/96.0);
	#e_int = scipy.integrate.cumtrapz( local_f1__initial, dx = (4.0*np.pi*10.0/96.0) )
	e_int = scipy.integrate.cumtrapz( init_f1_density, dx=(10.0/96.0)) #x=x_data) # = (4.0*np.pi*10.0/96.0) )
	e_int *= 4.0*np.pi
	
	clf()
	#plot( ( init_f1_density) )
	plot( e_int )
	plot( init_f1_density )
	title("ffff")
	local_f1__initial = Esolva.get_charge_density_for_species(0)
	local_f1__initial -= local_f0__initial
	plot( local_f1__initial[3:-3] )
	savefig( "f1__initial_int", dpi=200)
	
	print " ave of perurb: "
	print np.mean ( init_f1_density)
	print " ave of E: "
	print np.mean ( my_E )

	print local_f1__initial
	print init_f1_density[0:local_f1__initial.shape[0]]
	"""
	
	
	"""
	# let all nodes get local density
	local_f1__initial = Esolva.get_charge_density_for_species(0)
	print local_f1__initial.shape[0]
	
	local_f1__initial -= local_f0__initial
	num_boundry_cells_left  = STATE.E_field.num_boundry_cells[0][0]
	num_boundry_cells_right = STATE.E_field.num_boundry_cells[0][1]
	
	recieve_buffer = np.zeros (  STATE.spatial_axes__global.values.shape )
	mpi_info.comm.Alltoallv( local_f1__initial[ num_boundry_cells:-num_boundry_cells_right ], recieve_buffer )
	
	int_e = scipy.integrate.cumtrapz( local_f1__initial, dx=(STATE.spatial_axes__global[0].dx() ), initial = 0.0 )
	print int_e.shape[0]
	int_e *= 4.0*np.pi
	STATE.E_field.components[0][:] = int_e[:]
	if mpi_info.num_nodes == 1:
		update_boundries.calc(STATE, STATE.config)
	else:
		mpi_update.update(STATE, STATE.config)
	
	
	
	cumultive_factor = STATE.E_field.components[0][ num_boundry_cells - 1 ]
	STATE.E_field.components[0][ num_boundry_cells:-num_boundry_cells_right] +=
	#
	#
	#exit(-1)
	
	
	#STATE.E_field.components[0][:] = local_f1__initial[:]
	#STATE.E_field.components[0] *= -1.0
	#
	"""
	
	
	# write out basic simulation data...
	STATE.metadata.dt = dt
	STATE.metadata.output_interval__time_steps = STATE.num_timesteps_between_outputs
	STATE.metadata.set_dirty()
	
	initial_density = None
	initial_dist = None
	initial_dist_slice1 = None
	initial_dist_slice2 = None
	
	simple_E_buffer = None
	simple_E_buffer__index = 0
	
	output_interval = STATE.num_timesteps_between_outputs
	while( time < t_end):
		
		#if n >  6:
		#	exit(-1)
		
		# output!
		if time > 1.9:
			output_interval = STATE.num_timesteps_between_outputs
		if time > 2.1:
			output_interval = STATE.num_timesteps_between_outputs
		
		for i in range(0, num_species):
			if((n%output_interval) == 0):
				
				species_name = STATE.species[i].species_config.name
				species_name_prefix = ""
				if num_species > 1:
					species_name_prefix = species_name + "_"
				
				if True and n % 50 == 0:
					distribution_function = actual_dist_functions[i].get_full_lm( STATE.species[i].F, n, time)
					if(distribution_function != None):
						if initial_dist == None:
							initial_dist = np.copy ( distribution_function )
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
							savefig("./%sf_cheap_1_cross_section_%06d.png" % (species_name_prefix, n))
							
							clf()
							title("slices of F at 40 cheap at time=%f, step=%d" % (time, n) )
							plot(f_slice_axis, np.log(np.fabs(slut1)), color='green' )
							#plot(f_slice_axis, np.log(np.fabs(initial_dist[40,:])), color='black', alpha=.2)
							savefig("./%sf_cheap_2_cross_section_%06d.png" % (species_name_prefix, n))
							
							distribution_function -= initial_dist
							
							clf()
							spot = distribution_function.shape[1] / 5
							plot( distribution_function[:, spot-5] )
							plot( distribution_function[:, spot+5] )
							title("slices of F along x at cells +-5 at time=%f, step=%d" % (time, n) )
							savefig("./%sf_long_cheap_cross_section_%06d.png" % (species_name_prefix, n))
							"""
							pass
						
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
						if DaWaveDrivah != None:
							
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
						title("Dist function versus absolute value of Px\nstep %f, time %e" % (n, time))
						xlabel("Momentum along X")
						savefig("./%sf_cross_section_%06d.png" % (species_name_prefix, n) )
						
						clf()
						# cmap='PRGn'
						#imshow(np.log(np.fabs(distribution_function)), aspect='auto', cmap="copper", extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						#imshow(np.log(np.fabs(distribution_function)), cmap='Paired', aspect='auto', extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						imshow( distribution_function,  aspect='auto', cmap='Paired' , extent=[f_slice_axis[0], f_slice_axis[-1],0.0,STATE.spatial_axes__global[0].values[-1]]  )
						colorbar()
						#clim( (-4.5e-4, 4.5e4) )
						#xlim( (50,135) )
						"""
						X, Y = np.meshgrid(f_slice_axis, np.arange(0,distribution_function.shape[0]) )
						contourf( X,Y,distribution_function, 40, cmap=cm.Paired)
						"""
						# xlim( (-0.075, 0.075) )
						# xlim( (-0.5, 0.5) )
						# xlim( (f_slice_axis[0], f_slice_axis[-1]) )
						
						title("X versus Px\nStep %f, Time %e" % (n, time))
						savefig("./%sp1x1_%06d.png" % (species_name_prefix, n) )
						np.save("./numpy/p1x1_%06d.npy" % n, distribution_function)
						#xlim( (-2.0, 2.0) )
						#savefig("./%sp1x1_zoom_%06d.png" % (species_name_prefix, n) )
						
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
						#savefig("./harmonics_%06d.png" % n)
						subplot(212)
						bar( xlocations, np.log(actual_dist_functions[i].diag__harmonic_content_minus[:len(xlocations)]), width=width )
						xticks( xlocations + width/2.0, labels,  size="xx-small")
						xlim(0, xlocations[-1]+width*2.0)
						ylim([-10.0, 10.0])
						gca().get_yaxis().tick_left()
						savefig("./harmonics_%06d.png" % n)
						"""
				
				if i == 0:
					# This is the 'raw' E-field data from the simulation... so we need to remove the guard cells
					
					#if mpi_gather.gather( STATE.E_field__external.components[0], remove_guard_cells=True):
					if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
						if simple_E_buffer == None:
							simple_E_buffer = np.zeros( [100,mpi_gather.stiching_buffer.shape[0]] )
							simple_E_buffer__index = 0
						simple_E_buffer[simple_E_buffer__index, : ] = mpi_gather.stiching_buffer[:]
						simple_E_buffer__index  += 1
						if simple_E_buffer__index == simple_E_buffer.shape[0]:
							np.save("./numpy/E_%06d.npy" % (n, ) , simple_E_buffer)
							simple_E_buffer__index = 0
							print "SAVED 500 E-feild entries!!!!"
						
						#if n%25 == 0 and time > 240.0 and time < 300.0:
						if True:
							clf()
							subplot(211)
							plot(mpi_gather.stiching_buffer, marker="*", ms=.8)
							#ylim( [-0.0015, 0.0015] )
							#ylim( [-5e-5, 5e-5] )
							axhline( y = 0)
							title("E field at step %f, time %e" % (n, time));
							subplot(212)
							size = np.arange(0, 2.0, 2.0/mpi_gather.stiching_buffer.shape[0])
							size -= 1.0
							#fft_data = np.absolute(fft.fftshift(fft.fft(mpi_gather.stiching_buffer)))
							fft_data = np.absolute((fft.fftshift( fft.fft(mpi_gather.stiching_buffer ))))   		#, n=256)))
							
						
							
							plot(np.log10(fft_data))
							#ylim( [0.0, .04] )
							#ylim( [0.0, .005] )
							savefig("./png/E_%06d.png" % (n, ) )
							#np.save("./numpy/E_%06d.npy" % (n, ) , mpi_gather.stiching_buffer)
							#print "Signal to noise: %e" % scipy.stats.signaltonoise( mpi_gather.stiching_buffer)
							#print "length of fft buffer %d" % fft_data.shape[0]
						
				"""
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if mpi_gather.gather(  moments[i].get_temp_ev(), remove_guard_cells = True):
					tempp = moments[i].get_temp_ev()					
					clf(); plot( mpi_gather.stiching_buffer); 
					title("Temperatire (ev) at step %f, time %e" % (n, time)); draw();
					savefig("./png/%stemp_ev_%06d.png" % (species_name_prefix, n) )
					np.save("./numpy/%stemp_ev_%06d.npy" % (species_name_prefix,n ), mpi_gather.stiching_buffer)
				"""
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if False and n%50 == 0 and mpi_gather.gather(  moments[i].get_density(), remove_guard_cells=True ):
					if initial_density == None:
						initial_density = np.zeros ( mpi_gather.stiching_buffer.shape[0] )
						initial_density[:] = mpi_gather.stiching_buffer[:]
					else:
						mpi_gather.stiching_buffer -= initial_density
					
					clf(); 
					subplot(211)
					plot(mpi_gather.stiching_buffer); title("Density %.02f, time %.02e" % (n, time));
					#draw(); 
					subplot(212)
					size = np.arange(0, 2.0, 2.0/mpi_gather.stiching_buffer.shape[0]); size -= 1.0
					fft_data = np.absolute(fft.fftshift(fft.fft(mpi_gather.stiching_buffer)))
					fft_data += 1e-30
					plot(size, np.log10(fft_data))
					savefig("./png/%sdensity_%06d.png" % (species_name_prefix, n) )
					np.save("./numpy/%sdensity_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
					#print fft_data
					print "\t Total charge %e" % (np.sum(moments[i].get_density()), )
				
				"""
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if mpi_gather.gather(  moments[i].get_qx(), remove_guard_cells=True ):
					tempp = moments[i].get_qx()
					clf(); plot(mpi_gather.stiching_buffer); title("Qx at step %f, time %e" % (n, time));
					ylim( (-1e-8, 1e-8) ); axhline(y=0, color='black'); #lines_temp_qx[0].set_ydata( moments[i].get_qx())
					savefig("./png/%sqx_%06d.png" % (species_name_prefix, n) )
					np.save("./numpy/%sqx_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
				
				# the daignostics emmited from Moments already have the gaurd cells trimmed
				if mpi_gather.gather( moments[i].get_heat_conductivity_x(), remove_guard_cells=True ):
					clf(); plot(mpi_gather.stiching_buffer); title("Heat Conduction at step %f, time %e" % (n, time));
					axhline(y=0, color='black'); #lines_temp_kx[0].set_ydata(moments[i].get_heat_conductivity_x())
					ylim( (-4, +4) ); draw();
					savefig("./png/%skx_%06d.png" % (species_name_prefix, n) )
					np.save("./numpy/%skx_%06d.npy" % (species_name_prefix,n ) , mpi_gather.stiching_buffer)
				"""
				
				#output_f( STATE.species[i].F, n , time ,plot=plotme1.plot_ps)
				if mpi_info.rank == 0:
					print "Output written for n=%d" % n
				#exit(-1)
		
		# Do the actual simulation calculations!!!
		for i in range(0, num_species):
			if mpi_info.rank == 0:
				print "\t step %d for: %s" % ( n, STATE.species[i].species_config.name )
			solva.step( eval_functions[i], STATE.species[i], STATE, STATE.config )
			
			"""
			if DaWaveDrivah != None:
				if STATE.time <= (DaWaveDrivah.rise_time + DaWaveDrivah.sustain):
					print "doing Fucker-Plank"
					fokker_planks[i].calc_f00(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time)
					fokker_planks[i].calc_flm(STATE.species[i].F, STATE.species[i], STATE, STATE.config)
			"""
			
			#
			#
			#fake_data = np.ones(STATE.species[i].F.momentum_axis.num_points )
			#fake_data = np.arange(0, STATE.species[i].F.momentum_axis.num_points )
			#fokker_planks[i].reset_coeff( fake_data, 0.298462)
			#fokker_planks[i].implicit_advance( fake_data, 2)
			#exit(-1)
			#fokker_planks[i].reset_coeff( fake_data, 0.298462)
			#print "DONR!"
			#exit(-1)
		
		# wave driver
		if DaWaveDrivah != None:
			DaWaveDrivah.update(time, n, dt)
		
		if STATE.config.E_field_profile != None and time > STATE.config.E_field_profile.echo_time and True:
			# want to perturb system with output = current_state + amp*cos*current_state
			# let's go through the e field.... maybe need to goto F directly...
			if mpi_info.rank == 0:
				print "--------------------PERTURB PULSE!!!!!!!!!!!!!!!!!!!!!!!!!!=========="
			#buffah_accumulate = STATE.E_field.clone()
			#buffah_accumulate.scalar_mult(0.0)
			#print "BEFORE"
			#print buffah_accumulate
			try:
				second_pulse_amp = STATE.config.E_field_profile.second_pulse_amp
				STATE.config.E_field_profile.amp = second_pulse_amp
			except:
				pass
			try:
				second_pulse_k = STATE.config.E_field_profile.second_pulse_num_wavelengths
				STATE.config.E_field_profile.num_wavelengths = second_pulse_k
			except:
				pass
			
			"""
			STATE.config.E_field_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah = buffah_accumulate)
			#print "ADTERT"
			#print buffah_accumulate
			STATE.E_field.components[0] += buffah_accumulate.components[0]
			"""
			if( mpi_gather.gather(  moments[0].get_density(), remove_guard_cells=True ) ):
				echo_f0_density = np.copy ( mpi_gather.stiching_buffer )
			local_f0__echo = Esolva.get_charge_density_for_species(0)
			
			# hard code for one species for now...
			STATE.config.E_field_profile.apply_perturbation( STATE.species[0].F, STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset )
			# now try to update E to reflect this perturbation...
			# args are: Yh, species, state, CFG, Eout = None
			#STATE.E_field.components[0][:] = 0.0
			# HACK BUG BAD FIX ADAM
			
			moments[0].current_step -= 1
			if( mpi_gather.gather(  moments[0].get_density(), remove_guard_cells=True ) ):
				echo_f1_density = np.copy ( mpi_gather.stiching_buffer )
				echo_f1_density -= echo_f0_density
				
				clf()
				plot( echo_f0_density )
				savefig( "f0__initial", dpi=200)
				
				clf()
				plot( echo_f1_density )
				savefig( "f1__zecho", dpi=200)
			
			moments[0].current_step += 1
			
			local_f1__echo = Esolva.get_charge_density_for_species(0)
			local_f1__echo -= local_f0__echo
			STATE.E_field.components[0][:] -= local_f1__echo
			
			#eval_functions[0].current_solver.calc( STATE.species[0].F, None, None, None, None, Eout = STATE.E_field )
			
			# shut off the echo driver..
			STATE.config.E_field_profile.echo_time = 1e100
			
		# update boundry cells...
		# and update the super moment calculator.
		for i in range(0, num_species):
			if mpi_info.num_nodes == 1:
				update_boundries.calc(STATE, STATE.config)
			else:
				mpi_update.update(STATE, STATE.config)
			moments[i].notify_simulation_time(n, time)
		
		
		if mpi_info.rank == 0:
			TIMINGS.print_results()
			print "done with step %d (time=%f ) " % (n, time) 
		
		if STATE.metadata.is_dirty():
			STATE.metadata.save()
		
		
		time += dt
		n += 1
		STATE.n = n
		STATE.time = time
		
		
	pass
	if mpi_info.rank == 0:
		print "It is finished..."

def basic_e_on_dist_test1__init():
	# assume STATE is setup with something....
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# zero the DIST function...
	for xx in xrange(0, Y0.num_total_cells[0]):
		for l in xrange(0, Y0.num_l):
			Y0.get_harmonic_at_location(xx, l).momentum_projection[:] = 0.0
			#m_stop = l
			#if(m_stop > Y0.num_m):
			#	m_stop = Y0.num_m
			#for m in xrange(0, m_stop):
			#	Y0.get_harmonic_at_location(xx, l,m).momentum_projection[:] = 0.0
	# initalize a fake dist....
	# first up... spaital isotropic.. inital values only in Y(0,0)
	for xx in xrange(0, Y0.num_total_cells[0]):
		#p_vecter = Y0.get_harmonic_at_location(xx,0,0).momentum_projection
		p_vecter = Y0.get_harmonic_at_location(xx,0).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = pp*.1
	pass
	
	
def basic_e_on_dist_test1__loop():
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(0, 0).momentum_projection[:]
	print Yh.get_harmonic_at_location(0, 0).momentum_projection[:]
	
	evalf.e_on_dist.calc								(Y0, Yh, species, STATE, CFG)
	
	print Y0.get_harmonic_at_location(0, 1).momentum_projection[:]
	Y0.get_harmonic_at_location(0, 1).momentum_projection[0] = 666.0
	
	print Yh.get_harmonic_at_location(0, 1).momentum_projection[:]
	exit(-1)
	pass

def basic_adv_test1__init():
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# initalize a fake dist....
	# make spatital gradients for the advection codess....
	numx = 128
	period = 2*np.pi / (numx + 1)
	print "perios %f" % period
	
	for xx in xrange(0, Y0.num_total_cells[0]):
		cos_val = np.cos(period* float(xx));
		if cos_val < 0.0:
			cos_val *= -1.0
		p_vecter = Y0.get_harmonic_at_location(xx,0).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = .1*pp*cos_val
	pass

def basic_adv_test1__loop():
	evalf = eval_functions[0]
	print "I WIN!!!!"
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(2, 0).momentum_projection[:]
	
	evalf.spatial_advection_solver.calc								(Y0, Yh, species, STATE, CFG)
	
	
	print Yh.get_harmonic_at_location(2, 1).momentum_projection[:]
	exit(-1)
	pass


def testing1_init():
	
	# assume STATE is setup with something....
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# zero the DIST function...
	for xx in xrange(0, Y0.num_total_cells[0]):
		for l in xrange(0, Y0.num_l):
			m_stop = l
			if(m_stop > Y0.num_m):
				m_stop = Y0.num_m
			for m in xrange(0, m_stop):
				Y0.get_harmonic_at_location(xx, l).momentum_projection[:] = 0.0
	
	# fake a dumb dist function.... that has a current...
	for xx in xrange(0, Y0.num_total_cells[0]):
		p_vecter = Y0.get_harmonic_at_location(xx,1).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = pp*.1
	
	
	
def testing1_main(eval_functions):
	evalf = eval_functions[0]
	print "I WIN!!!!"
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(2, 1).momentum_projection[:]
	
	evalf.current_solver.calc								(Y0, Yh, species, STATE, CFG)
	
	
	print Yh.get_harmonic_at_location(2, 1).momentum_projection[:]
	exit(-1)
	pass
				
	
# set this externaly to configure things...
CONFIG_PROFILE = None
CONFIG_CALLBACK = None

STATE 	= None
FIGURE 	= None
FIGURE1 	= None
KEY		= None
GUI_LOOP	= None
filemode	= False
test_mode = False

def two_stream_setup():
	# the most basic structral info....
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						=  1
	num_cells__global[0] 		=  48 #576 #504  # 480 #528 # 512 = offical values
	system_size_min__global[0] = 0
	system_size_max__global[0] = 20.0
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	#STATE.init(-1243.7, 1243.7, 120, 1, 4, 1)
	#CFG = STATE.config
	CFG.num_species 				= 1
	CFG.num_l 						= 12
	CFG.num_m 						= 1
	CFG.dt 							= (.5/26.0) #0.0192308 #.005 #299383 #.01 #.002
	CFG.t_end 						= 50.0
	
	
	CFG.abs_p_max 					= 3.3 #0.15
	CFG.abs_p_num_cells 			= 96
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= .02  ##.05#0.001  # .1
	CFG.if_implicit1D 			= True
	
	CFG.E_field_profile 			=  None #EFieldSinusodialProfile( 1e-1, num_wavelengths = 1)
	
	# working (not crash but no waves) x= -.1,.1  vth = .024, e_amp=1e-4 dt = .01
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 2.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e23
		temp.ln_lambda 			= 10.2
		temp.zeta 					= 1.0
		# constant density gradient at 1.0
		
		# .024 is 300ev
		temp.pt_amb 				= 0.024
		temp.temp_profile 		= ConstantProfile( 0.024) #0.025 ) #.01 smalld t worked.
		temp.density_profile = ConstantProfile( 1.0 )
	return CFG
	
def ld_setup():
	
	# working: e amp: .0001
	#     x: 0-->.1:   96
	#    dt: 0.01
	# p_max: .1: 96
	#  num_l:  24
	#   temp:  0.0124
	#
	#
	
	# with 20 box, 30, 45 
	# the most basic structral info....
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						=  1
	num_cells__global[0] 		=  96 # 96 #96 #384 #192 #192 #504 #200  #96 # worked nice
	system_size_min__global[0] =  0.0
	system_size_max__global[0] =  10.0  # 1.0
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	# for echo attemp 55, 35 i used abs_p_num_cells=750 num_l=24 num_cells=192 and spacing was nice and good
	
	# i used math!... to derive another cool condition (almost all dealiased in k) 96 cells, pulses at 12,8 x=5.0
	
	#STATE.init(-1243.7, 1243.7, 120, 1, 4, 1)
	#CFG = STATE.config
	CFG.num_timesteps_between_outputs = 1
	CFG.colatate_factor 			= 500
	CFG.num_species 				= 1
	CFG.num_l 						= 16
	CFG.num_m 						= 1
	CFG.dt 							= .1 #2*np.pi/64.0 #0.1 		#0.005 #.02 #299383 #.01
	CFG.t_end 						= 2000.0 #3000.1  #3880.0
	
	CFG.abs_p_max 					= 1.0 			# .75 #.06 #.15 #0.15
	CFG.abs_p_num_cells 			= 192 #384 #192 #96  #150   #96
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= .002  ##.05#0.001  # .1
	CFG.if_implicit1D 			= True
	
	CFG.E_field_profile 			=  None# EFieldSinusodialProfile( 1e-2, num_wavelengths = 1)
	
	CFG.use_wave_driver			= True
	
	# working (not crash but no waves) x= -.1,.1  vth = .024, e_amp=1e-4 dt = .01
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 2.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e20
		temp.ln_lambda 			= 10.2
		temp.zeta 					= 1.0
		# constant density gradient at 1.0
		
		# .024 is 300ev
		temp.pt_amb 				= 0.024
		temp.temp_profile 		= ConstantProfile( 4.420e-2)  #.141 ) #4.420e-2)  #0.0124) #0.025 ) #.01 smalld t worked.
		temp.density_profile = ConstantProfile( 1.0 )
	return CFG

"""
	# abase line testing for now (fast version)
	#num_cells__global[0] 		=  192
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750		#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	
	CFG.dt 							= .4	# .1
	CFG.t_end 						= 1000.0
	
	#1.0e-18 is border line non-linear...
	
	CFG.E_field_profile 			= None
	CFG.E_field_profile 			= EFieldSinusodialProfile( 1.0e-18, num_wavelengths = 16)
	CFG.E_field_profile.echo_time 							= 250.0
	CFG.E_field_profile.second_pulse_amp 					= 1.0e-18
	CFG.E_field_profile.second_pulse_num_wavelengths 	= 32.0

"""
def ld_setup_3rd():
	
	# working: e amp: .0001
	#     x: 0-->.1:   96
	#    dt: 0.01
	# p_max: .1: 96
	#  num_l:  24
	#   temp:  0.0124
	#
	#
	
	# with 20 box, 30, 45 
	# the most basic structral info....
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						=  1
	num_cells__global[0] 		=  1152 # 96 #96 #384 #192 #192 #504 #200  #96 # worked nice
	system_size_min__global[0] =  0.0
	system_size_max__global[0] =  100.0  # 1.0
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	# for echo attemp 55, 35 i used abs_p_num_cells=750 num_l=24 num_cells=192 and spacing was nice and good
	
	# i used math!... to derive another cool condition (almost all dealiased in k) 96 cells, pulses at 12,8 x=5.0
	
	#STATE.init(-1243.7, 1243.7, 120, 1, 4, 1)
	#CFG = STATE.config
	CFG.num_timesteps_between_outputs = 1
	CFG.colatate_factor 			= 500
	CFG.num_species 				= 1
	CFG.num_l 						= 16
	CFG.num_m 						= 1
	CFG.dt 							= .1 #2*np.pi/64.0 #0.1 		#0.005 #.02 #299383 #.01
	CFG.t_end 						= 2000.0 #3000.1  #3880.0
	
	CFG.abs_p_max 					= 1.0 			# .75 #.06 #.15 #0.15
	CFG.abs_p_num_cells 			= 192 #384 #192 #96  #150   #96
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= .002  ##.05#0.001  # .1
	CFG.if_implicit1D 			= True
	
	CFG.E_field_profile 			=  None# EFieldSinusodialProfile( 1e-2, num_wavelengths = 1)
	
	CFG.use_wave_driver			= True
	
	# working (not crash but no waves) x= -.1,.1  vth = .024, e_amp=1e-4 dt = .01
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 2.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e20
		temp.ln_lambda 			= 10.2
		temp.zeta 					= 1.0
		# constant density gradient at 1.0
		
		# .024 is 300ev
		temp.pt_amb 				= 0.024
		temp.temp_profile 		= ConstantProfile( 4.420e-2)  #.141 ) #4.420e-2)  #0.0124) #0.025 ) #.01 smalld t worked.
		temp.density_profile = ConstantProfile( 1.0 )
	return CFG

def ld_setup__standing_wave__3rd_order():
	CFG = ld_setup_3rd()
	
	
	#num_cells__global[0] 		=  192
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	
	CFG.dt 							= .25	# .1
	CFG.t_end 						= 600.0
	
	#1.0e-18 is border line non-linear...
	
	CFG.E_field_profile 			= None
	CFG.E_field_profile 			= EFieldSinusodialProfile( 1.0e-8, num_wavelengths = 160)
	CFG.E_field_profile.echo_time 							= 200.0
	CFG.E_field_profile.second_pulse_amp 					= 1.0e-8
	CFG.E_field_profile.second_pulse_num_wavelengths 	= 160.0
	
	CFG.use_wave_driver			= False
	
	return CFG
	
def ld_setup__standing_wave():
	return ld_setup__standing_wave__3rd_order()
	
	CFG = ld_setup()
	
	
	#num_cells__global[0] 		=  192
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 1000	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	
	CFG.dt 							= .25	# .1
	CFG.t_end 						= 1050.0
	
	#1.0e-18 is border line non-linear...
	
	CFG.E_field_profile 			= None
	CFG.E_field_profile 			= EFieldSinusodialProfile( 1.0e-12, num_wavelengths = 16)
	CFG.E_field_profile.echo_time 							= 200.0
	CFG.E_field_profile.second_pulse_amp 					= 1.0e-12
	CFG.E_field_profile.second_pulse_num_wavelengths 	= 24
	
	CFG.use_wave_driver			= False
	
	return CFG
	
def spitzer_setup():

	# the most basic structral info....
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						=  1
	num_cells__global[0] 		=  120
	system_size_min__global[0] = -1243.7
	system_size_max__global[0] =  1243.7
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	
	#STATE.init(-1243.7, 1243.7, 120, 1, 4, 1)
	#CFG = STATE.config
	CFG.num_species 				= 1
	CFG.num_l 						= 5 #4
	CFG.num_m 						= 3 #1
	CFG.dt 							= .298462 #.299383 #.01
	CFG.t_end 						= 3880.0
	
	CFG.abs_p_max 					= 0.15
	CFG.abs_p_num_cells 			= 108
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 4
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= .1  ##.05#0.001  # .1
	CFG.if_implicit1D 			= True
	
	CFG.E_field_profile 			= None
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 2.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e23
		temp.ln_lambda 			= 10.2
		temp.zeta 					= 1.0
		# constant density gradient at 1.0
		
		# .024 is 300ev
		temp.pt_amb 				= 0.024
		temp.temp_profile 		= SinusodialProfile( 0.025 )
		temp.density_profile = ConstantProfile( 1.0 )
	return CFG
	
def gaussian_temp_setup():
	
	# the most basic structral info....
	#num_cells__global 		= np.zeros(3, dtype='int32')
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						= 1
	num_cells__global[0] 		= 144
	system_size_min__global[0] = 0.0
	system_size_max__global[0] = 10.0
	
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	# these all local values....
	#STATE.init(0.0, 10.0, 32, 1, 4, 1)
	CFG.num_timesteps_between_outputs = 1
	CFG.num_l 						= 12
	CFG.num_m 						= 1
	CFG.dt 							= .1
	CFG.t_end 						= 1000.0 #300.0 #.11    #1000.0
	
	CFG.abs_p_max 					= .4
	CFG.abs_p_num_cells 			= 96
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= 0.01
	CFG.if_implicit1D				= True
	
	
	CFG.E_field_profile 			= None
	
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 2.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e23
		temp.ln_lambda 			= 10.2
		temp.zeta 					= 1.0
		
		# constant density gradient at 1.0
		# .024 is 300ev
		temp.pt_amb 				= 0.044		
		#temp.temp_profile 		= GaussianProfile( .05, .8, 5.0 )	
		#temp.temp_profile 		= GaussianProfile( .3, 10.0, 5.0 )	
		temp.temp_profile 		= GaussianProfile( .35, 17.0, 5.0 )
		
		temp.density_profile = ConstantProfile( 1.0 )
	
	return CFG

def bland_const_temp():
	me = gaussian_temp_setup()
	me.num_cells__global[0] 		= 192 #60
	me.species_config[0].temp_profile = ConstantProfile( .024 ) #4.420e-2)
	return me
	

def e_p__plasma():
	
	# the most basic structral info....
	#num_cells__global 		= np.zeros(3, dtype='int32')
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						= 2
	num_cells__global[0] 		= 64
	system_size_min__global[0] = 0.0
	system_size_max__global[0] = 10.0
	
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	# these all local values....
	#STATE.init(0.0, 10.0, 32, 1, 4, 1)
	CFG.num_timesteps_between_outputs = 20
	
	CFG.num_l 						= 4
	CFG.num_m 						= 1
	CFG.dt 							= .01
	CFG.t_end 						= 100.0 
	
	CFG.abs_p_max 					= .4
	CFG.abs_p_num_cells 			= 96
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	CFG.if_tridiagonal 			= True
	
	CFG.small_dt 					= 0.01
	CFG.if_implicit1D				= True
	
	
	CFG.E_field_profile 			= None
	CFG.E_field_profile 			= EFieldSinusodialProfile( 1e-3, num_wavelengths = 1)
	
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
	
	temp = CFG.species_config[0]
	temp.name = "electrons"
	temp.charge 				= -1.0
	temp.mass 					= 1.0
	temp.vth 					= 2.0
	# the density in 1/cm^3
	temp.density_np 			= 1.0e21
	temp.ln_lambda 			= 10.2
	temp.zeta 					= 1.0
	# constant density gradient at 1.0
	# .024 is 300ev
	temp.pt_amb 				= 0.044		
	temp.temp_profile 		=  ConstantProfile( 0.05) #0.024)
	temp.density_profile = ConstantProfile( 1.0 )
	
	
	temp = CFG.species_config[1]
	temp.name = "positrons"
	temp.charge 				= 1.0
	temp.mass 					= 1.0
	temp.vth 					= 2.0
	# the density in 1/cm^3
	temp.density_np 			= 1.0e23
	temp.ln_lambda 			= 10.2
	temp.zeta 					= 1.0
	# constant density gradient at 1.0
	# .024 is 300ev
	temp.pt_amb 				= 0.044		
	temp.temp_profile 		=  ConstantProfile( 0.05) 
	temp.density_profile = ConstantProfile( 1.0 )
	
	
	return CFG
	
def setup():
	
	global STATE
	global filemode
	global moments
	global DaWaveDrivah
	
	filemode = True
	two_stream = False
	
		
	if two_stream:
		CFG = two_stream_setup()
	else:
		pass
		#print CFG.num_species
		#CFG = ld_setup()
		"""
		CFG = gaussian_temp_setup()
		CFG.species_config[0].pt_amb = .01
		CFG.species_config[0].temp_profile = ConstantProfile( 0.024, x_start=4.5, x_end=5.5)
		CFG.species_config[0].density_profile = ConstantProfile( 1.0, x_start=4.5, x_end=5.5)
		"""
		if CONFIG_PROFILE != None:
			CFG = CONFIG_PROFILE
		elif CONFIG_CALLBACK != None:
			CFG = CONFIG_CALLBACK()
		else:
			#CFG = bland_const_temp()
			# CFG = ld_setup()
			CFG = ld_setup__standing_wave()
			#CFG = e_p__plasma()
			#CFG = spitzer_setup()
		
		
	if test_mode:
		CFG = spitzer_setup()
	
	# create the official State Object!
	pickle.dump( CFG, open( "oshun_metadata__cfg.p", "wb" ) )
	
	
	
	STATE = DaState( CFG )
	CFG = STATE.config
	
	STATE.num_timesteps_between_outputs = CFG.num_timesteps_between_outputs
	
	#............ now use CFG  to make simulation objects...........
	
	# Setup Momentum |p| data... (this axis is not chopped up by MPI)
	abs_p_axis = Axis(CFG.abs_p_min, num_points=CFG.abs_p_num_cells, max=CFG.abs_p_max)
	
	# Configure the settings about boundry cells and solvas.
	# 		We require the number of spatial boundry cells (on each side) that 
	#		is the same as the of the solva with max order.
	num_boundry_cells = np.zeros( (3,2), dtype=int32 )
	
	# for 1d
	num_boundry_cells[0][0] = CFG.explicit_solver_order
	num_boundry_cells[0][1] = CFG.explicit_solver_order
	
	# now let's setup MPI sooner rather then later.
	# the main interface with MPI is via MPIAxis objects that chunk up global spatial axes
	#		into the appropiate bits or all the nodes.
	mpi_info.init()
	mpi_info.mpi_axes.append( AxisMpi(STATE.spatial_axes__global[0], num_boundry_cells, mpi_info) )
	# now use this MPI axis to fill in the State object about the node-related spatial divisions.
	# the STATE.spatial_axes will cop the mpi_axis's 'axis__local' object... this axis
	# starts has extent cells = [0..num_cells_on_node], x = [0.0, last_x_position_on_node]
	STATE.spatial_axes[0] = Axis( axis=mpi_info.mpi_axes[0].axis__local )
	
	
	STATE.E_field = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)
	STATE.B_field = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)
	
	# initialize any inital E-field, if requested.
	if CFG.E_field_profile != None:
		CFG.E_field_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset)
	
		vth = CFG.species_config[0].temp_profile.val
		w_wave = np.sqrt( 1 + 3*CFG.E_field_profile.wavenumber*CFG.E_field_profile.wavenumber*vth*vth)
		#print "freq (BG): %e " % w_wave
		#print "vphase: %e" % (w_wave/ CFG.E_field_profile.wavenumber)
		#print "K l_debye: %f" % (CFG.E_field_profile.wavenumber*vth)
		CFG.E_field_profile.v_phase = (w_wave/ CFG.E_field_profile.wavenumber)
	
	"""
	temp = zeros( STATE.E_field.components[0].shape[0] -8 )
	temp[:] =  STATE.E_field.components[0][4:-4]
	
	fft_temp = fft.fftshift( np.absolute(fft.fft( temp) ) )
	#fft_temp = (fft.fft( temp) ) 
	subplot(211)
	stem(np.arange(0, fft_temp.shape[0]), fft_temp)
	subplot(212)
	ttemp = STATE.E_field.components[0][4:-4]
	#stem( np.arange(0, STATE.E_field.components[0].shape[0]-8),ttemp )
	stem( STATE.spatial_axes__global[0].values,ttemp )
	print ttemp[1]
	print ttemp[-1]
	show()
	exit(-1)
	"""
	
	for i in range(0, CFG.num_species):
		# each species has it's own abs_p_axis to describe the |p| dimension...
		#	for now, they are all using the same one....
		species_config = CFG.species_config[i]
		temp_species = Species(	STATE.spatial_axes[0].num_points, CFG.num_l,
										momentum_axis = abs_p_axis, 
										species_config = species_config,
										num_boundry_cells = num_boundry_cells)
		#print "Just created offical species object for %s" % species_config.name
		STATE.species[i] = temp_species
		
		#	cheap hack. TODO: take out.
		if i == 0:
			mpi_info.mpi_axes[0].init_species_data( 	STATE.species[i],
																	STATE.E_field,			
																	STATE.species[i].F.cells[0].num_harmonics)
		
		
		# each species also has it's own J.. for now.. maybe not needed.
		temp_species.J = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)
		
				
		# set density profile...
		# set inital to 1.0 all places...num_cells
		density_profile = np.ones( temp_species.F.num_total_cells[0] ) # plain
		
		species_config.density_profile.apply_density( density_profile,STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local, .2, temp_species)
		density_profile *= 1.0
		
		# then you can mees with any values in these cells...
		temperature_grid = np.ones( temp_species.F.num_total_cells[0] )
		temp_amb = species_config.pt_amb
		
		"""
		test = HarmonicCentricView( np.array([48]), abs_p_axis)
		test_in = np.load('diff_dx_in.npy')
		test.data[:,:] = test_in[:,:]
		test.Dx()
		for ip in xrange(0, 96):
			print "\t",
			for ix in xrange(0, 48):
				print "%.9e" % test.data[ix, ip],
			print ""
		print ""
		exit(-1)
		"""
		
		if two_stream == True:
			temp_species.F.scalar_mult ( 0.0 )
			p_beam = 1.0
			p_beam_center = 1.0
			beam_temp = np.sqrt(.011) #.105    #0.011#0.105
			index_00 = temp_species.F.linear_index_for_harmonic(0)
			index_20 = temp_species.F.linear_index_for_harmonic(2)
			index_10 = temp_species.F.linear_index_for_harmonic(1)
			
			for ix in xrange(0, temp_species.F.num_plain_cells[0] ):
				#seed = np.cos(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix])
				"""
				if ix==0 and mpi_info.rank==0: 
					xval = 1.387778780781e-17
					seed = .1*np.sin(.20*np.pi*xval) # WORKS
				else:
					seed = .1*np.sin(.20*np.pi*(STATE.spatial_axes[0].values[ix])) # WORKS
				"""
				#if ix==0:
				#	STATE.spatial_axes[0].values[ix] += STATE.spatial_axes[0].dx()
				#if ix==temp_species.F.num_plain_cells[0]-1:
				#	STATE.spatial_axes[0].values[ix] -= STATE.spatial_axes[0].dx()
				#seed = .1*np.sin(.20*np.pi*(STATE.spatial_axes[0].values[ix])) # WORKS
				
				# niec deritivitce
				dx = STATE.spatial_axes[0].dx()
				seed = .1*np.sin(2*np.pi/(10.0+.5*dx)*(STATE.spatial_axes[0].values[ix]))
				
				
				#seed = .1*np.sin(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx())*(STATE.spatial_axes[0].values[ix] - STATE.spatial_axes[0].dx()))
				#seed = .1*np.sin(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix])
				#seed *= 10
				#seed *= 10.0 #check...
				#p_beam = 1.2*np.cos(2.0*np.pi/(20.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix]) + p_beam_center
				exp_arg_coeff =  -1.0/(2.0*np.power(beam_temp,2))
				factor_coeff = 1.0/( np.power(beam_temp, 3.0) *np.power(2.0*np.pi, 1.5) )
				
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					abs_p = STATE.species[i].momentum_axis.values[ip]
					exp_argument = exp_arg_coeff*(abs_p - p_beam)*(abs_p - p_beam)
					factor = .0104*factor_coeff
					factor *= np.exp(exp_argument)
					pedistal = .07*heaviside_step(p_beam-abs_p)
					
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_00].momentum_projection[ip] = pedistal + factor
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_20].momentum_projection[ip] = 2.0*factor
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_10].momentum_projection[ip] = factor*.2*seed #factor*(1.0/50.0)*seed
				
			
			#DEBUG( Y=temp_species.F,L=0, M=0)
			
			enable = False
			if mpi_info.rank == 0 and enable:
				x_points = [0,1,2,3]
				print ""
				for ix in x_points:
					print "\t%f" % STATE.spatial_axes[0].values[ix],
				print ""
				
				print "beam temp %f" % beam_temp
				print "EXP arg factor is: %f" % exp_arg_coeff
				print "FROT factor is: %f" % factor_coeff
				print "Y(0,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_00].momentum_projection[ip] ),
					print ""
				print ""
				print ""
				print "Y(1,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_10].momentum_projection[ip] ),
					print ""
				print ""
				print ""
				print "Y(2,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_20].momentum_projection[ip] ),
					print ""
				#exit(-1)
			
			
		else:
			#----------- Covert density and temperature profiles into the inital
			#------------	values for the |p| axes in the shperical harmonics.
			# As a baseline, set all |p| cells to 1.0 at all spatial locations 
			#		for the Y(0,0) harmonics..
			for x in xrange(0, temp_species.F.num_total_cells[0] ):
				temp_species.F.cells[x].harmonics[0].momentum_projection[:] = 1.0  
			
			# Then modulate Y(0,0)'s |p| entries to reflect the requested density profile
			# 		(we ignore the guard cells)
			for ix in xrange(0, temp_species.F.num_plain_cells[0] ):
				density_multiplier = density_profile[num_boundry_cells[0][0]+ ix]
				temp_species.F.cells[num_boundry_cells[0][0]+ ix].harmonics[0].momentum_projection[:] *= density_multiplier
				
			
			
			
			# -------------- Initalize the species's spatial temperature profile ---------------
			# 		Result is in temperature_grid
			#species_config.temp_profile.apply( temperature_grid, STATE.spatial_axes[0], temp_amb, temp_species)
			species_config.temp_profile.apply( temperature_grid, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local, temp_amb, temp_species)
			#for ixx in xrange(0, temp_species.F.num_total_cells[0]):
			#	print "%.04e %.04e" % ( temperature_grid[ixx], density_profile[ixx] ) 
			#exit(1)
			#-----------poke this temperature map into the Y0,0 harmonics as 
			#		a gaussian along the |p| axis.
			alpha = np.zeros( temp_species.F.num_total_cells[0]  )
			coeff = np.zeros( temp_species.F.num_total_cells[0]  )
			sh_linear_index 		= temp_species.F.linear_index_for_harmonic(0)
			sh_linear_index_10 	= temp_species.F.linear_index_for_harmonic(1)
			sh_linear_index_11 	= temp_species.F.linear_index_for_harmonic(1)
			
			for ix in xrange(0, temp_species.F.num_total_cells[0]):
				alpha[ix] = 1.0/(2.0*temperature_grid[ix]*temperature_grid[ix])
				coeff[ix] = 1.0/( math.sqrt(2.0*np.pi)*2.0*np.pi*temperature_grid[ix]*temperature_grid[ix]*temperature_grid[ix] )
			
			
			for p in xrange(0, abs_p_axis.num_points):
				p_squared = abs_p_axis.values[p] * abs_p_axis.values[p]
				for ix in xrange(0, temp_species.F.num_total_cells[0]):
					temp_species.F.cells[ix].harmonics[sh_linear_index].momentum_projection[p] *= math.exp(-1.0*(alpha[ix]*p_squared))
					temp_species.F.cells[ix].harmonics[sh_linear_index].momentum_projection[p] *= coeff[ix]
				pass
			
			
			
			
			"""
			if CFG.num_species > 1:
				print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %d" % i
				if i == 0:
					print "HI#############################"
					sh_linear_index_20 	= temp_species.F.linear_index_for_harmonic(2)
					for ix in xrange(20, 30): #temp_species.F.num_total_cells[0]):
						temp_species.F.cells[ix].harmonics[sh_linear_index_20].momentum_projection[1:30] += 25.0
				if i == 1:
					#sh_linear_index_20 	= temp_species.F.linear_index_for_harmonic(2)
					#for ix in xrange(0, temp_species.F.num_total_cells[0]):
					#	temp_species.F.cells[ix].harmonics[sh_linear_index_20].momentum_projection[5] += .05
					pass
			"""
				
	
	if test_mode:
		testing1_init()
	
	# setup the Output moments generator... (one per species)
	for i in range(0, CFG.num_species):
		moments.append( MomentCalculationsState(STATE.species[i], STATE) )
		
	# make directories...
	if mpi_info.rank==0:
		try:
			os.makedirs("./numpy")
		except:
			pass
		try:
			os.makedirs("./png")
		except:
			pass
	
	# C stuff...
	general_setup__c( CFG.num_species, STATE.E_field.components[0], STATE.E_field.num_boundry_cells, 
							STATE.spatial_axes[0])
	for i in xrange(0, CFG.num_species):
		general_species_setup( 	i, STATE.species[i].num_l, STATE.species[i].num_m, 
										STATE.species[i].F.momentum_axis,
										STATE.species[i].F.num_total_cells, STATE.species[i].F.num_plain_cells, 
										STATE.species[i].F.num_boundry_cells[i] )
		for ix in xrange(0, STATE.species[i].F.num_total_cells[0]):
			for ish in xrange(0, STATE.species[i].F.cells[ix].num_harmonics):
				transfer_species_F_buffers( i, ish, ix, STATE.species[i].F.cells[ix].harmonics[ish].momentum_projection )
	
	# ------------------------------- Done setting up inital temp/density profiles ------------------------
	#print "RANK %d is wating at block" % mpi_info.rank
	#mpi_info.comm.Barrier()	
	#print "RANK %d Done" % mpi_info.rank
	
	# this is stupid python stuff the will mostly go away....
	return True

def run_sim(config=None, config_callback=None):
	global CONFIG_PROFILE
	global CONFIG_CALLBACK
	CONFIG_PROFILE = config 
	CONFIG_CALLBACK = config_callback 
	go_on = setup()
	if go_on:
		main_loop()
	return

if __name__ == "__main__":
	go_on = setup()
	if go_on:
		main_loop()
	
	