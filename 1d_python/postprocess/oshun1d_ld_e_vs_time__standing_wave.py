import numpy as np
import numpy.fft
import glob
import os.path
import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy.optimize import curve_fit
import cPickle as pickle

import scipy.special
import scipy.optimize
import scipy.stats

from vqunits import Formulary 

import sys
sys.path.append("./trunk1d")

from Oshun1d import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile
#from oshun_hdf_lame import *

one_over_sqrt_2pi = 1.0 / np.sqrt(2.0*np.pi)
mone_over_sqrt_2pi = -1.0 / np.sqrt(2.0*np.pi)

class DoResults:
	def __init__(self):
		self.values = {}

def plot_echo_profile( k1__ve, k2__ve, k3__ve, ve):
	v_min = -2.0; v_max=2.0; num_steps = 400;
	v_vals = np.linspace( v_min, v_max, num=num_steps)
	t_vals = np.linspace( -10.0, 10.0, num=100)
	
	vals_at_t_real = np.zeros( t_vals.shape[0] )
	vals_at_t_imag = np.zeros( t_vals.shape[0] )
	
	for i in xrange(0, t_vals.shape[0]):
		time = t_vals[i]
		val_at_time = plot_echo_profile__get_val_at_t( k1__ve, k2__ve, k3__ve, ve, time, v_vals)
		vals_at_t_real[i] = val_at_time[0]
		vals_at_t_imag[i] = val_at_time[1]
	pass
	vals_at_t_real *= ve
	vals_at_t_imag *= ve
	mag = np.sqrt( np.power(vals_at_t_real, 2) + np.power(vals_at_t_imag, 2) )
	
	return (t_vals, vals_at_t_real, vals_at_t_imag, mag )

def plot_echo_profile__get_val_at_t( k1, k2, k3, ve, t, v_vals):
	temp_real = np.zeros( v_vals.shape[0] )
	temp_imag = np.zeros( v_vals.shape[0] )
	
	for i in xrange(0, v_vals.shape[0]):
		temp = plot_echo_profile__get_integrand_at_t_and_v( k1, k2, k3, ve, t, v_vals[i] )
		temp_real[i] = temp[0]
		temp_imag[i] = temp[1]
	
	dv = v_vals[1] - v_vals[0]
	val_real = np.sum( temp_real)
	val_real /= dv
	val_imag = np.sum( temp_imag)
	val_imag /= dv
	
	return ( val_real, val_imag )

	
def plot_echo_profile__get_integrand_at_t_and_v(k1, k2, k3, ve, t, v):
	f0_dv = mone_over_sqrt_2pi*v*np.exp(-1.0*v*v/(2.0*ve*ve))
	wave_factor_arg = k3*v*t
	wave_factor_real = np.cos( wave_factor_arg )
	wave_factor_imag = np.sin( wave_factor_arg )
	dielectric_k1 = plasma_epsilon_es( complex(0, k1*v), -1.0*k1, 1.0, ve)
	dielectric_k2 = plasma_epsilon_es( complex(0, -1.0*k2*v), k2, 1.0, ve)
	dielectric_k3 = plasma_epsilon_es( complex(0, -1.0*k3*v), k3, 1.0, ve)
	dielecric_factor = dielectric_k1*dielectric_k2*dielectric_k3
	
	integrand_real = f0_dv*wave_factor_real / dielecric_factor
	integrand_imag = f0_dv*wave_factor_imag / dielecric_factor
	integrand = integrand_real + complex(0,1)*integrand_imag
	return ( np.real(integrand), np.imag(integrand_imag) )
	
	
def plasma_epsilon_es( w, k, wp_e, vth__e, maxwellian_convention_factor = 2.0 ):
	chi_e = np.power( ( wp_e/(vth__e*k) ), 2.0) / maxwellian_convention_factor
	x = w / (np.sqrt(2)*vth__e*k)
	val = 1.0 - chi_e*plasma_dispersion_prime( x )
	return val

def plot_driver():
	spot = 50
	files = glob.glob( './m7__test7/f*.npy')
	files.sort()
	driver_data = np.zeros( len(files) )
	for i,file in enumerate(files):
		temp_data = np.load(file)
		print file
		d = temp_data[spot]
		
		driver_data[i] = d
	#plot( driver_data)
	#savefig('driver.png')
	return driver_data

def find_maxima( data, start_index=0, end_index = None, sort_by_freq=False, minima = False):
	if end_index == None:
		end_index = data.shape[0]
	
	"""
	#deritiv = np.gradient( data[start_index: end_index] )
	deritiv1 = np.copy ( data[start_index: end_index] )
	deritiv = np.zeros( deritiv1.shape[0] )
	deritiv[0:-1] = np.copy ( data[start_index+1: end_index] )
	deritiv -= deritiv1
	deritiv[0] = deritiv[1]; deritiv[-1] = deritiv[-2]
	
	fred =  (numpy.diff(numpy.sign(deritiv)) != 0)*1
	"""
	
	if minima:
		fred = (diff(sign(diff(data[start_index: end_index]))) > 0).nonzero()[0] + 1
	else:
		fred = (diff(sign(diff(data[start_index: end_index]))) < 0).nonzero()[0] + 1
	
	"""
	#a = diff(sign(diff(data))).nonzero()[0] + 1 # local min+max
	#b = (diff(sign(diff(data))) > 0).nonzero()[0] + 1 # local min
	"""
	
	maxima = [] 
	for i in xrange(0, fred.shape[0] ):
		#if fred[i] == 1:
		#maxima.append( ( start_index+fred[i],  data[start_index+i] ) )
		maxima.append( ( start_index+fred[i],  data[start_index+fred[i] ] ) )
	
	if sort_by_freq:
		# sort by freq
		maxima.sort(key=lambda x: x[0], reverse=False)
	else:
		# sort by value
		maxima.sort(key=lambda x: x[1], reverse=True)
	
	return maxima
	
def spectrum_test():
	
	box_length_x = 10.0
	num_x_cells = 96 #192
	cell_size_x = box_length_x / float( num_x_cells )
	vth = 4.420e-2
	wp_e = 1.0
	num_wavelengths = 16
	
	dt = .01
	sim_start_time = 0.0
	sim_end_time = 10.0
	
	dt_mode_fancy = False
	dt_mode_selfish = False
	dt_mode_fancy_num_cells = 1.0
	
	sample_spot = 50
	
	E_field = np.zeros( num_x_cells) 
	
	# create the wavelength/number for the wave we want to use...
	wavelength = box_length_x / float(num_wavelengths)
	wavenumber = 2.0*np.pi / wavelength
	
	# setup the x spitail grid spacing, cell size, axis, etc...
	cell_size__phase = 2.0*np.pi / ( float( num_x_cells ))
	k_phase = cell_size__phase*float(num_wavelengths) #* float( num_x_cells )/ float( num_x_cells +1)
	
	x_grids = np.arange(0.0, float( num_x_cells ) )
	x_grids_calculation = zeros( num_x_cells )
	x_grids_calculation[:] =  x_grids[:] 
	#x_grids_calculation += 0.5
	x_grids_with_gc = np.arange(-3.0, float( num_x_cells ) + 3.0 )
	
	#x_grids *= cell_size__phase
	k_grids = np.arange(0.0, float(num_x_cells) ); k_grids -= (float(num_x_cells)/2.0 )
	
	
	# take care of the plasma parameters....
	k__vth = wavenumber*vth
	complex_freq = landau_damping_frequency( wp_e, vth, wavenumber )
	w_wave = np.real( complex_freq )
	damping_rate = np.imag( complex_freq )
	
	# setup how we will timestep the wave...
	dt_quantum = cell_size__phase / w_wave 
	if dt_mode_fancy == False:
		# the most simple way of advancing the time-phase.
		#w_wave_phase = w_wave
		pass
	else:
		# more fancy....
		# choose a time step dt such that the phase from w_wave*dt = integral multiple of 'cell_size__phase'
		dt = dt_quantum*float(dt_mode_fancy_num_cells)
		#w_wave_phase = float(dt_mode_fancy_num_cells)*cell_size__phase

	current_amp = 1.0
	time = 0.0
	n = 0
	
	num_sim_steps = int(sim_end_time / dt)
	
	if num_sim_steps % 2 == 1:
		num_sim_steps -= 1
	
	dispersion = zeros ( [ num_sim_steps, num_x_cells] )
	single_spot = zeros ( num_sim_steps)
	
	while time < sim_end_time:
		# draw the sinewave...
		
		print w_wave
		print 2.0*np.pi / w_wave
		print time*w_wave
		print w_wave*time / cell_size__phase
		print dt_quantum
		
		E_field[:] = current_amp*np.sin( k_phase*x_grids_calculation - w_wave*time)
		
		if n < dispersion.shape[0]:
			dispersion[n, 0:num_x_cells] = E_field[:]
			single_spot[n] = E_field[ sample_spot ]
		"""
		# plot the E-field...
		clf()
		subplot(311)
		
		plot( x_grids_with_gc, np.zeros( x_grids_with_gc.shape[0] ), color='purple' )
		plot( x_grids, E_field, marker="*")
		plot( [-3,-2,-1,0], [E_field[-3],E_field[-2], E_field[-1],E_field[0] ] , color='purple', ms='x')
		plot( [float( num_x_cells ) - 1.0 ,float( num_x_cells )+0.0,float( num_x_cells )+1.0,float( num_x_cells )+2.0], [E_field[num_x_cells-1],E_field[0], E_field[1],E_field[2] ] , color='purple')
		#plot( [0,1,2,3], [.5,.5,.5,.5] )
		
		title("E field at time=%f (step=%d)" % ( time, n) )
		ylabel( "Electric Field")
		xlabel( "Position ( $x$ )")
		xlim( [ -3.0, float(num_x_cells)+3.0] )
		subplot(312)
		plot( k_grids, np.log10( np.absolute ( fftshift ( fft (  E_field ) ) ) ) )
		subplot(313)
		plot( k_grids, np.angle( fftshift ( fft (  E_field ) ) ))
		title('DFT Phase')
		xlabel('k (wavenumber)')
		ylabel(' Phase')
		#savefig("./driver_sim/png/E_%06d" % n)
		
		#print fft (  E_field ) 
		#np.save( "./driver_sim/numpy/E_%06d", E_field)
		"""
		print "Finished step %d" % n; time += dt; n += 1;
		
	fft_dispersion = np.log10( fftshift( np.real( ( fft2 ( dispersion ) ) ) ) )
	fft_phase =  fftshift( np.angle( ( fft2 ( dispersion ) ) ) )
	
	clf()
	title("dispersion")
	ylabel("freq")
	xlabel("wavenumber")
	imshow( fft_dispersion, origin="lower", cmap='Paired', aspect='auto')
	savefig("./driver_sim/png/e__w_vs_k_dispersion", dpi=200)
	
	title( " FFT phase")
	imshow( fft_phase, origin="lower", cmap='gray', aspect='auto')
	savefig("./driver_sim/png/e__phase_dispersion", dpi=200)
	
	title("Time versus psotion for E field")
	imshow( dispersion, origin="lower", cmap='Purples', aspect='auto')
	savefig("./driver_sim/png/e__time", dpi=200)
	
	clf()
	plot(single_spot)
	savefig("./driver_sim/png/e__single_spot", dpi=200)
	clf()
	plot( np.log10( np.real( fftshift( fft( single_spot) ) ) ) )
	savefig("./driver_sim/png/e__single_spot_fft", dpi=200)

def plasma_dispersion( value):
	a = scipy.special.wofz(value)
	#a *= 1j
	a *= np.sqrt(np.pi)*1j
	return a
def plasma_dispersion_prime( value):
	z = plasma_dispersion( value )
	return -2.0*(1.0+value*z)

def test_plasma_dispersion():
	#
	# define z = x + j*y  (j = sqrt(-1))
	# now reproduce figures 1-10 in Fried and Conte.. 
	#		This code should identically reproduce those figures...
	#
	
	# do figures 1,2
	( x_values, pd__real_values, pd_prime__real_values, pd__imag_values, pd_prime__imag_values )  = pdtest_plasma_dispersion_test_for_y( 0.0)
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Re(Z')$ versus $x$ for $y=0$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Re(Z')$ ")
	plot( x_values, pd__real_values, label="$\\Re(Z)$")
	plot( x_values, pd_prime__real_values, label="$\\Re(Z')$")
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-2.0, .8) )
	ax.set_yticks(numpy.arange(-2.0,1.2,0.4))
	grid()
	savefig("pd_test__fig_01_real_y_0.png")
	
	clf()
	ax = subplot(111)
	title("$\\Im(Z)$ and $\\Im(Z')$ versus $x$ for $y=0$")
	xlabel("$x$")
	ylabel("$\\Im(Z)$,$\\Im(Z')$ ")
	plot( x_values, pd__imag_values,label="$\\Im(Z)$")
	plot( x_values, pd_prime__imag_values, label="$\\Im(Z')$")
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-1.6, 1.8) )
	xlim( (0.0, 5.0) )
	ax.set_yticks(numpy.arange(-1.6,2.0,0.4))
	grid()
	savefig("pd_test__fig_02_imag_y_0.png")
	
	
	# do figs 3,4,5,6
	(  x_values, pd_r__01, pdp_r__01, pd_i__01, pdp_i__01) = pdtest_plasma_dispersion_test_for_y(1.0)
	(  x_values, pd_r__02, pdp_r__02, pd_i__02, pdp_i__02) = pdtest_plasma_dispersion_test_for_y(2.0)
	(  x_values, pd_r__05, pdp_r__05, pd_i__05, pdp_i__05) = pdtest_plasma_dispersion_test_for_y(5.0)
	(  x_values, pd_r__10, pdp_r__10, pd_i__10, pdp_i__10) = pdtest_plasma_dispersion_test_for_y(10.0)
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$")
	plot( x_values, pd_r__01, label='y=1'); plot( x_values, pd_r__02, label='y=2'); plot( x_values, pd_r__05, label='y=5'); plot( x_values, pd_r__10, label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.425, 0.0) )
	ax.set_yticks(numpy.arange(-0.4,0.05,0.05))
	grid()
	savefig("pd_test__fig_03_real_z_positive_y.png")
	
	
	clf()
	ax = subplot(111)
	title("$\\Im(Z)$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Im(Z)$")
	plot( x_values, pd_i__01, label='y=1'); plot( x_values, pd_i__02, label='y=2'); plot( x_values, pd_i__05, label='y=5'); plot( x_values, pd_i__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (0.0, 0.8) )
	ax.set_yticks(numpy.arange(0.0,0.9,0.1))
	grid()
	savefig("pd_test__fig_04_imag_z_positive_y.png")
	
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$")
	plot( x_values, pdp_r__01,label='y=1'); plot( x_values, pdp_r__02,label='y=2'); plot( x_values, pdp_r__05,label='y=5'); plot( x_values, pdp_r__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.5, 0.2) )
	ax.set_yticks(numpy.arange(-0.5,0.3,0.1))
	grid()
	savefig("pd_test__fig_05_real_zprime_positive_y.png")

	clf()
	ax = subplot(111)
	title("$\\Im(Z')$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Im(Z')$")
	plot( x_values, pdp_i__01,label='y=1'); plot( x_values, pdp_i__02,label='y=2'); plot( x_values, pdp_i__05,label='y=5'); plot( x_values, pdp_i__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.35, 0.0) )
	ax.set_yticks(numpy.arange(-0.35,0.05,0.05))
	grid()
	savefig("pd_test__fig_06_imag_zprime_positive_y.png")
	
	# do figs 7 and 8
	(  x_values, pd_r__m1, pdp_r__m1, pd_i__m1, pdp_i__m1) = pdtest_plasma_dispersion_test_for_y(-1.0)
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Im(Z)$ versus $x$ for $y=-1.0$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Im(Z)$ ")
	plot( x_values, pd_r__m1, label="$\\Re(Z)$"); plot( x_values, pd_i__m1, label="$\\Im(Z)$");
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (-7.0, 9.0) )
	ax.set_yticks(numpy.arange(-6.0,10.0,2.0))
	grid()
	savefig("pd_test__fig_07_real_and_imag_z__with_y_minus_1.png")	
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ and $\\Im(Z')$ versus $x$ for $y=-1.0$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$,$\\Im(Z')$ ")
	plot( x_values, pdp_r__m1, label="$\\Re(Z')$"); plot( x_values, pdp_i__m1, label="$\\Im(Z')$");
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-20.0, 10.0) )
	ax.set_yticks(numpy.arange(-20.0,12.0,4.0))
	grid()
	savefig("pd_test__fig_08_real_and_imag_zprime__with_y_minus_1.png")
	
	# do figs 9 and 10
	(  x_values, pd_r__m12, pdp_r__m12, pd_i__m12, pdp_i__m12) = pdtest_plasma_dispersion_test_for_y(-1.2)
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Im(Z)$ versus $x$ for $y=-1.2$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Im(Z)$ ")
	plot( x_values, pd_r__m12, label="$\\Re(Z)$"); plot( x_values, pd_i__m12, label="$\\Im(Z)$");
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (-12.0, 16.0) )
	ax.set_yticks(numpy.arange(-12.0,20.0,4.0))
	grid()
	savefig("pd_test__fig_09_real_and_imag_z__with_y_minus_1.2.png")	
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ and $\\Im(Z')$ versus $x$ for $y=-1.2$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$,$\\Im(Z')$ ")
	plot( x_values, pdp_r__m12, label="$\\Re(Z')$"); plot( x_values, pdp_i__m12, label="$\\Im(Z')$");
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-40.0, 20.0) )
	ax.set_yticks(numpy.arange(-40.0,24.0,8.0))
	grid()
	savefig("pd_test__fig_10_real_and_imag_zprime__with_y_minus_1.2.png")	
	
def pdtest_plasma_dispersion_test_for_y(y):

	# define z = x + j*y  (j = sqrt(-1))
	
	# plot of real(p.d.), real(p.d. prime) versus x for y =0
	x_values = np.arange(0, 10.0, 0.1)
	pd__real_values = np.zeros( x_values.shape[0] )
	pd_prime__real_values = np.zeros( x_values.shape[0] )
	pd__imag_values = np.zeros( x_values.shape[0] )
	pd_prime__imag_values = np.zeros( x_values.shape[0] )
	
	for i, x in enumerate( x_values ):
		value = complex(x,y)
		if y == 0.0:
			value = x
		
		pd__real_values[i] =  np.real( plasma_dispersion( value ) )
		pd__imag_values[i] = np.imag( plasma_dispersion( value ) )
		pd_prime__real_values[i] =  np.real ( plasma_dispersion_prime( value ) )
		pd_prime__imag_values[i] =  np.imag ( plasma_dispersion_prime( value ) )
	return (  x_values, pd__real_values, pd_prime__real_values, pd__imag_values, pd_prime__imag_values ) 
	
def plot_landau_dmapings( vth_e, wp_e = 1.0):
	# plot from .01 kvth to 1.0
	k_start = .01 / vth_e
	k_end = .9 / vth_e
	ks = np.arange(k_start, k_end, .5)
	
	reals = np.zeros( ks.shape[0] ); imags = np.zeros( ks.shape[0] );
	for i,k in enumerate(ks):
		crap = landau_damping_frequency( wp_e, vth_e, k )
		reals[i] = np.real( crap)
		imags[i] = np.imag( crap)
	
	print "done searching landau space"
	
	ks *= vth_e
	
	clf()
	title( "landau reals")
	stem( ks, reals)
	savefig("landau_reals_vth_%f.png" % vth_e, dpi=200)
	clf()
	title( "landau imags")
	stem( ks, imags)
	savefig("landau_imags_vth_%f.png" % vth_e, dpi=200)

def landau_dmapings_find_closest_k( vth_e, freq, damping, wp_e = 1.0, k_fundimental=1.0):
	k_start = .01 / vth_e
	k_end = .9 / vth_e
	ks = np.arange(k_start, k_end, .01)
	
	
	real_epsilon = 1e10
	imag_epsilon = 1e10
	closest_k_real = -1
	closest_freq = -1
	closest_k_imag = -1
	closest_damping = -1
	
	for i,k in enumerate(ks):
		crap = landau_damping_frequency( wp_e, vth_e, k )
		reals = np.real( crap)
		imags = np.imag( crap)
		delta_real = np.fabs( reals - freq)
		if delta_real < real_epsilon:
			real_epsilon = delta_real
			closest_k_real = k
			closest_freq = reals
		delta_imag = np.fabs( imags - damping)
		if delta_imag < imag_epsilon:
			imag_epsilon = delta_imag
			closest_k_imag = k
			closest_damping = imags
	
	closest_k_real /= k_fundimental
	closest_k_imag /= k_fundimental
	
	print "freq: %f \t closest was %f \t at k=%f (kvth=%f) (error=%e) " % ( freq,closest_freq,closest_k_real,  (closest_k_real*vth_e), real_epsilon)
	print "damp: %f \t closest was %f \t at k=%f (kvth=%f) (error=%e) " % ( damping,closest_damping,closest_k_imag, (closest_k_imag*vth_e),   imag_epsilon)	
	
		
def landau_damping_frequency(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0, inital_root_guess = None ):
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
	
		
	if inital_root_guess == None:
	#	# use the Bohm-Gross dispersion formulas to get an initial guess for w
		inital_root_guess = np.sqrt( wp_e*wp_e + 3*wavenumber*wavenumber*vth__e*vth__e)
	
	epsilon_root = scipy.optimize.newton( plasma_epsilon1, inital_root_guess )
	#epsilon_root = scipy.optimize.newton( plasma_epsilon1, 1.02)
	return epsilon_root * wavenumber * vth__e * np.sqrt(maxwellian_convention_factor)

# given a real w_wave, return the corresponding complex valued k.
def landau_spatial_damping( wp_e, vth__e, w_wave, maxwellian_convention_factor = 2.0, inital_root_guess = None ):
	exp_factor = w_wave / np.sqrt(maxwellian_convention_factor) / vth__e
	chi_factor = wp_e*wp_e/maxwellian_convention_factor / vth__e / vth__e
	
	def plasma_epsilon1( x ):
		pd = plasma_dispersion_prime( x )
		k_guess = exp_factor/ x
		chi_guess = chi_factor / pow( k_guess, 2)
		val = 1.0 - chi_guess*pd
		return val
	
	if inital_root_guess == None:
		# use the Bohm-Gross dispersion formulas to get an initial guess for w
		inital_root_guess = np.sqrt( np.power(w_wave,2.0) - np.power(wp_e,2.0) / ( 3.0*np.power(vth__e,2.0) ))
	
	epsilon_root = scipy.optimize.newton( plasma_epsilon1, inital_root_guess )	
	return exp_factor / epsilon_root

# from krall and trivvelpiece
def estimated_damping(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0):
	
	l_debye = vth__e / wp_e
	k_ldb = l_debye*wavenumber
	
	return np.sqrt( np.pi/8.0)*wp_e/np.power( k_ldb, 3.0)*np.exp( -1.0*(( 1.0/(2.0*k_ldb*k_ldb) ) + 1.5) )

def test__show_dispersion_relation( wp__e, vth__e, k_min, k_max, k_steps):
	(ks, w_real, w_imag) = test__calc_dispersion_relation(wp__e, vth__e, k_min, k_max, k_steps)
	
	k__vthe = np.copy( ks )
	k__vthe *= vth__e
	
	
	# also, calculate the Bohm Gross result/ approximate dmaping result, while we are here...
	# and keep track of the first value above lower_cutoff.. bcause my plasma-dispersion dolver goes off the rails there.
	cut_off_value = 1.05; cut_off_index = 0; 
	bohm_gross_freqs = np.zeros( k__vthe.shape[0] )
	krall_trivelpiece_damping = np.zeros( k__vthe.shape[0] )
	for i,kvthe in enumerate( k__vthe ):
		bohm_gross_freqs[i] = np.sqrt( wp__e*wp__e + 3*kvthe*kvthe)
		krall_trivelpiece_damping[i] = estimated_damping( wp__e, vth__e, ks[i])
		if kvthe > cut_off_value and cut_off_index == 0:
			cut_off_index = i
	
	clf()
	plot(k__vthe,w_real, label='Plasma Dispersion Calculation' )
	plot(k__vthe,bohm_gross_freqs, label='Bohm-Gross Relation' )
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.0, k__vthe[-1] ) ); ylim( ( 0.0, 2.0) );
	savefig('dispersion__real_freq__kinetic_electrostatic_dielectrix', dpi=200)
	
	clf()
	plot( k__vthe, w_imag*-1.0, label='Plasma Dispersion Calculation' )
	plot( k__vthe, krall_trivelpiece_damping, label='Krall-Trivelpiece approx. result' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '-1*Imag frequency (units of $w_p$)' )
	legend(loc=2, prop={'size':10})
	xlim( ( 0.0, k__vthe[-1] ) )
	axhline( y=0, color='black')
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix', dpi=200)
	
	# log versions...
	clf()
	plot( k__vthe, np.log( ( -1.0*w_imag ) ), label='Plasma Dispersion Calculation'  )
	plot( k__vthe, np.log( krall_trivelpiece_damping ), label='Krall-Trivelpiece approx. result' )
	axvline( x=.105, color='black', alpha=.5, label="Lower Limit Allowed in Code")
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '$log_{e}$  of -1*Imag frequency (units of $w_p$)' )
	xlim( ( 0.0, k__vthe[-1] ) )
	ylim( ( -40.0, 3.0) )
	axhline( y=0, color='black')
	legend(loc=4, prop={'size':10})
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__log', dpi=200)
	clf()
	plot( k__vthe, np.log10( ( -1.0*w_imag ) ), label='Plasma Dispersion Calculation' )
	plot( k__vthe, np.log10( krall_trivelpiece_damping ), label='Krall-Trivelpiece approx. result' )
	axvline( x=.105, color='black', alpha=.5, label="Lower Limit Allowed in Code")
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '$log_{10}$ of -1*Imag frequency (units of $w_p$)' )
	xlim( ( 0.0, k__vthe[-1] ) )
	ylim( ( -20.0, 1.0) )
	axhline( y=0, color='black')
	legend(loc=4, prop={'size':10})
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__logten', dpi=200)
	
	# another version.. I can match with come papers..
	clf()
	k__vthe2 = np.copy( k__vthe )
	np.power( k__vthe2 , 2.0, out = k__vthe2 )
	
	print k__vthe
	print k__vthe2
	plot(k__vthe2,w_real, label='Plasma Dispersion Calculation' )
	plot(k__vthe2,bohm_gross_freqs, label='Bohm-Gross Relation' )
	xlabel( '$V_{thermal}*Wavenumber$ squared ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.04, .12 ) ); ylim( ( 1.05, 1.25) );
	savefig('dispersion__real_freq__kinetic_electrostatic_dielectrix__squared_x', dpi=200)
	
	clf()
	plot( k__vthe2, np.log10( w_imag*-1.0 ), label='Plasma Dispersion Calculation' )
	plot( k__vthe2, np.log10(krall_trivelpiece_damping), label='Krall-Trivelpiece approx. result' )
	xlabel( '$V_{thermal}*Wavenumber$ squared ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.04, .12 ) ); ylim( ( -5.0, -1.0) );
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__squared_x', dpi=200)
	
def test__calc_dispersion_relation( wp__e, vth__e, v_the__k_min, v_the__k_max, k_steps):
	k_min = v_the__k_min/vth__e; k_max =  v_the__k_max/vth__e;
	ks = np.linspace(k_min, k_max, num=k_steps)
	w_reals = np.zeros( ks.shape[0]); w_imags = np.zeros( ks.shape[0]);
	print ks
	
	for i,k in enumerate(ks):
		print i
		print k
		w = landau_damping_frequency(wp__e, vth__e, k)
		w_reals[i] = np.real( w )
		w_imags[i] = np.imag( w )
	return (ks, w_reals, w_imags )
	
def plot_epilsion(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0):
	
	chi_e = np.power( ( wp_e/(vth__e*wavenumber) ), 2.0) / maxwellian_convention_factor
	def plasma_epsilon1( x ):
		val = 1.0 - chi_e*plasma_dispersion_prime( x )
		return val
	
	arg_const = 1.0 / (np.sqrt(2) * vth__e * wavenumber )
	bounds_real = np.linspace(2.547486, 2.547488, num=2000) 
	bounds_imag = np.linspace(-0.054607, -0.054609, num=2000)
	result = np.zeros(  (bounds_real.shape[0], bounds_imag.shape[0]) )
	reals = np.zeros( bounds_real.shape[0]*bounds_imag.shape[0] )
	imags = np.zeros( bounds_real.shape[0]*bounds_imag.shape[0] )
	
	
	#bounds_real *= arg_const
	#bounds_imag *= arg_const
	
	flat_idx = 0
	for ix in xrange( 0, bounds_real.shape[0] ):
		for iy in xrange( 0, bounds_imag.shape[0] ):
			#arg = arg_const * np.complex( bounds_real[ix], bounds_imag[iy] )
			arg =  np.complex( bounds_real[ix], bounds_imag[iy] )
			#print  arg
			rrr = plasma_epsilon1 ( arg )
			result[ ix, iy ] = np.absolute( rrr )
			reals[flat_idx] = np.real( rrr )
			imags[flat_idx] = np.imag( rrr )
			flat_idx += 1
	
	w = landau_damping_frequency(wp_e, vth__e, wavenumber)
	w *= arg_const
	wr = np.real(w); wi = np.imag(w);
	
	clf()
	imshow( np.log10( result) ,interpolation='nearest', extent=[bounds_real[0], bounds_real[-1],bounds_imag[0], bounds_imag[-1]], aspect='auto')
	axvline( x = wr); axhline( y = wi); 
	colorbar()
	savefig('aaaa', dpi=400)
	np.save( 'aaaaa', result)
	clf()
	plot( reals, imags)
	savefig('aaaa1')
	
	real_index = int(np.argmin(result) / bounds_imag.shape[0])
	imag_index = np.argmin(result) - real_index* bounds_imag.shape[0]  
	print imag_index
	print "PD min location:          %f + %fj    with value: %e" % (	wr,  	wi , np.absolute( plasma_epsilon1( w) ) )
	print "experiment min location:  %f + %fj    with value: %e" % ( bounds_real[real_index],  bounds_imag[imag_index] , result[real_index, imag_index ] )
	print "global min                                                 %e" % np.min(result)
	
def show_w_spect(freq_axis):
	try:
		data = np.load('e_vs_frquency.npy')
		axis = np.load('e_vs_frquency_freq_axis.npy')
		plot( freq_axis, data)
		xlim([-1.5, 1.5] )
		show()
	except:
		pass

def is_data_in_array_valid(data):
	test__for_invalids_values = np.isfinite( data )
	np.logical_not( test__for_invalids_values , out= test__for_invalids_values )
	if np.any( test__for_invalids_values):
		return False
	return True
	
	
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

def func_t_off(time):
	def func(x, a, b, c, phase, w):
		return a*np.cos(w*(x+time) + phase)*(np.exp(-b*(x+time))) + c
	return func

def func(x, a, b, c, phase, w):
	#return a*np.cos(w*x + phase)*(np.exp(-b*sqrt(2*np.pi)*1.5*x)) + c
	return a*np.cos(w*x + phase)*(np.exp(-b*x)) + c

def func_complex_exp(x, a, b, c, w):
	#return a*np.cos(w*x + phase)*(np.exp(-b*sqrt(2*np.pi)*1.5*x)) + c
	return a*np.cos(w*x)*(np.exp(-b*x)) + c

def func_cos(x, a, b, c, w):
	#return a*np.cos(w*x + phase)*(np.exp(-b*sqrt(2*np.pi)*1.5*x)) + c
	return a*np.cos(w*x + b) + c

def func_exp(x, a, b, c):
	return a*(np.exp(-b*x)) + c
	#return a*np.cos(w_wave*x + phase)*(np.exp(-b*x)) + c

def make_cos2_opt( inital_amp, inital_w):
	def func_decaying_cos_squared( x, b, w):
		cos2 = np.fabs( np.cos((inital_w+w)*x ) )
		return inital_amp*cos2*(np.exp(-b*x))
	return func_decaying_cos_squared

#plot_epilsion( 1.0, 0.042, .333 / 0.042 )
#exit( -1 )

#test__show_dispersion_relation(1.0, 0.014, 0.0, 0.9, 200)
#exit(-1)

#test_plasma_dispersion()
#exit(-1)

#spectrum_test()
#exit(-1)


"""
dirs = ['m3__l4']
files = glob.glob( './m3__l4/numpy/E_*.npy')
dirs = ['m3__l128']
files = glob.glob( './m3__l128/numpy/E_*.npy')
dirs = ['m3__l16']
files = glob.glob( './m3__l16/numpy/E_*.npy')
dirs = ['m3__test']
files = glob.glob( './m3__test/numpy/E_*.npy')
dirs = ['m4__test']
files = glob.glob( './m4__test/numpy/E_*.npy')

#dirs = ['m5__test_e_amp__1e40']
#files = glob.glob( './m5__test_e_amp__1e40/numpy/E_*.npy')

dirs = ['m5__test']
files = glob.glob( './m5__test/numpy/E_*.npy')

"""
#dirs = ['m6__test']
#files = glob.glob( './m6__test/numpy/E_*.npy')

#dirs = ['m7__test7']
#files = glob.glob( './m7__test7/numpy/E_*.npy')

#dirs = ['e_01']
#files = glob.glob( './e_01/numpy/E_*.npy')

#dirs = ['e_norm_04__1e12']
#files = glob.glob( './e_norm_04__1e12/numpy/E_*.npy') 

#dirs = ['e_norm_04__1e12__num_p_double']
#files = glob.glob( './e_norm_04__1e12__num_p_double/numpy/E_*.npy') 

dirs = ['.']
files = glob.glob( './numpy/E_*.npy') 

#dirs = ['mf1__withsolver_enabled_driving']
#files = glob.glob( './mf1__withsolver_enabled_driving/numpy/E_*.npy')

#dirs = ['mf1__withsolver_enabled_driving']
#dirs = ['m3__l128']
#dirs = ['m3__test']

#files = glob.glob( './mf1__withsolver_enabled_driving/numpy/E_*.npy')
#files = glob.glob( './m3__l128/numpy/E_*.npy')
#files = glob.glob( './m3__test/numpy/E_*.npy')

"""
#vals0 = plot_echo_profile( 1.0, 1.5,.5, 1.0)
vte = 1.0
vte1 = .0442
vals = plot_echo_profile( .5/vte, 1.5/vte,1.0/vte, vte)
clf()
plot( vals[0],vals[1] )
#plot( vals0[0], vals0[1] )
axvline(x=0, color='black'); axhline(y=0, color='black')
savefig('crap')
print vals
exit(-1)
"""

#print plasma_dispersion(5.0+.3j)
#print plasma_dispersion_prime(5.0+.3j)
#print "landau"
#print landau_damping_frequency(1.0, .0442, 2.9615)
#print estimated_damping( 1.0, .0442, 2.9615)
#exit(-1)



sim_metadata = pickle.load( open( "./%s/oshun_metadata.p" % (dirs[0],), "rb" ) )
sim_metadata_cfg = pickle.load( open( "./%s/oshun_metadata__cfg.p" % (dirs[0],), "rb" ) )


output_step 				= sim_metadata.output_interval__time_steps   #1
time_step_size 		 	= sim_metadata.dt   									#.106153

box_size__x 				= sim_metadata_cfg.system_size_max__global[0]
num_cells__x 				= sim_metadata_cfg.num_cells__global[0]
#mode_excited 				= 160


print 
if_echo_mode 							= False
if_echo_mode_2nd_order 				= False
if_echo_mode_3rd_order 				= False
wave_driver_enabled					= False
perturbation_driver_enabled		= False
spatial_echo_mode						= False

solve_fit_method_cos_times_exp		= False

driver_stop_time = 0.0
driver_stop_index = 0

if True:

	try:
		mode_excited				= sim_metadata_cfg.E_field_profile.num_wavelengths
		first_pulse_amp				= sim_metadata_cfg.E_field_profile.amp
		second_pulse_mode			= sim_metadata_cfg.E_field_profile.second_pulse_num_wavelengths
		second_pulse_time			= sim_metadata_cfg.E_field_profile.echo_time
		second_pulse_time_index = int( second_pulse_time/time_step_size/output_step )
		second_pulse_amp			= sim_metadata_cfg.E_field_profile.second_pulse_amp
		perturbation_driver_enabled	= True
		w_wave_driven				= 0.0
	except:
		mode_excited				= sim_metadata_cfg.wave_driver.num_wavelengths
		second_pulse_mode			= sim_metadata_cfg.wave_driver.echo_2nd_pulse_num_wavelengths
		second_pulse_time			= sim_metadata_cfg.wave_driver.echo_2nd_pulse_time
		second_pulse_time_index 	= int( second_pulse_time/time_step_size/output_step )
		first_pulse_amp				= sim_metadata_cfg.wave_driver.sustain_amp
		driver_stop_index			= 0
		try:
			w_wave_driven				= sim_metadata_cfg.wave_driver.w_wave
		except:
			w_wave_driven = 0.0
			pass
		
		driver_inital_pulse_amp_change_rate				= sim_metadata_cfg.wave_driver.amp_rate
		#driver_inital_pulse_rise_time					= sim_metadata_cfg.wave_driver.rise_time__time_units
		driver_inital_pulse_starting_amp				= sim_metadata_cfg.wave_driver.starting_amp
		driver_inital_pulse_final_amp					= sim_metadata_cfg.wave_driver.starting_amp
		driver_inital_pulse_sustain_time				= sim_metadata_cfg.wave_driver.sustain_time
		
		wave_driver_enabled		= True
		
	
	#second_pulse_mode			= mode_excited
	#second_pulse_time			= 200.0
	#second_pulse_amp			=	1e-4
	
	
	if second_pulse_mode == mode_excited:
		if_echo_mode_3rd_order = True
		if_echo_mode_2nd_order = False
	else:
		if_echo_mode_3rd_order = False
		if_echo_mode_2nd_order = True
	
	if if_echo_mode_2nd_order:
		echo_2nd_order_time			= float( second_pulse_time ) * float(second_pulse_mode) / ( float(second_pulse_mode) -float( mode_excited) )
		echo_2nd_order_start_index	= int( (echo_2nd_order_time*.9 ) / time_step_size / output_step )
		echo_2nd_order_end_index			= int( (echo_2nd_order_time * 1.2 ) / time_step_size / output_step )
		echo_pulse_time				= echo_2nd_order_time
		echo_pulse_start_index		= echo_2nd_order_start_index
		echo_pulse_end_index			= echo_2nd_order_end_index
	
	if if_echo_mode_3rd_order:
		echo_3rd_order_time			= float( second_pulse_time ) * 2.0*float(second_pulse_mode) / (  2.0*float(second_pulse_mode) -float( mode_excited) )	
		echo_3nd_order_start_index	= int( (echo_3rd_order_time*.9 ) / time_step_size / output_step )
		echo_3nd_order_end_index	= int( (echo_3rd_order_time * 1.2 ) / time_step_size / output_step )
		echo_pulse_time				= echo_3rd_order_time
		echo_pulse_start_index		= echo_3nd_order_start_index
		echo_pulse_end_index			= echo_3nd_order_end_index
		
	if_echo_mode				= True



try:
	if wave_driver_enabled:
		if sim_metadata_cfg.wave_driver.use_as_spatial_echo_driver:
			
			spatial_echo_x1 = float(sim_metadata_cfg.wave_driver.x_for_w1)
			spatial_echo_x2 = float(sim_metadata_cfg.wave_driver.x_for_w2)
			spatial_echo_w1 = float(sim_metadata_cfg.wave_driver.spatial_echo_w1)
			spatial_echo_w2 = float(sim_metadata_cfg.wave_driver.spatial_echo_w2)
			
			spatial_echo_mode = True
			
			x_global_index_for_w1 = int (sim_metadata.driver['x_global_index_for_w1'] )
			x_global_index_for_w2 = int ( sim_metadata.driver['x_global_index_for_w2'] )
			
			spatial_echo_factor = spatial_echo_w1 / (spatial_echo_w2 - spatial_echo_w1)
			spatial_echo_x_plus = spatial_echo_x2 + np.fabs( (spatial_echo_x2 - spatial_echo_x1) ) * spatial_echo_factor
			spatial_echo_x_minus = spatial_echo_x2 - np.fabs( (spatial_echo_x2 - spatial_echo_x1) ) * spatial_echo_factor
			

except:
	pass


density						= sim_metadata_cfg.species_config[0].density_profile.val
vth_e							= sim_metadata_cfg.species_config[0].temp_profile.val # * np.sqrt(1.5)
wp_e							= np.sqrt( density )


# TODO decimation needs to come first before all these calculations!!!!!

k_of_excited_mode 		= 2.0*np.pi / ( box_size__x / float( mode_excited)  )
k__vth_e						= k_of_excited_mode*vth_e
bounce_time				= 2.0*np.pi / np.sqrt( first_pulse_amp*k_of_excited_mode)

fit_start_index 			= int( 0.0/time_step_size/output_step )
fit_end_index 				= int( 100.0/time_step_size/output_step )
fit_offset_index 			= int( 00.0/time_step_size/output_step ) - fit_start_index

print perturbation_driver_enabled
print spatial_echo_mode
print wave_driver_enabled
	
if wave_driver_enabled:
	#sim_metadata_cfg.wave_driver.init( vth_e, time_step_size)
	rise_time = sim_metadata.driver['rise_time__time_units']
	driver_stop_time			= rise_time + sim_metadata_cfg.wave_driver.sustain_time
	driver_stop_index = int( driver_stop_time/time_step_size/output_step )
	#print driver_stop_index
	#print driver_stop_time
	

output_time_step = float(output_step) * time_step_size


f = Formulary()
density__cgs 			= 1e21
vth__cgs 				= f.convert_osi_vth_to_beta ( vth_e ) * 3e10
temp__ev 				= f.convert_osi_vth_to_ev( vth_e) 
wp_e__cgs 				= f.wpe( density__cgs )
skindepth__cgs 		= 3e10/wp_e__cgs
k__cgs 					= k_of_excited_mode / skindepth__cgs
debye_length 			= f.debye_length(temp__ev,  density__cgs)
k_vth__cgs 				= k__cgs*vth__cgs

print ""
print "w of driven wave                                   %f" % (w_wave_driven)
print "dt (time step size)                                %f" % (time_step_size)
print "L (Box Size)                                       %f" % (box_size__x)
print "number spatial cells                               %d" % (num_cells__x)
print "spatial cell size                                  %f" % (box_size__x / float(num_cells__x) )  # this is wrong-ish
print ""

print "v thermal, electrions, normized units from code:   %f" % ( vth_e )
print "k, normalized                                      %e" % k_of_excited_mode
print "bounce time (normaied to 1/wp)                     %e" % bounce_time
print "denisty, cgs                                       %f" % density__cgs
print "skindepth                                          %e cm" % skindepth__cgs
print "vth                                                %e cm/s" % vth__cgs
print "temp                                               %e ev" % temp__ev
print "debye length                                       %e" % debye_length
print "wp                                                 %e" % wp_e__cgs
print "k                                                  %e" % k__cgs
print "kvth                                               %e" % k_vth__cgs
print "k*l_debye                                          %e" % ( k__cgs * debye_length)
print "wp_e/(vth_e*k)                                     %e" % ( wp_e__cgs / ( vth__cgs * k__cgs) )
#ld = landau_damping_frequency( wp_e__cgs, vth__cgs, k__cgs , maxwellian_convention_factor=1.0)
#print ld
#print ld / wp_e__cgs


#plot_landau_dmapings( vth_e)  

index = 0
file_index = 0
current_step = 0
current_time = 0.0

maxs = []
fft_maxs = []
wave_k_vs_time__found_minima = []
snr =[]
files.sort()



filename = os.path.split(files[0])[1]
f1 = np.load( "./%s/numpy/%s" % ( dirs[0], filename) )
pack_factor = 1
data_entry_size = f1.shape[0]

if len(f1.shape) == 2:
	pack_factor = f1.shape[0]
	data_entry_size = f1.shape[1]
print "File packing factor: %d" % pack_factor
print "Dataset size: %d" % data_entry_size


#num_cells__x = data_entry_size
cell_size__x = box_size__x / float( num_cells__x )
sample_point_index 			= num_cells__x/2
max_rows_of_data 				= 20000
limit_data_to_before_time	= -1.0 #20.0
decimation_factor				= 1


# determine the amount of data...
# for now, only full data chunks are written out..
actual_number_of_timestep_datas = 0
if pack_factor == 1:
	actual_number_of_timestep_datas = len(files)
else:
	actual_number_of_timestep_datas = len(files)*pack_factor

sim_stop_time = float( actual_number_of_timestep_datas ) * output_time_step * float( pack_factor)
sim_num_time_steps = actual_number_of_timestep_datas

# is this too much?
if limit_data_to_before_time > 0.0:
	if limit_data_to_before_time < sim_stop_time:
		sim_stop_time = limit_data_to_before_time
		sim_num_time_steps = int( sim_stop_time / output_time_step ) - 1
		if sim_num_time_steps % 2 == 1:
			sim_num_time_steps -= 1

# enforce decimation
if sim_num_time_steps > max_rows_of_data:
	while True:
		decimation_factor += 1
		num_iterms_after_decimation = float(sim_num_time_steps) / float(decimation_factor)
		if num_iterms_after_decimation <= float(max_rows_of_data):
			# the current decimenation factor will due the job!
			break
		pass
	
	sim_num_time_steps = int( np.floor( float(sim_num_time_steps) / float(decimation_factor) ) )
	output_time_step *= float( decimation_factor )
	output_step *= int( decimation_factor )
	sim_stop_time = output_time_step * float( sim_num_time_steps )
	print "Data limited to  %d time entries... therefore, only every %dth files was used, giving you %d time points." % ( max_rows_of_data,decimation_factor, sim_num_time_steps ) 

dispersion_time_stop = sim_stop_time
dispersion_num_time_steps = sim_num_time_steps

print "steps between_output (output_step)                 %d" % (output_step)
print "Sim Time between output steps (output_time_step)   %e" % ( output_time_step ) 
	
#dispersion_time_stop = 300.0
#sim_stop_time = float( len( files ) ) * output_time_step * float( pack_factor)

#if sim_stop_time < dispersion_time_stop:
#	dispersion_time_stop = sim_stop_time
#dispersion_num_time_steps = dispersion_time_stop / output_time_step

#dispersion_num_time_steps = int(dispersion_num_time_steps) - 1
#if dispersion_num_time_steps % 2 == 1:
#	dispersion_num_time_steps -= 1

dispersion = np.zeros( (dispersion_num_time_steps, data_entry_size) )

for filen in files:
	#try:
	filename = os.path.split(filen)[1]
	
	
	#print filename
	f1 = np.load( "./%s/numpy/%s" % ( dirs[0], filename) )
	
	if pack_factor == 1:
		if file_index % decimation_factor > 0:
			#print "skipped"
			file_index += 1
			continue
		#print "used as %d (but really file %d)" % ( index, file_index)
		
		if index < dispersion.shape[0]:
			if is_data_in_array_valid( fl ):
				dispersion[index,0:f1.shape[0]] = f1[:]
			else:
				print "Invalid data found... stopping reading input data..."
				dispersion_num_time_steps = index
				if dispersion_num_time_steps%2 == 1:
					dispersion_num_time_steps -= 1
				dispersion_time_stop = float( dispersion_num_time_steps ) * output_time_step
				sim_stop_time = float( dispersion_num_time_steps ) * output_time_step
				
				dispersion_new = zeros( (dispersion_num_time_steps, data_entry_size) )
				dispersion_new[:,:] = dispersion[0:dispersion_num_time_steps,:]
				dispersion = None
				dispersion = dispersion_new
				break
		else:
			break
		#snr.append( scipy.stats.signaltonoise(f1) )
		fft_data = np.absolute(  numpy.fft.fft( numpy.fft.fftshift(f1 )) )
	
		maxs.append( f1[sample_point_index] )
		fft_maxs.append( fft_data[mode_excited]  ) 
	
		file_index += 1
		index += 1
		current_step += output_step
		current_time += output_time_step
	else:
		for ir in xrange(0, pack_factor):
			
			if file_index%decimation_factor != 0:
				#print "skipped"
				file_index += 1
				continue
			#print "used as %d (but really file %d)" % ( index, file_index)
			
			if index < dispersion.shape[0]:
				if is_data_in_array_valid( f1[ir,:] ):
					dispersion[index,0:data_entry_size] = f1[ir,:]
				else:
					print "Invalid data found... stopping reading input data..."
					dispersion_num_time_steps = index
					if dispersion_num_time_steps%2 == 1:
						dispersion_num_time_steps -= 1
					# rest the size of main array!
					dispersion_time_stop = float( dispersion_num_time_steps ) * output_time_step
					sim_num_time_steps	= dispersion_num_time_steps
					dispersion_new = zeros( (dispersion_num_time_steps, data_entry_size) )
					dispersion_new[:,:] = dispersion[0:dispersion_num_time_steps,:]
					dispersion = None
					dispersion = dispersion_new
					break
			else:
				break
			
			#snr.append( scipy.stats.signaltonoise(f1[ir,:]) )
			fft_data = np.absolute(  numpy.fft.fft( numpy.fft.fftshift(f1[ir,:] )) )
		
			maxs.append( f1[ir,sample_point_index] )
			fft_maxs.append( fft_data[mode_excited]  ) 
		
			#print "\t\t %d" % ir
			file_index += 1
			index += 1
			current_step += output_step
			current_time += output_time_step
			
	#except:
	#	break
	
	if len(maxs) > dispersion.shape[0]:
		del fft_maxs[-1]
		del maxs[-1]
		
print ""
print "Done reading files"
print ""
np.save( 'e_vs_time__efield_t_vs_x', dispersion) 


#fft2_dis = np.log( ( np.absolute( fftshift( fft ( dispersion,axis=1 ) ) ) ) )



#fft2_dis[:,5] = 1.0 #create a vettical stripe (along a row)
#fft2_dis[5,:] = 1.0 #create a horizontal stripe ( along a col)
# Therfore, dipsersion[row, col], so dispersion[index, :] = f1[:] fills in E at a certain time (along rows)
# ftt along axis 0 mean along the cols (y axis) so it converts t--> w
# ftt along axis 1 mean along the cols (x axis) so it converts x --> k
#fft2_dispersion = np.log10( fftshift( np.absolute( ( fft2 ( dispersion, s=(256,256) ) ) ) ) )
print dispersion.shape

fft2_dispersion = np.log10( fftshift( np.absolute( ( fft2 ( dispersion) ) ) ) )

fft2_time_vs_k = fft (dispersion, axis=1)

try:
	imaginary_k__vs_time = fft2_time_vs_k[:, mode_excited]
	imaginary_k__vs_time__echo = fft2_time_vs_k[:, second_pulse_mode]
except:
	imaginary_k__vs_time = None
	imaginary_k__vs_time__echo = None

fft2_time_vs_k = ( ( np.absolute( fftshift(fft2_time_vs_k, axes=1 ) ) ) )
#fft2_time_vs_k = ( ( np.absolute( fft2_time_vs_k ) ) )
#fft2_time_vs_k[:, fft2_time_vs_k.shape[1]/2 + 160] = .024


fft2_time_vs_k_phase = np.angle(  fftshift( fft ( dispersion, axis=1), axes=1 ) )
#fft2_time_vs_k_phase = np.angle( ( fft ( dispersion, axis=1)) )
fft2_w_vs_x = ( ( np.absolute( fftshift( fft ( dispersion, axis=0), axes=0 ) ) ) )

# create a normalized, centered axis for frequency. This axis has a constant cell size of delta_freq=2.0/(num_time_data_points).
# 		It is normalized such that it spans from -1.0 goes to 0.0 then continues to (1.0-delta_freq). We will use the 'fftshift' data to
#		reorder the fft data to nathc this centered arrangment. the value 1.0 is the normalized Niquist frequency.
#			For even number of time points (which we enforce..) buffer[0] has the largest 'negitive' frequency... the Niquist freq. 
#				Buffer[N/2] = average (w=0 mode)
#				Incoming data is real, so the FFT will be Hermetian. This implies that the w=0 mode and Niquist mode must be real. All other
#					vlaues will be complex in general.. but points summetic acoss the N/2 buffer location will always be complex cinjugates of eachother.
#centered_frequency_axis_norm = np.arange(0, 2.0, 2.0/dispersion.shape[0]); centered_frequency_axis_norm -= 1.0; centered_frequency_axis_norm *= f_nyquist
#centered_k_axis_norm = np.arange(0, dispersion.shape[1]); centered_k_axis_norm -= ( dispersion.shape[1]/2 )

# create axes that give real values (using Oshun's usual normalizations )
freq_axis = numpy.fft.fftshift( numpy.fft.fftfreq( fft2_dispersion.shape[0], d = (output_time_step/ (2*np.pi*wp_e)) ) ) 
k_axis = numpy.fft.fftshift( numpy.fft.fftfreq( fft2_dispersion.shape[1],  d = (cell_size__x/ (2*np.pi*wp_e)) ) )
t_axis = np.arange(0, float(dispersion_num_time_steps)*float(output_time_step), float(output_time_step))
x_axis = np.arange(0, 1.0, 1.0/dispersion.shape[1] )
x_axis *= box_size__x

# also create an E-field energy diagnostic..... by integrating E^2/2 over all x at each time step. 
E_squared = np.multiply( dispersion, dispersion)
E_squared /= 2.0
E_wave_energy = scipy.integrate.simps( E_squared, dx = cell_size__x , axis=1)

clf()
title("Wave Energy versus time."); xlabel( "Time"); ylabel("Wave Energy in E Field");
plot(t_axis[0:E_wave_energy.shape[0]], np.log10( E_wave_energy) ) 
savefig('!wave_energy_vs_time')

clf()
freq_axis_sq = np.copy( freq_axis )
freq_axis_sq /= 2.0
data_start_time = 2.0
data_end_time = 20.0
data_start_index 	= int(data_start_time / time_step_size/output_step)
data_end_index 		= int(data_end_time / time_step_size/output_step)
fft_of_energy = np.log10( numpy.fft.fftshift( np.absolute( numpy.fft.fft(E_wave_energy[data_start_index:data_end_index], n=E_wave_energy.shape[0]) ) ) ) 

maxima = find_maxima( fft_of_energy, start_index = fft_of_energy.shape[0]/2, sort_by_freq = True )

title("Log FFT of Wave Energy versus time."); xlabel( "Time"); ylabel("FFT of Wave Energy in E Field");
plot( freq_axis_sq, fft_of_energy)
#xlim([-5.0, 5.0] )
for i in xrange(0, len(maxima)):
	print "Wave energy freq hump at %f " % freq_axis_sq[maxima[i][0]]
	axvline( x = freq_axis_sq[maxima[i][0]] )
	if i > 2:
		break
savefig('!wave_energy_fft')

# also attempt a fit of the maxima of the wave energies...
maxima = find_maxima( E_wave_energy[data_start_index:data_end_index], start_index = 0, sort_by_freq = True )
# set what the t=0 zero point of the fit is...
if len(maxima) > 0:
	if True: #try:
		t_equals_zero_index =  maxima[0][0]
		#t_equals_zero_index =  driver_start_index
		#t_equals_zero_index =  data_start_index
		energy_fit_points__time = np.zeros( ( len(maxima) ) )
		energy_fit_points__value = np.zeros( ( len(maxima) ) )
		print energy_fit_points__time.shape
		print energy_fit_points__value.shape
		for i in xrange(0, len(maxima)):
			energy_fit_points__time[i] = t_axis[ maxima[i][0] - t_equals_zero_index]
			energy_fit_points__value[i] = maxima[i][1]
		popt, pcov = curve_fit( func_exp ,energy_fit_points__time,  energy_fit_points__value, p0= np.array( [energy_fit_points__value[0], .1, 0.0] ) )
		curve_fits = popt[0]*(np.exp(-popt[1]*energy_fit_points__time)) + popt[2]
		t_fit_data__put_back = np.copy( energy_fit_points__time )
		t_fit_data__put_back += t_axis[ t_equals_zero_index + data_start_index]
		print t_equals_zero_index
		print t_fit_data__put_back
		#exit(-1)
		
		clf()
		title("Wave Energy versus time."); xlabel( "Time"); ylabel("Wave Energy in E Field");
		plot(t_axis[0:E_wave_energy.shape[0]], np.log10( E_wave_energy) ) 
		plot(t_fit_data__put_back, np.log10( curve_fits) ) 
		xlim( (data_start_time, data_end_time))
		savefig('!wave_energy_vs_time__fit')
		print "Wave energy fit for the damping rate %e" % ( popt[1]/2.0)
		print ""
	#except:
	#	print ""
	#	print "Error doing the Wave energy fit."
	#	print ""

#print "pahse"
#print imaginary_k__vs_time__echo.shape[0]
#print imaginary_k__vs_time.shape[0]
#print second_pulse_time_index

try:
	clf()
	title('Excited K-mode in complex plane\n plotted over time')
	ylabel('Imaginary K')
	xlabel('Real K')
	plot ( np.imag( imaginary_k__vs_time[second_pulse_time_index:] ), np.real( imaginary_k__vs_time[second_pulse_time_index:] ) )
	savefig('!wave_phasor__excited__modes')
	
	index_jump = int (imaginary_k__vs_time.shape[0] / 100)
	if index_jump == 0:
		index_jump = 1
	index = index_jump; i = 0
	while index < imaginary_k__vs_time.shape[0]:
		clf()
		plot ( np.imag( real_k__vs_time[0:index] ), np.real( imaginary_k__vs_time[0:index] ) )
		plot ( np.imag( imaginary_k__vs_time[index-2:index] ), np.real( imaginary_k__vs_time[index-2:index] ), color='red' )
		savefig('./phase/phase__%05d'% i)
		index += index_jump
		i += 1
except:
	pass

if True:
	# try to constuct an nice phasor for the main Kmode and driver in same plot....
	if wave_driver_enabled:
		
		# construct the driver's phasor...
		fake_final_driver_amp = np.absolute( imaginary_k__vs_time[driver_stop_index-1] )
		fake_init_driver_amp = np.absolute( imaginary_k__vs_time[1] )
		
		fake_final_driver_amp_rate = fake_final_driver_amp / float(driver_stop_index) / output_time_step
		def driver_phasor(t, ampme):
			#amp = fake_init_driver_amp + fake_final_driver_amp_rate*t
			#amp = 
			#print amp
			#if amp > fake_final_driver_amp:
			#	amp = fake_final_driver_amp
			driver_phasor1 = np.complex( np.cos( w_wave_driven*t), np.sin( w_wave_driven*t) )
			return ampme*driver_phasor1
		"""
		mode_excited
		driver_stop_index
		first_pulse_amp
		w_wave_driven
		driver_inital_pulse_amp_change_rate
		driver_inital_pulse_rise_time
		driver_inital_pulse_sustain_time
		driver_inital_pulse_starting_amp
		driver_inital_pulse_final_amp
		"""
		
		drive_pahsor = np.zeros( driver_stop_index, dtype=complex )
		drive_pahsor_r = np.zeros( driver_stop_index )
		drive_pahsor_i = np.zeros( driver_stop_index )
		temp_time = 0
		for i in xrange( 0, driver_stop_index):
			temp = driver_phasor( temp_time, np.absolute( imaginary_k__vs_time[i]) )
			drive_pahsor[i] = temp
			drive_pahsor_r[i] = np.real( temp )
			drive_pahsor_i[i] = np.imag( temp )
			temp_time += output_time_step
		# first plot the K mode response from time=0 until time sriver is shut off...
		clf()
		subplot(211)
		title('Excited K-mode and Wave Driver in complex plane\n plotted over time')
		ylabel('Imaginary K')
		
		plot ( np.imag( imaginary_k__vs_time[0:driver_stop_index] ), np.real( imaginary_k__vs_time[0:driver_stop_index] ) )
		plot ( drive_pahsor_r, drive_pahsor_i, color='red') 
		subplot(212)
		temp1 = np.angle( imaginary_k__vs_time[0:driver_stop_index] )
		temp2 = np.angle( drive_pahsor )
		plot( t_axis[0:driver_stop_index], temp1, color='blue')
		plot( t_axis[0:driver_stop_index], temp2, color='red')
		title('Real K (final phase diff=%e' % (temp2[-1]-temp1[-1] ) )
		savefig('!wave_phasor__excited__modes_with_driver')
		clf()
		title(' phase difference betwwen driver and response vs time')
		phase_diff = temp2-temp1
		plot( t_axis[0:driver_stop_index], phase_diff )
		savefig('!wave_phasor__excited__modes_with_driver__diff')
		


#except:
#	pass


constant_postion_index =  sample_point_index
temp = np.copy( dispersion[0,:] )
if np.abs(temp[constant_postion_index+1]) > temp[constant_postion_index]:
	temp_max = temp[constant_postion_index]
	#current_index = constant_postion_index
	# first check that we are doing + numbers
	while temp[constant_postion_index+1] < 0.0:
		constant_postion_index+= 1
	
	while np.abs(temp[constant_postion_index+1]) > temp[constant_postion_index]:
		constant_postion_index += 1
		temp_max = temp[constant_postion_index+1]
		#print constant_postion_index
		
	#constant_postion_index = current_index

#show_w_spect(centered_frequency_axis)

"""
for i in xrange(0,20):
	clf()
	plot( fft2_time_vs_k_phase[i,:] )
	#print np.max ( fft2_time_vs_k_phase )
	#print fft2_time_vs_k_phase[i,:] 
	title( "Phase lineout at time index %05d" % i)
	savefig('crap%05d' % i)
"""
#imshow( fft2_dis, origin="lower", cmap='Purples', aspect='auto' )

clf()
imshow( fft2_dispersion, origin="lower", cmap='Paired', aspect='auto', extent = (k_axis[0], k_axis[-1], freq_axis[0],freq_axis[-1] ) )
title('Dispersion Relation')
ylabel( "Frequency (units of $w_p$)")
xlabel( "Wavenumber (k)")
colorbar()
savefig("e_vs_time__dispersion", dpi=400)
ylim( [-2.0,2.0] )
savefig("e_vs_time__dispersion_zoom2", dpi=400)

clf()
imshow( ( fft2_time_vs_k) , origin="lower", cmap='flag_r', aspect='auto', extent = (k_axis[0], k_axis[-1], t_axis[0],t_axis[-1] ) )
title('Time versus Wavenumber (K)')
ylabel( "Time (units of $w_p$)")
xlabel( "Wavenumber (k)")
colorbar()
savefig("e_vs_time__efield_time_vs_k", dpi=200)
clf()
imshow( np.log10(fft2_time_vs_k), origin="lower", cmap='Paired', aspect='auto', extent = (k_axis[0], k_axis[-1], t_axis[0],t_axis[-1] ) )
title('Time versus Wavenumber (K)')
ylabel( "Time (units of $w_p$)")
xlabel( "Wavenumber (k)")
colorbar()
savefig("e_vs_time__efield_time_vs_k__high_contrast", dpi=200)

clf()
imshow( fft2_time_vs_k_phase, origin="lower", cmap='gray', interpolation='nearest', aspect='auto', extent = (k_axis[0], k_axis[-1], t_axis[0],t_axis[-1] ) )
title('Time versus Wavenumber (K) PHASE')
ylabel( "Time (units of $w_p$)")
xlabel( "Wavenumber (k)")
colorbar()
savefig("e_vs_time__efield_time_vs_k__phase", dpi=400)

clf()
imshow( np.log(fft2_w_vs_x), origin="lower", cmap='Accent', aspect='auto', extent = (x_axis[0], x_axis[-1], freq_axis[0],freq_axis[-1] ) )
title('Frequency ( $w$ ) versus Position ( $x$ )\n k*vth=%f' % k__vth_e)
ylabel( "Frequency (w - normalized to Nyquist)")
xlabel( "Position ( $x$ normalized to box size)")
colorbar()
savefig("e_vs_time__efield_w_vs_x", dpi=200)

clf() # PuOr PRGn
imshow( ( dispersion) , origin="lower", cmap='Purples',interpolation='nearest', aspect='auto', extent = (x_axis[0], x_axis[-1], t_axis[0],t_axis[-1] ) )
title('Time versus Position ( $x$ )')
ylabel( "Time (units of $w_p$)")
xlabel( "Position ( units if $\\frac{c}{w_p}$ )")
colorbar()
#xlim( [0.0, 3.0 ] )
savefig("e_vs_time__efield_t_vs_x", dpi=400)
clf() # PuOr PRGn
imshow( ( dispersion) , origin="lower", cmap='flag',interpolation='nearest', aspect='auto', extent = (x_axis[0], x_axis[-1], t_axis[0],t_axis[-1] ) )
title('Time versus Position ( $x$ )')
ylabel( "Time (units of $w_p$)")
xlabel( "Position ( units if $\\frac{c}{w_p}$ )")
colorbar()
#xlim( [0.0, 3.0 ] )
savefig("e_vs_time__efield_t_vs_x__high_contrast", dpi=400)


clf()
time_start_index 	= int(20.0 / time_step_size/output_step)
time_start_index = 0
title('FFt of E field at mid-box versus time frequency'  )
fft_of_e_field_at_mid_box = fftshift( np.absolute(  fft( dispersion[time_start_index:, constant_postion_index], n=dispersion.shape[0] ) ) )
maxima = find_maxima( fft_of_e_field_at_mid_box, start_index = fft_of_e_field_at_mid_box.shape[0]/2, sort_by_freq = True )
plot( freq_axis, np.log10( fft_of_e_field_at_mid_box) )
#plot( freq_axis, np.log10(  fftshift( np.absolute(  fft( E_wave_energy ) ) ) ) )
xlim( (-5.0, 5.0 ) )

try:
	print ""
	print "The first 3 Temporal Frequency Maxima in the Fourier Spectrum of the E field at the center of the Box"
	print "------------"
	print "(Temporal) Frequency Mode Number: %d" % (maxima[0][0] - freq_axis.shape[0] / 2)
	print "(Temporal) Frequency              %f" % freq_axis[ maxima[0][0] ]
	print "Magnitude                         %e" % float( maxima[0][1] )
	print "------------"
	print "(Temporal) Frequency Mode Number: %d" % (maxima[1][0] - freq_axis.shape[0] / 2)
	print "(Temporal) Frequency              %f" % freq_axis[ maxima[1][0] ]
	print "Magnitude                         %e" % float( maxima[1][1] )
	
	print "(Temporal) Frequency Mode Number: %d" % (maxima[2][0] - freq_axis.shape[0] / 2)
	print "(Temporal) Frequency              %f" % freq_axis[ maxima[2][0] ]
	print "Magnitude                         %e" % float( maxima[2][1] )
	print "------------"
	print ""
	
	axvline( x = freq_axis[ maxima[0][0] ] )
	axvline( x = freq_axis[ maxima[1][0] ] )
	axvline( x = freq_axis[ maxima[2][0] ] )
except:
	pass
savefig("e_vs_time__central_box__fft")
temp = np.absolute(  fft( dispersion[:, constant_postion_index] ) )
fft_maximum_time_freq = np.argmax ( temp )

# attempt to get the minima in the fixed wavenumber versus time plot to use in determining
#	termporal frequencies in the data.
#try:
if True:
	
	wave_k_vs_time__time_between_found_minima = []
	wave_k_vs_time__time_between_found_minima__implied_freqs = None
	fft_maxs_array = np.array( fft_maxs)
	minima = find_maxima( fft_maxs_array , minima = True, sort_by_freq = True)
	
	
	for i in xrange(0, len(minima)):
		minima_time_location =  t_axis[ minima[i][0] ]
		wave_k_vs_time__found_minima.append( ( minima[i][0], minima_time_location,  minima[i][1] ) )
		if i > 0 and minima_time_location > driver_stop_time:
			#print wave_k_vs_time__found_minima[i][1]
			delta_time = wave_k_vs_time__found_minima[i][1] - wave_k_vs_time__found_minima[i-1][1]
			wave_k_vs_time__time_between_found_minima.append( delta_time )
			#print "delta time %e" % delta_time
	
	wave_k_vs_time__time_between_found_minima = np.array ( wave_k_vs_time__time_between_found_minima )
	if wave_k_vs_time__time_between_found_minima.shape[0] > 0:
		#wave_k_vs_time__time_between_found_minima__implied_freqs = np.copy( wave_k_vs_time__time_between_found_minima )
		wave_k_vs_time__time_between_found_minima__implied_freqs = np.unique(wave_k_vs_time__time_between_found_minima.round(decimals=4))
		for i in xrange(0,wave_k_vs_time__time_between_found_minima__implied_freqs.shape[0]):
			temp = 2.0*np.pi / wave_k_vs_time__time_between_found_minima__implied_freqs[i]
			temp *= .5
			wave_k_vs_time__time_between_found_minima__implied_freqs[i] = temp
	
	if wave_k_vs_time__time_between_found_minima__implied_freqs!= None and wave_k_vs_time__time_between_found_minima__implied_freqs.shape[0] > 0:
		print ""
		print "------------"
		print "Unique Time Frequencies implied time domain analyisis of fixed K (k of driver) versus Time data."
		for i in xrange(0, wave_k_vs_time__time_between_found_minima__implied_freqs.shape[0]):
			print "Frequency                     %f" % wave_k_vs_time__time_between_found_minima__implied_freqs[i]
		print ""
	
	if len( wave_k_vs_time__found_minima) > 0:
		print "As diagnostic of the driver, the first few minima in the fixed K-mode E-field verus time are: "
		for i in xrange(0, 30):
			time_occured = wave_k_vs_time__found_minima[i][1]
			e_field_value_at_minima = wave_k_vs_time__found_minima[i][2]
			if time_occured >= driver_stop_time:
				break
			print " minima %d:		time: %f				E-Field value: %e" % ( i, time_occured, e_field_value_at_minima)
	idx = np.argmax( fft_maxs_array )
	val_at_idx = fft_maxs_array[idx]
	print "The Absolut maximum occcured at time: %f (index=%d)			E-Field value: %e" % ( t_axis[idx], idx, val_at_idx)
	print ""
	
#except:
#	wave_k_vs_time__found_minima = []
#	wave_k_vs_time__time_between_found_minima = np.zeros( 0 )
#	wave_k_vs_time__time_between_found_minima__implied_freqs = np.zeros( 0 )
#	pass

clf()
title('Wave k versus time'); xlabel( "Time (units of $w_p$)"); ylabel('magnitude of k mode');
plot(t_axis[0:len(fft_maxs)], fft_maxs )
try:
	# plot the minima (if any.. and if the claculation was sucessfull...) on the fixed k versus time plot.
	for i in xrange(0, len(wave_k_vs_time__found_minima)):
		axvline( x=wave_k_vs_time__found_minima[i][1] )
except:
	pass
savefig("wave_k_vs_time")

print ";''''''''''''''''''''''''''''''''''''''''''''''"
if True:
	clf()
	title('E field at mid-box versus time Time')
	xlabel( "Time (units of $\\frac{1}{w_p}$)")
	ylabel( "E field")
	
	fucktard = dispersion[:,constant_postion_index]
	plot( t_axis[0:fucktard.shape[0]], fucktard)


	# data_fit.....
	# massage data to fit...
	if fit_end_index > fucktard.shape[0]:
		fit_end_index = fucktard.shape[0]
	
	
	if spatial_echo_mode: 
		#qucktest...
		index_of_space_echo = int ( spatial_echo_x_plus / cell_size__x + 1 )
		index_of_space_echo_minus = int ( spatial_echo_x_minus / cell_size__x + 1)
		clf()
		
		space_echo_data = dispersion[:, index_of_space_echo-1 ]; space_echo_data_minus = dispersion[:, index_of_space_echo_minus];
		space_echo_data_time = t_axis[0:space_echo_data.shape[0]]
		start_index = np.argmax(space_echo_data); start_index_minus = np.argmax(space_echo_data_minus);
		
		popt, pcov = curve_fit( func_cos ,space_echo_data_time[start_index:],  space_echo_data[start_index:], p0= np.array( [ np.max(space_echo_data), .1, 0.0, 1.0 ] ) )
		curve_fits =  popt[0]*np.cos(popt[3]*space_echo_data_time[start_index:]+ popt[1]) + popt[2]
		print "ECHO"
		print popt[3]
		
		#popt, pcov = curve_fit( func_cos ,space_echo_data_time[start_index_minus:],  space_echo_data_minus[start_index_minus:], p0= np.array( [ np.max(space_echo_data_minus), .1, 0.0, 1.0 ] ) )
		#curve_fits =  popt[0]*np.cos(popt[3]*space_echo_data_time[start_index_minus:]+ popt[1]) + popt[2]
		#print "ECHO MINUS"
		#print popt[3]
		
		clf(); title('Spatital echo at %f,%f (E versus time)' % (spatial_echo_x_plus, spatial_echo_x_minus) ); xlabel('Time');ylabel('Electric Field');
		plot(space_echo_data_time, space_echo_data ); plot(space_echo_data_time, space_echo_data_minus ); savefig( 'space_echo')
		
		clf(); title('Spatital echo data (E versus time) with fit overlay\n freq fit as %f' % popt[3]); xlabel('Time');ylabel('Electric Field');
		plot(space_echo_data_time, space_echo_data ); plot(space_echo_data_time[start_index:], curve_fits );savefig( 'space_echo_with_fit')	
		#plot( space_echo_data_time, dispersion[:, index_of_space_echo+1 ] )
		#plot( space_echo_data_time, dispersion[:, index_of_space_echo-1] )
		
		
	elif wave_driver_enabled:
		
		# get the index of when the driver stopped....
		time_index_offset = driver_stop_index
		#time_index_offset = int(900.0 / time_step_size/output_step)
		#fit_end_index = int(950.0 / time_step_size/output_step)
		data_to_fit = np.copy ( dispersion[time_index_offset:fit_end_index, constant_postion_index] )
		t_fit_data = np.copy( t_axis[0:fit_end_index-time_index_offset] )
		
		
		# now get the next max/min
		min = 0.0
		min_index = 0
	
		for ix in xrange(0, data_to_fit.shape[0] ):
			val = data_to_fit[ix]
			if np.fabs(val) > min:
				min = val
				min_index = ix
			else:
				break
		min_index += 0
		
		ext_is_min = True
		if data_to_fit[min_index] < 0.0:
			ext_is_min = False
		
		data_to_fit = np.copy( data_to_fit[min_index:] )
		t_fit_data = np.copy ( t_fit_data[0:-min_index] )
		
		
		print t_fit_data.shape
		print data_to_fit.shape
		print fit_end_index
		print driver_stop_index
		print min_index
		print ""
		print "------------"
		idx_of_max = np.argmax( dispersion[:, constant_postion_index] )
		max_val = dispersion[idx_of_max, constant_postion_index]
		time_at_max = t_axis[idx_of_max]
		print "The maximum E-field value (within spatial cell %d) is %e and occured at time %f (time index %d)" % (constant_postion_index, max_val, time_at_max, idx_of_max)
		print ""
		
		if solve_fit_method_cos_times_exp == False:
			
			if ext_is_min:
				# seek out minima of the curve and fit to expoential...
				minima = find_maxima( data_to_fit , minima = True, sort_by_freq = True)
				exp_fit_data__time = np.zeros( len( minima ) )
				exp_fit_data__amp  = np.zeros( len( minima ) )
				for i in xrange(0, len(minima)):
					exp_fit_data__time[i] = t_axis[ minima[i][0] ]
					exp_fit_data__amp[i]  = minima[i][1]
				
			else:
				# seek out maxima of the curve and fit to expoential...
				maxima = find_maxima( data_to_fit , minima = False, sort_by_freq = True)
				maxima.insert(0, (0, data_to_fit[0]) )
				exp_fit_data__time = np.zeros( len( maxima ) )
				exp_fit_data__amp  = np.zeros( len( maxima ) )
				for i in xrange(0, len(maxima)):
					exp_fit_data__time[i] = t_axis[ maxima[i][0] ]
					exp_fit_data__amp[i]  = maxima[i][1]
			
			popt, pcov = curve_fit( func_exp ,exp_fit_data__time,  exp_fit_data__amp, p0= np.array( [data_to_fit[0], .1, 0.0] ) )
			popt1 = np.zeros([ len(popt) + 1]); popt1[0:3] = popt[:]; popt1[-1] = 1.0; popt = popt1; 
			curve_fits = popt[0]*(np.exp(-popt[1]*exp_fit_data__time)) + popt[2]
			t_fit_data__put_back = np.copy( t_fit_data)
			t_fit_data__put_back_fit = exp_fit_data__time
		else:
			popt, pcov = curve_fit( func_complex_exp ,t_fit_data,  data_to_fit, p0= np.array( [data_to_fit[0], .1, 0.0, 1.0 ] ) )
			curve_fits =  popt[0]*np.cos(popt[3]*t_fit_data)*(np.exp(-popt[1]*t_fit_data)) + popt[2]
			t_fit_data__put_back = np.copy( t_fit_data)
			t_fit_data__put_back_fit = np.copy( t_fit_data)
		axis = subplot(111)
		print data_to_fit[0]
		print t_axis[min_index+time_index_offset]
		axvline( x=t_axis[min_index+time_index_offset], color='black' )
		#exit(-1)
		t_fit_data__put_back += ( t_axis[min_index+time_index_offset] )
		t_fit_data__put_back_fit += ( t_axis[min_index+time_index_offset] )
		plot(t_fit_data__put_back, data_to_fit)
		plot(t_fit_data__put_back_fit, (curve_fits))
		axhline( y=0, color='black')
		xlim([0.0, 25.0])
		#ylim([-1e-31, 1e-31])
		#axis.autoscale_view(True,True,True)
		savefig("e_vs_time__central_box")
	elif perturbation_driver_enabled:
		if True: #try:
			time_index_offset 	= int(2.0 / time_step_size/output_step)
			fit_end_index 		= int(10.0 / time_step_size/output_step)
			data_to_fit = np.copy ( dispersion[time_index_offset:fit_end_index, constant_postion_index] )
			t_fit_data = np.copy( t_axis[0:fit_end_index-time_index_offset] )
			#time_index_offset = 0
			print time_index_offset
			print fit_end_index
			
			print data_to_fit.shape[0]
			print ""
			# now let ignore the first max....
			min = 1e100
			min_index = -1
		
			for ix in xrange(0, data_to_fit.shape[0] ):
				val = data_to_fit[ix]
				if val < min:
					min = val
					min_index = ix
			min_index += 1
		
		
			data_to_fit = np.copy( data_to_fit[min_index:] )
			t_fit_data = np.copy ( t_fit_data[0:-min_index] )
			
			popt, pcov = curve_fit( func_complex_exp ,t_fit_data,  data_to_fit, p0= np.array( [data_to_fit[0], .1, 0.0, .5 ] ) )
			curve_fits =  popt[0]*np.cos(popt[3]*t_fit_data)*(np.exp(-popt[1]*t_fit_data)) + popt[2]
			
			t_fit_data__put_back = np.copy( t_fit_data)
			t_fit_data__put_back += ( t_axis[min_index + time_index_offset] )
			plot(t_fit_data__put_back, data_to_fit)
			plot(t_fit_data__put_back, (curve_fits))
			axhline( y=0, color='black')
			#xlim([0.0, 50.0])
			savefig("e_vs_time__central_box")
			
			# try to fit FFT data...
			#func_decaying_cos_squared( x, a, b, c, phase, w)
			
			fft_mode_fit_end_index = int ( 30.0 / output_time_step )
			fft_maxs__as_np_array = np.array( fft_maxs )
			data_to_fit_fft_mode = np.copy ( fft_maxs__as_np_array[min_index:fft_mode_fit_end_index ] )
			#data_to_fit_fft_mode now starts a first min.. for fitting the FFT coefficients, we want the 2nd maximum... so we start from this 1st min and find the next max..
			max_index = np.argmax( data_to_fit_fft_mode )
			data_fit_start_offset = max_index + min_index
			data_to_fit_fft_mode = np.copy ( fft_maxs__as_np_array[data_fit_start_offset: fft_mode_fit_end_index] )
			t_fit_data_fft_mode = np.copy( t_axis[0:data_to_fit_fft_mode.shape[0] ] )
			
			opt_funct = make_cos2_opt( data_to_fit_fft_mode[0],  1.0)		# opt_funct: b, w, c
			popt_cos2, pcov = curve_fit( opt_funct ,t_fit_data_fft_mode,  data_to_fit_fft_mode, p0= np.array( [ popt[1], 0.0] ) ) 		# popt[3]
			
			# create the fitted data to plot...
			
			popt_cos2[1] += 1.0
			temp2 = np.fabs( np.cos((popt_cos2[1])*t_fit_data_fft_mode) )
			curve_fits_fft_mode =   data_to_fit_fft_mode[0]*temp2*(np.exp(-1.0*popt_cos2[0]*t_fit_data_fft_mode)) 
			
			print "\n\n\|COS| FIT!!!!!!"
			print "\t amp: %e \t cos freq: %f \t damping: %e \t (gamma/w): %f  " % (data_to_fit_fft_mode[0],  popt_cos2[1], popt_cos2[0],  (popt_cos2[0]/popt_cos2[1]) )  
			
			t_fit_data_fft_mode += ( t_axis[data_fit_start_offset] )
			
			clf()
			title('Wave k versus time'); xlabel( "Time (units of $w_p$)"); ylabel('magnitude of k mode');
			plot(t_axis[0:len(fft_maxs)], fft_maxs )
			plot( t_fit_data_fft_mode, data_to_fit_fft_mode)
			plot( t_fit_data_fft_mode, curve_fits_fft_mode, color='black')
			savefig("wave_k_vs_time")
		#except:
		#	print "error fitting curves.. prolly because not enugh simulation steps have occured."
		
	if True:#try:
		print "\n\n\n"
		print "Fit Data:"
		print "\t amp: %e \t cos freq: %f \t damping: %e \t (gamma/w): %f \t y-offset: %e" % ( popt[0],  popt[3], popt[1], (popt[1]/popt[3]), popt[2])
		#pd_fit_paramters0 = landau_damping_frequency(wp_e, vth_e, k_of_excited_mode, inital_root_guess=popt[3], maxwellian_convention_factor=1.0)
		#print "PD says (factor=1) cos freq: %f \t damping: %f" % ( np.real(pd_fit_paramters0), np.imag(pd_fit_paramters0))
		pd_fit_paramters_expected = landau_damping_frequency(wp_e, vth_e, k_of_excited_mode) #, inital_root_guess=popt[3])
		print "\nPlasma Dispersion function gives (factor=2):"
		print	"\t                     \t cos freq: %f \t damping: %e \t (gamma/w): %f" % ( np.real(pd_fit_paramters_expected), np.imag(pd_fit_paramters_expected), (np.imag(pd_fit_paramters_expected)/ np.real(pd_fit_paramters_expected)) )
		print""
		
		#print  dispersion[-1, constant_postion_index]
		#Ev_phase =  np.real(pd_fit_paramters0) / k_of_excited_mode
		#print "vpahse %e" % v_phase
		v_phase =  np.real(pd_fit_paramters_expected) / k_of_excited_mode
		print "v_phase: %e \t v_thermal/v_phase: %e \t (v_thermal/v_phase)^2: %e" % ( v_phase, (vth_e/v_phase), ((vth_e/v_phase)*(vth_e/v_phase))  )
		print ""
		
		"""
		fundimental_k =  k_of_excited_mode / float( mode_excited )
		for mode in xrange (int( mode_excited ) - 5 , int( mode_excited ) + 5 ):
			if mode < 1:
				continue
			k_of_this_mode = fundimental_k * float( mode )
			kvth_of_this_mode = k_of_this_mode*vth_e
			if (kvth_of_this_mode > 1.8):
				continue
			
			pd_fit_paramters = landau_damping_frequency(wp_e, vth_e, k_of_this_mode, inital_root_guess=popt[3], maxwellian_convention_factor=2.0)
			print "MODE %d: \t k=%f \t k*vth: %f \t cos freq: %f \t damping: %f" % ( mode, k_of_this_mode,  kvth_of_this_mode, np.real(pd_fit_paramters), np.imag(pd_fit_paramters) )
		"""
		
		niquist_angular_freq = ( 2*np.pi / float(output_time_step) )
		num_time_samples = int ( dispersion.shape[0] ) # t_fit_data.shape[0] #int ( dispersion.shape[0] )
		delta_w = niquist_angular_freq / float( num_time_samples )
		
		print ""
		print "Frequency information:"
		print " \t Output interval (sampling period) %f" % ( output_time_step )
		print " \t Number output samples used        %d" % ( num_time_samples )
		print " \t Niqist angular frequency          %f" % ( niquist_angular_freq )
		print " \t Size of angular freq cells        %f" % ( delta_w )		
		print ""
		print " \t Fit freq/(angular freq cell size) %f" % ( popt[3] / delta_w )
		print " \t PD  freq/(angular freq cell size) %f" % ( np.real(pd_fit_paramters_expected) / delta_w )
		print ""
		print " \t Freq of largest E-field FFT peak  %f" % ( float(fft_maximum_time_freq) * delta_w )
		print " \t Freq of largest E-field FFT peak /(angular freq cell size) %f" % ( float(fft_maximum_time_freq) )
	#except:
	#	pass
#except:
#	print " Error during curve fit."
#	pass

#landau_dmapings_find_closest_k( vth_e, popt[3], -1.0* popt[1] , k_fundimental = 2*np.pi/box_size__x )
#print ""
#print "kvth %f" % k__vth_e



if if_echo_mode:
	data = dispersion[:, constant_postion_index]
	
	#if if_echo_mode_2nd_order and echo_pulse_start_index < data.shape[0]:
	if echo_pulse_start_index < data.shape[0]:
		if echo_pulse_end_index > data.shape[0]:
			echo_pulse_end_index = data.shape[0]
		
		echo_duration = 2.0 / ( k_of_excited_mode*vth_e)
		print "Echo duration theroy %e" % echo_duration
		
		#echo_pulse_time
		#echo_pulse_start_index
		#echo_pulse_end_index
		
		echo_duration = 15.0
		echo_pulse_start_index = int ( ( echo_pulse_time - echo_duration ) / time_step_size / output_step )
		echo_pulse_end_index =  int ( ( echo_pulse_time + echo_duration ) / time_step_size / output_step )
		
		#echo_pulse_start_index -= 1000
		#echo_pulse_end_index += 700
		data = dispersion[:, 17]
		#data = np.array( fft_maxs )
		#data = E_wave_energy
		try:
			np.save( 'shit_17__for_6_10__96cells__1ld', data[ echo_pulse_start_index:echo_pulse_end_index]  )
			clf()
			subplot(211)
			plot( t_axis[ echo_pulse_start_index:echo_pulse_end_index ],  data[ echo_pulse_start_index:echo_pulse_end_index]  )
			title('E-field wave centered around excepted 2nd order echo time.\n'); xlabel( "Time (units of $\\frac{1}{w_p}$)");ylabel( "E field");
			axvline( x = echo_pulse_time )
			subplot(212)
			data = E_wave_energy
			plot( t_axis[ echo_pulse_start_index:echo_pulse_end_index ],  data[ echo_pulse_start_index:echo_pulse_end_index]  )
			title('Integrated Total energy is E around  excpected echo time.')
			xlabel( "Time (units of $\\frac{1}{w_p}$)");ylabel( "E field Energy (frac{$E^2$}{2})");
			axvline( x = echo_pulse_time )
			savefig("!e_vs_time__2nd_order_echo_goodish")
			
			np.save('wave_energy', E_wave_energy)
		except:
			pass
		
		"""
		for i in xrange( dispersion.shape[1] ):
			print i
			clf()
			data = dispersion[:, i]
			plot( t_axis[ echo_pulse_start_index:echo_pulse_end_index ],  data[ echo_pulse_start_index:echo_pulse_end_index]  )
			title('E-field wave centered around excepted 2nd order echo time.\n'); xlabel( "Time (units of $\\frac{1}{w_p}$)");ylabel( "E field"); axvline( x = echo_pulse_time );
			savefig("./png/e_vs_time__central_box__lowest_order_echo_%03d" % i)
		"""
			
	if if_echo_mode_3rd_order and echo_3nd_order_start_index < data.shape[0]:
		if echo_3nd_order_end_index > data.shape[0]:
			echo_3nd_order_end_index = data.shape[0]
		
		clf()
		plot( t_axis[echo_3nd_order_start_index: echo_3nd_order_end_index], data[ echo_3nd_order_start_index:echo_3nd_order_end_index] )
		title('E centered 3rd order echo theo spot.\n'); xlabel( "Time (units of $\\frac{1}{w_p}$)");ylabel( "E field");
		savefig("e_vs_time__central_box__3rd_order_echo")
		
		clf()
		idxp = fft2_time_vs_k_phase.shape[1]/2 + mode_excited 
		idxm = fft2_time_vs_k_phase.shape[1]/2 + 2*second_pulse_mode
		print mode_excited
		print 2*second_pulse_mode
		
		subplot(211)
		plot( t_axis[0:fft2_time_vs_k_phase.shape[0]], fft2_time_vs_k_phase[:, idxp])
		subplot(212)
		plot( t_axis[0:fft2_time_vs_k_phase.shape[0]], fft2_time_vs_k_phase[:, idxm])
		savefig("3rd_order_echo__phase_k1_k2")
pass

# collect results......
results = DoResults()  
results.values['fixed_k_vs_time']			= fft_maxs
results.values['k_vth'] 	= k__vth_e
results.values['e_fit'] 	= {'freq': popt[3], 		'amp':popt[0], 			'gamma':popt[1], 	'y_offset':popt[2] }
try:
	results.values['fft_fit'] 	= {'freq': popt_cos2[1], 	'gamma': popt_cos2[0], 	'amp':data_to_fit_fft_mode[0] }
except:
	pass

results.values['theo_fit']	= {'freq': np.real(pd_fit_paramters_expected), 'gamma':np.imag(pd_fit_paramters_expected) }
pickle.dump( results, open( "run_results.p", "wb" ) )

print "bounce time			%e" % ( bounce_time )

lowest_excited_freq_at_all_positions = []
most_excited_freq_at_all_postions = []

for x in xrange(0, fft2_w_vs_x.shape[1]):
	maxima = find_maxima( fft2_w_vs_x[:, x], start_index = fft2_w_vs_x.shape[0]/2, sort_by_freq = True )
	lowest_excited_freq_at_all_positions.append( freq_axis[ maxima[0][0] ] )
	maxima = find_maxima( fft2_w_vs_x[:, x], start_index = fft2_w_vs_x.shape[0]/2, sort_by_freq = False )
	most_excited_freq_at_all_postions.append( freq_axis[ maxima[0][0] ] )

lowest_excited_freq_at_all_positions 	= np.array( lowest_excited_freq_at_all_positions )
most_excited_freq_at_all_postions 		= np.array( most_excited_freq_at_all_postions )
clf()
plot(x_axis, most_excited_freq_at_all_postions, color='red', label='Highest amplitude excited mode versus x')
plot(x_axis, lowest_excited_freq_at_all_positions, color='green', label='Lowest freq excited mode versus x')
legend()
savefig('!shit_bim')

print ""
print "------------"
print "List of all unique values of the lowest frequency excited temporal mode, across all spatial postions"
temp = np.unique(lowest_excited_freq_at_all_positions.round(decimals=4))
for i in xrange(0, temp.shape[0]):
	print "Frequency                     %f" % temp[i]

print "------------"
print "List of all unique values of the frequency of the highest amplitude temoral mode, across all spatial postions"
temp = np.unique(most_excited_freq_at_all_postions.round(decimals=4))
for i in xrange(0, temp.shape[0]):
	print "Frequency                     %f" % temp[i]
print ""

exit(-1)



real_data_start = int (driver_stop_time / output_time_step) + 1
data = np.array( maxs )

time_axis = np.linspace(0, (output_time_step*float(len( maxs))), len( maxs))
clf()

#print time_axis
#print pure_harmonic_dt
#exit(-1)


"""
driver_data = plot_driver() 
driver_axis = np.arange(0.0, float( len (  driver_data) ))
driver_axis *= time_step_size

popt, pcov = curve_fit(func, driver_axis, driver_data)
x_fit_data = np.copy ( driver_axis	)
curve_fit =  popt[0]*np.cos(popt[4]*x_fit_data + popt[3])*(np.exp(-popt[1]*x_fit_data)) + popt[2]


print "------------------------ DRIVER ================"
print popt
plot(driver_axis, driver_data)
savefig('drivem.png')
driver_dt =  numpy.fft.fftshift( numpy.fft.fftfreq( fft2_dispersion.shape[0], d = (.098/ (2*np.pi*wp_e)) ) ) 
print driver_dt
print fft2_dispersion.shape[0],
exit(-1)
"""

"""
clf()
plot( time_axis, np.array ( snr ) )
savefig("e_vs_time__snr")
"""

clf()
ax = subplot(111)
ylabel('Electric Field')
xlabel('Time (units of $w_p$)')

est_damp = estimated_damping(wp_e, vth__e, wavenumber)
max_value = 0.0

#try:
	
look = np.array ( data[fit_start_index:fit_end_index] )
look_fft = np.array ( fft_maxs[fit_start_index:fit_end_index] )
max_index = -1


for i in xrange(fit_offset_index, look.shape[0] ):
	
	if np.abs ( look[i] ) > np.abs (max_value):
		max_value = look[i]
		max_index = i



if max_index > -1:
	look /= max_value
	data /= max_value
	
	data_2_fit = look[max_index:]
	data_2_fit_fft = look_fft[max_index:]

	# try non-linear exponential fit....
	#x_fit_data = time_axis[0:(fit_end_index-fit_start_index - max_index)]
	fit_fun = None
	if fit_offset_index > 0.0:
		fit_fun = func_t_off( time_axis[max_index])
		#x_fit_data = time_axis[max_index:max_index+data_2_fit.shape[0] ]
		x_fit_data = time_axis[:(fit_end_index-fit_start_index - max_index)]
	else:
		fit_fun = func
		x_fit_data = time_axis[:(fit_end_index-fit_start_index - max_index)]
	
	#x_fit_data /= 2.0
	#x_fit_data /= np.pi()
	print x_fit_data
	#exit(-2)
	
	"""
	# pure exponential
	popt, pcov = curve_fit(func_exp, x_fit_data, data_2_fit)
	curve_fit = popt[0]*np.exp( -1.0*popt[1]*x_fit_data ) + popt[2]
	"""
	
	
	#exit(-2)
	
	# complex exponential
	popt, pcov = curve_fit(fit_fun, x_fit_data, data_2_fit)	
	curve_fits =  popt[0]*np.cos(popt[4]*x_fit_data + popt[3])*(np.exp(-popt[1]*x_fit_data)) + popt[2]
	
	print x_fit_data
	#popt, pcov = curve_fit(func, x_fit_data, data_2_fit_fft)	
	#curve_fits_fft =  popt[0]*np.cos(popt[4]*x_fit_data + popt[3])*(np.exp(-popt[1]*x_fit_data)) + popt[2]
	
	
	clf()
	plot( time_axis[0:len(curve_fits)], curve_fits)
	plot( time_axis[0:len(curve_fits)], data_2_fit)
	savefig("e_vs_time__fit_plot")
	clf()
	
	damping_rate_adjusted = ( popt[1]*np.sqrt(2.0*np.pi)*1.5)
	damping_rate_bohm_gross = np.sqrt( wp_e*wp_e + 3*wavenumber*wavenumber*vth__e*vth__e)
	new_damping = landau_damping_frequency(wp_e, vth__e, wavenumber, inital_root_guess=popt[4])
	
	print "wave number"
	print wavenumber
	print "vth"
	print vth__e
	print "k*vth"
	print vth__e*wavenumber
	print "vphase"
	print (w_wave / wavenumber)
	print "new vphase"
	print (popt[4] / wavenumber)
	
	print new_damping
	
	print ""
	print "\tClosest Root of kinetic, warm, elctrostatic plasma dielectric (k*vth=%e):" % k__vth
	print "\t                                    " + str(landau_damping_frequnecy__complex)
	print ""
	print "====================================================="
	print "\tfrequency-phase shift (from fit):\t%e" % popt[3]
	print "\t            amplitude (from fit):\t%e" % popt[0]
	print "\t             y-offset (from fit):\t%e" % popt[2]
	print ""
	print "\t            frequency (from fit):\t%e" % popt[4]
	print "%e" % ( popt[4]/ np.real ( landau_damping_frequnecy__complex ))
	print "%e" % (np.real ( landau_damping_frequnecy__complex ) /  popt[4])
	print "%e" % ( 2.0*np.pi*popt[4])
	print "\t                      (theory  ):\t%e" % np.real ( landau_damping_frequnecy__complex )
	print "\t                      (Bohm-Gro):\t%e" % np.real ( damping_rate_bohm_gross )
	print ""
	print "\t         damping rate (from fit):\t%e" % popt[1]
	print "\t                      (adjusted):\t%e" % damping_rate_adjusted
	print "\t                      ( Landau ):\t%e" % np.imag(landau_damping_frequnecy__complex)
	print "\t                    ( etimated ):\t%e" % est_damp
	
	print "======================================================"
	print ""
	
#except:
#	pass

title('Electric Field (at fixed spatial location) versus Time\n $kv_{th}=%f$ $E_{max}=%e$' % ( sim_metadata.driver['k__vth'], max_value ) )
plot(time_axis,  (data) )
#xlim([10.0,20.0])

np.save('e_vs_time', data)
print "fdgdfgdf"
exit(-1)

try:
	plot(time_axis[(fit_start_index+max_index):fit_end_index], curve_fits )
	#plot( time_axis[0:len(curve_fits)], curve_fits)
	
except:
	pass

axvline(x=driver_stop_time, color='black', ls='--')
#loc = ax.yaxis.get_ticklocs()[5]
#text(4.0,loc-10,'Driver Shut Off')
#ax.annotate('Driver Shut Off', xy=(4.0, -1.0),  xycoords='data',
#                xytext=(0.8, 0.05), textcoords='axes fraction',
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                horizontalalignment='right', verticalalignment='bottom',
#                )
#xlim([300,700])
#ylim( [-.1, .1])
savefig('e_vs_time')

clf()
ax = subplot(111)
title('Electric Field (at fixed spatial location) versus Time\n $kv_{th}=%f$' % sim_metadata.driver['k__vth'])
ylabel('Electric Field')
xlabel('Time (units of $w_p$)')

				
plot(time_axis,  data )
plot(time_axis,  np.fabs(data) )
axvline(x=driver_stop_time, color='black', ls='--')
"""
loc = ax.yaxis.get_ticklocs()[5]
text(4.0,loc-10,'Driver Shut Off')
ax.annotate('Driver Shut Off', xy=(4.0, -1.0),  xycoords='data',
                xytext=(0.8, 0.05), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='bottom',
                )
"""
savefig('e_vs_time_abs')


clf()
ax = subplot(111)
title('ass DFT of Electric Field (at fixed spatial location) along  the Time axis. \n $kv_{th}=%f$' % sim_metadata.driver['k__vth'])
ylabel('|Electric Field fourier component|')
xlabel('Frequency (units of $w_p$)')


freq_axis = numpy.fft.fftshift( numpy.fft.fftfreq( data.shape[0], d = (output_time_step/ (2*np.pi*wp_e)) ) ) #d = (1.0/(2.0*f_nyquist) ) ) )

fftfata = np.log10( np.absolute( numpy.fft.fftshift( numpy.fft.fft( data ) ) ) )
plot(freq_axis, fftfata )
#xlim([-1.5, 1.5])
#ax.set_xticks(numpy.arange(-1.5,1.5,0.1))
#for tick in ax.xaxis.get_major_ticks():
#	tick.label1.set_fontsize(6.0)
#
#grid()

np.save('e_vs_frquency', fftfata)
np.save('e_vs_frquency_freq_axis', freq_axis)
savefig('e_vs_time__fft')

half = int(data.shape[0] / 2)

"""
try:
	print " values of k vs fft of time:"
	for i in xrange(0, freq_axis.shape[0]/2):
		idx = half+2*i
		idx1 = idx + 1
		print "freq_bin %02d: \t freq: %f \t value: %f" % (2*i, freq_axis[idx], fftfata[idx]),
		print "\t\t freq_bin %02d: \t freq: %f \t value: %f" % (2*i + 1, freq_axis[idx1], fftfata[idx1])
	print ""
	print ""
	print ""
except:
	pass
"""
print freq_axis[0]
print freq_axis[-1]

clf()
ax = subplot(111)
title('Magnitude of Electric field''s spatial fouierr compoent (k of driver) versus Time \n $kv_{th}=%f$' % sim_metadata.driver['k__vth'])
ylabel('|Electric Field fourier component|')
xlabel('Time (units of $w_p$)')

						
plot(time_axis,  fft_maxs, marker="*", ms=.8 )
try:
	#plot(x_fit_data, curve_fits_fft)
	#print x_fit_data
	#exit(-1)
	pass
except:
	pass
axvline(x=driver_stop_time, color='black', ls='--')
"""
loc = ax.yaxis.get_ticklocs()[5]
text(4.0,loc-10,'Driver Shut Off')
ax.annotate('Driver Shut Off', xy=(4.0, -1.0),  xycoords='data',
                xytext=(0.8, 0.05), textcoords='axes fraction',
                arrowprops=dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='bottom',
                )
"""
savefig('e_vs_time__fft_mode')

try:
	print "dt was: %f (pure harmonic is %f) " % (time_step_size, pure_harmonic_dt)
	print "(output_time_step): each out timestep is: %f " % output_time_step
	print "w_wave was: %f or %f " % (w_wave, 2*np.pi*w_wave)
	print (popt[1] *2.0 *np.pi )
	#print popt[1] *2.0 *np.pi / .243
	#print .243 / popt[1] *2.0 *np.pi 
	#print time_axis
	print "Theoretical frequnecy from dipersion releation:"
	print landau_damping_frequnecy__complex	
	print "Ratio theoretical damping/ measured_damping % e" % ( -1.0*np.imag(landau_damping_frequnecy__complex) / popt[1] )
	print "Ratio theoretical damping/ measured_damping % e" % ( -1.0*popt[1]/np.imag(landau_damping_frequnecy__complex) )
	print "measured_damping*3/2*sqrt(2*pi) % e" % ( popt[1]*np.sqrt(2.0*np.pi)*1.5)
	print "wavenumber %e" % wavenumber
	print "wavelength %e" % (2.0*np.pi/wavenumber)
except:
	pass
print dispersion.shape

#k_0_org = 
k_0 =  k_axis[k_axis.shape[0] / 2 +1 ]
k_30 = k_axis[k_axis.shape[0] / 2 + 30 ]
k_45 = k_axis[k_axis.shape[0] / 2 + 45 ]
print k_45/ k_30
print k_45 - k_30
print  k_axis[k_axis.shape[0] / 2 + 15 ]
print k_axis[k_axis.shape[0] / 2 ]

num_k = k_axis.shape[0]
num_k2 = k_axis.shape[0] / 2


#for i in xrange( 0, num_k2):
#	print "Slot %d: k*vth: %f" % ( i, k_axis[num_k2+i]*vth__e)