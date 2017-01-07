import matplotlib
matplotlib.use('Agg')

import numpy as np
import numpy.fft
import oshun1d_utils_plasma_dispersion as plasma_dispersion

from pylab import *
import scipy.integrate
import scipy.stats as stats
from scipy.optimize import curve_fit
import cPickle as pickle
import glob as glob
import os
from vqunits import Formulary 
from vqscript import specrogram_k_vs_x

from Oshun1d import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile


class DoResults:
	def __init__(self):
		self.values = {}

#
# TODO: take this out....
# 
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

def func_exp(x, a, b, c):
	return a*(np.exp(-b*x)) + c

def load_packed_E_data(output_time_step = 1.0, limit_data_to_before_time = 0.0, wp_e=1.0, cell_size__x=1.0, box_size__x = 1.0):
	files = glob.glob( './numpy/E_*.npy')
	files.sort()
	
	# determine how many E steps are in a given E-file and the size of each E field entry.
	filename = os.path.split(files[0])[1]
	f1 = np.load( "./%s/numpy/%s" % ( '.', filename) )
	pack_factor = 1
	data_entry_size = f1.shape[0]
	
	if len(f1.shape) == 2:
		pack_factor = f1.shape[0]
		data_entry_size = f1.shape[1]
	
	number_of_timestep_datas_availiable = len(files)*pack_factor
	sim_stop_time = float( number_of_timestep_datas_availiable ) * output_time_step# * float( pack_factor)
	sim_num_time_steps = number_of_timestep_datas_availiable
	
	print sim_stop_time
	# Now we can modify how many E-steps we actually want to load in..
	
	# only load in data up to certain time.....
	if limit_data_to_before_time > 0.0:
		if limit_data_to_before_time < sim_stop_time:
			sim_stop_time = limit_data_to_before_time
			sim_num_time_steps = int( sim_stop_time / output_time_step ) - 1
			if sim_num_time_steps % 2 == 1:
				sim_num_time_steps -= 1
	
	
	# the data array!!!!!
	e_raw_data = np.zeros( (sim_num_time_steps, data_entry_size) )
	
	
	# load in data files.....
	current_index = 0
	for file_entry in files:
		filename = os.path.split(file_entry)[1]
		f1 = np.load( "./%s/numpy/%s" % ( '.', filename) )
		
		for ir in xrange(0, pack_factor):
			if current_index < e_raw_data.shape[0]:
				if is_data_in_array_valid( f1[ir,:] ):
					e_raw_data[current_index,0:data_entry_size] = f1[ir,:]
				else:
					print "Invalid data found... stopping reading input data..."
					sim_num_time_steps = index
					if sim_num_time_steps % 2 == 1:
						sim_num_time_steps -= 1
					sim_stop_time = float( sim_num_time_steps ) * output_time_step
					temp = np.zeros ( (sim_num_time_steps, data_entry_size) )
					temp[0:sim_num_time_steps, :] = e_raw_data[0:sim_num_time_steps,:]
					e_raw_data = None
					e_raw_data = temp
			current_index += 1
		pass
	
	print ""
	print "Done reading files"
	print ""
	np.save( 'e_vs_time__efield_t_vs_x', e_raw_data)
	
	# create axes that give real values (using Oshun's usual normalizations )
	freq_axis = numpy.fft.fftshift( numpy.fft.fftfreq( e_raw_data.shape[0], d = (output_time_step/ (2*np.pi*wp_e)) ) ) 
	k_axis = numpy.fft.fftshift( numpy.fft.fftfreq( e_raw_data.shape[1],  d = (cell_size__x/ (2*np.pi*wp_e)) ) )	
	t_axis = np.arange(0, float(sim_num_time_steps)*float(output_time_step), float(output_time_step))
	x_axis = np.arange(0, 1.0, 1.0/e_raw_data.shape[1] )
	x_axis *= box_size__x
	
	return (e_raw_data, freq_axis, k_axis, t_axis, x_axis)
	
def is_data_in_array_valid(data):
	test__for_invalids_values = np.isfinite( data )
	np.logical_not( test__for_invalids_values , out= test__for_invalids_values )
	if np.any( test__for_invalids_values):
		return False
	return True


def find_maxima( data, start_index=0, end_index = None, sort_by_location=False, minima = False):
	if end_index == None:
		end_index = data.shape[0]
	
	if minima:
		fred = (diff(sign(diff(data[start_index: end_index]))) > 0).nonzero()[0] + 1
	else:
		fred = (diff(sign(diff(data[start_index: end_index]))) < 0).nonzero()[0] + 1
	
	maxima = [] 
	for i in xrange(0, fred.shape[0] ):
		maxima.append( ( start_index+fred[i],  data[start_index+fred[i] ] ) )
	
	if sort_by_location:
		maxima.sort(key=lambda x: x[0], reverse=False)
	else:
		# sort by value
		maxima.sort(key=lambda x: x[1], reverse=True)
	
	return maxima
	

# global confugartion params that need to be moved...
limit_data_to_before_time = 25.0


#
# This gets run over a run (directory)...
#
def the_main():
	#sim_data["oshun_metadata__cfg"]
	# TODO: this needs to be formalized...
	# load information from the simulation...
	sim_metadata = pickle.load( open( "./%s/oshun_metadata.p" % ('.',), "rb" ) )
	sim_metadata_cfg = pickle.load( open( "./%s/oshun_metadata__cfg.p" % ('.',), "rb" ) )
	
	output_step 				= sim_metadata.output_interval__time_steps
	time_step_size 		 	= sim_metadata.dt   	
	output_time_step			= float(output_step) * time_step_size
	
	box_size__x 				= sim_metadata_cfg.system_size_max__global[0]
	num_cells__x 				= sim_metadata_cfg.num_cells__global[0]
	cell_size__x 				= box_size__x / float( num_cells__x )
	
	# also read in the plasma paramters...
	density						= sim_metadata_cfg.species_config[0].density_profile.val	
	vth_e							= sim_metadata_cfg.species_config[0].temp_profile.val
	wp_e							= np.sqrt( density )
	
	# read in wave excitation parameterers....
	mode_excited				= sim_metadata_cfg.E_field_profile.num_wavelengths
	
	# and related quantitites..
	k_of_excited_mode 		= 2.0*np.pi / ( box_size__x / float( mode_excited)  )
	k__vth_e						= k_of_excited_mode*vth_e
	
	print "K was caled as %e and kvt=%e" % ( k_of_excited_mode, k__vth_e)
	other_k = np.pi*2.0 / ( box_size__x + cell_size__x )
	other_k__vth_e = other_k * vth_e
	print "other K was caled as %e and kvt=%e" % ( other_k, other_k__vth_e)
	k_of_excited_mode = other_k; k__vth_e = other_k__vth_e;
	
	# and also calculate the expected complex frequeny from the plasma dipsersion function
	pd_fit_paramters_expected 	= plasma_dispersion.landau_damping_frequency(wp_e, vth_e, k_of_excited_mode)
	pd_real_freq					= np.real( pd_fit_paramters_expected )
	pd_imag_freq					= np.imag( pd_fit_paramters_expected )
	
	v_phase_pd = pd_real_freq/ k_of_excited_mode
	print v_phase_pd
	

	
	pp_string = "Fred"
	form=Formulary() 
	vth_in_ev = form.convert_osi_vth_to_ev( vth_e )
	try:
		if sim_metadata_cfg.collsions_enabled:
			densitity__1cm3 = sim_metadata_cfg.species_config[0].density_np
			
			pp = form.plama_paramter( vth_in_ev, densitity__1cm3 )
			print "Teperature in EV: %e" % vth_in_ev
			print "Plasma parameter %e (%0.0f ev and %0.1e $cm^{-3}$" % ( pp, vth_in_ev, densitity__1cm3)
			pp_string = "Plasma Parameter %0.1e (%0.0f ev and %0.1e $cm^{-3}$)" % ( pp, vth_in_ev, densitity__1cm3)
	except:
		pp_string = "Collisionless, %0.0f ev, $kv_{th}=%.2f$)" % (vth_in_ev, k__vth_e)
	
	# use this damping rate to calculate how much time is needed to damp to 1e^2
	#time_needed_to_damp_theoretical = np.fabs( 2.0 / pd_imag_freq )
	#limit_data_to_before_time = time_needed_to_damp_theoretical * 1.1
	limit_data_to_before_time = 1500.0
	if limit_data_to_before_time < 25.0:
		limit_data_to_before_time = 25.0
	
	
	print "\n\n"
	print "Reading until simulation time %f" % limit_data_to_before_time
	print "\n\n"
	
	# list all the E files for this run....
	files = glob.glob( './numpy/E_*.npy')
	files.sort()
	
	# determine how many E steps are in a given E-file and the size of each E field entry.
	filename = os.path.split(files[0])[1]
	f1 = np.load( "./%s/numpy/%s" % ( '.', filename) )
	pack_factor = 1
	data_entry_size = f1.shape[0]
	
	if len(f1.shape) == 2:
		pack_factor = f1.shape[0]
		data_entry_size = f1.shape[1]
	
	number_of_timestep_datas_availiable = len(files)*pack_factor
	sim_stop_time = float( number_of_timestep_datas_availiable ) * output_time_step# * float( pack_factor)
	sim_num_time_steps = number_of_timestep_datas_availiable
	
	print sim_stop_time
	# Now we can modify how many E-steps we actually want to load in..
	
	# only load in data up to certain time.....
	if limit_data_to_before_time > 0.0:
		if limit_data_to_before_time < sim_stop_time:
			sim_stop_time = limit_data_to_before_time
			sim_num_time_steps = int( sim_stop_time / output_time_step ) - 1
			if sim_num_time_steps % 2 == 1:
				sim_num_time_steps -= 1
	
	
	# the data array!!!!!
	e_raw_data = np.zeros( (sim_num_time_steps, data_entry_size) )
	
	
	# load in data files.....
	current_index = 0
	for file_entry in files:
		filename = os.path.split(file_entry)[1]
		f1 = np.load( "./%s/numpy/%s" % ( '.', filename) )
		
		for ir in xrange(0, pack_factor):
			if current_index < e_raw_data.shape[0]:
				if is_data_in_array_valid( f1[ir,:] ):
					e_raw_data[current_index,0:data_entry_size] = f1[ir,:]
				else:
					print "Invalid data found... stopping reading input data..."
					sim_num_time_steps = index
					if sim_num_time_steps % 2 == 1:
						sim_num_time_steps -= 1
					sim_stop_time = float( sim_num_time_steps ) * output_time_step
					temp = np.zeros ( (sim_num_time_steps, data_entry_size) )
					temp[0:sim_num_time_steps, :] = e_raw_data[0:sim_num_time_steps,:]
					e_raw_data = None
					e_raw_data = temp
			current_index += 1
		pass
	
	print ""
	print "Done reading files "
	np.save( 'e_vs_time__efield_t_vs_x', e_raw_data)
	
	# create axes that give real values (using Oshun's usual normalizations )
	freq_axis = numpy.fft.fftshift( numpy.fft.fftfreq( e_raw_data.shape[0], d = (output_time_step/ (2*np.pi*wp_e)) ) ) 
	k_axis = numpy.fft.fftshift( numpy.fft.fftfreq( e_raw_data.shape[1],  d = (cell_size__x/ (2*np.pi*wp_e)) ) )	
	t_axis = np.arange(0, float(sim_num_time_steps)*float(output_time_step), float(output_time_step))
	x_axis = np.arange(0, 1.0, 1.0/e_raw_data.shape[1] )
	x_axis *= box_size__x
	
	
	E_squared = np.multiply( e_raw_data, e_raw_data)
	E_squared /= 2.0
	E_wave_energy = scipy.integrate.simps( E_squared, dx = cell_size__x , axis=1)
	
	#clf()
	#plot( E_wave_energy)
	#ylim( [-30, -20] )
	#show()
	#exit(-1)
	
	# set any offset from the start of th data to begin the fit..
	data_start_time 		= 0.0
	data_end_time 			= 20.0
	data_start_index 		= int(data_start_time / time_step_size/output_step)
	data_end_index 		= int(data_end_time / time_step_size/output_step)
	
	maxima = find_maxima( E_wave_energy, start_index = 0, sort_by_location = True )
	minima = find_maxima( E_wave_energy, start_index = 0, minima=True, sort_by_location = True )
	
	default_steps_between_maxima = abs ( maxima[1][0] - maxima[0][0] )
	
	temp = np.zeros( len(maxima) )
	for i in range(0, len(maxima) ):
		temp[i] = maxima[i][1]
	maxima_of_maxima = find_maxima ( temp, start_index = 0, sort_by_location = True)
	
	
	if len(maxima) == 0:
		print "No maxima found!!! exiting."
		return
	if len(maxima) == 1:
		print "Only 1 maxima found!! can't fit to one point!!"
		return
	
	mylog = np.log
	if False:
		# take the first maxima value..
		# loop trough and see if se have enough maxima points such that we get to at least 1/e of this value.
		first_maxima_value = maxima[0][1]
		maxima_stop_threshold = np.exp(-2)*first_maxima_value
		index_of_last_maxima_to_use = 0
		for ix in xrange(1, len( maxima ) ):
			if maxima[ix][1] < maxima_stop_threshold:
				index_of_last_maxima_to_use = ix
				break
			print  (first_maxima_value / maxima[ix][1] ) / np.exp(2.0)
		if index_of_last_maxima_to_use < 2:
			print "Could not do a fit with Wave energy since there were not enough maxima above the 1/e threshold"
			return
	
		# now create the point to fit.... we will set the zero point at the first maxima...
		t_equals_zero_index =  maxima[0][0]
		number_of_fit_points = index_of_last_maxima_to_use - 1
		energy_fit_points__time = np.zeros( index_of_last_maxima_to_use )
		energy_fit_points__value = np.zeros( index_of_last_maxima_to_use )
		for ix in range( 0, index_of_last_maxima_to_use):
			energy_fit_points__time[ix] = t_axis[ maxima[ix][0] - t_equals_zero_index]
			energy_fit_points__value[ix] = maxima[ix][1]
		popt, pcov = curve_fit( func_exp ,energy_fit_points__time,  energy_fit_points__value, p0= np.array( [energy_fit_points__value[0], .1, 0.0] ) )
		curve_fits = popt[0]*(np.exp(-popt[1]*energy_fit_points__time)) + popt[2]
	
	
	final__a_s = -1.0
	final__b_s = -1.0
	energy_fit_points__time = None
	energy_fit_points__value = None
	if True:
		t_equals_zero_index =  maxima[0][0]
		current_diff = -1.0
		for ix in xrange(2, len( maxima ) ):
			# construct the fit data for this...
			energy_fit_points__time_temp = np.zeros( ix+1 )
			energy_fit_points__value_temp = np.zeros( ix+1 )
			for ixx in range( 0, ix+1):
				energy_fit_points__time_temp[ixx] 	= t_axis[ maxima[ixx][0] - t_equals_zero_index]
				energy_fit_points__value_temp[ixx]	= mylog( maxima[ixx][1] )
			
			(a_s,b_s,r,tt,stderr)=stats.linregress(energy_fit_points__time_temp,energy_fit_points__value_temp)
			
			print "\t r=%e using first %d maxima: %e with damping %e" % ( r, ix+1, r, a_s/2.0) ,
			print r,
			print maxima[0][0]
			
			
			"""
			# quick hack... make it so that fit occurs out to 100 always for now...
			if current_diff == -1.0:
				print ""
				current_diff = r
				final__a_s = a_s
				final__b_s = b_s
				energy_fit_points__time = energy_fit_points__time_temp
				energy_fit_points__value = energy_fit_points__value_temp
				final_r = r
			else:
				if t_axis[ maxima[ixx][0] - t_equals_zero_index ] > 100.0:
					break
				else:
					final__a_s = a_s
					final__b_s = b_s
					energy_fit_points__time = energy_fit_points__time_temp
					energy_fit_points__value = energy_fit_points__value_temp
					final_r = r
					print ""
			"""
			
			
			if current_diff == -1.0:
				current_diff = np.abs(r)
				final__a_s = a_s
				final__b_s = b_s
				energy_fit_points__time = energy_fit_points__time_temp
				energy_fit_points__value = energy_fit_points__value_temp
				final_r = r
			else:
				if  abs(r) < current_diff:
					print "... abs( r ) decreased.. so stopping search and using %d maxima" % ix
					break
				else:
					final__a_s = a_s
					final__b_s = b_s
					energy_fit_points__time = energy_fit_points__time_temp
					energy_fit_points__value = energy_fit_points__value_temp
					final_r = r
					current_diff = np.abs(r)
					print ""
			
		curve_fits = np.exp(b_s)*(np.exp(a_s*energy_fit_points__time))
		damping_estimate = final__a_s / 2.0;
	
	print "\n\n"
	t_fit_data__put_back = np.copy( energy_fit_points__time )
	t_fit_data__put_back += t_axis[ t_equals_zero_index + data_start_index]
	
	# temp.. only plot till 300
	#index_300 = int(300.0 / time_step_size/output_step)
	
	clf()
	#title("Log Wave Energy versus Time\n" + pp_string)
	xlabel( "Time ( $w_p^{-1}$ )"); ylabel("Log Wave Energy in E Field");
	plot(t_axis[0:E_wave_energy.shape[0]], np.log10( E_wave_energy), color='black', alpha=0.6 ) 
	#plot(t_axis[0:index_300], np.log10( E_wave_energy[0:index_300]), color='black', alpha=0.6, lw=4.0 ) 
	plot(t_fit_data__put_back, np.log10( curve_fits), label="Fit damping %0.2e" % damping_estimate, lw=4.0) 
	savefig('!!!super_cool_wave_energy_vs_time__autoscale')
	
	xlim( ( 0.0, 30.0) )
	ylim( [-11, -4] )
	savefig('!!!super_cool_wave_energy_vs_time__fit')
	ylim( [-9, -4] )
	#legend(prop={'size':10})
	#axvline( x=t_axis[maxima[0][0]] )
	savefig('!!!super_cool_wave_energy_vs_time__fit_zoomed')
	
	print "E wave energy size"
	print E_wave_energy.shape
	print ""
	
	# this is for a quick E field energy spectrogram.. y_limit=10.0,
	
	window_size = abs( 5*default_steps_between_maxima)
	overlap = 2*default_steps_between_maxima
	print "Spectrogram maxima spacing %d" % default_steps_between_maxima
	print "Spectrogram window spacing %d (overlap is %d) " % ( window_size, overlap)
	(bins, freqs, pxx, im, ax1, fig) = specrogram_k_vs_x( E_wave_energy, float(output_time_step), "!!e_energy_spectrogram.png",  args1={ 'NFFT':window_size, 'noverlap':overlap}, title="E Field Energy Specrogram", xlabel="frequency", ylabel="Time" )
	clf()
	imshow( np.log10( pxx), extent = [ freqs[0], freqs[-1], bins[0], bins[-1] ], aspect='auto', origin='bottom', cmap='gray')
	colorbar()
	xlabel('Frequency'); ylabel('Time')
	savefig("!!e_energy_spectrogram.png", dpi=200)
	#show()
	
	#
	# This is for the non-linear plot....
	#
	#E_field_for_bounce_time =  np.sqrt( maxima_of_maxima[0][1] ) * 0.9
	#bounce_time = 2.0*np.pi / np.sqrt(E_field_for_bounce_time*k_of_excited_mode)
	#time_of_first_maxima_of_maxima = t_axis[ maxima[maxima_of_maxima[0][0]][0] ]
	bounce_time = 1.0
	
	#t0 = [(time_of_first_maxima_of_maxima-bounce_time), time_of_first_maxima_of_maxima]
	#t1 = [np.log10(maxima_of_maxima[0][1]), np.log10(maxima_of_maxima[0][1])]
	
	clf()
	title("Log Wave Energy versus Time\n" + pp_string)
	xlabel( "Time ( $w_p^{-1}$ )"); ylabel("Log Wave Energy in E Field");
	figtext( .55, .85, "nominal bounce period %0.1f $w_p^{-1}$" % ( bounce_time ) )
	plot(t_axis[0:E_wave_energy.shape[0]], np.log10( E_wave_energy) )
	
	#plot( np.array(t0) , np.array(t1), color='purple')	
	#print "dEB"
	#print time_of_first_maxima_of_maxima
	#print E_field_for_bounce_time
	#print maxima_of_maxima
	#print k_of_excited_mode
	#axvline( x=time_of_first_maxima_of_maxima, color='purple')
	
	xlim( ( 0.0, 350.0) ) #350.0) )
	ylim(  [-20, -4] )
	#xlim( ( 0.0, 200.0) )
	#ylim(  [-8, -2] )

	savefig('!!!super_cool_wave_energy_vs_time')
	
	xlim( ( 0.0, 500.0) )
	savefig('!!!super_cool_wave_energy_vs_time_to_0500')
	xlim( ( 0.0, 1000.0) )
	savefig('!!!super_cool_wave_energy_vs_time_to_1000')
	xlim( ( 0.0, 1500.0) )
	savefig('!!!super_cool_wave_energy_vs_time_to_1500')
	# linear autoscale....
	clf()
	title("Log Wave Energy versus Time\n" + pp_string)
	xlabel( "Time ( $w_p^{-1}$ )"); ylabel("Log Wave Energy in E Field");
	figtext( .55, .85, "nominal bounce period %0.1f $w_p^{-1}$" % ( bounce_time ) )
	plot(t_axis[0:E_wave_energy.shape[0]],E_wave_energy )
	savefig('!!!super_cool_wave_energy_vs_time__linear_autoscale')
	
	# Also examine the FFT of the wave energy to find the dominate excited mode...
	fft_of_energy = np.log( numpy.fft.fftshift( np.absolute( numpy.fft.fft(E_wave_energy) ) ) )
	maxima = find_maxima( fft_of_energy, start_index = fft_of_energy.shape[0]/2, sort_by_location = False )
	freq_axis_sq = np.copy( freq_axis )
	freq_axis_sq /= 2.0
	
	dominate_freq_fft = -1.0
	if len(maxima) > 1:
		dominate_freq_fft = freq_axis_sq[ maxima[0][0] ]
	
	# a time domain analysis of the dominate freqsuency....
	freqs_from_time_domain_temp = []
	temp_prev = t_axis[ minima[0][0] ]
	for ix in xrange(0, energy_fit_points__time.shape[0]):
		try:
			temp =  t_axis[ minima[ix+1][0] ]
			freqs_from_time_domain_temp.append( np.pi/( temp - temp_prev ) )
			temp_prev = temp
		except:
			break
	freqs_from_time_domain = np.array( freqs_from_time_domain_temp )
	dominate_freq_time_domain = np.mean( freqs_from_time_domain )
	
	
	final_dominate_frequency = 0.0
	if len(maxima) < 2:
		print "No dominate frequency modes found in FFT.. using time domain."
		final_dominate_frequency = dominate_freq_time_domain
	else:
		final_dominate_frequency = freq_axis_sq[ maxima[0][0] ]
		
	
	frac_diff_real_from_time_domain = np.fabs( dominate_freq_time_domain - pd_real_freq)/pd_real_freq
	frac_diff_imag = np.fabs( damping_estimate - pd_imag_freq)/pd_imag_freq
	#frac_diff_imag = np.fabs( -1.0*popt[1]/2.0 - pd_imag_freq)/pd_imag_freq
	
	if dominate_freq_fft > 0.0:
		frac_diff_real_from_fft = np.fabs( dominate_freq_fft - pd_real_freq)/pd_real_freq
		print "Fit frequency (from fft)         is %f		(landau: %f)	(frac diff %f)" % (dominate_freq_fft, pd_real_freq, frac_diff_real_from_fft)
	
	print "For k*vth:                          %f" % k__vth_e
	print "---------------------------------------------------\n\n"
	print "Fit frequency (from time domain) is %f		(landau: %f)	(frac diff %f)" % (dominate_freq_time_domain, pd_real_freq, frac_diff_real_from_time_domain)
	print "Frequency selected for results   is %f" % final_dominate_frequency
	print "Fit Damping  Rate                is %e		(landau: %e)	(frac diff %f)" % (damping_estimate, pd_imag_freq, frac_diff_imag)
	
	print ""
	print "Wave energy maxima found at times: "
	print t_fit_data__put_back
	
	results 											= DoResults()
	results.values['theo_fit']					= {'freq': np.real(pd_real_freq), 'gamma':pd_imag_freq }
	results.values['data_fit']					= {'freq': np.real(final_dominate_frequency), 'gamma': damping_estimate }
	results.values['kvth']						= k__vth_e
	results.values['k_of_excited_mode']		= k_of_excited_mode
	results.values['vth_e']						= vth_e
	pickle.dump( results, open( "run_results.p", "wb" ) )

if __name__ == "__main__":
	the_main()
	print "\n\n done. \n\n"

	