import oshun_modules
import oshun_modules.OshunInputDeckHelper as Helper
import oshun_modules.setup.SetupTemplates as templates
import oshun_modules.Oshun1d as Oshun1d

import numpy as np
from  oshun_modules.Oshun1d import ConstantProfile, SinusodialProfile, EFieldSinusodialProfile, run_sim
from oshun_modules.InputDeckConfiguration import Config, SpeciesConfig

USE_C_VERSION_GLOBAL = True

from oshun_modules.profiles.SpatialProfileShapes import profiles
from oshun_modules.profiles.ArbitraryDensity import ArbitraryDensityProfile
from oshun_modules.profiles.ArbitraryTemperature import ArbitraryTemperatureProfile

CFG = None
def ld_setup_3rd():
	global CFG

	# OVERVIEW
	# The name of the game is to make a Config object and populate with simulation parameters.
	# then pass this Config object to Oshun1d.run_sim(config = CFG))

	# METHOD 1
	# The easiest way to configure is to use one of the predefined 'template' simulations.
	#	Calling one of these template functions passes back a Config object ready to go for some common
	#	use-case simulation. Then you can fongiure it further.
	# Example: Run the benchmark Spitzer-Harm simulation.
	#CFG = templates.spitzer_setup()
	#Oshun1d.run_sim( config = CFG)

	# METHOD 2
	# Another aletenative is to setup the object yourself (usually by copy-n-pasting) from another
	#	exampe or one of the 'template' configuration functions. Below is the copied-n-Pasted contents
	# 	of the benchmark Spitzer-Harm simulation
	# 	the most basic structral info....
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)

	num_species						=  1
	num_cells__global[0] 		=  1152/1
	system_size_min__global[0] = 0.0
	system_size_max__global[0] = 100.0 
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)


	CFG.num_species 				= 1
	CFG.num_l 						= 16 #4 #4
	CFG.num_m 						= 1 #3 #1
	CFG.dt 							= .1 #.299383 #.01
	CFG.t_end 						= 2000.0# 50.0*.298462 # 3880.0

	CFG.abs_p_max 					= 1.0
	CFG.abs_p_num_cells 			= 192
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 		= 3
	CFG.do_not_abort_on_abs_p_max_warning = False		# defaults to False. Oshun checks t see if you abs_p_max if large enough for the initalt temperature. If not it will print a warning and exit. To ignore the warning and still run th simulation, set this to True 

	# Currently, the support values are CFG.PERODIC and CFG.MIRROR
	CFG.boundryTypeEMF 				= CFG.PERODIC

	# feature (module) configuration
	CFG.use_wave_driver = False

	# paramteters controlling collisons..
	CFG.collsions_enabled			= False
	CFG.if_tridiagonal 				= True
	CFG.small_dt 					= .002  ##.05#0.001  # .1
	CFG.if_implicit1D 				= True

	# output configuration
	# For the following options, set the _interval to be non-zero in order for output to show.
	#	Note: Density is given as normalized density which is the density Normalized to the density_np input deck parameter (which is given in [cm^-3]).. so for example, a denisty of 1.0 means denisty_np in [cm^-3]
	CFG.num_timesteps_between_outputs 	= 100
	CFG.temperature_interval 			= 100
	CFG.Q_interval 						= 100
	CFG.heat_conduction_interval 		= 100
	CFG.density_interval				= 100
	CFG.v2_interval						= 0
	CFG.E_interval 						= 100


	# options when outputting the Distribution Function.
	CFG.F_output_interval 				= 0	# output X verus V... ie.e phasepsace.. i.e. the distirubtion function
	# set x limits on plot
	CFG.F_output_velocity_interval		= None 	# None means auto, to specify specific interval, set to something like [-0.1,0.1]
	# set vertical range (can do log, too)
	CFG.F_output_range					= None 	# None means auto, to specify specific interval, set to something like [-0.1,0.1]
	CFG.F_output_log_scale				= False	# display output on log scale
	# Normally, the output is 2D plot of colored squares (color bar on right shows value of distribution function)
	#	By setting the following options, we can also make an addtion plot where we fix space at certain cells, and look at velocity distrubtion
	CFG.F_output_lineouts				= [10]	# set spatial cells indexs at which to do velcity disoution lineouts. None means no lineouts. Otherwise set to something like [0] or something like [0,5,20]

	CFG.E_field_profile 			=  profiles.sin( amp=.1, num_wavelengths=1.0, wavelength=None, k=None, phase=0.0)

	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		
		# the density in 1/cm^3
		temp.density_np 			= 1.0e20
		temp.zeta 					= 1.0
		
		# set the boundry for the species distrubution function to the same as used by the EMF Fields.
		# (This is the only support mode for now.. you will be able to specifcy species boundries indepent of the EMF boundries soon)
		temp.boundryType 			= CFG.boundryTypeEMF
		
		# temp.min_allowed_temperature = 0.0001 by default
		#	This paramter ensures no zero temperatures are initalized. Any tmeperatures < min_allowed_temperature
		#	will be set to min_allowed_temperature.

		# To use an abritrary temperature profile, load in (or create) an array that is of size equal
		#	to num_cells__global[0] (i.e. the total number of spatial cells in the simulation). Say this array is named 'inital_profile'
		#	'inital_profile' can be specified in electron volts or in normalized proper momentum (i,e beta*gamma)
		#	Set the units='ev' if the profile is in electron volts or units='vth' if if is in normalized proper momentum
		#	An example that makes a boring 500 ev constant profile:
		ev_profile = np.ones( num_cells__global[0]  )
		ev_profile *= 1000
		temp.temp_profile = ArbitraryTemperatureProfile( profile=ev_profile, units='ev')
		
		# .024 is 300ev
		#temp.temp_profile  = SinusodialProfile( 0.025, pt_offset = 0.024 )
		# Note: another, newer way to use the sinusodial profile in terms of electron volts
		#temp.temp_profile  = SinusodialProfile( 20, pt_offset = 300,units='ev' )
		#temp.temp_profile = ConstantProfile( 4.420e-2)

		# If not specified, the density profile defaults to a constant profile of value 1.0
		# 	Here, for demonstration, we set a constant profile using the ArbitraryDensityProfile option.
		#	We don't specify a profile, but nstead only the number of cells. This returns an array of the
		#	appropiate size it it's .profile propery...which we can subsequently set to desired values.
		#	This also works for the temperature profile.
		#
		#	We could also simply use:
		#	temp.density_profile = ConstantProfile(1.0)
		#temp.density_profile = ArbitraryDensityProfile ( num_cells = num_cells__global[0] )
		#temp.density_profile.profile[:] = 1.0
		temp.density_profile = ConstantProfile(1.0)
	return CFG

def ld_setup__standing_wave__3rd_order():
	global CFG
	CFG = ld_setup_3rd()
	
	
	#num_cells__global[0] 		=  192
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	
	CFG.dt 							= .25	# .1
	CFG.t_end 						= 2000.0
	
	#1.0e-18 is border line non-linear...
	#k_min=( (2.0*np.pi)/(CFG.system_size_max__global[0] - CFG.system_size_min__global[0]))
	CFG.E_field_profile 			= None
	CFG.E_field_profile = 			profiles.sin( amp=1.0e-8, num_wavelengths=2, wavelength=None, k=None, phase=0.0, isCos=True)
	#CFG.E_field_profile = 			profiles.sin( amp=(-1.0e-1/k_min ), num_wavelengths=2, wavelength=None, k=None, phase=0.0, isCos=True)
	#CFG.E_field_profile 			= EFieldSinusodialProfile( 1.0e-8, num_wavelengths = 160)
	#CFG.E_field_profile.echo_time 							= 200.0
	#CFG.E_field_profile.second_pulse_amp 					= 1.0e-8
	#CFG.E_field_profile.second_pulse_num_wavelengths 	= 160.0
	CFG.use_wave_driver			= False
	
	return CFG

def ld_setup__standing_wave():
	global CFG 
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

# this worked very well and have leanr results all the way to t=2000 with half the spatial cells as above
def attempt2():
	CFG = ld_setup__standing_wave__3rd_order()
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)

	CFG.t_end 						= 2000

	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=2, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# this worked very well and have leanr results all the way to t=4000 with half the spatial cells as attempt2
def attempt3():
	# WITH num_cells__global[0] 		=  1152/2
	CFG = ld_setup__standing_wave__3rd_order()
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)

	CFG.t_end 						= 4000

	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=2, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# this worked very well and have leanr results all the way to t=4000 with 1/8 the spatial cells as attempt2
def attempt4():
	# WITH num_cells__global[0] 		=  1152/8 = 144
	CFG = ld_setup__standing_wave__3rd_order()
	
	CFG.num_l 						= 16	#24
	CFG.abs_p_max 					= .8
	CFG.abs_p_num_cells 			= 750	#600
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)

	CFG.t_end 						= 4000

	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=2, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# try a much higher kvth... set excited mode to 25 which is about k*vth=6.955615e-02
# works..
def attempt5():
	CFG = attempt4()
	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=25, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# k*vth = .105 works
def attempt6():
	CFG = attempt4()
	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=38, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# k*vth = .105 works
def attempt7():
	CFG = attempt4()
	CFG.E_field_profile = 			profiles.sin( amp=1.0e-18, num_wavelengths=75, wavelength=None, k=None, phase=0.0, isCos=True)
	return CFG

# ld_setup__standing_wave__3rd_order
#	work ok.. 
CFG = ld_setup__standing_wave__3rd_order()
#CFG = attempt2()
#CFG = attempt3()
#CFG = attempt4()
#CFG = attempt5()
#CFG = attempt6()
# CFG = attempt7()



Oshun1d.run_sim( config = CFG, profile=False, profileIncludeFilters=['Oshun1d'], profileExcludeFilters=None)


