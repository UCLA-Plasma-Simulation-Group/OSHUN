"""
	Temorpary deposit place for various routines that setup various classic simulation cases
		TODO: will move into a proper samples directory.
"""

from ..oshun1ed_ld import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile, run_sim, WaveDrivah

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
	CFG.num_l 						= 5 #4 #4
	CFG.num_m 						= 1 #3 #1
	CFG.dt 							= .298462 #.299383 #.01
	CFG.t_end 						= 3880.0 #3880.0# 50.0*.298462 # 3880.0
	
	CFG.abs_p_max 					= 0.15
	CFG.abs_p_num_cells 			= 108
	CFG.abs_p_min					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 4
	
	# feature (module) configuration
	CFG.use_wave_driver = False
	
	# paramteters controlling collisons..
	CFG.collsions_enabled		= True
	CFG.if_tridiagonal 			= True
	CFG.small_dt 					= .1  ##.05#0.001  # .1
	CFG.if_implicit1D 			= True
	
	# output configuration
	CFG.num_timesteps_between_outputs = 100
	CFG.F_output_interval = 0
	CFG.temperature_interval = 100
	CFG.Q_interval = 100
	CFG.heat_conduction_interval = 100
	
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
		temp.ln_lambda 				= 10.2
		temp.zeta 					= 1.0
		# constant density gradient at 1.0
		
		# .024 is 300ev
		temp.pt_amb 				= 0.024
		temp.temp_profile 		= SinusodialProfile( 0.025 )
		#temp.temp_profile  = ConstantProfile( 4.420e-2) 
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
	
	CFG.num_timesteps_between_outputs = 1
	CFG.F_output_interval = 0
	CFG.temperature_interval = 1
	CFG.Q_interval = 1
	CFG.heat_conduction_interval = 1
	CFG.density_interval = 1
	CFG.E_interval = 1
	
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
		CFG.E_field_profile 			= None
		CFG.use_wave_driver			= False
		
	return CFG

def setup_ben_chappman_profile():
	num_cells__global 		= np.zeros(3)
	system_size_min__global = np.zeros(3)
	system_size_max__global = np.zeros(3)
	
	num_species						= 1
	num_cells__global[0] 		= 120 #4000
	system_size_min__global[0] = -1243.7
	system_size_max__global[0] = 1243.7 #3580.0
	
	CFG = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	
	CFG.num_timesteps_between_outputs = 1
	CFG.F_output_interval = 0
	CFG.F_output_range = [-.6,.6]
	CFG.temperature_interval = 1
	CFG.Q_interval = 1
	CFG.heat_conduction_interval = 1
	CFG.density_interval = 1
	CFG.E_interval = 1
	
	CFG.num_l 						= 6
	CFG.num_m 						= 1
	CFG.dt 							= .01
	CFG.t_end 						= 1000.0 #300.0 #.11    #1000.0
	
	CFG.abs_p_max 					= 3.0
	CFG.abs_p_num_cells 			= 512 
	CFG.abs_p_min 					= CFG.abs_p_max/(2.0 * CFG.abs_p_num_cells)
	CFG.explicit_solver_order 	= 3
	CFG.if_tridiagonal 			= True
	
	CFG.collsions_enabled		= False
	CFG.small_dt 					= 0.001
	CFG.if_implicit1D				= True
	
	CFG.E_field_profile 			= None
	CFG.use_wave_driver			= False 
	
	for i in range(0, CFG.num_species):
		temp = SpeciesConfig()
		CFG.species_config[i] = temp
		
		temp.name = "Species %d" % i
		temp.charge 				= -1.0
		temp.mass 					= 1.0
		temp.vth 					= 0.0001

		# the density in 1/cm^3
		temp.density_np 			= 4.5e20
		temp.ln_lambda 			= 100.2
		temp.zeta 					= 10.0
		
		# .024 is 300ev
		temp.pt_amb 				= 1.396e-4 #4.020e-2 # 0.0024	
		temp.temp_profile 		=  ConstantProfile( 1.088e-1 ) #2.493e-1) #4.420e-2) #1.088e-1 ) #GaussianProfile( 0.33166, 178.0, 100.0 )
		
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