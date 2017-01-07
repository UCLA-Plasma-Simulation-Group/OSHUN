import numpy as np
import os.path

from vqunits import Formulary 

import oshun1d_ld
from oshun1d_ld import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile, run_sim

import sys
sys.path.append("../trunk1d")



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

CFG = ld_setup__standing_wave__3rd_order()


# if this is the actual script run, then start Oshun... if this script is imported by another script,
# then that script needs to call run_sim
if __name__ == "__main__":
	run_sim( config = CFG)
