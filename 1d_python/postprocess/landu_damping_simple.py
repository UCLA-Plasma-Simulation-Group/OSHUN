import numpy as np
import glob
from pylab import *
import cPickle as pickle
import sys

sys.path.insert(0, "../source/python")
import oshun_modules.OshunGlobalState as GLOBAL
from oshun_modules.formulas import formulary

GLOBAL.USE_C_VERSION = False

def getInitalTemperature():
	# now, read in the inital temperture
	temperatureProfileFilenames = glob.glob("./output/numpy/temp_ev_*.npy")
	if len(temperatureProfileFilenames) < 1:
		print "Error. Could not read in the inital temperature profile. Exiting"
		exit(-1)
	intialTemperatureProfile = np.load(temperatureProfileFilenames[0])
	temperatureMin = np.min(intialTemperatureProfile)
	temperatureMax = np.max(intialTemperatureProfile)
	if np.abs(temperatureMax-temperatureMin) > 1e-6:
		print "Error. The intial Temperature profile is not costant. difference between amx and min is %e We expect it to be so.. Exiting" % ( np.abs(temperatureMax-temperatureMin),  )
		exit(-1)
	return temperatureMin


# read in the inputdeck
sim_metadata     = pickle.load( open( "./output/misc/oshun_metadata.p", "rb" ) )
sim_metadata_cfg = pickle.load( open( "./output/misc/oshun_metadata__cfg.p", "rb" ) )


intialTemperature__ev 		= getInitalTemperature()
density						= 1.0 # sim_metadata_cfg.species_config[0].density_profile.val

output_step 				= sim_metadata.output_interval__time_steps   #1
time_step_size 				= sim_metadata.dt   									#.106153

box_size__x 				= sim_metadata_cfg.system_size_max__global[0]
num_cells__x 				= sim_metadata_cfg.num_cells__global[0]



vth_e						= formulary.convert_ev_to_vth(intialTemperature__ev, 1.0) # * np.sqrt(1.5)
wp_e						= np.sqrt( density )


#k_of_excited_mode 		= 2.0*np.pi / ( box_size__x / float( mode_excited)  )
k_of_excited_mode			= sim_metadata_cfg.E_field_profile.wavenumber
amplitude_of_excited_mode 	= sim_metadata_cfg.E_field_profile.amp
k__vth_e					= k_of_excited_mode*vth_e
bounce_time					= 2.0*np.pi / np.sqrt( amplitude_of_excited_mode*k_of_excited_mode)

output_time_step = float(output_step) * time_step_size

f = formulary
density__cgs 			= sim_metadata_cfg.species_config[0].density_np
vth__cgs 				= f.convert_osi_vth_to_beta ( vth_e ) * 3e10
wp_e__cgs 				= f.wpe( density__cgs )
skindepth__cgs 			= 3e10/wp_e__cgs
k__cgs 					= k_of_excited_mode / skindepth__cgs
debye_length 			= f.debye_length(intialTemperature__ev,  density__cgs)
k_vth__cgs 				= k__cgs*vth__cgs


print ""
#print "w of driven wave                                   %f" % (w_wave_driven)
print "dt (time step size)                                %f" % (time_step_size)
print "L (Box Size)                                       %f" % (box_size__x)
print "number spatial cells                               %d" % (num_cells__x)
print "spatial cell size                                  %f" % (box_size__x / float(num_cells__x) )  # this is wrong-ish
print ""

print "v thermal, electrions, normized units from code:   %f" % ( vth_e )
print "k, normalized                                      %e" % k_of_excited_mode
print "bounce time (normaied to 1/wp)                     %e" % bounce_time
print "denisty, cgs                                       %e" % density__cgs
print "skindepth                                          %e cm" % skindepth__cgs
print "vth                                                %e cm/s" % vth__cgs
print "temp                                               %e ev" % intialTemperature__ev
print "debye length                                       %e" % debye_length
print "wp                                                 %e" % wp_e__cgs
print "k                                                  %e" % k__cgs
print "kvth                                               %e" % k_vth__cgs
print "k*l_debye                                          %e" % ( k__cgs * debye_length)
print "wp_e/(vth_e*k)                                     %e" % ( wp_e__cgs / ( vth__cgs * k__cgs) )
print "                                                   %e" % (k__vth_e, )


EfieldVsTimeFilenames = glob.glob("./output/numpy/E_*.npy")
data = np.zeros( len( EfieldVsTimeFilenames) )

for i,fname in enumerate(EfieldVsTimeFilenames):
	dataAtTimestep = np.load(fname)
	dataAtTimestep *= dataAtTimestep
	maxAtTimestep = np.max(dataAtTimestep)
	data[i] = maxAtTimestep
plot(data)
title("E Field Energy Vs. Time")
xlabel("Time")
ylabel("E Field Energy")
savefig('test.png')

clf()
plot(np.log(data))
title("Log E Field Energy Vs. Time")
xlabel("Time")
ylabel("Log E Field Energy")
savefig('test_log.png')
