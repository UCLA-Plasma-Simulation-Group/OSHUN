# hi
import ProfileCommonUtilities
from ..OshunMpi import mpi_info
from ..formulas.formulary import formulary
import numpy as np

"""
	This class expects the global profile (i.e. the full profile )
"""
class ArbitraryTemperatureProfile:
	def __init__(self, profile=None, num_cells=None, units=None ):
		
		if(units == None):
			if( mpi_info.is_I_root()):
				print "Error. The units paramter must be set when using ArbitraryTemperatureProfile."
				print "		Valid values are 'vth' and 'ev'"

		self.units = units
		if( profile == None and num_cells == None ):
			print "Error setting up ArbitraryTemperatureProfile. You must specify either the profile or num_cells parameter."
		if(profile == None):
			self.profile = np.zeros( num_cells )
		else:
			self.profile = profile
	
	def apply( self, local_temperature_profile, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		if( spatial_axis_global.num_points != self.profile.shape[0] ):
			print "Error: The inital temperature profile was specified in a buffer that has %d elements. \n       But the simulation has %d spatial cells. These two sizes must be the same." % ( self.profile.shape[0], spatial_axis_global.num_points)
		
		this_nodes_starting_index =  mpi_axes[0].global_cell_offset[0]
		num_boundry_cells = mpi_axes[0].num_boundry_cells

		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			global_profile_index = this_nodes_starting_index + n
			local_profile_index = num_boundry_cells[0][0] + n
			temperature_profile_value = self.profile[ global_profile_index ]
			if( self.units == 'ev'):
				temperature_profile_value = formulary.convert_ev_to_vth( temperature_profile_value, species.species_config.mass)
			local_temperature_profile[local_profile_index] *= temperature_profile_value
		pass
