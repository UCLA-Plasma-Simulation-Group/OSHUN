# hi
import ProfileCommonUtilities
import numpy as np

"""
	This class expects the global profile (i.e. the full profile )
"""
class ArbitraryDensityProfile:
	def __init__(self, profile=None, num_cells=None ):
		if( profile == None and num_cells == None ):
			print "Error setting up ArbitraryDensityProfile. You must specify either the profile or num_cells parameter."
		if(profile == None):
			self.profile = np.zeros( num_cells )
		else:
			self.profile = profile
		
		self.apply_density = self.apply

	def apply( self, local_density_profile, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		if( spatial_axis_global.num_points != self.profile.shape[0] ):
			print "Error: The inital density profile was specified in a buffer that has %d elements. \n       But the simulation has %d spatial cells. These two sizes must be the same." % ( self.profile.shape[0], spatial_axis_global.num_points)
		
		this_nodes_starting_index =  mpi_axes[0].global_cell_offset
		num_boundry_cells = mpi_axes[0].num_boundry_cells

		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			global_profile_index = this_nodes_starting_index[0] + n
			local_profile_index = int(num_boundry_cells[0][0] + n)
			local_density_profile[local_profile_index] *= self.profile[ global_profile_index ]
		pass

