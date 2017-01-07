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

class Sine( ):
	def __init__( self, amp=1.0, num_wavelengths=None, wavelength=None, k=None, phase=0.0, isCos=False):
		self.amp = amp
		self.num_wavelengths = num_wavelengths
		self.wavelength = wavelength
		self.k = k
		self.phase = phase
		self.isCos = isCos

		if num_wavelengths == None and wavelength==None and self.k == None:
			self.num_wavelengths = 1.0

	# pass in an appropiate buffer 'data'
	def apply(self, data, spatial_axis_global, mpi_axes, species=None):
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
			mpi_axes[0].num_boundry_cells
			indic = []
			#for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			#	global_cell_index = float( global_cell_offset[0] + n )
			#	data[ E_field.num_boundry_cells[0][0] + n] = self.amp * np.sin(discrete_factor*global_cell_index + phase_cell_offset + phase_me )
			spatial_axis__local_in_global = mpi_axes[0].axis__local
			num = len(spatial_axis__local_in_global.values) + mpi_axes[0].num_boundry_cells[0][0] + mpi_axes[0].num_boundry_cells[0][1]
			sincos = np.sin
			if self.isCos:
				sincos = np.cos
				#print "USE SINE@@"
			xs1 = []
			for n in xrange(0, num):
				global_cell_index = float( mpi_axes[0].global_cell_offset[0] + n - mpi_axes[0].num_boundry_cells[0][0])
				data[n] = self.amp * sincos(2.0*np.pi*float( self.num_wavelengths)/(box_size +dx )*global_cell_index*dx + self.phase) # + phase_cell_offset + phase_me )
				
				xs1.append( global_cell_index * dx )
				
				#print "DRIVAVVV: DIVIDING BY K WHICH IS: %e" % (2.0*np.pi*float( self.num_wavelengths)/(box_size +dx ))
				#global_cell_index = float( global_cell_offset[0] + n - E_field.num_boundry_cells[0][0])
				##indic.append( global_cell_index )
				#data[n] = self.amp * np.cos(discrete_factor*global_cell_index) # + phase_cell_offset + phase_me )
			

			self.wavelength = box_size / float( self.num_wavelengths)
			self.wavenumber = 2*np.pi / self.wavelength
			#self.v_phase = 1.0/self.wavenumber
			return
			# TODO: wied AXIS.. need to go trough and normailze this once and for all
			#	Derivitives too.
			
		elif self.wavelength != None:
			# use a specifc wavelength...
			self.wavenumber = 2.0*np.pi/(self.wavelength)
			#self.num_wavelengths = box_size/self.wavenumber 
		else:
			self.wavenumber = self.k
			self.wavelength = 2*np.pi / self.wavenumber
			self.num_wavelengths = box_size/self.wavelength 

		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			data[ mpi_axes[0].num_boundry_cells[0][0] + n] = self.amp * np.sin(self.wavenumber*x_val + self.phase )
		self.v_phase = 1.0/self.wavenumber

	# this routine assume that the guardcells are in synch correctly.
	def computeDensityPerturbation(self,data, spatial_axis_global, mpi_axes, species=None,E_field_perturbation = None):

		# TODO: make the API better
		if ( E_field_perturbation == None):
			return


class Profiles:
	def __init__(self):
		self.sin = Sine

profiles = Profiles()