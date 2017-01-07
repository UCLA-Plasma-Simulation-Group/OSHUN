
# NOTES TODOs BUGS
# "density_multiplier" not set?!??! this could cause lots of problems...
#
import cProfile

import OshunGlobalState as GLOBAL
from OshunGlobalState import PRINT as PRINT

USE_C_VERSION_GLOBAL = GLOBAL.USE_C_VERSION

# python system libraries
import os
import math
import traceback
import inspect
import cPickle as pickle
import re

# numerical libraries
import numpy as np
import scipy.linalg
import scipy.integrate
import scipy.fftpack as fft
import scipy.special
import scipy.integrate
from scipy.stats import norm
import scipy.special
import scipy.optimize
import scipy.stats

# oshun libraries
from InputDeckConfiguration import Config, SpeciesConfig

import oshun_modules.OshunEventQueue as OshunEventQueue
from .diagnostics.moments import MomentCalculationsState
from OshunMpi import *
import Oshun1dEvents
from FokkerPlankCollisionOperator import *
from .formulas.formulary import formulary
 
# If we wish to use the C-Python code, then include the following.
if USE_C_VERSION_GLOBAL:
	PRINT("Using the python version mixed with c optimized routines!")
	import numpy.ctypeslib
	from OshunCInterface import *
	if mpi_info.is_I_root:
		printCConfigurationHeader()
	
from OshunImplictESolver import RK2ForImplicit
from oshun_modules.utils.timing import timings,timeit,TIMINGS
from oshun_modules.utils.debugging import *

#
#globals
# global timer utility (rembers timer names and line numers etc..)
#TIMINGS = timings()


"""
Stupid class that will probably be replaced that stores metadata about the simualtion.
	Different modules can write into this metadata store. Usually they will write out
	their relevent configuration paramters/calcualted values so that these can be used by
	postprocessing scripts. This data is written to disk (by the MPI root node) as a 
	Python 'pickle' file.
"""
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
			pickle.dump( self, open( GLOBAL.OUTPUT_DIR("oshun_metadata.p", "metadata"), "wb" ) )
		self.dirty = False

"""
Yet another basic Axis abstract class. The world is filled with so many.
	But rather then being just an empty container class, there is some important
	stuff in here. Namely, the way we do indexing, which the was decided by the Greek.
"""
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
			#PRINT("AXIS MADE!!!")
			#PRINT(min)
			#PRINT(max)
			#PRINT(num_points)
			#PRINT(self.values)
			#(self.values, self.step) = np.linspace(min, max, num=num_points, endpoint=False, retstep=True)
			return
		if min != None and step != None and num_points != None:
			self.min = min
			self.num_points = num_points
			self.max = num_points*step + min
			self.values = [(x*step+min) for x in range(0, num_points)]
			return
		raise Exception("Invalid Axis specification..")
		self.num_cells = self.num_points
	# These are just syntax sugar
	def dx(self):
		return self.step
	def dp(self):
		return self.step


"""
A 1D (in space) field with 'num_components' at each each spatial cell.

The underlying data is arranged as a 2D array:
	Field[component_index][spatial_cell_index]
		where 	component_index=0..num_components-1 (number of componets in Field)
				spatial_cell_index=0..num_cells-1 (number spatial cell along x direction)

	TODO: At this point, the first dimension is implemented as a Python list. So the overall strucuture
			is a list of arrays.. each array holds the value of a specific field component at each 
			spatial cell (ie. [component 0 at all x,component 1 at all x,component 2 at all x ] ).
			This needs to change to be a proper 2D numpy array.
	Note: When possible, specific operations (functions) are called on the Field. You should perfer this
			rather then exposing the Field data and doing explicit array calcualtions in other parts of the code.
			 Basically abstract away the underlying data layout of the Fields.
	Note: Axes and distances in this class are all local to the node, with guard cells.
"""
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
		
		self.temp_data = np.zeros( self.num_total_cells[0] )
		self.temp_data1 = np.zeros( self.num_total_cells[0] )
		
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
	def zeroOut( self ):
		for i in xrange(0, self.num_components):
			self.components[i] *= 0.0
	def X(self):
		return self.components[0]
	def Y(self):
		return self.components[1]
	def Z(self):
		return self.components[2]
	def comp(self, component):
		return self.components[component]
	def copyFrom(self, source, component=-1):
		if( component == -1):
			self.components[0][:] = source.components[0][:]
			self.components[1][:] = source.components[1][:]
			self.components[2][:] = source.components[2][:]
		else:
			self.components[component][:] = source.components[component][:]
	def copyTo(self, target, component=-1):
		if( component == -1):
			target.components[0][:] = self.components[0][:]
			target.components[1][:] = self.components[1][:]
			target.components[2][:] = self.components[2][:]
		else:
			target.components[component][:] = self.components[component][:]
	def Dx(self, component):
		
		data = self.components[component]
		if False:
			self.temp_data[0,:] = 0.0
			self.temp_data1[-1:0] = 0.0
			
			self.temp_data [1:,:] = self.data[:-1,:]; 
			self.temp_data1[:-1,:] = self.data[1:,:]; 
			self.temp_data1 -= self.temp_data
			self.temp_data1 *= 8.0
			
			self.temp_data1[2:,:] += self.data[:-2,:]; 
			self.temp_data1[:-2,:] -= self.data[2:,:];
			
			self.temp_data1[-1,:] = self.data[-1,:];
			
			self.data[:,:] = self.temp_data1[:,:]
			self.data /= 12.0
			self.data[:,0] = self.temp_data1[:,2]
			self.data[:,1] = self.temp_data1[:,2]
		
		if True:
			for x in range(0, self.num_total_cells[0]-2):
				data[x] -= data[x+2]
			for x in range(self.num_total_cells[0]-3, -1, -1):
				data[x+1] = data[x]
		
		# normalize
		# The factor of 2.0 is here becasue this is a central differnece derivitve.
		data /= 2.0*STATE.spatial_axes[component].dx()

"""
Class that represents a (velocity space) spherical harmonic component of the particle distribution function at
	a specific spatial cell. The core data this class contains is an array of length=number_abs_p_cells. Each value
	represents the specific harmonics's contribution to the ditribution function, at a speciifc spatial cell, 
	as it's radius is scaled to each of the values of the |p| cell grid.
"""	
class SHarmonic:
	def __init__(self, l, m, linear_index=-1, momentum_axis = None, distribution_function_data=None):
		self.my_l = l
		self.my_m = m
		self.my_linear_index = linear_index
		self.init(momentum_axis, distribution_function_data)
	def init(self, momentum_axis, distribution_function_data = None):
		self.momentum_axis = momentum_axis
		if distribution_function_data != None:
			self.momentum_projection = distribution_function_data
		elif momentum_axis != None:
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
				raise Exception("Cannot copy sphereical harmonic because 'momentum_axis' is Null is source.")
		
		self.momentum_axis = sh_to_copy.momentum_axis
		#np.copyto(self.momentum_projection, sh_to_copy.momentum_projection)
		#TODO: replace with faster
		self.momentum_projection[:] = sh_to_copy.momentum_projection[:]
		

"""
Class that represents a complete spatial cell in a particle ditribution function.
	Each spatial cell contains several (velocity space) sphereical harmonics (i.e. this class is
	a collection of SHarmonic classes)
"""	
class DistributionFunctionSpatialCell:
	def __init__(self, num_l, momentum_axis = None, distribution_function_data=None):
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
			self.harmonics.append( SHarmonic(l,0,len(self.harmonics),momentum_axis, distribution_function_data=distribution_function_data[l]) )
		pass
		
		# test quick optimzation when trying some big runs to speed for loops over assest...
		#	hold a list all of 'views' (basically pointers/references) to all this cell's SHarmonic momentum_projection arrays.
		#self.harmonic_data = [None]*num_l
		#for l in xrange(0,num_l):
		#	self.harmonic_data[l] = self.harmonics[l].momentum_projection
		
		
	def init(self):
		pass
	
	def copy(self, f_to_copy):
		if self.num_harmonics != f_to_copy.num_harmonics:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_harmonics' parameters")
		for i in xrange(0, self.num_harmonics):
			self.harmonics[i].copy( f_to_copy.harmonics[i] )
		pass


"""
Particle ditribution function (for a specific species of particles) which logically contains
	an array of DistributionFunctionSpatialCell objects.
"""
# cells and distances in this class are all local to the node...
class DistributionFunction:
	
	def __init__(self, num_cells, num_l, momentum_axis, num_boundry_cells):
		
		self.num_plain_cells 	= np.array( [num_cells,0,0] )
		self.num_boundry_cells 	= num_boundry_cells
		self.num_total_cells 	= np.zeros( (3,), np.int )
		self.num_total_cells[0] = self.num_plain_cells[0] + self.num_boundry_cells[0][0] + self.num_boundry_cells[0][1]
		#self.num_total_cells[1] = self.num_plain_cells[1] + self.num_boundry_cells[1][0] + self.num_boundry_cells[1][1]
		#self.num_total_cells[2] = self.num_plain_cells[2] + self.num_boundry_cells[2][0] + self.num_boundry_cells[2][1]
		#self.num_total_cells[1] = 0
		#self.num_total_cells[2] = 0
		
		
		self.dim_x = self.num_plain_cells[0]; 
		self.dim_y = self.num_plain_cells[1]; 
		self.dim_z = self.num_plain_cells[2];
		self.momentum_axis = momentum_axis

		#quick optimzation I put when doing a quick pass with a profilier.. 
		#	There were many scattered hot spots that involved the disjointed Dist function data layout
		#	(each spitial cell had its own set of buffers for the harmonic projection data.. so itterating
		#	over spatial axes was especialy expensive)
		#  So I changed things so that we allocate all distribution function memory as one chunk.. then
		# then we can use numpy indices to slice things (much faster then for-loops)
		# And we can now easily slice data over spaital indices/axes as well! 
		# Each level of the Dist function (spaital cells, SHarmonics, etc.. is given a pointer to the start
		#		of it's little slice of the single big buffer.
		#
		# the array (for now) is of shape:  [ self.num_total_cells[0], num_l, num_p]
		self.distribution_function_data = np.zeros( [ self.num_total_cells[0], num_l, momentum_axis.num_points] )
		
		# create spatial cell entries...
		self.cells= [DistributionFunctionSpatialCell(num_l, momentum_axis, distribution_function_data=self.distribution_function_data[x] ) for x in range(0, self.num_total_cells[0])]
		
		self.num_l = self.cells[0].num_l
		self.num_m = 1
		
		# create a mapping to quickly link (l,m) coords to linear..
		self.harmonic_ordinal_map = [None] * (self.num_l)
		for i in range(0, self.cells[0].num_harmonics):
			l = self.cells[0].harmonics[i].my_l
			m = 0
			linear_index = self.cells[0].harmonics[i].my_linear_index
			self.harmonic_ordinal_map[l] = linear_index;
	def clone_structure( self ):
		new_f = DistributionFunction( self.num_plain_cells[0], self.num_l, self.momentum_axis,  self.num_boundry_cells)
		return new_f
	
	def copy(self, f_to_copy):
		if f_to_copy.num_l != self.num_l:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_l' parameters")
		if f_to_copy.num_m != self.num_m:
			raise Exception("Cannot copy distrubtion function cause they gots unequal 'num_m' parameters")
		# SHIT LOOKIN THIS ADAM BUG
		for x in xrange(0, self.num_total_cells[0]):
			self.cells[x].copy( f_to_copy.cells[x] )
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
		#np.multiply( self.distribution_function_data , scalar, out = self.distribution_function_data )
		self.distribution_function_data *= scalar
		"""
		for x in xrange(0, self.num_total_cells[0]):
			for i in xrange(0, self.cells[x].num_harmonics):
				self.cells[x].harmonics[i].momentum_projection *= scalar
		pass
		"""
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
		hcv.data[:,:] = self.distribution_function_data[:, idx, :]
		#for x in xrange(0, self.num_total_cells[0]):
		#	hcv.data[x,:] = self.cells[x].harmonics[idx].momentum_projection[:]
	def get_num_harmonics_per_cell(self):
		return self.cells[0].num_harmonics
	def sub_from_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		for x in xrange(0, self.num_total_cells[0]):
			self.cells[x].harmonics[idx].momentum_projection[:] -= hcv.data[x,:]		
	def add_from_harmonic_centric_view(self, l, hcv):
		# change PROFILE
		idx = self.linear_index_for_harmonic(l)
		self.distribution_function_data[:, idx, :] += hcv.data[:,:]
		#idx = self.linear_index_for_harmonic(l)
		#for x in xrange(0, self.num_total_cells[0]):
		#	self.cells[x].harmonics[idx].momentum_projection[:] += hcv.data[x,:]
	def copy_from_harmonic_centric_view(self, l, hcv):
		idx = self.linear_index_for_harmonic(l)
		self.distribution_function_data[:, idx, :] = hcv.data[:,:]
		#for x in xrange(0, self.num_total_cells[0]):
		#	cells[x].harmonics[idx].momentum_projection[:] = hcv.data[x,:]
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

	# updates the spatial density profile to be pertrubed by 'pertrubation'. Perturbation should
	#	have the dimensions of the a Field1D object (boundry cells + number of spatial cells on this parrell partition node)
	#
	#	Pass in either the moment calculator to calcualte the density, or the density profile itself.
	#
	#	TODO: passing in the 'momentCalculator' is a bit akward. Change to have species know it's index.
	def applyDensityPerturbation( self, pertrubation, scratchSpatialCellBuffer=None, momentCalculator=None, presentDensityProfile=None ):

		if( pertrubation == None):
			return

		if scratchSpatialCellBuffer == None:
			scratchSpatialCellBuffer = np.copy( pertrubation )

		if momentCalculator != None:
			# get the present density profile for this Distribution function by doing hte 0th moment intergral.
			presentDensityProfile = momentCalculator.get_density(force_update = True)
			pass

		if presentDensityProfile == None:
			# TODO: throw error.
			return
				
		print presentDensityProfile
		scratchSpatialCellBuffer[:] += presentDensityProfile[:]

		# there are several viable operations here.
		#	Just alter to f00 at each spaital point:
		#		get_harmonic_at_location(x, 0)
		#	Renoamralize all harmonics at the cell:
		#		scalar_mult_spatial_cell(scalar, x)
		#	And you can also do (not recommended) to update the F distrubution directly (just f00 harmonic)
		#		F.distribution_function_data[ix,0,:] *= perturbattion_value
		#
		#	This loop also updates Guard cells.. so an additional MPI sync of guard cells is not needed.
		#		But additional Guard cells synch never hurt and to reduce interactions, we may do one here just in case if the number of module increases
		#		And it is assumed that this funtion is rarely called... maybe afew times a most during a simulation.
		for ix in xrange( 0, self.num_total_cells[0] ):
			perturbationScale = scratchSpatialCellBuffer[ix]
			self.scalar_mult_spatial_cell( perturbationScale, ix )

"""
This is a fundimental workhorse class that bundles all the data for a specific harmonic. It is a view that groups
	the information in a particle distribution function for imporved readablily and performance.

It is a 2D array that is of dimension:
	HarmonicCentricView[number_spatial_cells][number_abs_momentum_cells]

So this gives the velociy-space projection projection of a specific harmonic at each spatial cell.
"""
class HarmonicCentricView:
	
	def __init__(self, num_total_cells, momentum_axis):
		self.momentum_axis = momentum_axis
		self.num_total_cells = num_total_cells
		self.data = np.zeros( (self.num_total_cells[0], self.momentum_axis.num_points) )
		self.temp_data = np.zeros( (self.num_total_cells[0], self.momentum_axis.num_points) )
		self.temp_data1 = np.zeros( (self.num_total_cells[0], self.momentum_axis.num_points) )
	
	def scalar_mult(self, val):
		self.data *= val
	
	def Dp(self):
		abs_p_num_points = self.momentum_axis.num_points
		"""
		# This routine is a faster version using numpy broadcasting rather then 'for' loops
		# The original 'for' loop version is

		for x in xrange(0, self.num_total_cells[0]):
			for i in range(0, abs_p_num_points-2):
				self.data[x,i] -= self.data[x,i+2]
			for i in range(abs_p_num_points-3, -1, -1):
				self.data[x,i+1] = self.data[x, i]
		"""
		# do a derivitive along |p|... set zeroth element to 0
		#		for last element, pretend f(nump) = f(nump-1) and calc the deritivitive. 
		self.temp = np.zeros(self.num_total_cells[0] )
		self.temp[:] = self.data[:,-2]
		self.temp[:] -= self.data[:,-1]
		
		self.data[:, 1:-1] = ( self.data[:, :-2] - self.data[:, 2:] )

		self.data[:, 0] = 0.0
		self.data[:, -1] = self.temp[:]
	
	def __call__(self, abs_p, x):
		return self.data[x,abs_p]
	
	"""
		Repeat over all spatial cells:
			Element-wise multiply all values of |p| with the array 'mul_array' passed in
				(	'mul_array' is the same at all spatial cells and 
					'mul_array' must be the same size as momentum_axis.num_points )
	"""
	def mpaxis(self, mul_array):
		# This routine is made faster using numpy broadcasting rather then 'for' loop
		# Replace:
		#	for x in xrange(0, self.num_total_cells[0]):
		#		self.data[x,:] *= mul_array
		# With:
		self.data *= mul_array

	"""
		Repeat over all spatial cells.
			Element-wise multiply all values of |p| by a constant 'const'
				( 	'const' is different at each spatial cell
					'const' array must be same size as number of spatial cells)

		Input: values = 1d array of 'const' values for each spatial cell.
	"""
	def mpaxis__constant_at_each_spatial_cell(self, values):
		# values[x] evaluates to a constant.
		for x in xrange(0, self.num_total_cells[0]):
			self.data[x,:] *= values[x]
	
	#TODO: this can be done better using numpy magic...
	def Dx(self):
		if False:
			# faster version..
			# TODO: tme this and validate that it works..
			abs_p_num_points = self.momentum_axis.num_points
			num_total_cells = self.num_total_cells[0]
			
			
			self.temp_data[0,:] = 0.0
			self.temp_data1[-1:0] = 0.0
			
			self.temp_data [1:,:] = self.data[:-1,:]; 
			self.temp_data1[:-1,:] = self.data[1:,:]; 
			self.temp_data1 -= self.temp_data
			self.temp_data1 *= 8.0
			
			self.temp_data1[2:,:] += self.data[:-2,:]; 
			self.temp_data1[:-2,:] -= self.data[2:,:];
			self.temp_data1[-1,:] = self.data[-1,:];
			
			self.data[:,:] = self.temp_data1[:,:]
			self.data /= 6.0
			
			self.data[0,:] = self.data[2,:]
			self.data[1,:] = self.data[2,:]
			
			"""
			self.temp_data[:,0] = 0.0
			self.temp_data1[:,-1] = 0.0
			
			self.temp_data [:,1:] = self.data[:,:-1]; 
			self.temp_data1[:,:-1] = self.data[:,1:]; 
			self.temp_data1 -= self.temp_data
			self.temp_data1 *= 8.0
			
			#self.temp_data1[:,0:2] = 0.0
			#self.temp_data1[:,-1] = 0.0
			
			self.temp_data1[:,2:] += self.data[:,:-2]; 
			self.temp_data1[:,:-2] -= self.data[:,2:];
			self.data[:,:] = self.temp_data1[:,:]
			self.data /= 6.0  #(the other 2.0*h) will come later in the code
			"""
		
		# second order deritivite..
		if True:
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

"""
This class contains the Particle Distribution function and the properties of the particle type
	(i.e. mass, charge etc). This class actually conatians several versions of the PDisribution fucntion Collection of Particle Distribution functions and related particle quantities for similar particles in the simualtion.
"""
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

	def apply(self, temperature_grid, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			temperature_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
			#PRINT("Assigned (%d, %f) val %f %d" % ( species.F.num_boundry_cells[0][0] + n, x_val, temperature_grid[species.F.num_boundry_cells[0][0] + n], species.F.num_boundry_cells[0][0]))
		pass
	
	def apply_density( self, density_grid, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			density_grid[species.F.num_boundry_cells[0][0] + n] *= self.val
		pass	
		
class SinusodialProfile:
	def __init__(self, ampltitude, pt_offset = None, num_wavelengths=1.0, phase=0, average_value=0.0, units=None):
		self.ampltitude = ampltitude
		self.num_wavelengths = num_wavelengths
		self.phase = phase
		self.average_value = average_value

		if( pt_offset != None):
			# this is legacy.. i will remove it later..
			self.average_value = pt_offset

		self.units = units
		if(units == None):
			self.units = 'vth'

	def apply(self, temperature_grid, spatial_axis_global, mpi_axes, species=None):
		if(self.units=='ev'):
			spatial_axis__local_in_global = mpi_axes[0].axis__local
			box_size = (spatial_axis_global.values[-1] - spatial_axis_global.values[0] )/ self.num_wavelengths
			dx = spatial_axis_global.dx()
			for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
				temperature_grid[species.F.num_boundry_cells[0][0] + n] = self.average_value + self.ampltitude*np.sin(2.0*np.pi*x_val/(box_size+dx) + self.phase )
				temperature_grid[species.F.num_boundry_cells[0][0] + n] = formulary.convert_ev_to_vth( temperature_grid[species.F.num_boundry_cells[0][0] + n], species.species_config.mass) 
			pass			
		else:
			self.apply_vth(temperature_grid, spatial_axis_global, mpi_axes, species=species)

	def apply_vth(self, temperature_grid, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		self.umax = self.average_value
		self.Aver = self.ampltitude
		self.Ampl = 0.5*(self.umax*self.umax - self.Aver*self.Aver)
		self.Aver = 0.5*(self.umax*self.umax + self.Aver*self.Aver)
		
		box_size = (spatial_axis_global.values[-1] - spatial_axis_global.values[0] )/ self.num_wavelengths
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
			#PRINT("%d out of %d" % ( i, STATE.config.num_species))
			#PRINT(extant_charge_density[-1][:])
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
		
	# should prolly move this...
	def apply_perturbation(self, F, E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset ):
		# tkae incoming F and make it into (1+amp*cos(k*x))*F
		#PRINT("NO!!!")
		#exit(-1)
		# first make the cosine...
		box_size = spatial_axis_global.values[-1] - spatial_axis_global.values[0]
		
		
		if self.num_wavelengths != None:
			buffer = E_field.clone( )
			self.apply( E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset, buffah=buffer )
			
			buffer.components[0] += 1.0
			#PRINT(buffer.components[0])
			for ix in xrange( 0, F.num_total_cells[0] ):
				perturbattion_value = buffer.components[0][ix]
				#PRINT(perturbattion_value)
				F.distribution_function_data[ix,0,:] *= perturbattion_value
				#F.scalar_mult_spatial_cell( perturbattion_value, ix )
	
	def apply(self, E_field, spatial_axis_global, spatial_axis__local_in_global, global_cell_offset, buffah=None, phase_me = 0.0, sine=False):
		
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
			
			#PRINT("IN DRIVER..... %f" % box_size #spatial_axis_global.values.shape[0])
			#PRINT(spatial_axis_global.values.shape[0])
			#PRINT(dx)
			
			
			if buffah != None:
				data = buffah.components[0]
			else:
				data = E_field.components[0]
			
			indic = []
			#for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			#	global_cell_index = float( global_cell_offset[0] + n )
			#	data[ E_field.num_boundry_cells[0][0] + n] = self.amp * np.sin(discrete_factor*global_cell_index + phase_cell_offset + phase_me )
			
			num = len(spatial_axis__local_in_global.values) + E_field.num_boundry_cells[0][0] + E_field.num_boundry_cells[0][1]
			sincos = np.cos
			if sine:
				sincos = np.sin
				#PRINT("USE SINE@@")
			xs1 = []
			for n in xrange(0, num):
				global_cell_index = float( global_cell_offset[0] + n - E_field.num_boundry_cells[0][0])
				data[n] = self.amp * sincos(2.0*np.pi*float( self.num_wavelengths)/(box_size +dx )*global_cell_index*dx + phase_me) # + phase_cell_offset + phase_me )
				
				xs1.append( global_cell_index * dx )
				
				#PRINT("DRIVAVVV: DIVIDING BY K WHICH IS: %e" % (2.0*np.pi*float( self.num_wavelengths)/(box_size +dx )))
				#global_cell_index = float( global_cell_offset[0] + n - E_field.num_boundry_cells[0][0])
				##indic.append( global_cell_index )
				#data[n] = self.amp * np.cos(discrete_factor*global_cell_index) # + phase_cell_offset + phase_me )
						
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
		
		#PRINT("CCCCCCUNT")
		for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
			data[ E_field.num_boundry_cells[0][0] + n] = self.amp * np.sin(self.wavenumber*x_val + self.phase )
		self.v_phase = 1.0/self.wavenumber
		#PRINT("if wp=1, v_phase = wp/k, vphase=%f" % ( self.v_phase ))
	
class GaussianProfile:
	def __init__(self, ptx, stdd, center, constant_vtx_offset ):
		self.ptx = ptx
		self.stdd = stdd
		self.center = center
		self.constant_vtx_offset = constant_vtx_offset
	
	# This is based off the imporved Gaussian rpfile Ben and I worked on for the simualtions that Ben did with Rosmus and Chappman
	#
	# Added by ADAM: A quick change to output temperatures closer to the relativistic proper velocity we use in Osiris.
	# 				Mototvation: 	In the input deck, we specify temperatures via 'thermal momenta' (pth)
	#								(in our units, pth has the same value as vth since mass of the electron is set to 1.0)
	#								The code (with the exception of the collision operator at this point) is relativistic,
	#								so the velocities and momenta are all 'proper'. So the input deck values we specify are
	#								all 'proper' momenta/velocities.. and it is these 'proper' quantites that we distrubute
	#								according to a Maxwellian distrubution when initalizing tempuratures. The original temperture
	#								output diagnostic computed the moments of the non-relativistic velocity^2... 
	#								I added this daignositic to compute the moments of the 'proper' moemtum/velcoity then
	#								after having such moments, convert them into the non-proper velocity and then into electronVolts.
	#								This is closer to the concept of temperature we mean in OSHUN.. BUT this is an open question
	#								with no 'right' answer.. only a better answer.. becasue techincally we are using the non-relvatistic
	#								Maxwell-Blotzzman distrucition to represent temperatures in a relaitivc paramterizeed system...
	#								Technically we should be using a Junter but this has a bunch of other problems and the correction
	#								is genrally mild.
	#
	#	The semantics are that you set a Gassian width this.stdd and a center in space with self.center.. and you set the termperature at the spatial point
	#		at center of the Gaussian by seeting the paramter.ptx.. Let's say self.ptx = 6000ev.. Then you passin constant_vtx_offset = 2000ev for the
	#		ambianent temperature.. The means: At thte center of the Gaussian, the Terperature is ~= (6000ev+2000v).. then the Gaussian makles a profile
	#		jsut as if the ambaient temperature = 0... except that It's sitting a a 2000ev pedestal (and thus lowest temperature in system initially is the ambaient)
	def apply(self, temperature_grid, spatial_axis_global, mpi_axes, species=None):
		spatial_axis__local_in_global = mpi_axes[0].axis__local
		if( constant_vtx_offset< 0.0):
			constant_vtx_offset = 0.0;
		alpha = 1.0/(sqrt(2.0)*self.stdd)
		# temperature_grid is initalized to al 1.0
		for (n,x) in enumerate(spatial_axis__local_in_global.values):
			tmp = alpha * ( x -  self.center )
			gaussian_part__vtx2 = (self.ptx*self.ptx) * math.exp(-tmp*tmp);
			vtx2_total = gaussian_part__vtx2 + constant_vtx_offset*constant_vtx_offset;
			
			ptx2_total = vtx2_total / (1.0 - vtx2_total);
			ptx2_total = sqrt( ptx2_total );
			temperature_grid[n + species.F.num_boundry_cells[0][0] ] *= ptx2_total
			
	
	def apply_old(self, temperature_grid, spatial_axis_global, spatial_axis__local_in_global, pt_amb, species):
		# temperature_grid is initalized to al 1.0
		for (n,x) in enumerate(spatial_axis__local_in_global.values):
			arg = (x - self.center) / self.stdd
			val = self.ptx*self.ptx*math.exp(-1.0*arg*arg)
			#PRINT(temperature_grid[n])
			temperature_grid[n + species.F.num_boundry_cells[0][0] ] *= val
		pass
		
	

# TODO: Removes these after asummer when Andrewis finsihed.. but not until then so as no to detrub his input decks.		
# class SpeciesConfig:
# 	def __init__(self):
		
# 		self.name = "No name"
# 		self.charge = -1.0
# 		self.mass = 1.0
# 		self.temperature_profile = None
# 		self.min_allowed_temperature = .0001

# class Config:
	
# 	def __init__(self, num_cells__global, system_size_min__global, system_size_max__global, num_species):
				
# 		self.dim = 1
# 		self.num_spatial_dims = 1
# 		self.num_cells__global 			= np.array( num_cells__global , dtype='int32')
# 		self.system_size_min__global 	= np.array(system_size_min__global)
# 		self.system_size_max__global 	= np.array(system_size_max__global)
# 		self.system_size__global 		= self.system_size_max__global - self.system_size_min__global
		
		
# 		#if len(num_cells__global) == 1:
# 		self.total_num_cells__global 	= self.num_cells__global[0]
# 		self.cell_size 					= np.divide(self.system_size__global, self.num_cells__global[0])
# 		self.cell_volume 					= self.cell_size[0]
		
# 		self.num_l 					= -1
# 		self.num_m 					= 1
# 		self.num_species 			= num_species
		
# 		self.dt 						= 0.0
# 		self.t_end 					= 0.0
		
# 		self.E_num_guard_cells = np.array((3,2))
# 		self.B_num_guard_cells = np.array((3,2))
# 		self.J_num_guard_cells = np.array((3,2))
		
# 		self.species_config = [None] * num_species
		

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
		#	most calculations should make use of this... any explicit knowledge of gloabl positions should
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
		self.yess = 0
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
		PRINT("p0p1_sq: %f" % p0p1_sq)
		PRINT("inv_mp0p1_sq: %f" % inv_mp0p1_sq)
		PRINT("g_r: %f" % g_r)
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
		import warnings
		# more numpy freidnly version...
		operator_data = _operator.data
		#PRINT("CUNT")
		#PRINT( E_field.components[0][:, None])
		
		# TODO: generalize this...
		with warnings.catch_warnings(record=True) as w:
			operator_data *= E_field.components[0][:, None]  		# the ,None makes a column vector
			operator_data *= scalar
			target_distribution_function.distribution_function_data[:, target_l, :] += operator_data
			if len(w) > 0:
				PRINT(target_l)
				PRINT( E_field.components[0][:, None])
				PRINT("Warning in E-field advance! (Usually about numerical exception conditions)")
				PRINT(w)
				mpi_info.abort()
				exit(-1)
		
		"""
		# original version...
		num_abs_p_cells = self.momentum_axis.num_points
		operator_data = _operator.data
		
		for x in xrange(0, target_distribution_function.num_total_cells[0]):
			momentum_projection = target_distribution_function.cells[x].harmonics[sh_linear_index].momentum_projection
			value_Ex = E_field.read(x, 0) * self.species_charge
			multiplier = value_Ex * scalar
			operator_data[x,:] *= multiplier
			np.add( momentum_projection[:], operator_data[x,:], momentum_projection[:])
		pass
		"""
	
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
		
		#LOOK LATER: Added Adam 11/2013:
		if species.num_l < 3:
			return
		
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	2, Yh, self.A1[1], Ein, self.G)
		
		# m = 0, 1 < l < l0
		for l in xrange(2, species.num_l-1):
			self.MakeGH(l,Yin)
			# Ex *= self.A2[l,0]; Yh.SH(l-1,0) += H.mxy_matrix(Ex)
			self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(l-1, Yh, self.A2[l], Ein, self.H)
			# Ex *= self.A1[l,0]; Yh.SH(l+1,0) += G.mxy_matrix(Ex)
			# check this fix then delte the following commented out line.
			#self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(l+1, Yh, self.A1[1], Ein, self.G)
			self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(l+1, Yh, self.A1[l], Ein, self.G)
		
		#m = 0,  l = l0-1
		self.MakeGH(species.num_l-1, Yin)
		# Ex *= self.A2[num_l-1,0]; Yh.SH(l-2,0) += H.mxy_matrix(Ex)
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(species.num_l-2, Yh, self.A2[species.num_l-1], Ein, self.H)
		
		return
	def Implicit_Ex1d_calc(self, Yin, Yh, species, state, CFG, Ein=None):

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
		#LOOK LATER: Added Adam 11/2013:
		if species.num_l < 3:
			return
		# Ex *= A1(1,0); Yh.SH(2,0) += G.mxy_matrix(Ex);
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	2, Yh, self.A1[1], Ein, self.G)
		
		# m = 0, l = 2
		self.MakeGH(2,Yin)
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	1, Yh, self.A2[2], Ein, self.H)		

		# m = 0, l = 3
		self.MakeGH(3,Yin)
		self.state_harmonic_plus_Ex_x_scalar_mxy_matrix(	2, Yh, self.A2[3], Ein, self.H)		

		return
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
		
		#PRINT("-------")
		#PRINT("----%d \n %s" %( mpi_info.rank, str(jayx) ))
		
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
		#	PRINT(u_field.components[0])
		#	PRINT(u_field.components[0].shape)
			
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
		
		#PRINT("adam1")
		#PRINT("------------------------")
		#PRINT(np.max( self.G.data ))
		#exit(-1)
		Yh.add_from_harmonic_centric_view(0, self.G )
		
		
		for l in xrange(1, species.num_l):
			self.MakeGH( l, Yin )
			#PRINT(np.max( self.G.data ))
			#PRINT(np.max( self.H.data ))
			
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
			PRINT("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~AVECT THIS MUTHA FUCKA~~~~~~~~~~~~~~~~~~~")
		
		Yin.copy_into_harmonic_centric_view(0,self.fd1)
		self.fd1.Dx()
		
		if debug:
			PRINT("TOP m=0")
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
					PRINT("In L loop L=%d before 1st derivitive" % ( l,))
					self.dump_hcv( self.fd1, "")
			
			self.fd1.Dx()
			if debug:
				PRINT("In L loop L=%d" % ( l,))
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
		#PRINT(STATE.use_external_E_field)
		if STATE.use_external_E_field:
			#PRINT("------")
			#PRINT(STATE.E_field.components[0])
			self.E0.add( STATE.E_field__external)
			pass
			
		
		if self.debug_me():
			PRINT("Step1 start: Y0"); DEBUG(YYY=Y0, L=LL, M=MM);
		
		# step 1
		# Yh = F(Y0)
		F.sub_time_step = 0
		F(Y0, Yh, spec, state, CFG, Ein=self.E0, Eout=self.Eh)
		if self.debug_me():
			PRINT("Step1 end rigth after: Yh"); DEBUG(YYY=Yh, L=LL, M=MM);
			PRINT("Step1 end rigth after: Y0"); DEBUG(YYY=Y0, L=LL, M=MM);
			PRINT("Step1 end rigth after: Eh"); DEBUG_E(self.Eh);
			
		# Y0 = Y0 + h*Yh
		Yh.scalar_mult(self.h);											self.Eh.scalar_mult(self.h);	
		Y0.add(Yh);															self.E0.add(self.Eh);
		
		if self.debug_me():
			PRINT("Step2 start: Yh"); DEBUG(YYY=Yh, L=LL, M=MM);
			PRINT("Step2 start: Y0"); DEBUG(YYY=Y0, L=LL, M=MM);
			PRINT("Step2 start: Eh"); DEBUG_E(self.Eh);
        
		#step 2
		# Yh = F(Y0)
		F.sub_time_step = 1
		F(Y0, Yh, spec, state, CFG, Ein=self.E0, Eout=self.Eh);
		# Y0 = 1/4 * ( 3*Y + (Y0 + h*Yh) )
		Yh.scalar_mult(self.h);											self.Eh.scalar_mult(self.h);	
		Y0.add(Yh);															self.E0.add(self.Eh);
		Y0.madd(3.0, Y);													self.E0.madd(3.0, E)
		Y0.scalar_mult(0.25);											self.E0.scalar_mult(0.25);
		
		if self.debug_me():
			PRINT("Step3		start: Y0"); DEBUG(YYY=Y0, L=LL, M=MM);
			PRINT("Step3		start: E0");  DEBUG_E(self.E0);
			
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
			PRINT("AT END OF RK3, DIST is (1,0):"); DEBUG(YYY=Y, L=LL, M=MM);
			PRINT("AT END OF RK3, Ex is:"); DEBUG_E(E);
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
			#PRINT("%d,= %d" % (l, threshold))
			if(threshold > 0):
				self.trim_list.append( (l, 0, Yh.linear_index_for_harmonic( l), threshold) )
		
		if self.do_filter:
			#if mpi_info.rank == 0:
			#	PRINT("Trimming filter will occur!!!")
			pass
	
	def do_db(self):
		return False
		if self.time_step == 0 and self.sub_time_step == 0:
			return True
		return False

	def __call__(self, Y0, Yh, species, state, CFG, Ein=None, Eout=None):
		if self.do_db():
			PRINT("\t\t****************START EVAL Timestep %d,%d" % (self.time_step, self.sub_time_step))
		
		Yh.scalar_mult(0.0)
		Eout.scalar_mult(0.0)
		LL = 1
		MM = 0
		
		if self.do_db():
			PRINT("\t\tIN EVAL: Yin In at start!!"); DEBUG(YYY=Y0, L=LL, M=MM, indent=True);
			PRINT("\t\tIN EVAL:  Ein IN at start!!"); DEBUG_E(Ein, indent=True);
		
		self.e_on_dist.calc								(Y0, Yh, species, state, CFG, Ein = Ein)
		
		if self.do_db():
			PRINT("\t\tIN EVAL: Yh AFTER E effect!!"); DEBUG(YYY=Yh, L=LL, M=MM, indent=True);
		
		"""
		if DaWaveDrivah != None:
			if STATE.time >= (DaWaveDrivah.rise_time + DaWaveDrivah.sustain):
				self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
			else:
				if mpi_info.rank == 0:
					PRINT("Field solver skipped.")
				
		else:
			self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
		"""
		#if DaWaveDrivah.driver_on == False:
		self.current_solver.calc						(Y0, Yh, species, state, CFG, Eout = Eout)
		
		debug_avection = False
		if self.do_db():
			PRINT("\t\tIN EVAL: Yh EX AFTERJ!!"); DEBUG_E(Eout, indent=True);
			debug_avection = False
		
		
		self.spatial_advection_solver.calc			(Y0, Yh, species, state, CFG, debug = debug_avection)
		
		#self.moving_frame.calc			(Y0, Yh, species, state, CFG, debug = debug_avection)
		
		if self.do_db():
			PRINT("\t\tIN EVAL: Yh AFTER Avectection!!"); DEBUG(YYY=Yh, L=LL, M=MM, indent=True);
		
		if self.do_filter:
			for item in self.trim_list:
				idx = item[2]
				level_to_trim = item[3]
				for ix in xrange(0, Yh.num_total_cells[0]):
					Yh.cells[ix].harmonics[idx].momentum_projection[:level_to_trim] = 0.0
				#PRINT("zeros %d cells for %d,%d" % (  level_to_trim, item[0], item[1] ))
			pass
		
		
		
		
		if self.do_db():
			PRINT("\t\t****************END EVAL Timestep %d,%d" % (self.time_step, self.sub_time_step))
		
		#if self.time_step == 1 and self.sub_time_step == 2:
		#	DEBUG_BAIL()
		


plotme = None
plotme1 = None
moments = []
DaWaveDrivah = None
update_boundries = None


def main_loop():
	global DaWaveDrivah
	global update_boundries
	
	# setup the calculation..
	# TODO: move this, of course.
	time = 0.0
	n = 0
	STATE.n = n
	STATE.time = time
	dt = STATE.config.dt;
	
	callback_mode = False
	
	#setup the simulation notes output file....
	# This is used for debugging.. i will leave this fragment in case needed sometime.
	# if mpi_info.rank == 0:
	#	f = open('!sim_info.txt','w')
	#	f.write("\n\n")
	#	f.close()
	
	DaWaveDrivah = None
	if True:
		if STATE.config.use_wave_driver:
			#DaWaveDrivah = WaveDrivah(STATE.config.dt, STATE.species[0], 1.0, num_cells_wave_travels_in_dt=1.0)
			DaWaveDrivah = STATE.config.wave_driver
			DaWaveDrivah.init( STATE.species[0].species_config.temp_profile.val, dt )
	#except:
	#	pass
	
	
	# reconcile the boundry cells for the T=0 evaluatation.
	update_boundries = UpdateBoundryCells(STATE, STATE.config)
	if mpi_info.num_nodes == 1:
		# in the mirror case, the mpi code is idealy suited to also do the single node calculation.
		if STATE.config.boundryTypeEMF == STATE.config.MIRROR:
			mpi_update.update(STATE, STATE.config)	
		else:
			update_boundries.calc(STATE, STATE.config)
	else:
		mpi_update.update(STATE, STATE.config)
	moments[0].notify_simulation_time(n, time)

	
	# TODO: split this off into an event module.
	if STATE.config.E_field_profile != None and True and DaWaveDrivah == None:
		# calculate the E field profile, if needed.
		scratchField = STATE.E_field.clone()
		STATE.config.E_field_profile.apply( scratchField.components[0], STATE.spatial_axes__global[0],  mpi_info.mpi_axes)
		scratchField.Dx(0)
		scratchField.components[0] /= STATE.species[0].species_config.charge

		# For now, the perturbation gets applied to species[0]..
		# TODO: make this configfielurable.
		initialDensityProfile = moments[0].get_density(force_update=True)
		print "Inital denisty profile:  MIN: %e MAX: %e" %( np.min(initialDensityProfile),np.max(initialDensityProfile))
		STATE.species[0].F.applyDensityPerturbation( scratchField.components[0], presentDensityProfile=initialDensityProfile )

		print "This is the perturbed Ditribution function!"
		initialDensityProfile = moments[0].get_density(force_update=True)
		print "Perturbed denisty profile:  MIN: %e MAX: %e" %( np.min(initialDensityProfile),np.max(initialDensityProfile))
		print initialDensityProfile
		print initialDensityProfile.shape

		# now copy the initial field over the the Efield for the entire simulation.
		STATE.config.E_field_profile.apply( STATE.E_field.components[0], STATE.spatial_axes__global[0],  mpi_info.mpi_axes)

		# santity test:
		#	TODO: move to a unit test.
		# 	density/E (for sinusodial) should give k/q
		e_max = np.max( STATE.E_field.components[0][ (STATE.E_field.num_boundry_cells[0][0]):(-1*STATE.E_field.num_boundry_cells[0][1]) ] )
		density_max =  np.max(scratchField.components[0][ (STATE.E_field.num_boundry_cells[0][0]):(-1*STATE.E_field.num_boundry_cells[0][1]) ] )
		k_min =  2*np.pi/( STATE.spatial_axes__global[0].max - STATE.spatial_axes__global[0].min)

		print "Density max: %e  e max: %e  ratio n/e %e  k_min %e (ratio n/e)/k_min should be ~integer/charge %e" % ( density_max, e_max, (density_max/e_max) , k_min, (density_max/e_max)/k_min)

	if DaWaveDrivah != None:
		if DaWaveDrivah.choose_dt_from_driver:
			dt = DaWaveDrivah.suggested_dt
	
	PRINT("Timestep size: %f" % dt)
	
	t_end = STATE.config.t_end;
	num_species = STATE.config.num_species;
	
	# TODO: this needs to be confugrable you moron.
	solva = RK3(dt)
	
	eval_functions = []
	fokker_planks = []
	
	debug_fokker_plank_Fs = []
	debug_fokker_planks = []
	
	
	# take care of output intervals..
	STATE.out_diags = Oshun1dEvents.OutputDiagnostics ( STATE )
	STATE.output_diags = Oshun1dEvents.DiagDoer()
	STATE.simple_timer = Oshun1dEvents.SimpleTimer()
	STATE.simple_perf_timer = Oshun1dEvents.GeneralPerfMonitor(STATE, STATE.simple_timer  )
	
	
	for i in range(0, num_species):
		if STATE.config.implicitESolver:
			eval_functions.append( RK2ForImplicit(STATE.species[i], STATE, STATE.config,dt) )		
		else:
			eval_functions.append(F_explicit(STATE.species[i], STATE, STATE.config) )

		fokker_planks.append( FokkerPlankExplicit	(STATE.species[i], STATE, STATE.config, use_c_version = USE_C_VERSION_GLOBAL) )			

	# write out basic simulation data...
	STATE.metadata.dt = dt
	STATE.metadata.output_interval__time_steps = STATE.num_timesteps_between_outputs
	STATE.metadata.set_dirty()
	
	tt_crap_vuf = None; simple_hhh_buffer = None; simple_hhh_buffer__index = 0; 
	
	initial_density = None
	initial_dist = None
	initial_dist_slice1 = None
	initial_dist_slice2 = None
	
	simple_E_buffer = None
	simple_E_buffer__index = 0
	
	output_interval = STATE.num_timesteps_between_outputs
	
	# create a record of the simulation configuration (mainly for records/providance and for postprocessing)
	# 	Before the main loop starts.
	if mpi_info.rank==0:
		try:
			pickle.dump( STATE.config, open( GLOBAL.OUTPUT_DIR("oshun_metadata__cfg.p","metadata"), "wb" ) )
		except:
			pass

	while( time < t_end):
		STATE.simple_timer.start_timer("main_loop_timing")
		OshunEventQueue.execute_event("main_loop_tick_start", locals(), STATE)


		# time the diagnostic step (todo: make a decorator)
		#TIME_IO_START = time_utils.time()
		#STATE.simple_timer.start_timer("diags_timing")
		
		# do the main diagnostic output
		STATE.output_diags.do_it_for_this_timestep( STATE, STATE.n, STATE.time, moments)
		
		#STATE.simple_timer.end_timer("diags_timing")
		#tima = STATE.simple_timer.getTimer("diags_timing")
		#PRINT("Daigs: %f  total: %f " % ( tima["interval"][-1],tima["accumlated_time"] ) )

		#TIMER__CHECKING_FOR_GRADUAL_SLOWDOWNS___TIMER_FOR_IO_ACCUMULATOR += ( time_utils.time() - TIME_IO_START )
		if STATE.config.implicitESolver:
			# Do the implicit calculations.
			for i in range(0, num_species):
				eval_functions[i].advance_p1(STATE.species[i], STATE, STATE.config)

				# need some communications.
				# but not for Elextrostatic, since there is no B field.
				
				# do the e-e collsions for f00 (explicit)
				fokker_planks[i].calc_f00(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
				# do the main implicit E solve
				eval_functions[i].advanceE(STATE.species[i], STATE, STATE.config, fokker_planks[i], n, time)
				# do last part of Explicit advacne
				eval_functions[i].advance_p2(STATE.species[i], STATE, STATE.config)				

				# do the implicit part of the collsion operator.
				if STATE.out_diags.collsions_enabled:
					fokker_planks[i].calc_flm(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)

		else:
			# Use the explicit operator.
			for i in range(0, num_species):
				#PRINT("\t step %d for: %s (collisions enabled: %s but NOT f00 (EXPLICT)))" % ( n, STATE.species[i].species_config.name, str(collsions_enabled) )			)
				#PRINT("\t step %d for: %s (collisions enabled: %s but NOT flm))" % ( n, STATE.species[i].species_config.name, str(collsions_enabled) ))
				PRINT("\t step %d for: %s (collisions enabled: %s))" % ( n, STATE.species[i].species_config.name, str(STATE.out_diags.collsions_enabled) ))			

				
				solva.step( eval_functions[i], STATE.species[i], STATE, STATE.config )

				if STATE.out_diags.collsions_enabled:
					fokker_planks[i].calc_f00(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
					fokker_planks[i].calc_flm(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
				if(n%10==0):
					fokker_planks[i].printStatus()

		
		# wave driver
		if DaWaveDrivah != None:
			DaWaveDrivah.update(time, n, dt)
		
		# update boundry cells...
		# and update the super moment calculator.
		# TODO: unify this. Aloow differnet boundries for different objects.
		for i in range(0, num_species):
			if mpi_info.num_nodes == 1:
				# in the mirror case, the mpi code is idealy suited to also do the single node calculation.
				if STATE.config.boundryTypeEMF == STATE.config.MIRROR:
					mpi_update.update(STATE, STATE.config)	
				else:
					update_boundries.calc(STATE, STATE.config)
			else:
				STATE.simple_timer.start_timer("mpi_boundry_cells")
				mpi_update.update(STATE, STATE.config)
				STATE.simple_timer.end_timer("mpi_boundry_cells")
			moments[i].notify_simulation_time(n, time)
			
		
		#if mpi_info.rank == 0:
		#TIMINGS.print_results()
		PRINT("done with step %d (time=%f ) " % (n, time) )
		
		# Disabled for now..
		# TODO: re-enable
		if STATE.metadata.is_dirty():
			PRINT("SAVING METADATA")
			STATE.metadata.save()
		
		# Disabled for now
		# TODO: Unify all the timers and re-enable this.
		# The idea for this timer was to look for long-term slow downs (for memoery leaks, for exmaple) as the simulation runs.
		#STATE.simple_perf_timer.take_measurement(STATE, STATE.n)
		
		
		time += dt
		n += 1
		STATE.n = n
		STATE.time = time
		
		if USE_C_VERSION_GLOBAL:
			tick_cpp_sim_clock( STATE.n, time, dt)
		
		STATE.simple_timer.end_timer("main_loop_timing")
		OshunEventQueue.execute_event("main_loop_tick_end", locals(), STATE)
	pass

	tima = STATE.simple_timer.getTimer("main_loop_timing")
	PRINT("")
	PRINT("main_loop_timing: %f  total: %f " % ( tima["interval"][-1],tima["accumlated_time"] ) )
	
	if mpi_info.rank == 0:
		TIMINGS.print_results()
		PRINT("It is finished...")

STATE 	= None
KEY		= None


"""
	Key routine called to setup simulation.. this routine assumes that all configuration
		is finished. This routine then setup up the various proper structures for a given configuration.
		(e.g. it actually setps up a temperature is specified.)
	TODO: Setup will become an 'ibject' with each setup step a distinct member functions
		This will allow any particualr setup state to be proxied/overidden.

	The two optional named parameters specify how the simulation will be configured. One (but not both)
		of these parameters must be specified.
		 CONFIG_PROFILE: set this to a Configuration object
		 CONFIG_CALLBACK: set this to a function to be called. This function must return a Configuration object.

"""	
def setup(CONFIG_PROFILE = None, CONFIG_CALLBACK=None):
	
	global STATE
	global filemode
	global moments
	global DaWaveDrivah

	if CONFIG_PROFILE != None:
		CFG = CONFIG_PROFILE
	elif CONFIG_CALLBACK != None:
		CFG = CONFIG_CALLBACK()
	else:
		raise Exception("OSHUN failed in the setup(..) method because neither a Configuration object (CONFIG_PROFILE) or callback (CONFIG_CALLBACK) was passed in.") 
	

	STATE = DaState( CFG )
	CFG = STATE.config
	GLOBAL.STATE = STATE
	
	STATE.num_timesteps_between_outputs = CFG.num_timesteps_between_outputs
	
	#............ now use CFG  to make simulation objects...........
	
	# Setup Momentum |p| data... (this axis is not chopped up by MPI)
	abs_p_axis = Axis(CFG.abs_p_min, num_points=CFG.abs_p_num_cells, max=CFG.abs_p_max)
	
	# Configure the settings about boundry cells and solvas.
	# 		We require the number of spatial boundry cells (on each side) that 
	#		is the same as the of the solva with max order.
	# First make the numpy array
	num_boundry_cells = np.zeros( (3,2), dtype=np.int32 )
	# now fill it for 1d
	# TODO: loft this operation to a central function becasue its very important and needs to be very explicit.
	num_boundry_cells[0][0] = CFG.explicit_solver_order
	num_boundry_cells[0][1] = CFG.explicit_solver_order
	
	# now let's setup MPI sooner rather then later.
	# the main interface with MPI is via MPIAxis objects that chunk up global spatial axes
	#		into the appropiate bits or all the nodes.
	mpi_info.init()

	# report to the user telling the observed MPI confugration.

	PRINT("----------------------------------")
	if mpi_info.num_nodes == 1:
		PRINT("Configuring OSHUN to run in serial mode on one CPU.")
	else:
		PRINT("Configuring OSHUN to run in Parallel across %d MPI processes." % (mpi_info.num_nodes,) )
	PRINT("----------------------------------")

	mpi_info.mpi_axes.append( AxisMpi(STATE.spatial_axes__global[0], num_boundry_cells, mpi_info) )
	# also, tell MPI the type of boundries we have.
	if( CFG.boundryTypeEMF == CFG.MIRROR ):
		mpi_update.setBoundryType(isMirror=True)
	else:
		mpi_update.setBoundryType(isPerodic=True)

	# now use this MPI axis to fill in the State object about the node-related spatial divisions.
	# the STATE.spatial_axes will cop the mpi_axis's 'axis__local' object... this axis
	# starts has extent cells = [0..num_cells_on_node], x = [0.0, last_x_position_on_node]
	STATE.spatial_axes[0] = Axis( axis=mpi_info.mpi_axes[0].axis__local )
	
	# create the Field Data structures.
	STATE.E_field = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)
	STATE.B_field = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)

	# Keep track of the implied extra charge density added to simultion ereinitally to enforce
	#	charge neutrality (for example, most simulations specify just 1 species - the elctrons. 
	#	Protons are automatically assumed by the code in order to enforce charge neutrality and for the collsion operator.
	#	'assumedNetrualizingChargeDensity' keeps track of the exact assumed denosty profile, whcih is usefull in som perturbation calcuations/diagnostics.
	STATE.assumedNetrualizingChargeDensity = STATE.E_field.clone()
	STATE.assumedNetrualizingChargeDensity.zeroOut()

	
	# initialize any inital E-field, if requested.
	"""
	if False and CFG.Field_profile != None:
		CFG.E_field_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset)
		
		vth = CFG.species_config[0].temp_profile.val
		w_wave = np.sqrt( 1 + 3*CFG.E_field_profile.wavenumber*CFG.E_field_profile.wavenumber*vth*vth)
		#PRINT("freq (BG): %e " % w_wave)
		#PRINT("vphase: %e" % (w_wave/ CFG.E_field_profile.wavenumber))
		#PRINT("K l_debye: %f" % (CFG.E_field_profile.wavenumber*vth))
		CFG.E_field_profile.v_phase = (w_wave/ CFG.E_field_profile.wavenumber)
	"""
	
	for i in range(0, CFG.num_species):
		# each species has it's own abs_p_axis to describe the |p| dimension...
		#	for now, this axis is setup with identical paramters for each of the species.
		#	This can be changed later if needed (but is frought with other issues in the collision operator so needs to be done carefully)
		species_config = CFG.species_config[i]
		temp_species = Species(	STATE.spatial_axes[0].num_points, CFG.num_l,
										momentum_axis = abs_p_axis, 
										species_config = species_config,
										num_boundry_cells = num_boundry_cells)
		temp_species.index = i
		STATE.species[i] = temp_species
		
		# cheap hack. TODO: take out. LATER: I forgot why this is a hack.. Looks OK. Not great. 
		#	But sometimes shit go down like this in the ghetto. oh well.. look again someday.
		if i == 0:
			mpi_info.mpi_axes[0].init_species_data( 	STATE.species[i],
																	STATE.E_field,			
																	STATE.species[i].F.cells[0].num_harmonics)
		
		
		# each species also has it's own J.. for now.. maybe not needed but no real harm.
		temp_species.J = Field1d(STATE.spatial_axes[0].num_points, 1, num_boundry_cells)
		
		
		# The denisty/temperature initiak profiles are setup in two stages.
		#	First a simple array of numbers that corresponds to the denisty (temperature) at each spatial is created
		#		from the input simulation configuration (arrays 'density_profile' and 'temperature_grid')
		#	The second step: Given the array from step 1, alter the Distribution function ( ie.  set appropiate Sphereical Harmonic projection values)
		#		so that the requeted denisty (temperature) profile is setut/represented by the inital Distrbution function.
		
		# set density profile...
		# set inital to 1.0 all places... (num_total_cells means all cells including guard cells)
		density_profile = np.ones( temp_species.F.num_total_cells[0] )
		
		# now apply any the density profile specified in the simupation input deck.
		#	These profiles make thier changes by multiplying the input array (or in rare cases outright replacing it)
		species_config.density_profile.apply_density( density_profile,STATE.spatial_axes__global[0], mpi_info.mpi_axes, species=temp_species)
		
		# create a charge density profile that is the same as this species, but with opposite charge. This is the charge density implitly added
		#	to the simulation when charge neutrality is enforced. Add this offsetting denisty in a global accumlator since we only care about
		#	overall net summed over all species. 
		scratchData = np.copy  ( density_profile )
		scratchData *= ( -1.0*species_config.zeta)
		STATE.assumedNetrualizingChargeDensity.components[0] +=  scratchData


		#PRINT(dir(mpi_info.mpi_axes[0].axis__local))
		#PRINT(mpi_info.mpi_axes[0].axis__local.min)
		#PRINT(mpi_info.mpi_axes[0].axis__local.num_points)
		#exit(-1)

		#----------- Covert density and temperature profiles into the inital
		#------------	values for the |p| axes in the shperical harmonics.
		# As a baseline, set all |p| cells to 1.0 at all spatial locations 
		#		for the Y(0,0)  harmonics. ( Y(0,0) is the lowest order/l=0/S/isotropic harmonic)
		for x in xrange(0, temp_species.F.num_total_cells[0] ):
			temp_species.F.cells[x].harmonics[0].momentum_projection[:] = 1.0  
		
		# Then multiplicativly modulate the Y(0,0)'s ( |p| entries to reflect the requested density profile
		# 		(we ignore the guard cells)
		#	Note: 'num_plain_cells' is number of cells EXCLUDING guard cells. So this density calcualtion
		#		is a usefull fragment of code for writing to cell/guard cells. The pattern of these lines is reapeated many times.
		for ix in xrange(0, temp_species.F.num_plain_cells[0] ):
			density_multiplier = density_profile[num_boundry_cells[0][0]+ ix]
			temp_species.F.cells[num_boundry_cells[0][0]+ ix].harmonics[0].momentum_projection[:] *= density_multiplier
		

		# -------------- Initalize the species's spatial temperature profile ---------------
		# 		Result is in temperature_grid
		#		The 'temperature_grid' is a grid of temperatures give the temperature at ach girb location. The temperatures
		#			are specified as "relativistic proper momentum/velocities".. Especially note: the number in the temperature map
		#			are "p_thermal" NOT "p_thremal^2" (we explictly take the square root in the Profile creation functions if needed)
		
		# setup the temperature profile.
		temperature_grid = np.ones( temp_species.F.num_total_cells[0] )
		min_allowed_temperature = species_config.min_allowed_temperature
		#species_config.temp_profile.apply( temperature_grid, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local, temp_amb, temp_species)
		species_config.temp_profile.apply( temperature_grid, STATE.spatial_axes__global[0], mpi_info.mpi_axes, species=temp_species)

		# Ensure that a global (non-negitive minimum for temperature is enforced at all gird locations.. A <=0 temperature is catatrophic
		for ix in xrange(0, temp_species.F.num_total_cells[0]):
				if( temperature_grid[ix] < (min_allowed_temperature) ):
					# This is Mikhail's version..
					#temperature_grid[ix] += temp_amb*temp_amb
					# this is Adam's version which is better. Its a reflection that English is Adam native language
					#	and this is what ambient means.. not that other Greek defintion.
					temperature_grid[ix] = min_allowed_temperature
		
		#-----------poke this temperature map into the Y0,0 harmonics as 
		#		The temperature will be used to set a gaussian in p^2 along the |p| axis for the Y(0,0) harmonic.
		#		This Gaussian mulitplies the resuls from the density calculation (as it sould be.. as the integral 
		#			under the velocity ditbutions should give the density at a given spatial cell.)
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

		if( mpi_info.is_I_root()):
			#TODO: gather the lists from all MPI nodes.
			# check to make sure that the maxwell-boltzmann distribution is qdepqautly represented in velocity space
			flagged_spatial_points = []
			flagged_spatial_points_highest_p = []

			for ix in xrange(0, temp_species.F.num_plain_cells[0]):
				highest_p_value = temp_species.F.cells[num_boundry_cells[0][0]+ ix].harmonics[sh_linear_index].momentum_projection[-1]
				if( np.abs(highest_p_value) > 1e-4 ):
					flagged_spatial_points.append( ix )
					flagged_spatial_points_highest_p.append(  highest_p_value )
			if(len(flagged_spatial_points) > 0):
				PRINT("")
				PRINT("")
				PRINT("---------------- WARNING ------------------")
				PRINT("The Inital velocity distribution may not be adequately represented in the selected region of velocty space in %d (out of %d) spatial cells." % (len(flagged_spatial_points), temp_species.F.num_plain_cells[0]))
				PRINT("It is recommended that you increase the 'abs_p_max' input parameter from %.4e to a higher value." % (STATE.config.abs_p_max, ))
				PRINT("The following is the value of inital velocity distribution at the highest defined velocity for the suspicious points:")
				PRINT("----------------------")
				for i in xrange(0,len(flagged_spatial_points) ):
					PRINT("Spatial cell %d: %.4e" %(flagged_spatial_points[i], flagged_spatial_points_highest_p[i] ))
				PRINT("")
				PRINT("This warning occurs for any points where the inital velocity distribution is >1e-4 at the higest defined velocity.")
				PRINT("By default, the simulation will now abort. To force the simulation to run, set the input parameter: do_not_abort_on_abs_p_max_warning = True")
				PRINT("-------------------------")

				# TODO: make configuration object return 'None' for undefined paramters rather then throwing exceptions.
				try:
					if STATE.config.do_not_abort_on_abs_p_max_warning == False :
						exit(-1)
				except:
					# this is executed if 'do_not_abort_on_abs_p_max_warning' is not defined.. in that case, do the default which is to exit.
					exit(-1)

	# setup the Output moments generator... (one per species)
	for i in range(0, CFG.num_species):
		moments.append( MomentCalculationsState(STATE.species[i], STATE) )
		
	# make directories...
	if mpi_info.rank==0:
		try:
			os.makedirs(GLOBAL.OUTPUT_DIR("./numpy"))
		except:
			pass
		try:
			os.makedirs(GLOBAL.OUTPUT_DIR("./png"))
		except:
			pass
		try:
			os.makedirs(GLOBAL.OUTPUT_DIR("./png/dist"))
		except:
			pass
		try:
			os.makedirs(GLOBAL.OUTPUT_DIR("./dist"))
		except:
			pass
		try:
			os.makedirs(GLOBAL.OUTPUT_DIR("./misc"))
		except:
			pass


	if USE_C_VERSION_GLOBAL:
		# C stuff...
		general_setup__c( CFG.num_species, STATE.E_field.components[0], STATE.E_field.num_boundry_cells, 
								STATE.spatial_axes[0])

		for i in xrange(0, CFG.num_species):
			general_species_setup( 	i, STATE.species[i].num_l, STATE.species[i].num_m, 
											STATE.species[i].F.momentum_axis,
											STATE.species[i].F.num_total_cells, STATE.species[i].F.num_plain_cells, 
											STATE.species[i].F.num_boundry_cells[i], 
											STATE.species[i].F.distribution_function_data )
			for ix in xrange(0, STATE.species[i].F.num_total_cells[0]):
				for ish in xrange(0, STATE.species[i].F.cells[ix].num_harmonics):
					transfer_species_F_buffers( i, ish, ix, STATE.species[i].F.cells[ix].harmonics[ish].momentum_projection )
		pass
	pass

	# ------------------------------- Done setting up inital temp/density profiles ------------------------
	# I will leave the following code fragment in.. it can be usefull for debugging (and for MPI debugging)
	#	It will make all MPI nodes block here until the all nodes finish initialization.
	#PRINT("RANK %d is wating at block" % mpi_info.rank)
	#mpi_info.comm.Barrier()	
	#PRINT("RANK %d Done" % mpi_info.rank)
	
	# this 'true' means that setup completed sucessfully.
	return True

# setupNotNeeded: skip running setup. Defaults to False so that setup is run (once)
#						This is the correct behavior is almost all cases. But if you do you own setup
#						elsewhere the, pass setupNotNeeded=True to prevent setup from being called (and overwriting any of your custom setup) 
def run_sim(config=None, config_callback=None, profile=False, profileIncludeFilters=None, profileExcludeFilters=None, setupNotNeeded = False, resourceProfile=False):
	go_on = True
	if setupNotNeeded == False:
		go_on = setup(CONFIG_PROFILE = config, CONFIG_CALLBACK=config_callback)
	
	# if profiling is requested, then add the profiling event to the event queue.
	if profile:
		import oshun_modules.profile_main_loop as profiler
		profiler.register(profileIncludeFilters=profileIncludeFilters,profileExcludeFilters=profileExcludeFilters)
	if resourceProfile:
		import oshun_modules.utils.ResourceUsage as resourceProfiler
		resourceProfiler.register()

	
	if go_on:
		OshunEventQueue.execute_event("before_main_loop", locals(), STATE)
		try:
			main_loop()
		finally:
			OshunEventQueue.execute_event("after_main_loop", locals(), STATE)
			pass
	return



# This bit specifies what is executed if this python script is run directly from the command line 
#	(i.e. python this-file.py). It is not exeuted if tihis file is being run becasue it was included with an 'import'
#	statement from another .py file.
# TODO: Add some command line parameters here so this can actually work.
if __name__ == "__main__":
	# this will always generate an error to the user.. its not setup yet.
	run_sim(config=None, config_callback=None )

	

	