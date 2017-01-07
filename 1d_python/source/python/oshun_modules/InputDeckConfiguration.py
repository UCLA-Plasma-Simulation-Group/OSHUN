import numpy as np

# TODO: Removes these after asummer when Andrewis finsihed.. but not until then so as no to detrub his input decks.		
class SpeciesConfig:
	def __init__(self):
		self.MIRROR = 1
		self.PERODIC = 0

		self.name = "No name"
		self.charge = -1.0
		self.mass = 1.0
		self.temperature_profile = None
		self.min_allowed_temperature = .0001
		self.boundryType = self.PERODIC

class Config:
	
	def __init__(self, num_cells__global, system_size_min__global, system_size_max__global, num_species):
		self.MIRROR = 1
		self.PERODIC = 0

		self.dim = 1
		self.num_spatial_dims = 1
		self.num_cells__global 			= np.array( num_cells__global , dtype='int32')
		self.system_size_min__global 	= np.array(system_size_min__global)
		self.system_size_max__global 	= np.array(system_size_max__global)
		self.system_size__global 		= self.system_size_max__global - self.system_size_min__global
		self.boundryTypeEMF 			= self.PERODIC
		
		#if len(num_cells__global) == 1:
		self.total_num_cells__global 	= self.num_cells__global[0]
		self.cell_size 					= np.divide(self.system_size__global, self.num_cells__global[0])
		self.cell_volume 					= self.cell_size[0]
		
		self.num_l 					= -1
		self.num_m 					= 1
		self.num_species 			= num_species
		
		self.dt 						= 0.0
		self.t_end 					= 0.0
		
		self.E_num_guard_cells = np.array((3,2))
		self.B_num_guard_cells = np.array((3,2))
		self.J_num_guard_cells = np.array((3,2))
		
		self.species_config = [None] * num_species
