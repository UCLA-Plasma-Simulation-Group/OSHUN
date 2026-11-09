import numpy as np
import os.path
import sys
#sys.path.append("/u/home/m/muffinlo/utils/trunk1d")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE\\run_utils")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE\\utils")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE\\postprocess")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE\\c_version")
#sys.path.append("C:\\Users\\lobby\\Desktop\\vlasov\\MERGE_c_build")
#sys.path.append("C:\\lab_code\\python_modules\\adam_python_modules\\os_tools\\parser")
#C:\Users\lobby\Desktop\vlasov\MERGE_c_build\spitzer_tests

#from vqunits import Formulary
sys.path.append("../")
#sys.path.append("../oshun_modules")
sys.path.append("../c_version")
import Oshun1d
from Oshun1d import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile, run_sim, WaveDrivah



CFG = Oshun1d.spitzer_setup()


# if this is the actual script run, then start Oshun... if this script is imported by another script,
# then that script needs to call run_sim
if __name__ == "__main__":
	run_sim( config = CFG)

# if a template, others can call this function...
def run(config = CFG):
	run_sim( config )

# helper to make user's like easier..
def process_paramter(val, type_is_int=False, type_is_str=False, type_is_float=False):
	pass

# these are very hack-ish temporary functions to help configging util I make a real system..
def set_basic_paramters(config=CFG, box_size = None, num_cells = None, num_species = None):
	
	system_size_min__global  	= config.system_size_min__global
	if box_size == None:
		system_size_max__global = config.system_size_max__global
	if num_species == None:
		num_species 				= config.num_species
	if num_cells == None:
		num_cells 					=  config.num_cells
	
	
	value = None
	try:
		value = box_size[0]
		temp =  float( box_size[0] )
		system_size_max__global = np.array( [ temp, 0.0, 0.0 ] )
	except:
		# not an array.. maybe a value..
		try:
			value = float( box_size )
			system_size_max__global = np.array( [ value, 0.0, 0.0 ] )
		except:
			print "A unknown data type was passed into 'set_basic_paramters' for the paramter 'box_size' "
			exit(-1)
	
	value = None
	try:
		value = num_cells[0]
		temp =  int( num_cells[0] )
		num_cells__global = np.array( [ temp, 0, 0 ] )
	except:
		# not an array.. maybe a value..
		try:
			value = int( num_cells )
			num_cells__global = np.array( [ value, 0, 0 ] )
		except:
			print "A unknown data type was passed into 'set_basic_paramters' for the paramter 'num_cells' "
			exit(-1)
	
	value = None
	try:
		value = num_species[0]
		temp =  int( num_species[0] )
		num_species = temp
	except:
		# not an array.. maybe a value..
		try:
			value = int( num_species )
			num_species = value
		except:
			print "A unknown data type was passed into 'set_basic_paramters' for the paramter 'num_species' "
			exit(-1)
	
	#print "IN BASC ISETUP TWHERE WE GOT"
	#print num_cells__global
	#print system_size_min__global
	#print system_size_max__global
	#print num_species
	
	new_cfg = Config( num_cells__global, system_size_min__global, system_size_max__global, num_species)
	new_cfg = ld_setup__standing_wave__3rd_order(cfg = new_cfg )
	
	#print new_cfg
	print CFG
	
	global CFG
	CFG = new_cfg
	
	#print CFG
	return CFG
	
