import numpy as np
import ctypes as np__c
from numpy.ctypeslib import ndpointer
import os.path
import os

import OshunGlobalState as GLOBAL
from OshunGlobalState import PRINT

#print "----------------------------------------------"
#print os.getcwd()
#exit()
#try:
#directoryOfScript = os.path.dirname( __file__ )
#sharedLibrary = directoryOfScript + "/" + 'oshun_c_routines_for_speedup'
#sharedLibrary = os.path.normpath( sharedLibrary )
_oshun_lib = np.ctypeslib.load_library('oshun_c_routines_for_speedup', GLOBAL.GET_SHARED_LIBRARY_PATH('oshun_c_routines_for_speedup') )
#_oshun_lib = np.ctypeslib.load_library('oshun_c_routines_for_speedup', '.')

#print "OK"
#exit(-1)
#-----general setup---------
#print np__c.POINTER(np__c.c_double)
#print np__c.c_int
#print _oshun_lib.general_setup.argtypes

#utility functions.. there is probably a better way to do this.. but this works for now.
_oshun_lib.get_pointer_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_pointer_size_in_bytes.argtypes = []
_oshun_lib.get_int_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_int_size_in_bytes.argtypes = []
_oshun_lib.get_long_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_long_size_in_bytes.argtypes = []
_oshun_lib.get_longlong_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_longlong_size_in_bytes.argtypes = []
_oshun_lib.get_float_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_float_size_in_bytes.argtypes = []
_oshun_lib.get_double_size_in_bytes.restype = np__c.c_uint
_oshun_lib.get_double_size_in_bytes.argtypes = []


def getPointerSizeInBytes():
	return _oshun_lib.get_pointer_size_in_bytes()
def getIntSizeInBytes():
	return _oshun_lib.get_int_size_in_bytes()
def getLongSizeInBytes():
	return _oshun_lib.get_long_size_in_bytes()
def getLongLongSizeInBytes():
	return _oshun_lib.get_longlong_size_in_bytes()
def getFloatSizeInBytes():
	return _oshun_lib.get_float_size_in_bytes()
def getDoubleSizeInBytes():
	return _oshun_lib.get_double_size_in_bytes()

def printCConfigurationHeader():
	PRINT("")
	PRINT("-------C/cpu confiugration info--------" )
	PRINT("Pointer  size (in bytes): %d " % ( getPointerSizeInBytes () ) )
	PRINT("Int      size (in bytes): %d " % ( getIntSizeInBytes () ) )
	PRINT("Long     size (in bytes): %d " % ( getLongSizeInBytes () ) )
	PRINT("LongLong size (in bytes): %d " % ( getLongLongSizeInBytes () ) )
	PRINT("Float    size (in bytes): %d " % ( getFloatSizeInBytes () ) )
	PRINT("Double   size (in bytes): %d " % ( getDoubleSizeInBytes () ) )
	PRINT("---------------------------------------")
	PRINT("")

# for 64 bit windows, I have to change following call definition to use (3rd arg: was np__c.POINTER(np__c.c_long) --> np__c.POINTER(np__c.c_longlong) )
_oshun_lib.general_setup.restype = None
_oshun_lib.general_setup.argtypes= [	np__c.c_int, np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_longlong), np__c.c_int, np__c.POINTER(np__c.c_int),
												np__c.POINTER(np__c.c_double), np__c.c_int]
												

_oshun_lib.general_species_setup.restype = None
_oshun_lib.general_species_setup.argtypes = [	np__c.c_int, np__c.c_int, np__c.c_int, 
																np__c.POINTER(np__c.c_double), np__c.c_int,
																np__c.POINTER(np__c.c_int),  np__c.POINTER(np__c.c_int),
																np__c.POINTER(np__c.c_int), np__c.POINTER(np__c.c_double) ] 
															
_oshun_lib.transfer_species_F_buffers.restype = None
_oshun_lib.transfer_species_F_buffers.argtypes=  [np__c.c_int, np__c.c_int, np__c.c_int, np__c.POINTER(np__c.c_double)] 
							
							
#----- Effect of E ------
_oshun_lib.init_effect_of_E_on_distribution.restype = None
_oshun_lib.init_effect_of_E_on_distribution.argtypes = [	np__c.c_int, np__c.c_int, np__c.c_int,
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double),
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double),
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double),
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double),
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double),
																			np__c.POINTER(np__c.c_double),np__c.POINTER(np__c.c_double)]
#---------------------																	
#----- collisons.....
#---------------------

# callback defintion for using numpy matrix solving libraries.
# returns int (for no reason)
# expects double* right_hand_side_data, int current_x_index, int current_l_index, double dt
matrix_solver_callback_template = np__c.CFUNCTYPE( np__c.c_int, np__c.POINTER(np__c.c_double), np__c.c_int, np__c.c_int, np__c.c_double)
debug_callback_template = np__c.CFUNCTYPE( np__c.c_int, np__c.c_int )

#  for 64 bit windows had tp change 2nd arg  (2nd arg: was np__c.POINTER(np__c.c_long) --> np__c.POINTER(np__c.c_longlong) )
_oshun_lib.init__fokker_plank__implicit.restype = None
_oshun_lib.init__fokker_plank__implicit.argtypes = [	np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_longlong),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double),
																		#np__c.POINTER(np__c.c_double),
																		#np__c.POINTER(np__c.c_double),
																		np__c.c_double, np__c.c_double, np__c.c_double,
																		np__c.c_int, np__c.c_int, np__c.c_int,
																		np__c.POINTER(np__c.c_double),  np__c.POINTER(np__c.c_double),
																		matrix_solver_callback_template, 
																		matrix_solver_callback_template,
																		debug_callback_template]

_oshun_lib.fokker_plank__implicit__calc_flm.restype = None
_oshun_lib.fokker_plank__implicit__calc_flm.argtypes = [ np__c.c_int,  np__c.c_double,  np__c.c_double]

_oshun_lib.fokker_plank__implicit__advance_1.restype = None
_oshun_lib.fokker_plank__implicit__advance_1.argtypes = [ np__c.c_int,  np__c.c_double,  np__c.c_double]

_oshun_lib.fokker_plank__implicit__advance_flm.restype = None
_oshun_lib.fokker_plank__implicit__advance_flm.argtypes = [ np__c.c_int,  np__c.c_double,  np__c.c_double]


_oshun_lib.fokker_plank__implicit__advance.restype = None
_oshun_lib.fokker_plank__implicit__advance.argtypes = [ 	np__c.POINTER(np__c.c_double), 
																			np__c.c_int, np__c.c_double, np__c.c_double, np__c.c_int]

_oshun_lib.fokker_plank__implicit__reset_coeff.___logEE_ctype = np__c.c_double()
_oshun_lib.fokker_plank__implicit__reset_coeff.restype = None
_oshun_lib.fokker_plank__implicit__reset_coeff.argtypes = [ np__c.POINTER(np__c.c_double), np__c.c_double]



#  for 64 bit windows had tp change 2nd arg  (2nd arg: was np__c.POINTER(np__c.c_long) --> np__c.POINTER(np__c.c_longlong) )
_oshun_lib.init__fokker_plank__explicit.restype = None
_oshun_lib.init__fokker_plank__explicit.argtypes = [	np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_longlong),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.c_double, np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double),
																		np__c.POINTER(np__c.c_double), np__c.c_double,
																		np__c.c_int, np__c.c_double,  np__c.c_int, np__c.c_int]

_oshun_lib.fokker_plank_explicit__calc_f00.restype = None
_oshun_lib.fokker_plank_explicit__calc_f00.argtypes = [ np__c.c_int,  np__c.c_double,  np__c.c_double]

_oshun_lib.fokker_plank_explicit_slope_calc.restype = None
_oshun_lib.fokker_plank_explicit_slope_calc.argtypes = [ np__c.POINTER(np__c.c_double), np__c.POINTER(np__c.c_double) ]


# 
#-------- Setup.................					
def general_setup__c( num_species, E_field_data, E_file_boundry_cells, spatial_axis__x):
	_oshun_lib.general_setup( 	num_species, 
										E_field_data.ctypes.data_as(np__c.POINTER(np__c.c_double)), E_field_data.ctypes.shape, E_field_data.ndim,
										E_file_boundry_cells.ctypes.data_as(np__c.POINTER(np__c.c_int)),
										spatial_axis__x.values.ctypes.data_as(np__c.POINTER(np__c.c_double)), spatial_axis__x.num_points)
	
#	_oshun_lib.general_setup( 	num_species, 
#										E_field_data.ctypes.data_as(np__c.POINTER(np__c.c_double)), E_field_data.ctypes.shape, E_field_data.ndim,
#										E_file_boundry_cells.ctypes.data_as(np__c.POINTER(np__c.c_int)),
#										spatial_axis__x.values.ctypes.data_as(np__c.POINTER(np__c.c_double)), spatial_axis__x.num_points)
def general_species_setup( species_idx,  num_l, num_m, momentum_axis, F_num_total_cells, num_plain_cells, num_boundry_cells, F_data):
	_oshun_lib.general_species_setup(	species_idx,num_l, num_m, 
													momentum_axis.values.ctypes.data_as(np__c.POINTER(np__c.c_double)), momentum_axis.num_points,
													F_num_total_cells.ctypes.data_as(np__c.POINTER(np__c.c_int)),
													num_plain_cells.ctypes.data_as(np__c.POINTER(np__c.c_int)),
													num_boundry_cells.ctypes.data_as(np__c.POINTER(np__c.c_int)),
													F_data.ctypes.data_as(np__c.POINTER(np__c.c_double)) )
def transfer_species_F_buffers( species_idx, harmonic_idx, ix, data):
	_oshun_lib.transfer_species_F_buffers( species_idx, harmonic_idx, ix, data.ctypes.data_as(np__c.POINTER(np__c.c_double)) )


#----- Effect of E ------
def init_effect_of_E_on_distribution( species_index, num_l, num_m, A1, A2, B1, B2, C1, C2, C3, C4, Hp0, G, H, invpr):
	_oshun_lib.init_effect_of_E_on_distribution(		species_index, num_l, num_m,
																	A1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	A2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	B1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	B2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	C1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	C2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	C3.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	C4.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	Hp0.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	G.data.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	H.data.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																	invpr.values.ctypes.data_as(np__c.POINTER(np__c.c_double)) )

def effect_of_E_on_distribution_calc():
	pass 

	
#----------- collisons.....
def init__fokker_plank_implicit__cpp(	vr, vr3, 
													df0, ddf0, 
													Scattering_Term, Alpha, AlphaTri, 
													J1m, TriI1, TriI2, 
													IvDnDm1, IvDnDp1, Ivsq2Dn, 
													I0, kpre, zeta, f00_factor,
													is_tridiagonal, if_implicit1D, species_index,
													_LOGee, _ZLOGei,
													tridiagonal_matrix_solver, general_matrix_sovler, debug_callback,
													solverStats):
	_oshun_lib.init__fokker_plank__implicit( 	vr.ctypes.data_as(np__c.POINTER(np__c.c_double)), vr.ctypes.shape,
															vr3.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															df0.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															ddf0.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															Scattering_Term.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															Alpha.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															AlphaTri.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															J1m.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															TriI1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															TriI2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															IvDnDm1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															IvDnDp1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															Ivsq2Dn.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															I0.ctypes.data_as(np__c.POINTER(np__c.c_double)),
															kpre, zeta, f00_factor,
															is_tridiagonal, if_implicit1D, species_index,
															np__c.byref( np__c.c_double ( _LOGee ))  ,
															np__c.pointer( np__c.c_double ( _ZLOGei )),
															tridiagonal_matrix_solver, general_matrix_sovler, debug_callback,
															solverStats.ctypes.data_as(np__c.POINTER(np__c.c_double))
															);

#															np__c.POINTER( np__c.c_double ( _LOGee ))  ,
#															np__c.pointer( np__c.c_double ( _ZLOGei )),
#															tridiagonal_matrix_solver, general_matrix_sovler);
															#_LOGee.ctypes.data_as(np__c.POINTER(np__c.c_double)), _ZLOGei.ctypes.data_as(np__c.POINTER(np__c.c_double)),)

def fokker_plank__calc_flm__cpp( Yin, species, state, CFG, n_step, time, dt ):	
#def fokker_plank__calc_flm__cpp( n_step,time, dt ):
	_oshun_lib.fokker_plank__implicit__calc_flm( n_step,time, dt );
def fokker_plank__advance_1__cpp( Yin, species, state, CFG, n_step, time, dt ):	
	_oshun_lib.fokker_plank__implicit__advance_1( n_step,time, dt );
def fokker_plank__advance_flm__cpp( Yin, species, state, CFG, n_step, time, dt ):	
	_oshun_lib.fokker_plank__implicit__advance_flm( n_step,time, dt );

# note: 'implcit_FP_object' is assigned in the FokkerPlank __init__ functoion
# TODO: make this less random and magical
def fokker_plank__implicit__rest_coeff__cpp(fin, Dt) :
	# make a refence for passing back the logEE
	_oshun_lib.fokker_plank__implicit__reset_coeff.___logEE_ctype.value = -1.0
	_oshun_lib.fokker_plank__implicit__reset_coeff (  fin.ctypes.data_as(np__c.POINTER(np__c.c_double)), Dt, np__c.byref(shun_lib.fokker_plank__implicit__reset_coeff.___logEE_ctype ) )
	fokker_plank__implicit__rest_coeff__cpp.implcit_FP_object._LOGee = _oshun_lib.fokker_plank__implicit__reset_coeff.___logEE_ctype.value

def fokker_plank__implicit__advance__cpp(fin, LL, LogEE, dt, current_x_index):
	#if is_tridiagonal:
	_oshun_lib.fokker_plank__implicit__advance( 	fin.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																LL, LogEE, dt, current_x_index)
	#else:
	#	_oshun_lib.fokker_plank__implicit__advance( 	fin.ctypes.data_as(np__c.POINTER(np__c.c_double)),
	#																LL, LogEE, Dt, is_tridiagonal)
	
"""
# function will be called FROM c/c++ in order t a hunk of python code (a callback)
#		that will, as a first pass, use numpy Matrix Solving routines.
def fokker_plank__matrix_solver_callback_setup__cpp(matrix, ):

	# setup the callback... the 1st argument (a c_int  in this case) is the return type.... all others are the paramets deilvered to
	#		to calllback function upon invokation.
	matrix_solver_callback_template = np__c.CFUNCTYPE( np__c.c_int, POINTER(c_double) )
	
"""


def init__fokker_plank_explicit__cpp( vr, U4, U4m1, U2, U2m1, U1, U1m1, J1, I4, U3, Qn, Pn, I2, 
												G_constant_1, G_constant_vr3, G_constant_vr5, G_constant_vr5_2, G_constant_vr7,
												data_density_np, NB, c_kpre, num_subcycling_steps, species_index):
	_oshun_lib.init__fokker_plank__explicit(		vr.ctypes.data_as(np__c.POINTER(np__c.c_double)), vr.ctypes.shape,
																U4.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U4m1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U2m1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U1m1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																J1.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																I4.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																U3.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																Qn.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																Pn.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																I2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																G_constant_1,
																G_constant_vr3.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																G_constant_vr5.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																G_constant_vr5_2.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																G_constant_vr7.ctypes.data_as(np__c.POINTER(np__c.c_double)),
																data_density_np, NB, c_kpre, num_subcycling_steps, species_index)

def fokker_plank__calc_f00__cpp( Yin, species, state, CFG, n_step, time, dt ):	
#def fokker_plank__calc_f00__cpp( n_step,time, dt ):
	_oshun_lib.fokker_plank_explicit__calc_f00( n_step,time, dt );

def fokker_plank__explicit__advance__cpp(fin, fh):
	_oshun_lib.fokker_plank_explicit_slope_calc( fin.ctypes.data_as(np__c.POINTER(np__c.c_double)), 
																fh.ctypes.data_as(np__c.POINTER(np__c.c_double)) )


_oshun_lib.testa_cpp_passing.restype = None
_oshun_lib.testa_cpp_passing.argtypes = [ np__c.POINTER(np__c.c_double) ]

def testa_cpp( pointer_to_a_float):
	_oshun_lib.testa_cpp_passing(  pointer( np__c.c_double ( pointer_to_a_float ))  )
	

_oshun_lib.testb_cpp_passing.restype = None
_oshun_lib.testb_cpp_passing.argtypes = [ np__c.POINTER(np__c.c_double) ]
def testb_cpp( pointer_to_a_float):
	_oshun_lib.testa_cpp_passing(  pointer( np__c.c_double ( pointer_to_a_float ))  )


_oshun_lib.advance_cpp_sim_clock.restype = None
_oshun_lib.advance_cpp_sim_clock.argtypes = [np__c.c_int, np__c.c_double, np__c.c_double  ]	
def tick_cpp_sim_clock( time_step, time, dt):
	_oshun_lib.advance_cpp_sim_clock( time_step, time, dt )