import os.path
import os
from utils.logging import *
import numpy as np
import matplotlib.pyplot as pyplot
from  OshunMpi import mpi_info

USE_C_VERSION = True
SHARED_LIBRARY_PATH_BASE_PATH = None

OUTPUT_ROOT_DIR =  os.path.normpath('output')
OUTPUT_CATEGORIES = {
	'metadata': 'misc',
	"numpy": "numpy",
	"png": "png"
}

STATE = None


def PRINT( message ):
	if mpi_info.is_I_root():
		print message

def NP_SAVE( filename, data ):
	np.save( OUTPUT_DIR( filename, category='numpy'), data )

def SAVEFIG( filename ):
	pyplot.savefig( OUTPUT_DIR( filename, category='png') )

def OUTPUT_DIR( filename, category=None ):
	# make sure only the root node actually writes to disk
	if(mpi_info.is_I_root() == False):
		return

	if category in OUTPUT_CATEGORIES:
		category = OUTPUT_CATEGORIES[category]
	if category:
		filename = os.path.join(category, filename)

	filename = os.path.normpath( filename )
	return os.path.join(OUTPUT_ROOT_DIR, filename)

def GET_SHARED_LIBRARY_PATH( libname ):
	if( SHARED_LIBRARY_PATH_BASE_PATH ) == None:
		FatalError( "You attempted to load a sharad library but the proper initalization has not yet occured.")
	return SHARED_LIBRARY_PATH_BASE_PATH

def GET_SHARED_LIBRARY_FULLPATH( libname, errorIfNonexistant = True ):
	if( SHARED_LIBRARY_PATH_BASE_PATH ) == None:
		FatalError( "You attempted to load a sharad library but the proper initalization has not yet occured.")
	libFilename = os.path.normpath( SHARED_LIBRARY_PATH_BASE_PATH + "/" + libname )
	
	if( errorIfNonexistant ):
		if( os.path.exists( libFilename ) == False):
			FatalError( "The shared library %s cannot be found at %s" % ( libFilename,  SHARED_LIBRARY_PATH_BASE_PATH ))
		pass
	return libFilename

def set_SHARED_LIBRARY_PATH_BASE_PATH( value ):
	global SHARED_LIBRARY_PATH_BASE_PATH
	if( SHARED_LIBRARY_PATH_BASE_PATH != None ):
		StrongWarning( "The cirutcal configuration variable SHARED_LIBRARY_PATH_BASE_PATH was set more then once. This is fishy.")
	SHARED_LIBRARY_PATH_BASE_PATH = value
