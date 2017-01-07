"""
import os
import glob
import pkgutil
import os.path

import shared_libraries as sl
import OshunGlobalState as GLOBAL
from utils.logging import *
from  OshunMpi import mpi_info


def setupSharedLibraries( ):
	
	sharedLibrariesInsideOshun = sl.sharedLibraries

	destinationPathForSharedLibraries = pathWhereScriptIsRun + "/.shared_libraries"
	destinationPathForSharedLibraries = os.path.normpath( destinationPathForSharedLibraries )

	# Tell the global state where the shared libraries are located. 
	GLOBAL.set_SHARED_LIBRARY_PATH_BASE_PATH( destinationPathForSharedLibraries )

	# create the .shared_libraries directory
	try:
		os.makedirs( destinationPathForSharedLibraries )
	except:
		pass

	# check if it was actually created..
	if( os.path.exists( destinationPathForSharedLibraries ) == False ):
		FatalError( "Could not create the .shared_libraries folder in the current directory.")

	# delete any existing files...
	for filename in  glob.glob( os.path.normpath( destinationPathForSharedLibraries + "/*") ):
		os.remove(  os.path.normpath(filename))

	print sharedLibrariesInsideOshun
	# loop over shared library names.. copy each to destination.
	for libFilename in sharedLibrariesInsideOshun:
		# get the file contents from the oshun module (this cool function gets the data no matter if
		#	the module is just a directory or is stored in a .zip file etc..)
		libFile = pkgutil.get_data('oshun_modules', 'shared_libraries/' + libFilename )

		# copy the file over.
		destinationFilename = os.path.normpath( destinationPathForSharedLibraries + "/" + libFilename )
		f = open(destinationFilename, 'wb')
		f.write(  libFile )
		f.close()

	# check that the shared libraries were copied over OK
	for libFilename in sharedLibrariesInsideOshun:
		 # GET_SHARED_LIBRARY_PATH will throw a ftal error if the library does not exist.
		 GLOBAL.GET_SHARED_LIBRARY_PATH( libFilename )


pathWhereScriptIsRun = os.getcwd()

# Do any initazation that should only be done for 1 (the root) MPI node
#	and have all other nodes wait until the root is done with its task.
if( mpi_info.is_I_root() ):
	setupSharedLibraries()
mpi_info.barrierWait()

"""