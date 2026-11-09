# Bootstrap
# This the prime Bootstrapping Code.. This is the start of the execution path
# there are lots of try-n-catches in this so we can report erors to the user..


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#	Imports for the Bootstrapping
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#print "----Oshun is starting the Bootstrap Process"
import traceback

# TODO: configure so that this is not imported if not needed (minor)
# TODO: Only import this on the root node (to stop file contention during import process)
# This command needs to be executed BEFORE any other matplotlib commands in the current execution environment.
#	If a specific non-interactive (like 'agg') output method is not chosen, a default one will be slected that often invokes X-Windows.
import matplotlib
matplotlib.use('Agg')

# This surpresses python warnings... numpy spits out a bunch so disable for now.
# TODO: decide if warnings should be surpressed
import warnings
warnings.filterwarnings("ignore")

import os
import glob
import pkgutil
import os.path
import sys

def reportError( message ):
	e = sys.exc_info()[0]
	print "\t %s" % message
	print "--------------------------------------------------"
	e = sys.exc_info()[0]
	print "%s" % str(e)
	traceback.print_exc()
	print ""
	exit(-1)



import OshunGlobalState as GLOBAL
from utils.logging import *

try:
	from  OshunMpi import mpi_info
except:
	reportError( "Error in initial import of the OshunMpi module." )

try:
	import shared_libraries as sl
except:
	reportError( "Error in initial import of the Oshun shared_libraries module." )

try:
	import OshunGlobalState as GLOBAL
except:
	reportError( "Error in initial import of the OshunGlobalState module." )

try:
	from utils.logging import *
except:
	reportError( "Error in initial import of the Oshun utils.logging module." )



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#	Functions used in the Bootstrapping
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def setupSharedLibraries( ):
	
	sharedLibrariesInsideOshun = sl.sharedLibraries
	print sl.sharedLibraries
	print "sl.sharedLibraries"

	destinationPathForSharedLibraries = pathWhereScriptIsRun + "/.shared_libraries"
	destinationPathForSharedLibraries = os.path.normpath( destinationPathForSharedLibraries )

	# Tell the global state where the shared libraries are located. 
	GLOBAL.set_SHARED_LIBRARY_PATH_BASE_PATH( destinationPathForSharedLibraries )

	if( mpi_info.is_I_root() ):
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

		#print sharedLibrariesInsideOshun
		# loop over shared library names.. copy each to destination.
		for libFilename in sharedLibrariesInsideOshun:
			# check to see that the filename is non-trival
			if(libFilename.strip() == ""):
				continue
			
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



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#	Main Code Executed in the Bootstrapping
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
pathWhereScriptIsRun = os.getcwd()

# Do any initazation that should only be done for 1 (the root) MPI node
#	and have all other nodes wait until the root is done with its task.
setupSharedLibraries()
mpi_info.barrierWait()


#print "----Bootstrapping completed"
