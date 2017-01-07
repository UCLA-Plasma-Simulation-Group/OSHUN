import glob
import os
from os.path import *
import zipfile
import sys

# gets the directory of the shard_libraries directory..
sharedLibraries = []
directoryOfSharedLibaries = os.path.dirname( __file__ )

# Test if the module we are in is in a .zip file or a simple directory.
# TODO: There is probably a better way to this this. MAybe just have build process inject all the names.
#		I have worries about 1000s of processes trying to itterate though this.

zipRefererence = None
try:
	# look two levels up in the directory struture. Tey to load that path as a zip file.
	zipPath = os.path.split( __file__ )[0]
	zipPath = os.path.split( zipPath )[0]
	zipPath = os.path.split( zipPath )[0]
	zipRefererence = zipfile.ZipFile(zipPath, 'r')
except:
	pass

if( zipRefererence ):
	filesInZip = zipRefererence.namelist()
	for name in filesInZip:
		if( name.find("__init__") > 0 ):
			continue
		if( name.find("oshun_modules/shared_libraries") < 0 ):
			continue
		sharedLibraries.append( os.path.split(name)[1] )
	zipRefererence.close()
else:
	globPattern = directoryOfSharedLibaries + "/*"
	globPattern = normpath( globPattern )
	sharedLibraries = [x for x in glob.glob( globPattern ) if x.find("__init__") < 0] 
	sharedLibraries = [os.path.basename( x ) for x in sharedLibraries]