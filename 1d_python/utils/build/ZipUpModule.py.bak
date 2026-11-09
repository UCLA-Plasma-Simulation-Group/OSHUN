import shutil
import fnmatch
from glob import *
import os
import os.path
import sys
import zipfile

# example
# 		python build_utils\zipUpModule.py . oshun.zip oshun_modules oshun_modules
def zipUpModule( path, outputFilename, prefixInArchive = None, mainInitDir = None ):

	# erase all the .pyc files, if they exist.
	path = os.path.normpath( path )
	outputFilename = os.path.normpath( outputFilename )

	for root, dirnames, filenames in os.walk(path):
		for filename in fnmatch.filter(filenames, '*.pyc'):
	  		try:
	  			os.remove( os.path.join(root, filename) )
			except:
				pass
		pass

	# now zip up the files..
	if prefixInArchive:
		shutil.make_archive( outputFilename, 'zip', path, prefixInArchive)
	else:
		shutil.make_archive( outputFilename, 'zip', path)

	if mainInitDir != None:
		if( os.path.exists(  os.path.join(mainInitDir,"__main__.py")  ) == False):
			print "You requested that a __main__.py be used for this zip copied from the directory %s. But it could not be found" % (  os.path.join(mainInitDir,"__main__.py") , )
		addToZip( outputFilename+".zip", os.path.join(mainInitDir,"__main__.py"), destFilename = "__main__.py" )
	pass

def addToZip( zipFilename, filepathToAdd, destDirPrefix = None, destFilename = None):

	filepathToAdd = os.path.normpath( filepathToAdd )

	zipFilename = os.path.normpath( zipFilename )
	zipFile = zipfile.ZipFile(zipFilename, "a")
	if destDirPrefix:
		destDirPrefix = os.path.normpath(destDirPrefix)

	if( os.path.isdir(filepathToAdd)  ):
		for root, dirnames, filenames in os.walk(filepathToAdd):
			for filename in fnmatch.filter(filenames, '*.*'):
		  		try:
		  			outputFilename = filename
		  			if( destDirPrefix != None ):
		  				outputFilename = os.path.join(destDirPrefix, filename)

		  			zipFile.write( os.path.join(root, filename),   outputFilename )
				except:
					pass
			pass
		pass
	elif( os.path.isfile(filepathToAdd) ):
		outputFilename = filepathToAdd
		if( destDirPrefix != None and destFilename != None ):
			outputFilename = os.path.join( destDirPrefix, destFilename)
		elif( destDirPrefix != None ):
			outputFilename = os.path.join( destDirPrefix, filepathToAdd)
		elif( destFilename != None ):
			outputFilename = destFilename		
		
		zipFile.write( filepathToAdd, outputFilename)

	zipFile.close()


if __name__ == "__main__":
	if  len(sys.argv) < 3:
		print "Wrong number of arguments provided.. the format of the command is: python zipUpModule directory zip_file_name"
		exit(-1)

	prefix = None
	mainInitDir = None
	if( len( sys.argv) > 3):
		prefix = sys.argv[3]
	if( len( sys.argv) > 4):
		mainInitDir = sys.argv[4]

	zipUpModule( sys.argv[1], sys.argv[2], prefix, mainInitDir)
 