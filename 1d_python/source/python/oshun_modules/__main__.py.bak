# This file is used when the Oshun Module is called as 'main'
#	That is, when Oshun is invoked from the command line like 'python oshun'
#	This works if 'oshun' is a directory or is a zip file ('python oshun.zip')
#
#	We expect that an input deck be specified.
#		python oshun input_deck.py

#	If no input deck is specified then default to the name 'oshun-input.py'
import sys
import os
import imp
import traceback

def reportException( message ):
	e = sys.exc_info()[0]
	print "\t %s" % message
	print "--------------------------------------------------"
	e = sys.exc_info()[0]
	print "%s" % str(e)
	traceback.print_exc()
	print ""
	exit(-1)

def loadInputDeck(  inputDeckNameSpecified ):
	
	inputDeckName = inputDeckNameSpecified

	if( inputDeckNameSpecified == None):
		print ""
		print "No Input Deck was specified, so the default 'oshun-input.py' is being used."
		print ""
		inputDeckName = 'oshun-input.py'

	inputDeckName = os.path.normpath( inputDeckName )

	
	if( os.path.exists(inputDeckName ) == False ):
		newInputDeckName = None
		mod_name,file_ext = os.path.splitext(os.path.split(inputDeckName)[-1])
		# be a bit charitiable here. Try using the adding or removeing the .py extension"
		if file_ext.lower() == '.py':
			newInputDeckName = mod_name
		else:
			newInputDeckName = mod_name + ".py"

		if(  newInputDeckName != None and os.path.exists(newInputDeckName ) == True ):
			print "The Input Deck '%s' was not found. Using Input Deck '%s' instead." % (inputDeckName, newInputDeckName )
			inputDeckName = newInputDeckName
		else:
			print ""
			print "The Input Deck %s could not be found. Aborting." % ( inputDeckName, )
			print ""
			exit(-1)

	mod_name,file_ext = os.path.splitext(os.path.split(inputDeckName)[-1])
	py_mod = imp.load_source(mod_name, inputDeckName)
	
	mainFunction = None
	try:
		mainFunction = getattr(py_mod, mainFunction)
	except:
		pass


	if( mainFunction != None ):
		py_mod.mainFunction()


# TODO:
#	Come back to this. It works.. but the except clause is always thrown.. I think it should
#	not be thrown in the .zip file case but it is.. come back and understand why.
#	It works as-is, but I would like ot know for my own edification.

# try to find the oshun module. This will suceed if oshun is run from a .zip
#	file but will fail if run directly from the oshun_modules directory.
#	So if we cant find the module, we construct the parect directory of this __main__.py
#	and add it to the python path.. this will make things work when oshun is run from the directory.

try:
	import oshun_modules
except:
	parentDirectoryOfThisScript = os.path.split(__file__)[0] 
	try:
		parentDirectoryOfThisScript = os.path.split( os.path.split(__file__)[0] )[0]
	except:
		pass

	if( parentDirectoryOfThisScript != '' ):
		sys.path.insert(0,parentDirectoryOfThisScript  )
	pass

# now try to import again.. if it fails then bail
try:
	import oshun_modules
except:
	reportException( "\nFatal Error:\n" + "\t Unable to import the ohsun_modules python module (i.e. the main OSHUN code).\n")


# now get the filename of the input deck, import it in, and run it
inputDeckNameSpecified = None
if( len(sys.argv) > 1):
	inputDeckNameSpecified = sys.argv[1]

loadInputDeck( inputDeckNameSpecified )