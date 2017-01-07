# Stupid utility class to encapsulation Error handling
# TODO: maybe put real 3rd party logging solution here

def FatalError( message ):
	print ""
	print "A Fatal Error Occured. The error was:"
	print "---------------------------------------------"
	print message
	print ""
	exit(-1)

def StrongWarning( message ):
	print ""
	print "A Strong Warning was issued:"
	print "---------------------------------------------"
	print message
	print ""
