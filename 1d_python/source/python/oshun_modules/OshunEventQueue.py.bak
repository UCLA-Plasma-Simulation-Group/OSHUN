# This is a basic event queue for OSHUN
# TODO: 	switch to a 3rd party event system for more robust handeling and 
#			one that uses Python features like function decoration and intercepting.
#
#
# For now implement at loose binding system based on a hash where various modules can
#	register via strings. This can easily replaced later, if wanted.
#
# There may be for complicated configurations and dependencies needed like something called before somehting else..
#		but implement that when needed
#

event_queue = {}

class Event:
	def __init__(self, listerner_name, function_to_call, implicit_state=None, event_properties=None ):
		self.function_to_call = function_to_call
		self.listerner_name = listerner_name
		self.event_properties = event_properties
		self.implicit_state = implicit_state
		if( event_properties == None):
			self.event_properties = {}
		if( self.listerner_name == None ):
			raise Error("Error creating Event. The listerner name must be set.")
		if( self.function_to_call == None ):
			raise Error("Error creating Event. The function_to_call upon event must be set.")

def add_event(event_name, event):

	if event_name not in event_queue:
		event_list = []
		event_queue[event_name] = event_list
	else:
		event_list = event_queue[event_name]

	event_list.append( event )

def execute_event( event_name, local_variables, simulation_state):
	# TODO: use the dictionarydefault object in the collection package to always return []
	if event_name not in event_queue:
		return
	listeners = event_queue[event_name]
	for i in xrange(0, len( listeners )):
		listeners[i].function_to_call( local_variables, simulation_state)



# Simple testing routines..
def event_queue_unit_test():
	def dumbtest_printme( val1, val2 ):
		print val1
	add_event("Hi", Event("printme 1", dumbtest_printme ) )
	add_event("Hi", Event("printme 2", dumbtest_printme ) )
	add_event("Hi There", Event("printme 1", dumbtest_printme ) )

	execute_event("Hi", "Nose", "Ear")

#event_queue_unit_test()