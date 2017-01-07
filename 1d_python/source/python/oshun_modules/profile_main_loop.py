# Simple module that profiles the main loop using the 'statprof'
#	statisical profiler. Note this this profiler requires the installation
#	of statprof (pip install statprof) and does NOT work under windows.

# for now, we explictly reisgter this conponent.. we can make it cooler later
import OshunEventQueue
import platform
from  OshunMpi import mpi_info
from OshunGlobalState import PRINT as PRINT

import re

if platform.system().lower() != 'windows':
	try:
		import statprof
	except:
		PRINT("Note: Profiling disabled. Unable to load the 'statprof' Python module. Install 'statprof'.")
		statprof = None
else:
	PRINT("Note: Profiling disabled when running under Windows because 'statprof' does not work.")
	PRINT("Note: Profiling disabled. Unable to load the 'statprof' Python module.")
	statprof = None

filtersForFiltersInclude = None
filtersForFiltersExclude = None

def register(profileIncludeFilters=None,profileExcludeFilters=None):
	global filtersForFiltersInclude
	global filtersForFiltersExclude
	filtersForFiltersInclude = profileIncludeFilters
	filtersForFiltersExclude = profileExcludeFilters

	if(statprof == None ):
		return

	start_event = OshunEventQueue.Event( "main_loop_profiler", profile_main_loop_start)
	end_event =   OshunEventQueue.Event( "main_loop_profiler", profile_main_loop_end)
	OshunEventQueue.add_event("before_main_loop", start_event)
	OshunEventQueue.add_event("after_main_loop", end_event)

def profile_main_loop_start( local_state, simulation_state):

	if platform.system().lower() != 'windows':
		if(statprof == None ):
			PRINT("Note: Profiling disabled. Unable to load the 'statprof' Python module. Install 'statprof'.")
			return
		statprof.start()
	else:
		PRINT("Note: Profiling disabled when running under Windows.")

def profile_main_loop_end( local_state, simulation_state):

	if platform.system().lower() != 'windows':
		if(statprof == None ):
			PRINT("Note: Profiling disabled. Unable to load the 'statprof' Python module. Install 'statprof'.")
			return
		statprof.stop()
		filter_statprof (statprof)
	else:
		PRINT("Note: Profiling disabled when running under Windows.")

# helper function to make the output of statprof nicer...
# Filter is a list of filenames. The output is filtered to include only these files.
#	If None (or not provided, then all filename are included)
def filter_statprof( statprof ):
	# only report the profiling for the Root node.. that's almost always what is needed.
	if mpi_info.is_I_root == False:
		return

	if platform.system().lower() != 'windows':
		if(statprof == None ):
			PRINT("Note: Profiling disabled. Unable to load the 'statprof' Python module. Install 'statprof'.")
			return
	else:
		PRINT("Note: Profiling disabled when running under Windows because 'statprof' does not work.")
		return
	pass


	PRINT("Final Profiling Info")
	PRINT("-----------------------------------")
	
	PRINT("\t Include filters:")
	filtersInclude = None
	if filtersForFiltersInclude != None and len(filtersForFiltersInclude) > 0:
		filtersInclude = []
		for filter in filtersForFiltersInclude:
			PRINT("\t\t %s" %(str(filter),))
			filtersInclude.append( str(filter) )
	else:
		PRINT("\t\t None")

	PRINT("\t Exclude filters:")
	filtersExclude = None
	if filtersForFiltersExclude != None and len(filtersForFiltersExclude) > 0:
		filtersExclude = []
		for filter in filtersForFiltersExclude:
			PRINT("\t\t %s" %(str(filter),))
			filtersExclude.append( str(filter) )
	else:
		PRINT("\t\t None")

	PRINT(""); PRINT("");

	outp = open('profile_data.txt', 'w')
	statprof.display(fp = outp)
	outp.close()
	
	inp = open('profile_data.txt', 'r')
	lines = inp.readlines()
	inp.close()
	
	error_lines = []
	profile_info = []
	line_number = -1


	for line in lines:
		line_number += 1
		# filter (include filter)
		if filtersInclude != None:
			found = False
			for name in filtersInclude:
				if line.find (name) > -1:
					found = True
					break
			if found == False:
				continue
		pass

		# filter (exclude filter)
		if filtersExclude != None:
			found = False
			for name in filtersExclude:
				if line.find (name) != -1:
					found = True
			if found == True:
				continue
		pass

		if line_number > 1:
			try:
				numbers = re.findall(r"([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", line)
				percent_time = float(numbers[0][0])
				cum_seconds = float(numbers[1][0])
				self_seconds = float(numbers[2][0])
					
				text = re.search(r"([a-zA-Z][\.\w\:]+)", line)
				info = text.group()
				#sort(key=lambda tup: tup[1])
				profile_info.append ( ( percent_time, cum_seconds, self_seconds, info) )
			except:
				error_lines.append( line )
	

	profile_info.sort( key = lambda sme: sme[0] )
	
	for line in profile_info:
		PRINT("%f \t %f	\t %f \t %s" % ( line[0],line[1],line[2],line[3]))
	if len( error_lines ) > 0:
		PRINT("---------------------")
		PRINT("Error lines!")
		PRINT("")
		for line in error_lines:
			PRINT("\t %s" % line)
	PRINT("")
	PRINT("Raw output in file 'profile_data.txt' and CSV data in 'profile.csv'")
	PRINT("")
	
	outp = open('profile_data.csv', 'w')
	for line in profile_info:
		outp.write("%f,%f,%f,%s\n" % ( line[0],line[1],line[2],line[3]))
	outp.close()
