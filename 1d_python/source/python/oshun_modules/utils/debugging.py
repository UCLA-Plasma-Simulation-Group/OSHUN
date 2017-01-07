
import numpy as np

def dump_F_at_xlm(F, x_points, l, m, title=None):
	try:
		temp = len(x_points)
	except:
		x_points = [x_points]
	
	num_boundry_cells = F.num_boundry_cells
	sh_index = F.linear_index_for_harmonic(l,m)
	
	print ""
	if title != None:
		print title
		print ""
	for ip in xrange(0, F.momentum_axis.num_points):
		print "dump_F_at_xlm: %d: "%( ip,) ,
		for ix in x_points:
			print "\t%f"%( F.cells[num_boundry_cells[0][0]+ix].harmonics[sh_index].momentum_projection[ip] ),
		print ""
	pass	


def DEBUG_E(E, bail=False, indent=False):
	if mpi_info.rank == 0:
		limit = 50
		if E.num_plain_cells[0] < limit:
			limit = E.num_plain_cells[0]
		num_boundry_cells = E.num_boundry_cells
		
		for ix in xrange(20, 30):
			#print "%d: %e" %( ix, E.components[0][num_boundry_cells[0][0] + ix])
			if indent:
				print "\t\t",
			print "%d: %.8e" %( ix,  E.components[0][ix])
		if bail:
			exit(-1)

def DEBUG_BAIL():
	if mpi_info.rank == 0:
		exit(-1)
	pass
def DEBUG(YYY=None, L=0, M=0, bail=False, indent=False):
	if mpi_info.rank == 0:
		temp_species = STATE.species[0]
		if YYY == None:
			print "Setting DEBUG Y to defualt for species"
			YYY = temp_species.F
		
		index_ml = YYY.linear_index_for_harmonic(L,M)
		num_boundry_cells = YYY.num_boundry_cells
		#print "INDEX OF HARM:" + str(index_ml)
		x_points = [0,1,2,3,4,5,6]
		#print id(YYY)
		if indent:
			print "\t\t",
		
		for ix in x_points:
			print "\t%.9e" % STATE.spatial_axes[0].values[ix],
		print ""
		if indent:
			print "\t\t",	
		print "--------------------------------------------------------"
		
		#print "beam temp %f" % beam_temp
		#print "EXP arg factor is: %f" % exp_arg_coeff
		#print "FROT factor is: %f" % factor_coeff
		#for ip in xrange(0, temp_species.momentum_axis.num_points):
		
		for ip in xrange(20, 30):
			if indent:
				print "\t\t",
			
			print "%d: "%( ip,) ,
			for ix in x_points:
				#print "\t%e"%( YYY.cells[num_boundry_cells[0][0]+ix].harmonics[index_ml].momentum_projection[ip] ),
				if( ix == num_boundry_cells[0][0]):
					print "",
				print "\t%.8e"%( YYY.cells[ix].harmonics[index_ml].momentum_projection[ip] ),
			print ""
		if bail:
			exit(-1)

def DEBUG_HCV(HCV, L=0, M=0, bail=False, indent=False):
	if mpi_info.rank == 0:
		x_points = [0,1,2,3,4,5,6]
		if indent:
			print "\t\t",
		for ix in x_points:
			print "\t%.8e" % STATE.spatial_axes[0].values[ix],
		print ""
		temp_species = STATE.species[0]
		nump = temp_species.momentum_axis.num_points
		num_boundry_cells = temp_species.F.num_boundry_cells
		#for ip in xrange(0, nump):
		for ip in xrange(20, 30):
			if indent:
				print "\t\t",
			
			print "%d: "%( ip,) ,
			for ix in x_points:
				if( ix == num_boundry_cells[0][0]):
					print "",
				#print "\t%e" % ( HCV.data[num_boundry_cells[0][0]+ix,ip],) ,
				print "\t%.8e" % ( HCV.data[ix,ip],) ,
			print ""
		if bail:
			exit(-1)


"""
In the main loop there was code for debuggin the Fokker Plank

Near top of main_loop 3 real lines.. then commented out FP debuggin code...
for i in range(0, num_species):
	eval_functions.append( F_explicit(STATE.species[i], STATE, STATE.config) )
	
	fokker_planks.append( FokkerPlankExplicit	(STATE.species[i], STATE, STATE.config, use_c_version = USE_C_VERSION_GLOBAL) )		
	# this is for denugging fokker plank (repeating calcuations...)
	#debug_fokker_plank_F = STATE.species[i].F.clone_structure()
	#debug_fokker_plank_Fs.append( debug_fokker_plank_F )
	#debug_fokker_planks.append( FokkerPlankExplicit ( STATE.species[i],  STATE, STATE.config, use_c_version = False) )

After the line:
	if STATE.out_diags.collsions_enabled:
		
		# Debugging: do the calcuation using python for reference
		#debug_fokker_plank_Fs[i].copy( STATE.species[i].F )
		#debug_fokker_planks[i].calc_f00(debug_fokker_plank_Fs[i], STATE.species[i], STATE, STATE.config, n, time, dt)
		#debug_fokker_planks[i].calc_flm(debug_fokker_plank_Fs[i], STATE.species[i], STATE, STATE.config, n, time, dt)
		
		# then 2 linws are real code then commented out:
		fokker_planks[i].calc_f00(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
		fokker_planks[i].calc_flm(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
		
		
		#temp = STATE.species[i].F.distribution_function_data[0,0,0]
		#STATE.species[i].F.distribution_function_data[0,0,0] = 666.0
		#temp1 = debug_fokker_plank_Fs[i].distribution_function_data[0,0,0]
		#print "compring FP... first elemtn test:"
		#print "\t %e %e" % ( STATE.species[i].F.distribution_function_data[0,0,0], temp1)
		#STATE.species[i].F.distribution_function_data[0,0,0] = temp
		#
		#if not np.allclose(STATE.species[i].F.distribution_function_data, debug_fokker_plank_Fs[i].distribution_function_data):
		#	print STATE.species[i].F.distribution_function_data.shape
		#	print STATE.species[i].F.distribution_function_data.shape[2]
		#	print STATE.species[i].F.distribution_function_data.shape[1]
		#	print STATE.species[i].F.distribution_function_data.shape[0]
		#	for x in xrange(0, STATE.species[i].F.distribution_function_data.shape[0] ):
		#		for L in xrange(0, STATE.species[i].F.distribution_function_data.shape[1] ):
		#			for P in xrange(0, STATE.species[i].F.distribution_function_data.shape[2] ):
		#				data1 = STATE.species[i].F.distribution_function_data[x,L,P]
		#				data2 = debug_fokker_plank_Fs[i].distribution_function_data[x,L,P]
		#				delta = np.fabs( data1 - data2)
		#				if delta > 1e-5:
		#					print "FAIL!!!"
		#					print "%d, %d, %d %e, %e" % ( x, L, P, data1,data2)
		#					exit(-1)
		#	
		#	print "The python and c++ FP operators did not match!"
		#	exit(-1)
		#else:
		#	print "\t FP's match!"
		
"""