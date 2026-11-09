from mpi4py import MPI
import numpy as np


class UpdateBoundriesMPI:
	def __init__(self):
		pass
	def init(self):
		pass

	def setBoundryType(self, isMirror=False, isPerodic=True):
		global mpi_info

		mpi_axis = mpi_info.mpi_axes[0]
		
		self.isMirror = isMirror
		self.isPerodic = isPerodic

		if( self.isMirror):
			self.update = self.update_mirror

			# a hack but one that i want to employ over th whole file...
			# 	if this an edge node, then alter the 'prepare_send_buffer' function of mpi_axis.
			#	The alteration consists od calling the original one, then applying the signs to the send buffer for mirror boundries.
			if mpi_info.is_I_root() or mpi_info.is_I_last_node():
				mpi_axis.prepare_send_buffer__core = mpi_axis.prepare_send_buffer 
				mpi_axis.prepare_send_buffer = mpi_axis.prepare_send_buffer_mirror__wrapper
			else:
				pass

		else:
			self.update = self.update_perodic
			self.isPerodic = True
		pass

	def verifiy_distibution(self, F):
		pass
	
	def test_distribution(self, state):
		F = state.species[0].F
		num_cells = len(F.cells)
		num_harmonics = F.cells[0].num_harmonics
		num_p = F.momentum_axis.num_points
		
		for ip in xrange(0, num_p):
			for sh in xrange(0, num_harmonics):
				for ix in xrange(0,num_cells):
					F.cells[ix].harmonics[sh].momentum_projection[ip] = mpi_info.rank
		return
	
	# TODO: These two update funtion should become OSHUN events.
	# TODO: stupid way to write this but I am so tired.. should be subclassed object and just take over for the non _mirror function
	def update_mirror(self, state, CFG ):
		mpi_axis = mpi_info.mpi_axes[0]
		# fix up species distributions..
		for i in xrange(0, len(state.species)):
			mpi_axis.reconcile_boundries_F_mirror( state.species[i].F )
		
		# fix up field distributions..
		data = state.E_field.components[0]
		mpi_axis.reconcile_boundries_field_mirror( data )
		# in 1-D, the Ex field does not need any sign changes when using mirror boundries.

		if state.use_external_E_field:
			
			data = state.E_field__external.components[0]
			mpi_axis.reconcile_boundries_field_mirror( data )
			#print "DONE DIING FILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		
		# fix up field distributions..
		data = state.B_field.components[0]
		mpi_axis.reconcile_boundries_field_mirror( data )
		# maybe Jay too..

	def update_perodic(self, state, CFG):
		
		mpi_axis = mpi_info.mpi_axes[0]
		
		# this writes a test pattern inot the distirbution function
		#self.test_distribution(state)
		
		# fix up species distributions..
		for i in xrange(0, len(state.species)):
			mpi_axis.reconcile_boundries_F( state.species[i].F )
		
		# fix up field distributions..
		data = state.E_field.components[0]
		mpi_axis.reconcile_boundries_field( data )
		
		if state.use_external_E_field:
			
			data = state.E_field__external.components[0]
			mpi_axis.reconcile_boundries_field( data )
			#print "DONE DIING FILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		
		# fix up field distributions..
		data = state.B_field.components[0]
		mpi_axis.reconcile_boundries_field( data )
		# maybe Jay too..

	def update_diagnostic_data(self, buffer):
		if len(buffer.shape) == 3:
			pass
		if len(buffer.shape) == 2:
			pass
		if len(buffer.shape) == 1:
			mpi_axis.reconcile_boundries_field( data )

# try a few ways!!
class GatherOutputMPI:
	
	def __init__(self):
		pass
	
	def init(self, mpi_axis = None):
		if mpi_axis == None:
			mpi_axis = mpi_info.mpi_axes[0]
		self.stiching_buffer_backup = np.zeros( mpi_axis.global_num_cells[0] )
		num_cell_with_gc = mpi_axis.global_num_cells[0] + mpi_info.num_nodes* (mpi_axis.num_boundry_cells[0][0] + mpi_axis.num_boundry_cells[0][1])
		self.stiching_buffer_backup_with_gc = np.zeros( num_cell_with_gc )
		self.buffer_cache = {}
		pass
	
	#TODO: think... make a consistent API.... I actually dont know what axis is gonna be used
	# do really cant do the buffa scitching unless that info is given!
	def gather(self, data, remove_guard_cells = False, data_good_as_is = False ):
		
		mpi_axis = mpi_info.mpi_axes[0]
		self.stiching_buffer = None
		
		if mpi_info.num_nodes == 1:
			
			if data_good_as_is:
				self.stiching_buffer = data
				return True
			
			if remove_guard_cells:
				# assume this is subject to the mpi_axis[0] and proceed..
				lower_limit = mpi_axis.num_boundry_cells[0][0]
				upper_limit = mpi_axis.num_total_cells__local[0] -  mpi_axis.num_boundry_cells[0][1]
				self.stiching_buffer = self.stiching_buffer_backup
				self.stiching_buffer = data[lower_limit:upper_limit]
			else:
				self.stiching_buffer = self.stiching_buffer_backup_with_gc
				self.stiching_buffer = data[:]
			return True
		
		if data_good_as_is:
			self.stiching_buffer = self.stiching_buffer_backup
			mpi_info.comm.Gather(data, self.stiching_buffer)
			return mpi_info.rank == 0
		
		if remove_guard_cells:
			lower_limit = mpi_axis.num_boundry_cells[0][0]
			upper_limit = mpi_axis.num_total_cells__local[0] -  mpi_axis.num_boundry_cells[0][1]
			self.stiching_buffer = self.stiching_buffer_backup
			mpi_info.comm.Gather(data[lower_limit:upper_limit], self.stiching_buffer)
		else:
			self.stiching_buffer = self.stiching_buffer_backup_with_gc
			mpi_info.comm.Gather(data, self.stiching_buffer)
		return mpi_info.rank == 0
	
	def get_2d_data(self, name):
		return self.buffer_cache[name]
	
	def gather_2d(self, name, full_x, full_y, data):
		if mpi_info.num_nodes == 1:
			self.buffer_cache[name] = data
			return True
		
		if mpi_info.rank == 0:
			hash_string = "(%d_%d)" % (full_x, full_y )
			# TODO: look up better search in hash
			if hash_string in self.buffer_cache.keys():
				buffa = self.buffer_cache[hash_string]
			else:
				buffa = np.zeros( (full_x, full_y) )
				self.buffer_cache[hash_string] = buffa
				self.buffer_cache[name] = buffa
			mpi_info.comm.Gather(data, buffa)
		else:
			# other ranks don't need a recieve Buffa
			mpi_info.comm.Gather(data, None)
		return mpi_info.rank == 0
	
	def gather_full_F_with_gc(self, F, F_with_gc = None, verify =True):
		from oshun import  DistributionFunction 
		
		num_cells = len(F.cells)
		num_harmonics = F.cells[0].num_harmonics
		num_p = F.momentum_axis.num_points
		data_along_x = np.zeros( num_cells )
		num_cells_dest = self.stiching_buffer_backup_with_gc.size
		
		if F_with_gc == None and mpi_info.rank==0:
			F_with_gc  = DistributionFunction(num_cells_dest, F.num_l, F.num_m, F.momentum_axis, F.num_boundry_cells)
		
		# get a x-slice, use MPI to get a such slices from alll nodes.. keep guard cells.
		# this is for debugging MPI stuff.
		for ip in xrange(0, num_p):
			for sh in xrange(0, num_harmonics):
				# get data for a specific harmonic, |p| along the x diection.
				for ix in xrange(0,num_cells):
					data_along_x[ix] = F.cells[ix].harmonics[sh].momentum_projection[ip]
					pass
				pass
				# save data for a specific harmonic, |p| along the x diection (with all guard cells) to destination.
				if self.gather( data_along_x, keep_guard_cells=True):
					for ix in xrange(0,num_cells_dest):
						F_with_gc.cells[ix].harmonics[sh].momentum_projection[ip] = self.stiching_buffer[ix]
					pass
				pass
			pass
			
		if mpi_info.rank > 0:
			return False
		if verify == False:
			return False
		
		# now  verify that the boundry cells are perodic...
		num_brdy_cells_L = mpi_info.mpi_axes[0].num_boundry_cells[0][0]
		num_brdy_cells_U = mpi_info.mpi_axes[0].num_boundry_cells[0][1]
		num_total_cells_per_node = F.num_total_cells[0]
		count_elements_processed = 0
		
		for ip in xrange(0, num_p):
			for sh in xrange(0, num_harmonics):
				# get data for a specific harmonic, |p| along the x diection.
				for ix in xrange(0,num_cells_dest):
					self.stiching_buffer[ix] = F_with_gc.cells[ix].harmonics[sh].momentum_projection[ip]
				
				for node in xrange(0, mpi_info.num_nodes):
					#check lower boundry cells...
					node_L_bndry_cells_base_idx = num_total_cells_per_node*node
					node_U_bndry_cells_base_idx = num_total_cells_per_node*(node+1) - num_brdy_cells_U
					lower_node_data_idx = node_L_bndry_cells_base_idx - num_brdy_cells_U - num_brdy_cells_L
					upper_node_data_idx = num_total_cells_per_node*(node+1) + num_brdy_cells_L
					if( node == (mpi_info.num_nodes - 1) ):
						upper_node_data_idx = 0 + num_brdy_cells_L
					
					for n in xrange(0, num_brdy_cells_L):
						if self.stiching_buffer[node_L_bndry_cells_base_idx + n] != self.stiching_buffer[lower_node_data_idx + n]:
							print " mismatch at lower b-cell %d for node %d at sh=%d, p=%d" % (n, node, sh, ip )
							exit(-1)
					for n in xrange(0, num_brdy_cells_U):
						if self.stiching_buffer[node_U_bndry_cells_base_idx + n] != self.stiching_buffer[upper_node_data_idx + n]:
							print " mismatch at uppper b-cell %d for node %d at sh=%d, p=%d" % (n, node, sh, ip )
							print "\t node %d has %e, upper node has %e" % (n, self.stiching_buffer[node_U_bndry_cells_base_idx + n], self.stiching_buffer[upper_node_data_idx + n])
							print " node idex is %d, upper is %d" % (node_U_bndry_cells_base_idx + n, upper_node_data_idx + n)
							exit(-1)
					count_elements_processed += num_brdy_cells_L + num_brdy_cells_U
				pass
			pass
			
		print "Distribution function has valid boundry cells (%d elements processed.)" % count_elements_processed
		return
		
	
# copies a nornal axis, partions it.. and assigns.
# TODO: HACK
# HACK: take out species realteds stuff from here.. logiacally thos axis ahoudl not be tied to species (maybe?)
#			will sortra assume 1 species for now...
class AxisMpi:
	def __init__(self, axis, num_boundry_cells, mpi_info):
		# TODO: bad form.. get rid of this (need to get rid of cicrular depencies)
		from Oshun1d import Axis
		
		self.axis = axis
		self.min = axis.min
		self.max = axis.max
		self.num_points = axis.num_points
		
		# for now, simple splitting for the work...
		chunk_size  = int(self.num_points / mpi_info.num_nodes)
		chunk_start = mpi_info.rank*chunk_size
		chunk_end   = chunk_start + chunk_size
		
		# last node takes on any extra leftover work
		if	mpi_info.is_I_last_node():
			if chunk_size*mpi_info.num_nodes != self.num_points:
				extra_work = (self.num_points - chunk_size* mpi_info.num_nodes)
				chunk_size += extra_work
				chunk_end += extra_work
				print "EXTRA!"
		
		self.axis__local = Axis( axis.values[chunk_start], max=axis.values[chunk_end-1], num_points=chunk_size)
		self.num_boundry_cells 				= num_boundry_cells
		self.num_plain_cells__local 		= np.zeros ( 3 )
		self.num_plain_cells__local[0] 	= chunk_size 
		self.num_total_cells__local 		= np.zeros( 3 )
		self.num_total_cells__local[:] 	= num_boundry_cells[:,0]+self.num_plain_cells__local[:]+num_boundry_cells[:,1] 

		self.global_offset = np.zeros( 3)
		self.global_cell_offset = np.zeros( 3)
		self.global_cell_offset[0] = chunk_start
		self.global_offset[0] = axis.values[chunk_start]
		self.global_num_cells = np.zeros( 3)
		self.global_num_cells[0] = self.num_points
		
		# attach boundry consition types here for now..
		self.boundry__is_perodic = True
		
		if mpi_info.is_I_root():
			self.lower_neighbor = mpi_info.num_nodes - 1
			self.upper_neighbor = mpi_info.rank      + 1
		elif mpi_info.is_I_last_node():
			self.lower_neighbor = mpi_info.rank      - 1
			self.upper_neighbor = 							0
		else:
			self.lower_neighbor = mpi_info.rank      -1
			self.upper_neighbor = mpi_info.rank      +1
		
		mpi_gather.init(mpi_axis=self)
	
	def init_species_data(self, species, E_field, num_harmonics):
		
		p_num_points = species.momentum_axis.num_points
		
		# prepare buffers for distribution...
		self.buffer_send_lower = np.empty( (species.F.num_boundry_cells[0][1]*num_harmonics*p_num_points) )
		self.buffer_send_upper = np.empty( (species.F.num_boundry_cells[0][0]*num_harmonics*p_num_points) )
		self.buffer_recv_lower = self.buffer_send_upper
		self.buffer_recv_upper = self.buffer_send_lower
		
		# prepare buffers for fields...
		self.field_send_lower = np.empty( (E_field.num_boundry_cells[0][1],) )
		self.field_send_upper = np.empty( (E_field.num_boundry_cells[0][0],) )
		
	
	def prepare_send_buffer_mirror__wrapper(self,buffer, F, lower=True, mirror=False ):
		# call the normal prepare_send_buffer
		self.prepare_send_buffer__core( buffer, F, lower=lower, mirror=True)

		if mpi_info.is_I_last_node() and mpi_info.is_I_root():
			pass
		elif mpi_info.is_I_root() and lower==False:
			return
		elif mpi_info.is_I_last_node() and lower==True:
			return
		elif mpi_info.is_I_last_node() == False and mpi_info.is_I_root()==False:
			return

		#if lower:
		#	print "NODE %d doing LOWER" % ( mpi_info.rank, )
		#else:
		#	print "NODE %d doing UPPER" % ( mpi_info.rank, )
		# now apply the needed signs
		num_harmonics = F.cells[0].num_harmonics; p_num_points = F.momentum_axis.num_points;
		num_bdry_cells_L = self.num_boundry_cells[0,0]; num_bdry_cells_U = self.num_boundry_cells[0,1];
		num_total_cells = F.num_total_cells[0]
		dest_start = 0; dest_end = p_num_points;
		if lower:
			num_cells = num_bdry_cells_U
		else: 
			num_cells = num_bdry_cells_L
			
		for n in xrange(0, num_cells):
			for sh in xrange(0, num_harmonics):
				# the sign for reflection is (-1)^(m+n)
				sign = 1.0 - 2.0*( (F.cells[n].harmonics[sh].my_l+F.cells[n].harmonics[sh].my_m) % 2)
				buffer[dest_start:dest_end] *= sign
				dest_start += p_num_points
				dest_end   += p_num_points
		pass
	def prepare_send_buffer(self, buffer, F, lower=True, mirror=False):
		num_harmonics = F.cells[0].num_harmonics; p_num_points = F.momentum_axis.num_points;
		num_bdry_cells_L = self.num_boundry_cells[0,0]; num_bdry_cells_U = self.num_boundry_cells[0,1];
		num_harmonics = F.cells[0].num_harmonics; num_total_cells = F.num_total_cells[0]
		if lower:
			dest_start = 0; dest_end = p_num_points;
			src_cell_start = num_bdry_cells_L;
			num_cells = num_bdry_cells_U
		else: 
			dest_start = 0; dest_end = p_num_points;
			src_cell_start = (num_total_cells-num_bdry_cells_L-num_bdry_cells_U);
			num_cells = num_bdry_cells_L
				
		src_cell_end=src_cell_start+num_cells; src_cell_step = 1;
		do_mirror = False
		if  mpi_info.is_I_root() and lower:
			self.do_mirror = mirror
		if  mpi_info.is_I_last_node() and lower==False:
			self.do_mirror = mirror

		if do_mirror:
			# swap src_start and src_end (and subtract 1 for both) and set inteeration step to -1
			src_cell_step = src_cell_end -1; src_cell_end = src_cell_start -1; src_cell_start = src_cell_step;
			src_cell_step = -1;
			#print "start %d to %d " % ( src_cell_start, src_cell_end)
		#for n in xrange(0, num_cells):
		for n in xrange(src_cell_start,src_cell_end,src_cell_step):
			for sh in xrange(0, num_harmonics):
				buffer[dest_start:dest_end] = F.cells[n].harmonics[sh].momentum_projection[:]
				dest_start += p_num_points
				dest_end   += p_num_points
			#src_cell += 1	
	def process_recv_buffer(self, buffer, F, lower=True):
		num_harmonics = F.cells[0].num_harmonics; p_num_points = F.momentum_axis.num_points;
		num_bdry_cells_L = self.num_boundry_cells[0,0]; num_bdry_cells_U = self.num_boundry_cells[0,1];
		num_harmonics = F.cells[0].num_harmonics; num_total_cells = F.num_total_cells[0]
		if lower:
			src_start = 0; src_end = p_num_points; dest_cell = 0; 
			num_cells = num_bdry_cells_L
		else: 
			src_start = 0; src_end = p_num_points; dest_cell = (num_total_cells-num_bdry_cells_U);  
			num_cells = num_bdry_cells_U
		for n in xrange(0, num_cells):
			for sh in xrange(0, num_harmonics):
				#if mpi_info.rank==0:
				#	print "cell=%d, sh=%d" % (dest_cell,sh )
				F.cells[dest_cell].harmonics[sh].momentum_projection[:] = buffer[src_start:src_end]
				src_start += p_num_points
				src_end   += p_num_points
			dest_cell += 1
	

	def reconcile_boundries_field(self, field):
		num_bdry_cells_L = self.num_boundry_cells[0,0]; num_bdry_cells_U = self.num_boundry_cells[0,1];
		num_total_cells = field.shape[0]
		
		if(mpi_info.rank %2 == 0):
			src_start = num_bdry_cells_L; src_end = num_bdry_cells_L+num_bdry_cells_U;
			mpi_info.comm.Send(field[src_start:src_end] ,   dest = self.lower_neighbor)
			mpi_info.comm.Recv(field[0:num_bdry_cells_L], source = self.lower_neighbor)
			src_start = (num_total_cells-num_bdry_cells_L-num_bdry_cells_U); src_end = num_total_cells-num_bdry_cells_U
			mpi_info.comm.Send(field[src_start:src_end],   dest = self.upper_neighbor)
			dest_start = (num_total_cells-num_bdry_cells_U); dest_end = num_total_cells
			mpi_info.comm.Recv(field[dest_start:dest_end], source = self.upper_neighbor)
		else:
			dest_start = (num_total_cells-num_bdry_cells_U); dest_end = num_total_cells
			mpi_info.comm.Recv(field[dest_start:dest_end], source = self.upper_neighbor)
			src_start = (num_total_cells-num_bdry_cells_L-num_bdry_cells_U); src_end = num_total_cells-num_bdry_cells_U
			mpi_info.comm.Send(field[src_start:src_end],   dest = self.upper_neighbor)
			mpi_info.comm.Recv(field[0:num_bdry_cells_L], source = self.lower_neighbor)
			src_start = num_bdry_cells_L; src_end = num_bdry_cells_L+num_bdry_cells_U;
			mpi_info.comm.Send(field[src_start:src_end] ,   dest = self.lower_neighbor)			
	
	#def reconcile_boundries_F(self, F):
	#	self.reconcile_boundries(F, self.prepare_send_buffer_F, self.process_recv_buffer_F)
	
	#def reconcile_boundries(self, F, prepare_send_buffer, process_recv_buffer):
	def reconcile_boundries_F(self, F):
		# TODO: move this..
		# do blocking for now since too easy in 1D (then try non-blocking)
		if(mpi_info.rank %2 == 0):
			self.prepare_send_buffer(self.buffer_send_lower, F, lower=True)
			#np.save('send_lower_0', self.buffer_send_lower)
			mpi_info.comm.Send(self.buffer_send_lower,   dest = self.lower_neighbor)
			mpi_info.comm.Recv(self.buffer_recv_lower, source = self.lower_neighbor)
			#np.save('recv_lower_0', self.buffer_recv_lower)
			self.process_recv_buffer(self.buffer_recv_lower, F, lower=True)
			
			self.prepare_send_buffer(self.buffer_send_upper, F, lower=False)
			#np.save('send_upper_0', self.buffer_send_upper)
			mpi_info.comm.Send(self.buffer_send_upper,   dest = self.upper_neighbor)
			mpi_info.comm.Recv(self.buffer_recv_upper, source = self.upper_neighbor)
			#np.save('recv_upper_0', self.buffer_recv_upper)
			self.process_recv_buffer(self.buffer_recv_upper, F, lower=False)
		else:
			mpi_info.comm.Recv(self.buffer_recv_upper, source = self.upper_neighbor)
			#np.save('recv_upper_1', self.buffer_recv_upper)
			self.process_recv_buffer(self.buffer_recv_upper, F, lower=False)
			self.prepare_send_buffer(self.buffer_send_upper, F, lower=False)
			#np.save('send_upper_1', self.buffer_send_upper)
			mpi_info.comm.Send(self.buffer_send_upper,   dest = self.upper_neighbor)
			
			mpi_info.comm.Recv(self.buffer_recv_lower, source = self.lower_neighbor)
			#np.save('recv_lower_1', self.buffer_recv_lower)
			self.process_recv_buffer(self.buffer_recv_lower, F, lower=True)
			self.prepare_send_buffer(self.buffer_send_lower, F, lower=True)
			mpi_info.comm.Send(self.buffer_send_lower,   dest = self.lower_neighbor)
			#np.save('send_lower_1', self.buffer_send_lower)

	def reconcile_boundries_field_mirror(self, field):
		num_bdry_cells_L = self.num_boundry_cells[0,0]; num_bdry_cells_U = self.num_boundry_cells[0,1];
		num_total_cells = field.shape[0]
		#print "there last 8 now at BEFORE"
		#print field[-8:]
		# even nodes do this.
		if(mpi_info.rank %2 == 0):
			src_start = num_bdry_cells_L; src_end = num_bdry_cells_L+num_bdry_cells_U;
			if  mpi_info.is_I_root():
				source_start = src_end-1;source_end = src_start-1
				field[0:num_bdry_cells_L] = field[source_start:source_end:-1]
				"""print "aaa to copy START"
				print field[source_start:source_end:-1]
				print "there first 8 now at START"
				print field[:8]"""

			else:
				mpi_info.comm.Send(field[src_start:src_end] ,   dest = self.lower_neighbor)
				mpi_info.comm.Recv(field[0:num_bdry_cells_L], source = self.lower_neighbor)


			src_start = (num_total_cells-num_bdry_cells_L-num_bdry_cells_U); src_end = num_total_cells-num_bdry_cells_U
			dest_start = (num_total_cells-num_bdry_cells_U); dest_end = num_total_cells
			if mpi_info.is_I_last_node():
				src_start -= 1;src_end -= 1;
				field[dest_start:dest_end] = field[src_end:src_start:-1]
				"""print "aaa to copy END"
				print field[src_end:src_start:-1]
				print "there last 8 now at END"
				print field[-8:]"""
			else:
				mpi_info.comm.Send(field[src_start:src_end],   dest = self.upper_neighbor)				
				mpi_info.comm.Recv(field[dest_start:dest_end], source = self.upper_neighbor)

		# odd nodes do this.
		else:
			src_start = (num_total_cells-num_bdry_cells_L-num_bdry_cells_U); src_end = num_total_cells-num_bdry_cells_U
			dest_start = (num_total_cells-num_bdry_cells_U); dest_end = num_total_cells
			if mpi_info.is_I_last_node():
				src_start -= 1;src_end -= 1;
				field[dest_start:dest_end] = field[src_end:src_start:-1]
			else:
				mpi_info.comm.Recv(field[dest_start:dest_end], source = self.upper_neighbor)			
				mpi_info.comm.Send(field[src_start:src_end],   dest = self.upper_neighbor)

			src_start = num_bdry_cells_L; src_end = num_bdry_cells_L+num_bdry_cells_U;
			if  mpi_info.is_I_root():
				source_start = 2*num_bdry_cells_L-1;source_end = num_bdry_cells_L-1
				field[0:num_bdry_cells_L] = field[source_start:source_end:-1]
			else:
				mpi_info.comm.Recv(field[0:num_bdry_cells_L], source = self.lower_neighbor)
				mpi_info.comm.Send(field[src_start:src_end] ,   dest = self.lower_neighbor)			
	
	#def reconcile_boundries_F(self, F):
	#	self.reconcile_boundries(F, self.prepare_send_buffer_F, self.process_recv_buffer_F)
	
	def reconcile_boundries_F_mirror(self, F):
		# TODO: move this..
		# do blocking for now since too easy in 1D (then try non-blocking)
		
		# even nodes
		if(mpi_info.rank %2 == 0):
			if  mpi_info.is_I_root():
				self.prepare_send_buffer(self.buffer_send_lower, F, lower=True, mirror=True)
				self.process_recv_buffer(self.buffer_send_lower, F, lower=True)
			else:
				self.prepare_send_buffer(self.buffer_send_lower, F, lower=True)
				#np.save('send_lower_0', self.buffer_send_lower)
				mpi_info.comm.Send(self.buffer_send_lower,   dest = self.lower_neighbor)
				mpi_info.comm.Recv(self.buffer_recv_lower, source = self.lower_neighbor)
				#np.save('recv_lower_0', self.buffer_recv_lower)
				self.process_recv_buffer(self.buffer_recv_lower, F, lower=True)
			
			if mpi_info.is_I_last_node():
				self.prepare_send_buffer(self.buffer_send_upper, F, lower=False, mirror=True)
				self.process_recv_buffer(self.buffer_send_upper, F, lower=False)
			else:				
				self.prepare_send_buffer(self.buffer_send_upper, F, lower=False)
				#np.save('send_upper_0', self.buffer_send_upper)
				mpi_info.comm.Send(self.buffer_send_upper,   dest = self.upper_neighbor)
				mpi_info.comm.Recv(self.buffer_recv_upper, source = self.upper_neighbor)
				#np.save('recv_upper_0', self.buffer_recv_upper)
				self.process_recv_buffer(self.buffer_recv_upper, F, lower=False)

		# odd nodes
		else:
			if mpi_info.is_I_last_node():
				self.prepare_send_buffer(self.buffer_send_upper, F, lower=False, mirror=True)
				self.process_recv_buffer(self.buffer_send_upper, F, lower=False)
			else:
				mpi_info.comm.Recv(self.buffer_recv_upper, source = self.upper_neighbor)
				#np.save('recv_upper_1', self.buffer_recv_upper)
				self.process_recv_buffer(self.buffer_recv_upper, F, lower=False)
				self.prepare_send_buffer(self.buffer_send_upper, F, lower=False)
				#np.save('send_upper_1', self.buffer_send_upper)
				mpi_info.comm.Send(self.buffer_send_upper,   dest = self.upper_neighbor)
			
			if  mpi_info.is_I_root():
				self.prepare_send_buffer(self.buffer_send_lower, F, lower=True, mirror=True)
				self.process_recv_buffer(self.buffer_send_lower, F, lower=True)
			else:			
				mpi_info.comm.Recv(self.buffer_recv_lower, source = self.lower_neighbor)
				#np.save('recv_lower_1', self.buffer_recv_lower)
				self.process_recv_buffer(self.buffer_recv_lower, F, lower=True)
				self.prepare_send_buffer(self.buffer_send_lower, F, lower=True)
				mpi_info.comm.Send(self.buffer_send_lower,   dest = self.lower_neighbor)
				#np.save('send_lower_1', self.buffer_send_lower)
		
		
		
# simple MPI implementation....
class MPIInfo:

	def __init__(self):
		self.comm = MPI.COMM_WORLD
		self.rank = self.comm.Get_rank()
		self.num_nodes = self.comm.size
		self.mpi_axes = []

	def init(self):
		# basic MPI startup...
		#self.comm = MPI.COMM_WORLD
		#self.rank = self.comm.Get_rank()
		#self.num_nodes = self.comm.size
		#print "I am rank %d!!!" % self.rank
		#
		#self.mpi_axes = []
		pass
	def abort(self):
		self.comm.Abort()
	def is_I_root(self):
		return self.rank == 0
	def is_I_last_node(self):
		return ( self.rank == (self.num_nodes - 1) )

	# TODO: move these to more logical place one day
	def barrierWait( self ):
		mpi_info.comm.Barrier()
	def bail( self ):
		mpi_info.comm.MPI_Finalize()
	

mpi_info = MPIInfo()
mpi_update = UpdateBoundriesMPI()
mpi_gather = GatherOutputMPI()





def non_blocking_test():
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	if rank == 0:
		print "Hi! I am root! I will do by best to serve honorably in this most esteemed position."
		data = numpy.arange(10000, dtype='i')
		info = comm.Isend([data, MPI.INT], dest=1)
		print info
	else:
		data = numpy.empty(10000, dtype='i')
		print "Hi i am node %d" % rank
		info = comm.Irecv([data, MPI.INT], source=0)
		print "waiting for data!"
		while not info.Test():
			print "........"
		print "Node %d says: I just got:" % rank
		print data
			
		pass
	

#non_blocking_test()