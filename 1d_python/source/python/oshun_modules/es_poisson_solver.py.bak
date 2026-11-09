
class PerterbationESolver:
	def __init__(self):
		self.solve_mode_moments 			= True
		self.solve_mode_local 				= False
		self.solve_mode_fft 				= False
		
		self.init_f0_density  				= None
		self.init_f1_density  				= None
		
		self.local_f0__initial 				= None
		self.local_f1__initial			 	= None
		
		self.E_field_perturbation		 	= STATE.E_field.clone()
		self.full_E_field_unperturbed		= None
		self.full_E_field_perturbation	= None
		pass
	
	def set_F0_state(self, species_index, moments):
		
		# use moments to get density..
		if self.solve_mode_moments:
			if( mpi_gather.gather(  moments[species_index].get_density(force_update = True), remove_guard_cells=True ) ):
				self.init_f0_density = np.copy ( mpi_gather.stiching_buffer )
				#print "initial"
				#print self.init_f0_density
		
		if self.solve_mode_local:
			# save inital F
			self.local_f0__initial = Esolva.get_charge_density_for_species(0)
	
	def set_F1_state(self, species_index, moments):
		
		if self.solve_mode_moments:
			if(mpi_gather.gather(  moments[species_index].get_density(force_update = True), remove_guard_cells=True ) ):
				self.init_f1_density = np.copy ( mpi_gather.stiching_buffer )
				#print "F1 initial"
				#print self.init_f1_density
				
				self.init_f1_density -= self.init_f0_density
				
				#self.init_f1_density *= -1.0
				#print "F1 after"
				#print self.init_f1_density
				
				#clf()
				#plot( self.init_f1_density); savefig( 'cunt_f1.png')
				#clf()
				#plot( self.init_f1_density); savefig( 'cunt_f0.png')
		
		# record new F and subtract form unperturbed F to get 'f1'
		if self.solve_mode_local:
			self.local_f1__initial = Esolva.get_charge_density_for_species(0)
			self.local_f1__initial -= self.local_f0__initial
			self.local_f1__initial *= -1.0
	
	def solve_for_E(self, moments, species_index, E_out=None, overwrite_E=True, add_to_E =False, swap_guard_cells=True ):
		
		buffer = STATE.E_field.clone()
		STATE.config.E_field_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=buffer, sine=True )
		buffer.components[0] *= STATE.species[species_index].species_config.charge
		
		# E is the integral of the charge density... so we also need to divide by a factor of 1/k
		box_size = 	STATE.spatial_axes__global[0].max
		num_wavelengths = STATE.config.E_field_profile.num_wavelengths
		wavelength = box_size / float(STATE.config.E_field_profile.num_wavelengths)
		dx = STATE.spatial_axes__global[0].dx()
		k_for_mode_plain = 2.0*np.pi / wavelength
		k_for_mode = 2.0*np.pi*float( num_wavelengths)/(box_size +dx )
		amp = STATE.config.E_field_profile.amp
		
		#k_for_mode *= .8495
		#print "DIVIDING BY K WHICH IS k_for_mode: %e to give norminal %e or actual %e" % ( k_for_mode, amp / k_for_mode,  amp / k_for_mode_plain)
		#print "ratio from plain %e %e" % ( (k_for_mode_plain/k_for_mode), (k_for_mode/k_for_mode_plain) )
		
		buffer.components[0] /= k_for_mode_plain #k_for_mode
		
		self.full_E_field_perturbation
		
		if overwrite_E:
			STATE.E_field.components[0][:] = buffer.components[0][:]
			if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
				self.full_E_field_perturbation = np.copy( mpi_gather.stiching_buffer )
		if add_to_E:
			STATE.E_field.components[0][:] += buffer.components[0][:]
			if mpi_gather.gather(  buffer.components[0][:], remove_guard_cells=True):
				self.full_E_field_perturbation = np.copy( mpi_gather.stiching_buffer )
		
		"""
		self.solve_for_E1( moments, species_index, E_out=buffer, overwrite_E=True, add_to_E =False, swap_guard_cells=True )
		if mpi_info.rank == 0:
			print  buffer.components[0] 
		print "MAX of integrated E %e" % np.max( buffer.components[0] )
		print "MAX div 4pi %e" %  np.max( buffer.components[0] /4.0 / np.pi )
		print "MAX div k %e" %  np.max( buffer.components[0] /k_for_mode )
		"""
	
	def solve_for_E1(self, moments,species_index, E_out=None, overwrite_E=True, add_to_E =False, swap_guard_cells=True ):
		# do perturbation
		#STATE.config.E_field_profile.apply_perturbation( STATE.species[0].F, STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset )
		
		
		number_field_data_cells = STATE.E_field.num_plain_cells[0]
		num_boundry_cells_left  = STATE.E_field.num_boundry_cells[0][0]
		num_boundry_cells_right = STATE.E_field.num_boundry_cells[0][1]
		num_cells = STATE.spatial_axes__global[0].num_points
		box_size = 	STATE.spatial_axes__global[0].max
		#print "box size: \t\t %f" % float(box_size)
		#print "num_cells: \t\t %f" % float(num_cells)
		#print "delta_x: \t\t %f (offical from axis: %f)" % ( (float(box_size) / float(num_cells)), STATE.spatial_axes__global[0].dx()  )
		delta_x = float(box_size) / ( float(num_cells) + 1.0 )
		delta_x_2 = ( float(box_size) + STATE.spatial_axes__global[0].dx() ) / float ( num_cells )
		#print "adjusted dx  \t\t %f" % delta_x
		#print "adjusted dx2 \t\t %f" % delta_x_2
		self.recieve_buffer_for_E_field = np.zeros (  STATE.spatial_axes__global[0].values.shape )
		
		recieve_buffer_for_E_field = self.recieve_buffer_for_E_field
		init_f1_density = self.init_f1_density
		
		
		# use Moments to solve for E
		if(self.solve_mode_moments and mpi_gather.gather(  moments[species_index].get_density(force_update=True), remove_guard_cells=True ) ):
			
			
			init_f1_density__integrate_buffer = np.zeros(  (STATE.spatial_axes__global[0].num_points + 2 ) )
			#init_f1_density__integrate_buffer[1:] =  init_f1_density[:]
			#init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-1]
			init_f1_density__integrate_buffer[0:-2] =  init_f1_density[:]
			init_f1_density__integrate_buffer[-2] =  init_f1_density__integrate_buffer[0]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			
			#init_f1_density__integrate_buffer[1:-1] =  init_f1_density[:]
			init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-2]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			#print init_f1_density__integrate_buffer
			
			# we fed in a vector of length N+1 and will get out a vector of length N
			temp = scipy.integrate.cumtrapz( init_f1_density__integrate_buffer, dx=( delta_x_2 ) )
			
			recieve_buffer_for_E_field[0] = temp[-2]
			recieve_buffer_for_E_field[1:] = temp[0:-2]
			print "MEAN!!!!!!! %e" % ( np.mean( recieve_buffer_for_E_field ) )
			
			recieve_buffer_for_E_field *= (4.0*np.pi)
			# clf(); plot( recieve_buffer_for_E_field ); savefig('initial_pulse__E_field_perturbation')
			
			"""
			# basiclly works.. just need to hunt down normalization factors!
			# and try an FFT!!!
			dx = STATE.spatial_axes__global[0].dx()
			k_axis = np.arange(0.0, float(num_cells) ) * dx
			k_axis *= np.pi
			k_axis *= ( 2.0*np.pi / ( float(box_size) +  dx ) )
			fft_n = fft.fft( init_f1_density )
			
			symetic_k_axis = np.arange( float(-1*num_cells/2), float(num_cells/2) )
			symetic_k_axis *= np.pi
			symetic_k_axis *= dx
			symetic_k_axis /= ( 2.0*np.pi / ( float(box_size) +  dx ) )
			
			
			fft_n = fftshift( fft_n )
			
			for ik in xrange( 1 , num_cells ):
				if ik == (num_cells/2):
					print " skipped %f" % symetic_k_axis[ ik ]
					continue
				
				fft_n[ik] /= np.complex( 0, -1.0* symetic_k_axis[ ik ] )
				#fft_n[ik] /= np.complex( 0, -1.0* k_axis[ ik ] )
			
			E = fft.ifft( ifftshift( fft_n ) )
			E1 = E / dx
			E2 = E / np.sqrt(dx)
			print E1
			print E2
			"""
		
			
		
		# This gets E from local densoty...
		if( self.solve_mode_local and  mpi_gather.gather(  local_f1__initial, remove_guard_cells=True ) ):
			# only rank==0 can get in here.. and we have all 'local_f1__initial' buffers stiched togther!
			init_f1_density__integrate_buffer = np.zeros(  (STATE.spatial_axes__global[0].num_points + 2 ) )
			#init_f1_density__integrate_buffer[1:] =  mpi_gather.stiching_buffer[:]
			#init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-1]
			init_f1_density__integrate_buffer[0:-2] =  mpi_gather.stiching_buffer[:]
			init_f1_density__integrate_buffer[-2] =  init_f1_density__integrate_buffer[0]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			
			#init_f1_density__integrate_buffer[1:-1] =  mpi_gather.stiching_buffer[:]
			init_f1_density__integrate_buffer[0] = init_f1_density__integrate_buffer[-2]
			init_f1_density__integrate_buffer[-1] = init_f1_density__integrate_buffer[1]
			#print init_f1_density__integrate_buffer
			
			# we fed in a vector of length N+1 and will get out a vector of length N
			temp = scipy.integrate.cumtrapz( init_f1_density__integrate_buffer, dx=( delta_x_2 ) )
			
			
			recieve_buffer_for_E_field[0] = temp[-2]
			recieve_buffer_for_E_field[1:] = temp[0:-2]
			#print "MEAN!!!!!!! %f" % ( np.mean( recieve_buffer_for_E_field ) )
			
			recieve_buffer_for_E_field *= (4.0*np.pi)
			clf(); plot( recieve_buffer_for_E_field ); savefig('ass_cunt_E')
			#print recieve_buffer_for_E_field
			#print "MEAN final!!!!!!! %f" % ( np.mean( recieve_buffer_for_E_field ) )
			#exit(-1)
		
		self.full_E_field_perturbation = recieve_buffer_for_E_field
		data_start_index = mpi_info.mpi_axes[0].global_cell_offset[0]
		number_valid_data_cells = STATE.E_field.num_plain_cells[0]
		data_end_index = data_start_index + number_valid_data_cells
		
		#print "broadcast the integrated E field!"
		recieve_buffer_for_E_field = mpi_info.comm.bcast( recieve_buffer_for_E_field, root=0)
		#print "USE data from others to get what E is in our Local neighborhood"		
		appropiate_data = recieve_buffer_for_E_field[data_start_index:data_end_index ]
		
		# save this nodes portion of the E filed perurbation in case it's needed for something else.
		self.E_field_perturbation.components[0][ num_boundry_cells_left:num_boundry_cells_left+number_valid_data_cells] = appropiate_data[:]
		
		
		if E_out == None:
			target_E = STATE.E_field
		else:
			target_E = E_out
		
		if add_to_E:
			target_E.components[0][ num_boundry_cells_left:num_boundry_cells_left+number_valid_data_cells] += appropiate_data[:]
		elif overwrite_E:
			target_E.components[0][ num_boundry_cells_left:num_boundry_cells_left+number_valid_data_cells] = appropiate_data[:]
		else:
			return appropiate_data
		
		
		if swap_guard_cells:
			#print "NOw updates guard cells."
			# update everyones's gaurd cells...
			if mpi_info.num_nodes == 1:
				update_boundries.calc(STATE, STATE.config)
			else:
				mpi_update.update(STATE, STATE.config)		
			#print "Everyone's guard cells are good now!"
