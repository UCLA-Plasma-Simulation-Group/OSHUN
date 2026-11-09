"""
	Module for optional ES wave driver. It makes an external electric field to stimulate 
		a moncromatic (well, as monochomatic as possible depending on the situation)


"""

"""
	This module relies on the following external resources:
		The global simulation STATE function
			It writes all of it's configuration info into the "STATE.metadata.driver" variable
				so that these paramters can be used by any postprocessing (it is saved via Python's Pickling)
		The 'EFieldSinusodialProfile' function

"""

import oshun_modules.formulas.plasma_dispersion
import OshunGlobalState as GLOBAL
from OshunGlobalState import SAVEFIG

class WaveDrivah:
	def __init__( self, num_wavelengths ):
		pass
		
		self.driver_on = False 
		
		# the number of wavelengths this driver is...
		self.num_wavelengths = num_wavelengths
		
		self.starting_amp 	=  1e-8
		self.sustain_amp 		=	1e-8
		self.final_amp			=	None
		self.amp_rate 			= 	0.0
		self.fall_rate 		= 	0.0
		
		# these are in units of driver period.
		#		fall_time = 0 means no fall_time
		#		rise_time = 0 means no rise_time
		self.rise_time_in_time_steps = False
		self.rise_time 		= 2.5
		self.sustain_time		= 0.0
		self.fall_time 		= 0.0
		
		# should dt be chosen by the driver?
		self.choose_dt_from_driver = False
		self.choose_dt_from_driver__factor = 50.0
		
		# time when echo pulse should be enacted...
		self.enable_2nd_pulse		= False
		self.echo_2nd_pulse_time 	= 0.0
		self.echo_2nd_pulse_num_wavelengths = 1.0
		self.has_2nd_pulse_not_fired	= True
		self.echo_2nd_pulse_time_rise_time = -1.0
		
		# this is a stic (Standing) wave?
		self.is_standing_wave = True
		
		# these are settings for spatial echo drivering(for finding sptial_echos)
		# this mode will still use the rise time settings...
		self.use_as_spatial_echo_driver = False
		
		self.x_global_index_for_w1 = None
		self.x_global_index_for_w2 = None
		self.x_for_w1 = None
		self.x_for_w2 = None

		self.spatial_echo_w1 = 0.0
		self.spatial_echo_w2 = 0.0
		
		self.w_wave		= None
		
	def init(self, vth, dt, wp=1.0):
		
		#TODO: replace...
		self.vth = vth
		self.wp = wp
		
		# internal variables...
		self.current_amp = self.starting_amp
		spatial_axis_global = STATE.spatial_axes__global[0]
		dx = spatial_axis_global.dx()
		box_size = (spatial_axis_global.values[-1] - spatial_axis_global.values[0] )
		
		self.num_cells_in_box = spatial_axis_global.values.shape[0]
		#self.wavelength = ( box_size + dx ) / (float( self.num_wavelengths))
		self.wavelength = ( box_size ) / (float( self.num_wavelengths))
		self.wavenumber = 2.0*np.pi/( self.wavelength )
		
		self.buffah = STATE.E_field.clone()
		self.buffah_accumulate = STATE.E_field.clone()
		
		self.e_profiles = []
		
		#for labmdaa in xrange(int(self.num_wavelengths),int(self.num_wavelengths)+1): #(25, 26):
		#	self.e_profiles.append( EFieldSinusodialProfile(self.starting_amp, num_wavelengths=labmdaa) )
		
		self.e_profiles.append( EFieldSinusodialProfile(self.starting_amp, num_wavelengths=self.num_wavelengths) )
		
		self.w_wave__complex = plasma_dispersion.landau_damping_frequency(self.wp, self.vth, self.wavenumber)
		self.w_wave__theory = np.real( self.w_wave__complex )
		self.w_damping__theory = np.imag( self.w_wave__complex )
		
		#static wave
		if self.is_standing_wave:
			self.w_wave = 0
		elif self.w_wave == None:
			self.w_wave = np.abs( self.w_wave__theory )
		
		
		#self.w_wave = np.sqrt( 1 + 3*self.wavenumber*self.wavenumber*self.vth*self.vth)
		self.v_phase =  (self.w_wave/ self.wavenumber)
		
		
		if self.use_as_spatial_echo_driver:
			# determine the indices of the spatital locations...
			if self.x_for_w1 != None and self.x_for_w2:
				# determine the indices for these locations (and update the locations with exact grid values)
				self.x_global_index_for_w1 = int( self.x_for_w1 /  STATE.spatial_axes__global[0].dx() )
				self.x_global_index_for_w2 = int( self.x_for_w2 /  STATE.spatial_axes__global[0].dx() )
				self.x_for_w1 = STATE.spatial_axes__global[0].values[ self.x_global_index_for_w1 ]
				self.x_for_w2 = STATE.spatial_axes__global[0].values[ self.x_global_index_for_w2 ]				 
			else:
				self.x_for_w1 = STATE.spatial_axes__global[0].values[ self.x_global_index_for_w1 ]
				self.x_for_w2 = STATE.spatial_axes__global[0].values[ self.x_global_index_for_w2 ]
			self.update = self.update_spatial_echo_driver
			
		
		# interperate the rise_time as the number of driver wave cycles....
		if self.w_wave > 0.0:
			self.rise_time__time_units = (np.pi * 2.0 * self.wp / self.w_wave) * self.rise_time
		else:
			self.rise_time__time_units = (np.pi * 2.0 * self.wp /  np.abs( self.w_wave__theory ) ) * self.rise_time
		
		
		
		self.driver_off = False
		
		if  self.w_wave > 0.0:
			self.suggested_dt = (2.0*np.pi/ self.choose_dt_from_driver__factor ) * self.w_wave
		else:
			self.suggested_dt = (2.0*np.pi/ self.choose_dt_from_driver__factor ) *  np.abs( self.w_wave__theory ) 
		
		if self.use_as_spatial_echo_driver:
			self.suggested_dt = (2.0*np.pi/ self.choose_dt_from_driver__factor ) * self.x_for_w1
		
		if self.choose_dt_from_driver:
			self.dt = self.suggested_dt
		else:
			self.dt = dt
		dt = self.dt
		
		# actually, the rise time was meant to be in usints if timesteps if this flag set..
		if self.rise_time_in_time_steps:
			self.rise_time__time_units = self.rise_time * dt
		
		
		if self.rise_time__time_units > 0.0:
			self.amp_rate = (self.sustain_amp - self.starting_amp) / (self.rise_time__time_units / dt )
		else:
			self.amp_rate = 0.0;
			
		
		try:
			self.amp_fall_rate = 0.0
		except:
			self.amp_fall_rate = 0.0
		
		
		

			
		"""
		f = open('!sim_info.txt','a')
		f.write("\n\n")
		f.write("Thermal velocity (%s) (ev): %f \r\n" % (self.species.species_config.name, vth_in_ev) )
		f.write( "expect wavelength (imposed by driver such that %d spatial periods fit in the box of size %e):  %f \r\n" % (self.num_wavelengths, box_size, self.wavelength) )
		f.write( "assigned wavenumber (2*pi/( wavelength)):  %f \r\n" % (self.wavenumber,) )
		f.write( "freq (from Bohm-gross dispersion, t=%f (normalized) wavenumber %f, wp~=1 (normalized) ): %e \r\n" % (self.vth, self.wavenumber, self.w_wave) )
		f.write( "vphase ((bohm gross w)/(wavenumber)): %e  \r\n" % self.v_phase)
		f.write( "l_debye ( (thermal velocity (normalized)) / (wp (normalized)) ): %f c/wp \r\n" % (self.vth / self.wp) )
		f.write( "k*vth: %f \r\n" % ( (self.wavenumber*self.vth) ) )
		#f.write( "T_bounce: %f \r\n" % ( (  ) ) )
		f.write("\n\n")
		f.close()
	
		self.start_of_2nd_pulse = True
		self.export_extern_field	= True
		
		STATE.metadata.driver = {}
		STATE.metadata.driver['driver_stop_time'] = self.sustain + self.rise_time
			
		STATE.metadata.driver['starting_amp'] =  self.starting_amp
		STATE.metadata.driver['ending_amp_after_rise'] =  self.ending_amp
		STATE.metadata.driver['amp_rate'] = self.amp_rate
		STATE.metadata.driver['final_amp_after_fall'] =  self.final_amp
		STATE.metadata.driver['rise_time'] = self.rise_time
		STATE.metadata.driver['sustain']  = self.sustain
		STATE.metadata.driver['fall_time'] = self.fall_time
		
		
		
		STATE.metadata.driver['vth__ev'] = vth_in_ev
		STATE.metadata.driver['num_wavelengths'] = self.num_wavelengths
		STATE.metadata.driver['wavelength'] = self.wavelength
		STATE.metadata.driver['wavenumber'] = self.wavenumber
		STATE.metadata.driver['w_wave'] = self.w_wave
		STATE.metadata.driver['v_phase'] = self.v_phase
		STATE.metadata.driver['lambda_debye'] = (self.vth / self.wp)
		STATE.metadata.driver['k__vth'] = (self.wavenumber*self.vth)
		
		STATE.metadata.driver['num_cells_wave_travels_in_dt'] = self.num_cells_wave_travels_in_dt
		STATE.metadata.driver['phase_step_per_dt'] = self.phase_step_per_dt
		STATE.metadata.driver['pure_harmonic_dt'] = self.pure_harmonic_dt
		STATE.metadata.set_dirty()
		"""
		
		STATE.metadata.driver = {}
		STATE.metadata.driver['rise_time__time_units'] = self.rise_time__time_units
		STATE.metadata.driver['num_wavelengths'] = self.num_wavelengths
		STATE.metadata.driver['lambda_debye'] = (self.vth / self.wp)
		STATE.metadata.driver['k__vth'] = (self.wavenumber*self.vth)
		STATE.metadata.driver['w_wave__theory'] = self.w_wave__theory
		STATE.metadata.driver['w_damping__theory'] = self.w_damping__theory
		STATE.metadata.driver['w_wave'] = self.w_wave
		
		STATE.metadata.driver['x_global_index_for_w1'] = self.x_global_index_for_w1
		STATE.metadata.driver['x_global_index_for_w2'] = self.x_global_index_for_w2
		STATE.metadata.driver['x_for_w1'] = self.x_for_w1
		STATE.metadata.driver['x_for_w2'] = self.x_for_w2
		
		
			
		
		STATE.metadata.set_dirty()
		
		STATE.use_external_E_field = True
		self.export_extern_field	= True
		STATE.E_field__external = STATE.E_field.clone()
		STATE.E_field__external.scalar_mult(0.0)
	
	"""
	def normalized_vth_to_ev(self, normalized_vth):
		vth = normalized_vth 
		ev1 = vth/np.sqrt(vth*vth +1)*3e10/4.19e7
		return ev1*ev1*self.species.species_config.mass
	"""
	
	# working (but Bulgarian) full box traveling wave...
	def update_spatial_echo_driver(self, time, n, dt):
		
		
		# with this driver, we don't really ever turn it off... we have it run contineously.
		# with maybe rise times for w1,w2
		self.driver_on = False
		STATE.use_external_E_field = False
		
		if time > ( self.rise_time__time_units + self.sustain_time + self.fall_time):
			return
		
		self.buffah_accumulate.scalar_mult(0.0)
		
		num_cells_in_this_node = len(mpi_info.mpi_axes[0].axis__local.values) # + STATE.E_field.num_boundry_cells[0][0] + STATE.E_field.num_boundry_cells[0][1]
		starting_global_cell = mpi_info.mpi_axes[0].global_cell_offset[0]
		ending_global_cell = starting_global_cell + num_cells_in_this_node
		
		if mpi_info.rank == 0:
			print " In Spatial echo Driver!!!!! doing CONTINeous exciting!!!!!"
		
		if self.x_global_index_for_w1 >= starting_global_cell and self.x_global_index_for_w1 <= ending_global_cell:
			# we are the node that nneds to cycle w1...
			local_spatial_cell_index_w1 = self.x_global_index_for_w1 - starting_global_cell +  STATE.E_field.num_boundry_cells[0][0]
			self.buffah_accumulate.components[0][local_spatial_cell_index_w1] = self.current_amp*np.sin( self.spatial_echo_w1 * time )
			
		if self.x_global_index_for_w2 >= starting_global_cell and self.x_global_index_for_w2 <= ending_global_cell:
			# we are the node that nneds to cycle w2...
			local_spatial_cell_index_w2 = self.x_global_index_for_w2 - starting_global_cell +  STATE.E_field.num_boundry_cells[0][0]
			self.buffah_accumulate.components[0][local_spatial_cell_index_w2] = self.current_amp*np.sin( self.spatial_echo_w2 * time )
		
		self.driver_on = True
		STATE.use_external_E_field = True
		
		 
		#update the external E field...
		STATE.E_field__external.copy( self.buffah_accumulate )
		
		if time < self.rise_time__time_units:
			self.current_amp += self.amp_rate
			pass
		
	def update(self, time, n, dt):
		
		
		if mpi_info.rank == 0:
			print " In Driver!!!!! with w %f and rise time %f" % ( self.w_wave, self.rise_time__time_units )
		
		self.driver_on = False
		STATE.use_external_E_field = False
		
		go_away 	= False
		if time > ( self.rise_time__time_units + self.sustain_time + self.fall_time):
			go_away = True
		if self.enable_2nd_pulse:
			if time >= self.echo_2nd_pulse_time and time < (self.echo_2nd_pulse_time + self.rise_time__time_units + self.sustain_time+ self.fall_time):
				go_away = False
				time = time - self.echo_2nd_pulse_time
				if self.has_2nd_pulse_not_fired:
					self.current_amp = self.starting_amp
					k2 = self.echo_2nd_pulse_num_wavelengths # now try 48 /4 #45 #self.num_wavelengths+25
					self.e_profiles[0] = EFieldSinusodialProfile(self.starting_amp, num_wavelengths=k2)	
					self.export_extern_field = True
				self.has_2nd_pulse_not_fired = False
				if mpi_info.rank == 0:
					print "ECHO PUSLEEEE!!!!!!!!!!!!!! with %d" % self.e_profiles[0].num_wavelengths
				pass
			pass
		pass
		
		if mpi_info.rank == 0:
			print " In Driver!!!!! goways is"
			print go_away
		
		if go_away:			
			return
		
		self.driver_on = True
		STATE.use_external_E_field = True
		
		
		wt = -1.0*time*self.w_wave
		self.buffah_accumulate.scalar_mult(0.0)
		
		for i in xrange( 0,len(self.e_profiles) ):
			self.e_profiles[i].amp = self.current_amp
			self.e_profiles[i].apply( STATE.E_field__external, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah, phase_me=wt )
			self.buffah_accumulate.add( self.buffah)
		
		STATE.E_field__external.copy( self.buffah_accumulate )
		
		if self.export_extern_field:
			if( mpi_gather.gather(  STATE.E_field__external.components[0], remove_guard_cells=True ) ):
				clf(); plot( mpi_gather.stiching_buffer ); SAVEFIG( "E__external_driver_%06d" % n);
				np.save( GLOBAL.OUTPUT_DIR( "E__external_driver_%06d" % n, 'numpy' ) , mpi_gather.stiching_buffer)
		self.export_extern_field = False
		
		#STATE.E_field.copy( self.buffah_accumulate )
		
		if time < self.rise_time__time_units:
			self.current_amp += self.amp_rate
			pass
		
		if time > ( self.rise_time__time_units+ self.sustain_time):
			#self.current_amp += self.amp_fall_rate
			if mpi_info.rank == 0:
				print "lowered driver to drivah (amp=%e) (rate=%e) " % (self.current_amp,  self.amp_fall_rate)
		else:
			if mpi_info.rank == 0:
				print "Done with drivah (amp=%e) (rate=%e) " % (self.current_amp,  self.amp_rate)
	
	
	def update1(self, time, n, dt):
		# attempt at simple static wave amping...
		if n % self.update_interval_steps != 0:
			if mpi_info.rank == 0:
				if self.driver_off:
					print "E field driver is off."
				else:
					print "%d steps until E field push. (cuurent amp=%e) " % ( (self.update_interval_steps - (n % self.update_interval_steps)), self.current_amp)
			return
		
		"""
		multiple = time / self.update_interval 
		if multiple <= self.last_seen_multiple:
			return
		self.last_seen_multiple = multiple
		"""
		
		self.current_amp *= self.amp_rate
		if self.current_amp >= self.ending_amp:
			if mpi_info.rank == 0:
				print "wave driver upped wave to %f already and is OFF." % self.current_amp
				self.driver_off = True
			return
		
		self.e_profile.amp = -1.0*self.current_amp
		self.e_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah=self.buffah )
		STATE.E_field.add( self.buffah)
		
		if mpi_info.rank == 0:
			print "wave driver upped wave to %f" % self.current_amp
		
		#vth = CFG.species_config[0].temp_profile.val
		#w_wave = np.sqrt( 1 + 3*CFG.E_field_profile.wavenumber*CFG.E_field_profile.wavenumber*vth*vth)
		#print "freq (BG): %e " % w_wave
		#print "vphase: %e" % (w_wave/ CFG.E_field_profile.wavenumber)
		#print "K l_debye: %f" % (CFG.E_field_profile.wavenumber*vth)
		#CFG.E_field_profile.v_phase = (w_wave/ CFG.E_field_profile.wavenumber)
