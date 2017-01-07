

	# this was baisc unit testing code for the deritivitives..
	#TODO:		move this to a real unit test.
"""
	crap_org = np.array([1,34,12,54,-7,-9,5,7,3,5,8,2,41,23,53])
	crap = np.copy( crap_org )
	cunt = crap[-2] - crap[-1]
	print crap
	outme = np.zeros( crap.shape )
	print "Grreeek"
	print test_central_difference_derivitive_micheal( crap,1.0)
	print "numpy"
	
	
	crap = np.array([1,34,12,54,-7,-9,5,7,3,5,8,2,41,23,53])
	on_me = HarmonicCentricView( np.array( [1,0,0] ), Axis( 0, num_points=crap.shape[0],max=100.0) )
	on_me.data[0,:] = crap_org[:]
	on_me.Dp()
	print "HARMONICS!!!!!!!"
	print on_me.data

	crap = zeros( (crap_org.shape[0]) )
	crap[:] = crap_org[:]
	#crap = np.array([1,34,12,54,-7,-9,5,7,3,5,8,2,41,23,53,cunt])
	outme = np.zeros(crap.shape[0] )
	#ass = np.gradient( crap)
	outme[1:-1] = ( crap[:-2] - crap[2:]   )
	outme[0] = 0.0
	outme[-1] = crap_org[-2] - crap_org[-1];
	print outme
	
	#ass *= -2.0
	#ass[0] = 0.0; ass[-1] = crap_org[-2] - crap_org[-1];
	#print ass
	exit(-1)
"""

	# also other unit testing type code for basic numerics
"""
		if False:
		data = np.array([ [1.,34.,12.,54.,-7.,-9.,5.,7.,3.,5.,8.,2.,41.,23.,53.] ])
		temp_data = zeros( data.shape)
		temp_data1 = zeros( data.shape)
		
		temp_data [:,1:] = data[:,:-1]
		print temp_data
		temp_data1[:,:-1] = data[:,1:]
		print temp_data1
		temp_data1 -= temp_data
		print temp_data1
		temp_data1 *= 8.0
		temp_data1[:,2:] += data[:,:-2]
		temp_data1[:,:-2] -= data[:,2:]
		data[:,:] = temp_data1[:,:]
		data /= 12.0  #(the other 2.0*h) will come later in the code
		
		print data
		exit(-1)
"""


"""
# In depth unit testing code for comparing agains the c++ version

# These was a globl flag 'test_mode = False'
Which was called in 3 places.. one in main loop
	if test_mode:
		testing1_init()

Another right after basic State setup:
	if test_mode:
		CFG = spitzer_setup()
And another in main
if 	test_mode:
		testing1_main(eval_functions)

# and here are the funcitons themselves..
def basic_e_on_dist_test1__init():
	# assume STATE is setup with something....
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# zero the DIST function...
	for xx in xrange(0, Y0.num_total_cells[0]):
		for l in xrange(0, Y0.num_l):
			Y0.get_harmonic_at_location(xx, l).momentum_projection[:] = 0.0
			#m_stop = l
			#if(m_stop > Y0.num_m):
			#	m_stop = Y0.num_m
			#for m in xrange(0, m_stop):
			#	Y0.get_harmonic_at_location(xx, l,m).momentum_projection[:] = 0.0
	# initalize a fake dist....
	# first up... spaital isotropic.. inital values only in Y(0,0)
	for xx in xrange(0, Y0.num_total_cells[0]):
		#p_vecter = Y0.get_harmonic_at_location(xx,0,0).momentum_projection
		p_vecter = Y0.get_harmonic_at_location(xx,0).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = pp*.1
	pass
	
	
def basic_e_on_dist_test1__loop():
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(0, 0).momentum_projection[:]
	print Yh.get_harmonic_at_location(0, 0).momentum_projection[:]
	
	evalf.e_on_dist.calc								(Y0, Yh, species, STATE, CFG)
	
	print Y0.get_harmonic_at_location(0, 1).momentum_projection[:]
	Y0.get_harmonic_at_location(0, 1).momentum_projection[0] = 666.0
	
	print Yh.get_harmonic_at_location(0, 1).momentum_projection[:]
	exit(-1)
	pass

def basic_adv_test1__init():
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# initalize a fake dist....
	# make spatital gradients for the advection codess....
	numx = 128
	period = 2*np.pi / (numx + 1)
	print "perios %f" % period
	
	for xx in xrange(0, Y0.num_total_cells[0]):
		cos_val = np.cos(period* float(xx));
		if cos_val < 0.0:
			cos_val *= -1.0
		p_vecter = Y0.get_harmonic_at_location(xx,0).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = .1*pp*cos_val
	pass

def basic_adv_test1__loop():
	evalf = eval_functions[0]
	print "I WIN!!!!"
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(2, 0).momentum_projection[:]
	
	evalf.spatial_advection_solver.calc								(Y0, Yh, species, STATE, CFG)
	
	
	print Yh.get_harmonic_at_location(2, 1).momentum_projection[:]
	exit(-1)
	pass


def testing1_init():
	
	# assume STATE is setup with something....
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	
	# zero the DIST function...
	for xx in xrange(0, Y0.num_total_cells[0]):
		for l in xrange(0, Y0.num_l):
			m_stop = l
			if(m_stop > Y0.num_m):
				m_stop = Y0.num_m
			for m in xrange(0, m_stop):
				Y0.get_harmonic_at_location(xx, l).momentum_projection[:] = 0.0
	
	# fake a dumb dist function.... that has a current...
	for xx in xrange(0, Y0.num_total_cells[0]):
		p_vecter = Y0.get_harmonic_at_location(xx,1).momentum_projection
		for pp in xrange(0, Y0.momentum_axis.num_points):
			p_vecter[pp] = pp*.1
	
	
	
def testing1_main(eval_functions):
	evalf = eval_functions[0]
	print "I WIN!!!!"
	CFG = STATE.config
	species_config = CFG.species_config[0]
	species = STATE.species[0]
	Y0 = species.Y0
	Yh = species.Yh
	Yh.copy( Y0 )
	
	# fake an E field...
	STATE.E_field.components[0][:] = 0.4
	
	print Y0.get_harmonic_at_location(2, 1).momentum_projection[:]
	
	evalf.current_solver.calc								(Y0, Yh, species, STATE, CFG)
	
	
	print Yh.get_harmonic_at_location(2, 1).momentum_projection[:]
	exit(-1)
	pass
"""



""" IN EFieldSinusodialProfile:apply

			# smaple of adressing scheme
			clf(); 
			ddatt = np.zeros( ( num ) )
			xx = np.array( xs1)
			k = 2.0*np.pi*float( self.num_wavelengths)/(box_size +dx )
			print "nox length = "
			print (2.0*np.pi / k * self.num_wavelengths)
			print "xx"
			print xx
			for i in xrange(0, num):
				ddatt[i] = self.amp * sincos(k *xx[i] + phase_me)
			datt = ddatt[3:-3]
			the_xs = xx[3:-3]
			print "datt"
			print datt
			print datt.shape
			print the_xs.shape
			print "k (adj): %e" % k
			print "k (old): %e" % ( (2.0*np.pi / box_size), )
			old_xs = np.copy( the_xs )
			old_xs *= (2.0*np.pi / box_size)
			old_data = np.cos ( old_xs )
			print "ols data"
			print old_data
			ddd = k - (2.0*np.pi / box_size)
			print (2.0*pi /ddd)
			print (2.0*pi /ddd) * box_size
			clf()
			stem( the_xs, datt)
			axhline(y=0)
			savefig('shsvh', dpi=200)
			clf()
			ffft = np.abs(fftshift(fft.fft( datt )) ) + 1e-30
			plot( np.log10( ffft ) )
			savefig('shsvh1', dpi=200)
			#self.v_phase = 1.0/self.wavenumber
			#print phase_cell_offset
			#print phase_me
			#print "FOR %d (%d): \n %s" % ( mpi_info.rank, len( indic), np.array( indic ) )
			exit(-1)
"""

""" FRom EFieldSinusodialProfile:appy_frequency

data = E_field.components[0]
if buffah != None:
	data = buffah.components[0]
	dx = spatial_axis_global.dx()
	self.wavelength = box_size/ (float( self.num_wavelengths))
	self.wavenumber = 2.0*np.pi/(self.wavelength + dx)
	discrete_factor = np.pi*2*float( self.num_wavelengths)/spatial_axis_global.values.shape[0]
	
	data = E_field.components[0]
	if buffah != None:
		data = buffah.components[0]
	
	#for (n,x_val) in enumerate(spatial_axis__local_in_global.values):
pass

"""

"""
HARMONICCENTRICVIEW:dx
	#These was another commented out implementations that had all the indices reverse

	self.temp_data[:,0] = 0.0
	self.temp_data1[:,-1] = 0.0

	self.temp_data [:,1:] = self.data[:,:-1]; 
	self.temp_data1[:,:-1] = self.data[:,1:]; 
	self.temp_data1 -= self.temp_data
	self.temp_data1 *= 8.0

	#self.temp_data1[:,0:2] = 0.0
	#self.temp_data1[:,-1] = 0.0

	self.temp_data1[:,2:] += self.data[:,:-2]; 
	self.temp_data1[:,:-2] -= self.data[:,2:];
	self.data[:,:] = self.temp_data1[:,:]
	self.data /= 6.0  #(the other 2.0*h) will come later in the code





"""



"""
Stupid stuff to make interactive running simulations..

import threading
#import msvcrt
from multiprocessing import Pipe
import PyQt4.QtCore as QtCore
from PyQt4 import QtGui

def launch_plots():
	global pipe1main, pipe1upd	
	pipe1main, pipe1upd = Pipe()
	
def pause_sim(interval = .1):
	while( KEY.pause_sim ):
		show(block=False)
		draw()
		#time.sleep(interval)
	return

class GuiThread(threading.Thread):
	def __init__(self):
		super(GuiThread, self).__init__()
	def run(self):
		print("Starting GUI Loop")
		show()
		print("ending GUI Loop")

# this function was commented out
class KeyboardInterrupt(threading.Thread):
	def __init__(self):
		super(KeyboardInterrupt, self).__init__()
		self._stop = False
		self.quit_sim = False
		self.pause_sim = False
	def run(self):
		#set initial 'ch', alternatively check for ch != None in the while loop
		ch = '1'	
 
		#loop until 
		while (not self._stop):
			#msvcrt.getch is blocking so wrap it in a hit test loop
			while (msvcrt.kbhit()):
				ch = msvcrt.getch()
				if (ch in ['q']):
					print "Keyboard interrupt triggered, exiting program"
					self.quit_sim = True
				if(ch in [' ']):
					if self.pause_sim == False:
						self.pause_sim = True
					else:
						self.pause_sim = False
			
			time.sleep(0.030)
 
	def stop (self):
		self._stop = True

"""


"""
# Scraps of unit testing code about deritivitives

#matplotlib.interactive(True)
# standard central difference....
def central_difference_derivitive(vals, dx):
	result = np.ones( len(vals) )
	for i in range(1,len(vals)-1):
		result[i] = vals[i+1] - vals[i-1]
		result[i] /= 2.0*dx
	return result

# Micheal's central difference from the C++ code:
# 		is -1 times the normal one
#		does not do the division by "2h"
#		defines values that the normal central differnece does not:
#			element at 0: set to be the same as element 1
#			element at last position: unchanged.
def test_central_difference_derivitive_micheal(vals, dx):
	for i in range(0, len(vals)-2):
		vals[i] -= vals[i+2]
	for i in range(len(vals)-3, -1, -1):
		vals[i+1] = vals[i]
	print vals

# convert standrd central differnce results to the ones from Micheal just
#		as a test that I is doing it right.
def test_central_difference_boring_against_micheal(vals, dx):
	results = central_difference_derivitive(vals, dx)
	results *= (2*dx)
	results *= -1.0
	results[0] = results[1]
	results[len(vals)-1] = vals[len(vals)-1]
	return results

"""


"""
Scraps realted to testing the Two Stream simulations

		if two_stream == True:
			temp_species.F.scalar_mult ( 0.0 )
			p_beam = 1.0
			p_beam_center = 1.0
			beam_temp = np.sqrt(.011) #.105    #0.011#0.105
			index_00 = temp_species.F.linear_index_for_harmonic(0)
			index_20 = temp_species.F.linear_index_for_harmonic(2)
			index_10 = temp_species.F.linear_index_for_harmonic(1)
			
			for ix in xrange(0, temp_species.F.num_plain_cells[0] ):
				#seed = np.cos(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix])
				"""
				if ix==0 and mpi_info.rank==0: 
					xval = 1.387778780781e-17
					seed = .1*np.sin(.20*np.pi*xval) # WORKS
				else:
					seed = .1*np.sin(.20*np.pi*(STATE.spatial_axes[0].values[ix])) # WORKS
				"""
				#if ix==0:
				#	STATE.spatial_axes[0].values[ix] += STATE.spatial_axes[0].dx()
				#if ix==temp_species.F.num_plain_cells[0]-1:
				#	STATE.spatial_axes[0].values[ix] -= STATE.spatial_axes[0].dx()
				#seed = .1*np.sin(.20*np.pi*(STATE.spatial_axes[0].values[ix])) # WORKS
				
				# niec deritivitce
				dx = STATE.spatial_axes[0].dx()
				seed = .1*np.sin(2*np.pi/(10.0+.5*dx)*(STATE.spatial_axes[0].values[ix]))
				
				
				#seed = .1*np.sin(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx())*(STATE.spatial_axes[0].values[ix] - STATE.spatial_axes[0].dx()))
				#seed = .1*np.sin(2.0*np.pi/(10.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix])
				#seed *= 10
				#seed *= 10.0 #check...
				#p_beam = 1.2*np.cos(2.0*np.pi/(20.0+STATE.spatial_axes[0].dx() )*STATE.spatial_axes[0].values[ix]) + p_beam_center
				exp_arg_coeff =  -1.0/(2.0*np.power(beam_temp,2))
				factor_coeff = 1.0/( np.power(beam_temp, 3.0) *np.power(2.0*np.pi, 1.5) )
				
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					abs_p = STATE.species[i].momentum_axis.values[ip]
					exp_argument = exp_arg_coeff*(abs_p - p_beam)*(abs_p - p_beam)
					factor = .0104*factor_coeff
					factor *= np.exp(exp_argument)
					pedistal = .07*heaviside_step(p_beam-abs_p)
					
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_00].momentum_projection[ip] = pedistal + factor
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_20].momentum_projection[ip] = 2.0*factor
					temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_10].momentum_projection[ip] = factor*.2*seed #factor*(1.0/50.0)*seed
				
			
			#DEBUG( Y=temp_species.F,L=0, M=0)
			
			enable = False
			if mpi_info.rank == 0 and enable:
				x_points = [0,1,2,3]
				print ""
				for ix in x_points:
					print "\t%f" % STATE.spatial_axes[0].values[ix],
				print ""
				
				print "beam temp %f" % beam_temp
				print "EXP arg factor is: %f" % exp_arg_coeff
				print "FROT factor is: %f" % factor_coeff
				print "Y(0,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_00].momentum_projection[ip] ),
					print ""
				print ""
				print ""
				print "Y(1,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_10].momentum_projection[ip] ),
					print ""
				print ""
				print ""
				print "Y(2,0)"
				print "----------------------------------------------"
				for ip in xrange(0, STATE.species[i].momentum_axis.num_points):
					print "%d: "%( ip,) ,
					for ix in x_points:
						print "\t%.10e"%( temp_species.F.cells[num_boundry_cells[0][0]+ix].harmonics[index_20].momentum_projection[ip] ),
					print ""
				#exit(-1)
			
"""


"""
Fragemnts for E sovler inital profile and wave driver stuff that was in main loop


"""
# this is the frist field solver.. really.. given and E, produces a perturbation in F that is self consistent
#	this E..
if STATE.config.E_field_profile != None and False:
	
	solva = EFieldSolve1d()
	print "Gonna solve!"
	exit(-1)
	solva.solve( STATE.E_field )
	print "Done solvaing"
	# redo boundry conditions...
	if mpi_info.num_nodes == 1:
		update_boundries.calc(STATE, STATE.config)
	else:
		mpi_update.update(STATE, STATE.config)		
"""


#------------------ inital F perturbation and making E consistent...
init_f0_density = None
init_f1_density = None
Esolva = EFieldSolve1d()
F_from_External_E_Solver = None

if STATE.config.E_field_profile != None and True and DaWaveDrivah == None:
	F_from_External_E_Solver = PerterbationESolver()
	species_index = 0
	#print "FUCK"
	
	F_from_External_E_Solver.set_F0_state( species_index, moments)
	
	"""
	if mpi_gather.gather( moments[i].get_density( force_update = True ), remove_guard_cells=True ):
		#shitface = np.copy(  mpi_gather.stiching_buffer )
		#shitface += F_from_External_E_Solver.init_f1_density
		clf();subplot(211); title('IInit Density from Momet'); plot( mpi_gather.stiching_buffer  )
		da_fft = fftshift( np.absolute( fft.fft( mpi_gather.stiching_buffer  ) ) )
		da_fft += 1e-30
		subplot(212);  title('FFt Init Density from Momet ');  plot( np.log10(da_fft) ); savefig( "bopVVimomet__density_pertub_fft_%06d"% ( n, ), dpi=200)
	"""
	
	# do perturbation
	# THIS IS THE FUCKED UP SHITFACE HARMONIC MAKER
	STATE.config.E_field_profile.apply_perturbation( STATE.species[species_index].F, STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset )
	
	"""
	if mpi_gather.gather( moments[i].get_density( force_update = True ), remove_guard_cells=True ):
		#shitface = np.copy(  mpi_gather.stiching_buffer )
		#shitface += F_from_External_E_Solver.init_f1_density
		clf();subplot(211); title('IInit Density from Momet'); plot( mpi_gather.stiching_buffer  )
		da_fft = fftshift( np.absolute( fft.fft( mpi_gather.stiching_buffer  ) ) )
		da_fft += 1e-30
		subplot(212);  title('FFt Init Density from Momet ');  plot( np.log10(da_fft) ); savefig( "cunt__density_pertub_fft_%06d"% ( n, ), dpi=200)
	"""
	
	# mark down this changed state for comparison...
	F_from_External_E_Solver.set_F1_state(species_index, moments)		
	# now get the E field that corresponds to that perterbation..
	# for now, kill the E filed and replace it with the calculated field...
	
	print STATE.E_field.components[0][:]
	STATE.E_field.components[0][:] = 0.0
	## just to test perturbation onlyy... no field..
	F_from_External_E_Solver.solve_for_E(moments, species_index, E_out=None, overwrite_E=True, add_to_E =False, swap_guard_cells=True )
	
	#plot( F_from_External_E_Solver.self.init_f1_density
	
	# some quick diagnostics...
	if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
		#clf(); subplot(311); title('Unperturbed F Field'); plot(  E_field_unperturbed )
		charge_density = F_from_External_E_Solver.init_f1_density[:] * STATE.species[species_index].species_config.charge
		clf();subplot(211); title('Initial Pulse: charge Density Perturbation'); plot( charge_density );
		subplot(212); title('Current E-field (Sum of Original and Perturbation)'); plot(  mpi_gather.stiching_buffer  ); savefig('initial_pulse__e_fields', dpi=200);
		clf();subplot(211); title('Initial Pulse: Unperturbed Plamsa Density'); plot( F_from_External_E_Solver.init_f0_density )
		subplot(212);  title('Initial Pulse: Plamsa Density Perturbation');  plot( F_from_External_E_Solver.init_f1_density ); savefig( "initial_pulse__f0_and_f1", dpi=200)
		
		"""
		shitface = np.copy( F_from_External_E_Solver.init_f0_density)
		shitface += F_from_External_E_Solver.init_f1_density
		clf();subplot(211); title('Initial Pulse: init + preturb Density '); plot( shitface )
		da_fft = fftshift( np.absolute( fft.fft( shitface ) ) )
		subplot(212);  title('FFt Init + pertub Density ');  plot(da_fft ); savefig( "initial_pulse__density_pertub_fft", dpi=200)
		
		clf();subplot(211); title('Initial charge den '); plot( charge_density )
		da_fft = fftshift( np.absolute( fft.fft( charge_density ) ) )
		subplot(212);  title('FFt Init + inital charage Density ');  plot(da_fft ); savefig( "initial_charge_denity_fft", dpi=200)
		
		clf();subplot(211); title('Initial Pulse: Density Perturbation'); plot( F_from_External_E_Solver.init_f1_density )
		da_fft = fftshift( np.absolute( fft.fft( F_from_External_E_Solver.init_f1_density ) ) )
		subplot(212);  title('FFt Density Perturbation');  plot(da_fft ); savefig( "initial_pulse__density_pertub_fft", dpi=200)
		"""
	


"""





"""
moving ion code that was in main loop


right after the lines:

for i in range(0, num_species):
	eval_functions.append( F_explicit(STATE.species[i], STATE, STATE.config) )

I had these lines commented out

# play with frame...
#TODO: TAKE OUT
#STATE.species[i].moving_frame_u_field[:] = DaWaveDrivah.v_phase
#STATE.species[i].moving_frame_u_field[:] = DaWaveDrivah.v_phase
#STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset)


#fred = EFieldSinusodialProfile( 0.02, 	num_wavelengths = 3.0 )
#fred.apply(STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset,  buffah= STATE.species[i].moving_frame_u_field)

#mpi_info.mpi_axes[0].reconcile_boundries_field( STATE.species[i].moving_frame_u_field.components[0] )
#print STATE.species[i].moving_frame_u_field.components[0]

# TESTING FLUID IONS...
#if DaWaveDrivah != None:
#	STATE.species[i].moving_frame_u_field.components[0][:] = DaWaveDrivah.v_phase

#if mpi_info.rank > 0:
#	print "------------ Fludid velocity---------------------"
#	print STATE.species[i].moving_frame_u_field.components[0][:]


"""

"""
#wave driver code that was in main loop

# After the following real line
solva.step( eval_functions[i], STATE.species[i], STATE, STATE.config )

			
			
			#if DaWaveDrivah != None:
			#	if STATE.time <= (DaWaveDrivah.rise_time + DaWaveDrivah.sustain):
			#		print "doing Fucker-Plank"
			#		fokker_planks[i].calc_f00(STATE.species[i].F, STATE.species[i], STATE, STATE.config, n, time, dt)
			#		fokker_planks[i].calc_flm(STATE.species[i].F, STATE.species[i], STATE, STATE.config)

Right after FP part of main loop

if STATE.config.E_field_profile != None and time > STATE.config.E_field_profile.echo_time and STATE.config.E_field_profile.echo_time > 0.0 and True and DaWaveDrivah == None:
	if mpi_info.rank == 0:
		print "--------------------PERTURB PULSE!!!!!!!!!!!!!!!!!!!!!!!!!!=========="
	
	abort_echo_pulse = False
	try:
		second_pulse_amp = STATE.config.E_field_profile.second_pulse_amp
		STATE.config.E_field_profile.amp = second_pulse_amp
	except:
		pass
		abort_echo_pulse = True
	try:
		second_pulse_k = STATE.config.E_field_profile.second_pulse_num_wavelengths
		STATE.config.E_field_profile.num_wavelengths = second_pulse_k
	except:
		pass
		abort_echo_pulse = True
	
	if STATE.config.E_field_profile.echo_time < 0.0:
		abort_echo_pulse = True
	
	if abort_echo_pulse:
		if mpi_info.rank == 0:
			print "Echo pulse was not launched due to confugation error(s)"
			
	else:
		
		# save the unpertubed E-field...
		if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
			E_field_unperturbed = np.copy( mpi_gather.stiching_buffer )
		
		
		#F_from_External_E_Solver = PerterbationESolver()
		species_index = 0
		F_from_External_E_Solver.set_F0_state( species_index, moments)
		
		# do perturbation
		# want to perturb system with output = current_state + amp*cos*current_state
		STATE.config.E_field_profile.apply_perturbation( STATE.species[species_index].F, STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset )
		
		# mark down this changed state for comparison...
		F_from_External_E_Solver.set_F1_state(species_index, moments)

		
		
		
		# now get the E field that corresponds to that perterbation..
		# for now, kill the E filed and replace it with the calculated field...
		F_from_External_E_Solver.solve_for_E(moments, species_index, E_out=None, overwrite_E = False, add_to_E = True, swap_guard_cells=True )
		
		# some quick diagnostics...
		if mpi_gather.gather( STATE.E_field.components[0], remove_guard_cells=True):
			clf(); subplot(311); title('Unperturbed F Field'); plot(  E_field_unperturbed )
			print "Preturabtion BEfore:"
			print  F_from_External_E_Solver.full_E_field_perturbation
			subplot(312); title('Perturbation to the E field'); plot( F_from_External_E_Solver.full_E_field_perturbation)
			print "Preturabtion After:"
			print  F_from_External_E_Solver.full_E_field_perturbation
			subplot(313); title('Current E-field (Sum of Original and Perturbation)'); plot(  mpi_gather.stiching_buffer  )
			savefig('echo_pulse__e_fields', dpi=200)
			clf(); subplot(211); title('Echo Pulse: Unperturbed Plamsa Density'); plot( F_from_External_E_Solver.init_f0_density )
			subplot(212);  title('Echo Pulse: Plamsa Density Perturbation');  plot( F_from_External_E_Solver.init_f1_density ); savefig( "echo_pulse__f0_and_f1", dpi=200)
		
		# shut off the echo driver..
		STATE.config.E_field_profile.echo_time = 1e100
		
		# normal field solver...
		#eval_functions[0].current_solver.calc( STATE.species[0].F, None, None, None, None, Eout = STATE.E_field )
		
		# this block as all commented out
		#STATE.config.E_field_profile.apply( STATE.E_field, STATE.spatial_axes__global[0], mpi_info.mpi_axes[0].axis__local,  mpi_info.mpi_axes[0].global_cell_offset, buffah = buffah_accumulate)
		#print "ADTERT"
		#print buffah_accumulate
		#STATE.E_field.components[0] += buffah_accumulate.components[0]
		
"""