import numpy as np
import scipy.special
import scipy.optimize

#
# import matplotlib, but don't trow an error if it is not found...
#		matplotlib is only used plots when testing the plasma dispersion code.
#
try:
	import matplotlib
	#import pylab
except:
	pass


def plasma_dispersion( value):
	a = scipy.special.wofz(value)
	#a *= 1j
	a *= np.sqrt(np.pi)*1j
	return a

def plasma_dispersion_prime( value):
	z = plasma_dispersion( value )
	return -2.0*(1.0+value*z)
	
def landau_damping_frequency(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0, inital_root_guess = None ):
	#
	# we are gonna solve for the roots of the electrostatic warm plasma dispersion relation
	#	e.g. 	solve:  	eplison(w,k) = 0
	#			where:		eplison(w,k) = 1 + sum_of_terms_for_all_species_of_form( chi_s(w,k) ) 
	#			here:				( for example, chi_e = electron susceptibility or chi_i = ion susceptibility etc....)
	#			now:						chi_s(x,s) = - ( wp_s/(2*vth_s*k_wave) )^2 * plasma_dispersion_prime( w_wave / (sqrt(2)*vth_s*k_wave) )
	#												for any species s.
	#														wp_s = plasma frew for species s
	#														vth_s = thermal speed f species s
	#														k_wave = the wavenumber of the wave traveling in the plasma.
	#														w_wave = the frequency of the wave traveling in the plasma.
	#
	#		We will use this function to, given a wave's wavenumber as input, solve the dipersion realtion to find the (complex) frequency (i.e. pole
	#			of the dispersion equation)
	#
	#		The default guess_root is 1.01 just because we expect that the roots we want have a real frequency that is always > 1 (you can see this from
	#			the Bohm-Gross relation, for example...) We'll see how good it works to see if we need to be any smarter...
	#
	# example:
	#					wp_e = 1.0						
	#					vth_e = sqrt( .1 )
	#					wavenumber = 1.0				# setup so that vth_e*wavenumber = sqrt(.1)  
	#														#		any other combinations will work.. only the product vth_e*wavenumber matters
	#														#		Note: vth_e*wavenumber is often refered to as k*lambda_debye... this is an apporximation
	#														#			that assumes that the wave frequency is very close in value to wp_e (w_wave ~= wp_e)
	#														#			A rule of thumb is that any value of k*lambda_debye greater then .33 or so is pushing it.. 
	#														#			so dont't go much higher.. [also notice that's not an appoximation when in the coefficent chi..
	#														#			that reduces to 1/ k*lambda_debye... but it is an appoximation inside the eponential. ]
	#
	#														# for vth_e*wavenumber = sqrt(.1), wp_e=1.0, the verifed answer is 1.179161-0.0184479j
	#
	#					print landau_damping_frequency( wp_e, vth_e, wavenumber)
	#														
	chi_e = np.power( ( wp_e/(vth__e*wavenumber) ), 2.0) / maxwellian_convention_factor
	def plasma_epsilon1( x ):
		val = 1.0 - chi_e*plasma_dispersion_prime( x )
		return val
	
		
	if inital_root_guess == None:
	#	# use the Bohm-Gross dispersion formulas to get an initial guess for w
		inital_root_guess = np.sqrt( wp_e*wp_e + 3*wavenumber*wavenumber*vth__e*vth__e)
	
	epsilon_root = scipy.optimize.newton( plasma_epsilon1, inital_root_guess )
	#epsilon_root = scipy.optimize.newton( plasma_epsilon1, 1.02)
	return epsilon_root * wavenumber * vth__e * np.sqrt(maxwellian_convention_factor)

# given a real w_wave, return the corresponding complex valued k.
def landau_spatial_damping( wp_e, vth__e, w_wave, maxwellian_convention_factor = 2.0, inital_root_guess = None ):
	exp_factor = w_wave / np.sqrt(maxwellian_convention_factor) / vth__e
	chi_factor = wp_e*wp_e/maxwellian_convention_factor / vth__e / vth__e
	
	def plasma_epsilon1( x ):
		pd = plasma_dispersion_prime( x )
		k_guess = exp_factor/ x
		chi_guess = chi_factor / pow( k_guess, 2)
		val = 1.0 - chi_guess*pd
		return val
	
	if inital_root_guess == None:
		# use the Bohm-Gross dispersion formulas to get an initial guess for w
		inital_root_guess = np.sqrt( np.power(w_wave,2.0) - np.power(wp_e,2.0) / ( 3.0*np.power(vth__e,2.0) ))
	
	epsilon_root = scipy.optimize.newton( plasma_epsilon1, inital_root_guess )	
	return exp_factor / epsilon_root

# from krall and trivvelpiece
def estimated_damping(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0):
	
	l_debye = vth__e / wp_e
	k_ldb = l_debye*wavenumber
	
	return np.sqrt( np.pi/8.0)*wp_e/np.power( k_ldb, 3.0)*np.exp( -1.0*(( 1.0/(2.0*k_ldb*k_ldb) ) + 1.5) )

def test__show_dispersion_relation( wp__e, vth__e, k_min, k_max, k_steps):
	(ks, w_real, w_imag) = test__calc_dispersion_relation(wp__e, vth__e, k_min, k_max, k_steps)
	
	k__vthe = np.copy( ks )
	k__vthe *= vth__e
	
	
	# also, calculate the Bohm Gross result/ approximate dmaping result, while we are here...
	# and keep track of the first value above lower_cutoff.. bcause my plasma-dispersion dolver goes off the rails there.
	cut_off_value = 1.05; cut_off_index = 0; 
	bohm_gross_freqs = np.zeros( k__vthe.shape[0] )
	krall_trivelpiece_damping = np.zeros( k__vthe.shape[0] )
	for i,kvthe in enumerate( k__vthe ):
		bohm_gross_freqs[i] = np.sqrt( wp__e*wp__e + 3*kvthe*kvthe)
		krall_trivelpiece_damping[i] = estimated_damping( wp__e, vth__e, ks[i])
		if kvthe > cut_off_value and cut_off_index == 0:
			cut_off_index = i
	
	clf()
	plot(k__vthe,w_real, label='Plasma Dispersion Calculation' )
	plot(k__vthe,bohm_gross_freqs, label='Bohm-Gross Relation' )
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.0, k__vthe[-1] ) ); ylim( ( 0.0, 2.0) );
	savefig('dispersion__real_freq__kinetic_electrostatic_dielectrix', dpi=200)
	
	clf()
	plot( k__vthe, w_imag*-1.0, label='Plasma Dispersion Calculation' )
	plot( k__vthe, krall_trivelpiece_damping, label='Krall-Trivelpiece approx. result' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '-1*Imag frequency (units of $w_p$)' )
	legend(loc=2, prop={'size':10})
	xlim( ( 0.0, k__vthe[-1] ) )
	axhline( y=0, color='black')
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix', dpi=200)
	
	# log versions...
	clf()
	plot( k__vthe, np.log( ( -1.0*w_imag ) ), label='Plasma Dispersion Calculation'  )
	plot( k__vthe, np.log( krall_trivelpiece_damping ), label='Krall-Trivelpiece approx. result' )
	axvline( x=.105, color='black', alpha=.5, label="Lower Limit Allowed in Code")
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '$log_{e}$  of -1*Imag frequency (units of $w_p$)' )
	xlim( ( 0.0, k__vthe[-1] ) )
	ylim( ( -40.0, 3.0) )
	axhline( y=0, color='black')
	legend(loc=4, prop={'size':10})
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__log', dpi=200)
	clf()
	plot( k__vthe, np.log10( ( -1.0*w_imag ) ), label='Plasma Dispersion Calculation' )
	plot( k__vthe, np.log10( krall_trivelpiece_damping ), label='Krall-Trivelpiece approx. result' )
	axvline( x=.105, color='black', alpha=.5, label="Lower Limit Allowed in Code")
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric (Damping)')
	xlabel( '$V_{thermal}*Wavenumber$ ( units of inverse skindepth)')
	ylabel( '$log_{10}$ of -1*Imag frequency (units of $w_p$)' )
	xlim( ( 0.0, k__vthe[-1] ) )
	ylim( ( -20.0, 1.0) )
	axhline( y=0, color='black')
	legend(loc=4, prop={'size':10})
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__logten', dpi=200)
	
	# another version.. I can match with come papers..
	clf()
	k__vthe2 = np.copy( k__vthe )
	np.power( k__vthe2 , 2.0, out = k__vthe2 )
	
	print k__vthe
	print k__vthe2
	plot(k__vthe2,w_real, label='Plasma Dispersion Calculation' )
	plot(k__vthe2,bohm_gross_freqs, label='Bohm-Gross Relation' )
	xlabel( '$V_{thermal}*Wavenumber$ squared ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.04, .12 ) ); ylim( ( 1.05, 1.25) );
	savefig('dispersion__real_freq__kinetic_electrostatic_dielectrix__squared_x', dpi=200)
	
	clf()
	plot( k__vthe2, np.log10( w_imag*-1.0 ), label='Plasma Dispersion Calculation' )
	plot( k__vthe2, np.log10(krall_trivelpiece_damping), label='Krall-Trivelpiece approx. result' )
	xlabel( '$V_{thermal}*Wavenumber$ squared ( units of inverse skindepth)')
	ylabel( 'Real frequency (units of $w_p$)' )
	title('Nominal Dispersion Relation\nElectrostatic (Kinetic) Dielectric')
	legend(loc=2, prop={'size':10})
	xlim( ( 0.04, .12 ) ); ylim( ( -5.0, -1.0) );
	savefig('dispersion__imag_freq__kinetic_electrostatic_dielectrix__squared_x', dpi=200)
	
def test__calc_dispersion_relation( wp__e, vth__e, v_the__k_min, v_the__k_max, k_steps):
	k_min = v_the__k_min/vth__e; k_max =  v_the__k_max/vth__e;
	ks = np.linspace(k_min, k_max, num=k_steps)
	w_reals = np.zeros( ks.shape[0]); w_imags = np.zeros( ks.shape[0]);
	print ks
	
	for i,k in enumerate(ks):
		print i
		print k
		w = landau_damping_frequency(wp__e, vth__e, k)
		w_reals[i] = np.real( w )
		w_imags[i] = np.imag( w )
	return (ks, w_reals, w_imags )
	
def plot_epilsion(wp_e, vth__e, wavenumber, maxwellian_convention_factor = 2.0):
	
	chi_e = np.power( ( wp_e/(vth__e*wavenumber) ), 2.0) / maxwellian_convention_factor
	def plasma_epsilon1( x ):
		val = 1.0 - chi_e*plasma_dispersion_prime( x )
		return val
	
	arg_const = 1.0 / (np.sqrt(2) * vth__e * wavenumber )
	bounds_real = np.linspace(2.547486, 2.547488, num=2000) 
	bounds_imag = np.linspace(-0.054607, -0.054609, num=2000)
	result = np.zeros(  (bounds_real.shape[0], bounds_imag.shape[0]) )
	reals = np.zeros( bounds_real.shape[0]*bounds_imag.shape[0] )
	imags = np.zeros( bounds_real.shape[0]*bounds_imag.shape[0] )
	
	
	#bounds_real *= arg_const
	#bounds_imag *= arg_const
	
	flat_idx = 0
	for ix in xrange( 0, bounds_real.shape[0] ):
		for iy in xrange( 0, bounds_imag.shape[0] ):
			#arg = arg_const * np.complex( bounds_real[ix], bounds_imag[iy] )
			arg =  np.complex( bounds_real[ix], bounds_imag[iy] )
			#print  arg
			rrr = plasma_epsilon1 ( arg )
			result[ ix, iy ] = np.absolute( rrr )
			reals[flat_idx] = np.real( rrr )
			imags[flat_idx] = np.imag( rrr )
			flat_idx += 1
	
	w = landau_damping_frequency(wp_e, vth__e, wavenumber)
	w *= arg_const
	wr = np.real(w); wi = np.imag(w);
	
	clf()
	imshow( np.log10( result) ,interpolation='nearest', extent=[bounds_real[0], bounds_real[-1],bounds_imag[0], bounds_imag[-1]], aspect='auto')
	axvline( x = wr); axhline( y = wi); 
	colorbar()
	savefig('aaaa', dpi=400)
	np.save( 'aaaaa', result)
	clf()
	plot( reals, imags)
	savefig('aaaa1')
	
	real_index = int(np.argmin(result) / bounds_imag.shape[0])
	imag_index = np.argmin(result) - real_index* bounds_imag.shape[0]  
	print imag_index
	print "PD min location:          %f + %fj    with value: %e" % (	wr,  	wi , np.absolute( plasma_epsilon1( w) ) )
	print "experiment min location:  %f + %fj    with value: %e" % ( bounds_real[real_index],  bounds_imag[imag_index] , result[real_index, imag_index ] )
	print "global min                                                 %e" % np.min(result)

#
# generate the famous figures from Freid and Conte (figures 1 through 8)
#		This routine should produce plots that are virtually identical to the ones in FandC
#		Note: requires matplotlib to be installed.
def test_plasma_dispersion():
	#
	# define z = x + j*y  (j = sqrt(-1))
	# now reproduce figures 1-10 in Fried and Conte.. 
	#		This code should identically reproduce those figures...
	#
	
	# do figures 1,2
	( x_values, pd__real_values, pd_prime__real_values, pd__imag_values, pd_prime__imag_values )  = test_plasma_dispersion__generate_values_for_a_given_y( 0.0)
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Re(Z')$ versus $x$ for $y=0$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Re(Z')$ ")
	plot( x_values, pd__real_values, label="$\\Re(Z)$")
	plot( x_values, pd_prime__real_values, label="$\\Re(Z')$")
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-2.0, .8) )
	ax.set_yticks(numpy.arange(-2.0,1.2,0.4))
	grid()
	savefig("pd_test__fig_01_real_y_0.png")
	
	clf()
	ax = subplot(111)
	title("$\\Im(Z)$ and $\\Im(Z')$ versus $x$ for $y=0$")
	xlabel("$x$")
	ylabel("$\\Im(Z)$,$\\Im(Z')$ ")
	plot( x_values, pd__imag_values,label="$\\Im(Z)$")
	plot( x_values, pd_prime__imag_values, label="$\\Im(Z')$")
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-1.6, 1.8) )
	xlim( (0.0, 5.0) )
	ax.set_yticks(numpy.arange(-1.6,2.0,0.4))
	grid()
	savefig("pd_test__fig_02_imag_y_0.png")
	
	
	# do figs 3,4,5,6
	(  x_values, pd_r__01, pdp_r__01, pd_i__01, pdp_i__01) = test_plasma_dispersion__generate_values_for_a_given_y(1.0)
	(  x_values, pd_r__02, pdp_r__02, pd_i__02, pdp_i__02) = test_plasma_dispersion__generate_values_for_a_given_y(2.0)
	(  x_values, pd_r__05, pdp_r__05, pd_i__05, pdp_i__05) = test_plasma_dispersion__generate_values_for_a_given_y(5.0)
	(  x_values, pd_r__10, pdp_r__10, pd_i__10, pdp_i__10) = test_plasma_dispersion__generate_values_for_a_given_y(10.0)
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$")
	plot( x_values, pd_r__01, label='y=1'); plot( x_values, pd_r__02, label='y=2'); plot( x_values, pd_r__05, label='y=5'); plot( x_values, pd_r__10, label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.425, 0.0) )
	ax.set_yticks(numpy.arange(-0.4,0.05,0.05))
	grid()
	savefig("pd_test__fig_03_real_z_positive_y.png")
	
	clf()
	ax = subplot(111)
	title("$\\Im(Z)$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Im(Z)$")
	plot( x_values, pd_i__01, label='y=1'); plot( x_values, pd_i__02, label='y=2'); plot( x_values, pd_i__05, label='y=5'); plot( x_values, pd_i__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (0.0, 0.8) )
	ax.set_yticks(numpy.arange(0.0,0.9,0.1))
	grid()
	savefig("pd_test__fig_04_imag_z_positive_y.png")
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$")
	plot( x_values, pdp_r__01,label='y=1'); plot( x_values, pdp_r__02,label='y=2'); plot( x_values, pdp_r__05,label='y=5'); plot( x_values, pdp_r__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.5, 0.2) )
	ax.set_yticks(numpy.arange(-0.5,0.3,0.1))
	grid()
	savefig("pd_test__fig_05_real_zprime_positive_y.png")

	clf()
	ax = subplot(111)
	title("$\\Im(Z')$ versus $x$ for $y=1,2,5,10$")
	xlabel("$x$")
	ylabel("$\\Im(Z')$")
	plot( x_values, pdp_i__01,label='y=1'); plot( x_values, pdp_i__02,label='y=2'); plot( x_values, pdp_i__05,label='y=5'); plot( x_values, pdp_i__10,label='y=10')
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-0.35, 0.0) )
	ax.set_yticks(numpy.arange(-0.35,0.05,0.05))
	grid()
	savefig("pd_test__fig_06_imag_zprime_positive_y.png")
	
	# do figs 7 and 8
	(  x_values, pd_r__m1, pdp_r__m1, pd_i__m1, pdp_i__m1) = test_plasma_dispersion__generate_values_for_a_given_y(-1.0)
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Im(Z)$ versus $x$ for $y=-1.0$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Im(Z)$ ")
	plot( x_values, pd_r__m1, label="$\\Re(Z)$"); plot( x_values, pd_i__m1, label="$\\Im(Z)$");
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (-7.0, 9.0) )
	ax.set_yticks(numpy.arange(-6.0,10.0,2.0))
	grid()
	savefig("pd_test__fig_07_real_and_imag_z__with_y_minus_1.png")	
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ and $\\Im(Z')$ versus $x$ for $y=-1.0$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$,$\\Im(Z')$ ")
	plot( x_values, pdp_r__m1, label="$\\Re(Z')$"); plot( x_values, pdp_i__m1, label="$\\Im(Z')$");
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-20.0, 10.0) )
	ax.set_yticks(numpy.arange(-20.0,12.0,4.0))
	grid()
	savefig("pd_test__fig_08_real_and_imag_zprime__with_y_minus_1.png")
	
	# do figs 9 and 10
	(  x_values, pd_r__m12, pdp_r__m12, pd_i__m12, pdp_i__m12) = test_plasma_dispersion__generate_values_for_a_given_y(-1.2)
	clf()
	ax = subplot(111)
	title("$\\Re(Z)$ and $\\Im(Z)$ versus $x$ for $y=-1.2$")
	xlabel("$x$")
	ylabel("$\\Re(Z)$,$\\Im(Z)$ ")
	plot( x_values, pd_r__m12, label="$\\Re(Z)$"); plot( x_values, pd_i__m12, label="$\\Im(Z)$");
	axhline( y=0, color='black')
	legend( loc=1)
	ylim( (-12.0, 16.0) )
	ax.set_yticks(numpy.arange(-12.0,20.0,4.0))
	grid()
	savefig("pd_test__fig_09_real_and_imag_z__with_y_minus_1.2.png")	
	
	clf()
	ax = subplot(111)
	title("$\\Re(Z')$ and $\\Im(Z')$ versus $x$ for $y=-1.2$")
	xlabel("$x$")
	ylabel("$\\Re(Z')$,$\\Im(Z')$ ")
	plot( x_values, pdp_r__m12, label="$\\Re(Z')$"); plot( x_values, pdp_i__m12, label="$\\Im(Z')$");
	axhline( y=0, color='black')
	legend( loc=4)
	ylim( (-40.0, 20.0) )
	ax.set_yticks(numpy.arange(-40.0,24.0,8.0))
	grid()
	savefig("pd_test__fig_10_real_and_imag_zprime__with_y_minus_1.2.png")	
	
def test_plasma_dispersion__generate_values_for_a_given_y(y):
	
	# define z = x + j*y  (j = sqrt(-1))
	
	# plot of real(p.d.), real(p.d. prime) versus x for y =0
	x_values = np.arange(0, 10.0, 0.1)
	pd__real_values = np.zeros( x_values.shape[0] )
	pd_prime__real_values = np.zeros( x_values.shape[0] )
	pd__imag_values = np.zeros( x_values.shape[0] )
	pd_prime__imag_values = np.zeros( x_values.shape[0] )
	
	for i, x in enumerate( x_values ):
		value = complex(x,y)
		if y == 0.0:
			value = x
		
		pd__real_values[i] =  np.real( plasma_dispersion( value ) )
		pd__imag_values[i] = np.imag( plasma_dispersion( value ) )
		pd_prime__real_values[i] =  np.real ( plasma_dispersion_prime( value ) )
		pd_prime__imag_values[i] =  np.imag ( plasma_dispersion_prime( value ) )
	return (  x_values, pd__real_values, pd_prime__real_values, pd__imag_values, pd_prime__imag_values ) 
	
def plot_landau_dmapings( vth_e, wp_e = 1.0):
	# plot from .01 kvth to 1.0
	k_start = .01 / vth_e
	k_end = .9 / vth_e
	ks = np.arange(k_start, k_end, .5)
	
	reals = np.zeros( ks.shape[0] ); imags = np.zeros( ks.shape[0] );
	for i,k in enumerate(ks):
		crap = landau_damping_frequency( wp_e, vth_e, k )
		reals[i] = np.real( crap)
		imags[i] = np.imag( crap)
	
	print "done searching landau space"
	
	ks *= vth_e
	
	clf()
	title( "landau reals")
	stem( ks, reals)
	savefig("landau_reals_vth_%f.png" % vth_e, dpi=200)
	clf()
	title( "landau imags")
	stem( ks, imags)
	savefig("landau_imags_vth_%f.png" % vth_e, dpi=200)

def landau_dmapings_find_closest_k( vth_e, freq, damping, wp_e = 1.0, k_fundimental=1.0):
	k_start = .01 / vth_e
	k_end = .9 / vth_e
	ks = np.arange(k_start, k_end, .01)
	
	
	real_epsilon = 1e10
	imag_epsilon = 1e10
	closest_k_real = -1
	closest_freq = -1
	closest_k_imag = -1
	closest_damping = -1
	
	for i,k in enumerate(ks):
		crap = landau_damping_frequency( wp_e, vth_e, k )
		reals = np.real( crap)
		imags = np.imag( crap)
		delta_real = np.fabs( reals - freq)
		if delta_real < real_epsilon:
			real_epsilon = delta_real
			closest_k_real = k
			closest_freq = reals
		delta_imag = np.fabs( imags - damping)
		if delta_imag < imag_epsilon:
			imag_epsilon = delta_imag
			closest_k_imag = k
			closest_damping = imags
	
	closest_k_real /= k_fundimental
	closest_k_imag /= k_fundimental
	
	print "freq: %f \t closest was %f \t at k=%f (kvth=%f) (error=%e) " % ( freq,closest_freq,closest_k_real,  (closest_k_real*vth_e), real_epsilon)
	print "damp: %f \t closest was %f \t at k=%f (kvth=%f) (error=%e) " % ( damping,closest_damping,closest_k_imag, (closest_k_imag*vth_e),   imag_epsilon)	
