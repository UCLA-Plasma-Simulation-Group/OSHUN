
import math


def heaviside_step(x):
	if x< 0:
		return 0.0
	return 1.0;

	# This is an ultra cheapo units implementation... atleast all these formulas
	# will be in place so that as I find mistakes, at least they are all fixed at this central
	# place.
	#
	# And I have centeralized place to put successes until I make a more unified thingy

use_higher_percision = False

pi 							= 3.1415926535897932

elementaryCharge__c			= 1.602176565e-19
elementaryCharge__statc		= 4.80320425e-10

electronMass__kg 			= 9.10938291e-31
electronMass__g 			= 9.10938291e-28

protonMass__kg				= 1.672621777e-27
protonMass__g				= 1.672621777e-24

protonElectronMassRatio		= protonMass__kg/electronMass__kg

c__m_s						= 299792458
c__cm_s						= 29979245800

Kb__j_k						= 1.3806488e-23
kb__j_k						= 1.3806488e-23
boltzmannConstant__j_k		= 1.3806488e-23
Kb__erg_k					= 1.3806488e-16
kb__erg_k					= 1.3806488e-16
boltzmannConstant__erg_k	= 1.3806488e-16

epsilon0__F_m				= 8.854187817e-12
permittivityFreespace__F_m 	= 8.854187817e-12
dielectricConstant__F_m 	= 8.854187817e-12

mu0__H_m					= 4.0*pi*1e-7
permeabilityConstant__H_m	= 4.0*pi*1e-7
magneticConstant__H_m		= 4.0*pi*1e-7

erg__to__j					= 1e-7
j__to__erg					= 1e7

def calc_log_lambda(temp, density):
	debye_length_e = 7.43e2*math.sqrt(temp)/math.sqrt(density)
	classical_closest_approach = 1.44e-7/math.sqrt(temp)
	deBrogilieLength = 2.76e-8/math.sqrt(temp)
	
	bmin = classical_closest_approach
	if classical_closest_approach < deBrogilieLength:
		bmin = deBrogilieLength
	
	llambda = debye_length_e/bmin
	log_lam = math.log(llambda)
	if log_lam < 2.0:
		log_lam = 2.0
	return log_lam
	
def electron_mfp__cm(temp_e_in_ev, density, z):
	log_lam = calc_log_lambda(temp_e_in_ev, density)
	mfp = (4.0/z) * (temp_e_in_ev / 1000.0)*(temp_e_in_ev / 1000.0)*(10.0/log_lam)*(1e22/density) * .37*1e-4
	return mfp
def electron_ion_collision_wo_laser(temp_e_in_ev, density, z):
	log_lam = calc_log_lambda(temp_e_in_ev, density)
	return 3e-6 * log_lam * density * z / (temp_e_in_ev * math.sqrt(temp_e_in_ev))
	
	
def ion_collision_freq(temp_i_in_ev,ion_density, Z, mu=1):
	return 4.80e-8*Z*Z*Z*Z*ion_density*calc_log_lambda(temp_i_in_ev,ion_density)/(temp_i_in_ev*math.sqrt(temp_i_in_ev)*math.sqrt(mu))
def electron_collision_freq(temp_e_in_ev,electron_density):
	return 2.91e-6*electron_density*calc_log_lambda(temp_e_in_ev,electron_density)/(temp_e_in_ev*math.sqrt(temp_e_in_ev))
def wp(density, mass):
	constant = 5.64e4
	if use_higher_percision:
		constant = math.sqrt( 4.0*pi*elementaryCharge__statc*elementaryCharge__statc/electronMass__g)	
	return constant*math.sqrt(density)/math.sqrt(mass)
def wpe(density):
	return wp(density, mass=1.0)
def wpi(density, number_protons, charge):
	if use_higher_percision:
		return wp(density, mass=number_protons*protonElectronMassRatio)*charge
	return 1.32e5*charge*sqrt(density)/sqrt(mass)
def wpp(density):
	if use_higher_percision:
		return wp(density, mass=protonElectronMassRatio)
	return wpi( density, 1.0, 1.0)
def classical_closest_approach(temp_ev):
	return 1.44e-7/math.sqrt(temp_ev)
def debye_length(temp_ev, density):
	constant = 7.43e2
	if use_higher_percision:
		constant = math.sqrt( (elementaryCharge__c*j__to__erg)/(4.0*pi*elementaryCharge__statc*elementaryCharge__statc) )
	return constant*math.sqrt(temp_ev)/math.sqrt(density)
def deBroglie_length(temp_ev):
	return 2.76e-8/math.sqrt(temp_ev)
# number of particles in Debye shpere (spmetimes called plasma parameter)
def plasma_parameter(temp_ev, density):
	# here, the constant is the one for the debye_length, taken to the 3rd power, and multiplied by 4*pi/3
	constant = 1.72e9
	if use_higher_percision:
		constant = (elementaryCharge__c*j__to__erg)/(4.0*pi*elementaryCharge__statc*elementaryCharge__statc)
		constant *= math.sqrt(constant)
		constant *= 4.0*pi/3.0
	return constant*temp_ev*math.sqrt(temp_ev)/math.sqrt(density)
	# could also do this in another way:
	#db_length = debye_length(temp_ev,  density)
	#return db_length*db_length*db_length*density*4*pi/3.0

"""
convert_ev_to_vth
	temp_in_ev: 		temperature in ev
	mass:				mass normalized to the electron mass (i.e. for electron this is 1.0; for a proton this is protonElectronMassRatio ~= 1860)
	proper_velocity:	If True, output is a relativistic proper velocity
	normalized:			If True, output is a unitless velocity that is is in units o the speed of light. 
						If False, the output is in units of cm/sec
	If both proper_velocity and normalized are True, then the output is a relativistic proper velocity normalized to speed of light (i.e. the unitless quantity beta*gamma)
"""
def convert_ev_to_vth(temp_in_ev, mass, proper_velocity=True, normalized=True):	
	# In this constant, 4.19e7 = sqrt(e/Me * 10e7)   (e: charge of electron in coulmbs, Me:mass electron in grams, 1e7 conversion of joules to ergs)
	# 	temp_in_ev*e gives energy in joules, j__to__erg converts this in energy in ergs (cgs), then dividing by electron mass in grams gives (velocity in cm/sec)^2
	constant = 4.19e7/3e10
	if( use_higher_percision ):
		constant = math.sqrt( (elementaryCharge__c/electronMass__g) * j__to__erg )
		constant /= c__cm_s
	beta = constant*math.sqrt(temp_in_ev)/math.sqrt(mass)
	gamma = 1.0/math.sqrt(1.0-beta*beta)
	v = beta
	if not normalized:
		if( use_higher_percision ):
			v *= c__cm_s
		else:
			v *= 3e10		
	if proper_velocity:
		v *= gamma
	return v

def convert_ev_to_vth_cgs(temp_in_ev, mass=1.0):
	vth_in_in_cgs = convert_ev_to_vth(temp_in_ev, mass)
	if use_higher_percision:
		vth_in_in_cgs *= c__cm_s
	else:
		vth_in_in_cgs *= 3e10
	return vth_in_in_cgs

def convert_osi_vth_to_beta(temp_in_vth, mass=1.0):
	gamma = math.sqrt(1.0+temp_in_vth*temp_in_vth)
	# wrong, I think
	#return temp_in_vth/math.sqrt((temp_in_vth*temp_in_vth +1)*mass)
	return temp_in_vth/gamma
	

"""
convert_vth_to_ev
	temp_in_vth:		thermal velocity
	mass:				mass normalized to the electron mass (i.e. for electron this is 1.0; for a proton this is protonElectronMassRatio ~= 1860)

"""
def convert_vth_to_ev(temp_in_vth, mass=1.0, proper_velocity=True, normalized=True):
	# see convert_ev_to_vth for dicussion of the conversion constant. The constant here is 1.0/(constant in convert_ev_to_vth).
	constant = 1.0/4.19e7
	if( use_higher_percision ):
		constant = 1.0/math.sqrt( (elementaryCharge__c/electronMass__g) * j__to__erg )
	if normalized:
		if( use_higher_percision ):
			constant *= c__cm_s
		else:
			constant *= 3e10
	if proper_velocity:
		gamma = math.sqrt(1.0+temp_in_vth*temp_in_vth)
		ev1 = (temp_in_vth/gamma)*constant
	else:
		ev1 = (temp_in_vth)*constant
	return ev1*ev1*mass

#---- basic laser formulas!!!!!
def laser_group_speed__sim(denisty__sim):
	return math.sqrt(1-denisty__sim)
def a0_change_from_vaccum(density__sim):
	return 2.0/(1+math.sqrt(1.0-density__sim))
def electron_ion_collisions_small_laser(temp_ev,density,ion_charge,w_laser):
	log_lamda = modified_log_lambda(temp_ev, density, w_laser)
	return 3e-6 * log_lamda * density * ion_charge / (temp_ev * math.sqrt(temp_ev))
# This is the modified log lambda as decribed by Kruer
def modified_log_lambda(temp_ev, density, w_laser):
	
	b_max = convert_ev_to_vth_cgs(temp_ev, 1.0)/w_laser
	closest_approach = classical_closest_approach(temp_ev)
	deBrogilieLength = deBroglie_length(temp_ev)
	b_min = closest_approach
	if(closest_approach < deBrogilieLength):
		b_min = deBrogilieLength
	log_lam = math.log(b_max/b_min)
	if (log_lam < 2.0):
		log_lam = 2.0
	return log_lam
	

# normalized stuffs

# normalized debye length:
#
#	v_thermal= l_debeye* w_p      (nice indentity - strictly in real units for easier arguemnt)
#	l_debye = v_thermal/w_p
#
#   now we want to take simension out of this equaiton.... we will normalize time to 1/w_0 (w_o some frequency
#   		of our choosing (usualy the laser frequency or the plasma frequency but can be anything).
#			then length is normed to c/w_0. So the the debye length l_debye' in this system of units is:
#	
#	l_debye' = l_debye/(c/w_0) = (1/(c/w_0))*vthermal/w_p = w_0/w_p*v_thermal/c   (so no units on wither side!)
#			we know w_0/w_p = sqrt(n_0/n_p) = sqrt(1/density__norm)
#	l_debye' = sqrt(n_norm/n_p)*beta_vth   (note ususally n_norm is 1.0 for w_0 laser.. but it is the endensity that goes with w_0.)
def debye_length_norm(vth_osiris, density__norm, mass=1.0):
	ev = convert_osi_vth_to_ev(vth_osiris)
	print "\t\tdebug: ev: %f" % ev
	cgs = convert_ev_to_vth_cgs(ev)
	print "\t\tdebug: cgs: %f" % cgs
	beta1 = cgs / 3e10
	print "\t\tdebug: beta1: %f" % beta1
	beta = convert_osi_vth_to_beta(vth_osiris, mass=mass)
	print "\t\tdebug: beta: %f" % beta
	return beta/math.sqrt(density__norm)

# container to hold these formula so that we won't pollute namespaces...
class Formulary:
	def __init__(self, use_higher_percision_constants = False):
		global use_higher_percision
		use_higher_percision = use_higher_percision_constants
		self.hi = 1
		self.wp = wp
		self.wpe = wpe
		self.wpp = wpp
		self.calc_log_lambda = calc_log_lambda
		self.electron_mfp__cm = electron_mfp__cm
		self.deBroglie_length = deBroglie_length
		self.ion_collision_freq = ion_collision_freq
		self.electron_collision_freq = electron_collision_freq
		self.debye_length = debye_length
		self.plasma_parameter = plasma_parameter
		self.convert_ev_to_vth = convert_ev_to_vth
		self.convert_ev_to_vth_cgs = convert_ev_to_vth_cgs
		#self.vth_into_ev = vth_into_ev
		self.convert_vth_to_ev = convert_vth_to_ev
		self.convert_osi_vth_to_beta = convert_osi_vth_to_beta
		self.laser_group_speed__sim = laser_group_speed__sim
		self.electron_ion_collisions_small_laser = electron_ion_collisions_small_laser
		self.modified_log_lambda = modified_log_lambda
		self.electron_ion_collision_wo_laser = electron_ion_collision_wo_laser
		
		
		# a temporary spot for normalized formulas....
		self.debye_length_norm = debye_length_norm

formulary = Formulary(use_higher_percision_constants=True)

		