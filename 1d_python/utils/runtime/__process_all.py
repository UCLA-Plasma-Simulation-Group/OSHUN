import sys
sys.path.append("/u/home/mori/muffinlo/utils/trunk1d")

import numpy as np
import glob
import os
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
from pylab import *
from Oshun1d import Config, SpeciesConfig, ConstantProfile, EFieldSinusodialProfile


class DoResults:
	def __init__(self):
		self.values = {}


files = glob.glob('*')
valid_dirs = []
for f in files:
	if os.path.isdir(f):
		sub_glob = glob.glob('%s/*run.py' % f)
		if len(sub_glob) > 0:
			valid_dirs.append( f )

if len(valid_dirs) == 0:
	print "No valid dirs found!"
	exit(-1)


valid_dirs.sort()

# setup lists for aggregate data...
fit1_freq = []; fit1_gamma = [];
fit2_freq = []; fit2_gamma = [];
theo_freq = []; theo_gamma = [];
kvths = []
bgs		= []

fixed_k_vs_time = []

for idx,dir in enumerate( valid_dirs ):
	
	try:
		sim_metadata 		= pickle.load( open( "./%s/oshun_metadata.p" 		% (dir,), "rb" ) )
		sim_metadata_cfg 	= pickle.load( open( "./%s/oshun_metadata__cfg.p" 	% (dir,), "rb" ) )
		sim_results 		= pickle.load( open( "./%s/run_results.p" 			% (dir,), "rb" ) )
	except:
		print "Error getting results for result in directory %s" % dir
		continue
	
	print "Success in running results from directory %s" % dir
	
	
	# now we have the data..... do stuff......
	kvths.append( sim_results.values['k_vth'] )
	fit1_freq.append( sim_results.values['e_fit']['freq'] )

	theo_freq.append( sim_results.values['theo_fit']['freq'] )
	fit1_gamma.append( sim_results.values['e_fit']['gamma'] )
	
	theo_gamma.append( -1.0*sim_results.values['theo_fit']['gamma'] )
	bgs.append( np.sqrt( 1.0 + kvths[-1]*kvths[-1]*3.0 ) )
	
	try:
		fit2_freq.append( sim_results.values['fft_fit']['freq'] )
		fit2_gamma.append( sim_results.values['fft_fit']['gamma'] )
	except:
		pass
	if True:
		fixed_k_vs_time.append( sim_results.values['fixed_k_vs_time'] )
	#except:
	#	pass
	
clf()
plot( kvths, fit1_freq, label='fit 1' )

if len(fit2_freq) > 0:
	plot( kvths, fit2_freq, label='fit 2' )

plot( kvths, theo_freq, label='Theoretical (plasma dispersion)' )
plot( kvths, bgs, label='Bohm-Gross result' )
title('Dispersion Relation (real freq)')
xlabel('k*vth')
ylabel('Real freq')
legend( prop={'size':10} )
savefig('dispersion_all', dpi=200)

fit1_freq = np.array( fit1_freq ); theo_freq = np.array ( theo_freq )
fit1_gamma = np.array( fit1_gamma ); theo_gamma = np.array( theo_gamma )


clf()
title('Dispersion Relation (Imag freq)')
xlabel('k*vth')
ylabel('Damping Rate ( Imag freq )')
plot( kvths, fit1_gamma, label='fit 1' )
plot( kvths, theo_gamma, label='Theoretical (plasma dispersion)' )
legend( prop={'size':10} )
savefig('dispersion_damping_all', dpi=200)



v_th__vphase__fit1				= np.zeros( fit1_freq.shape[0] )
v_th__vphase__theo				= np.zeros( fit1_freq.shape[0] )
vth_e							= sim_metadata_cfg.species_config[0].temp_profile.val

for i in xrange(0, fit1_freq.shape[0]):
	k = kvths[i] / vth_e
	vphase_fit1 = 				fit1_freq[i]/k
	val_fit1 =					vth_e / vphase_fit1
	val_fit1 *= 				val_fit1
	
	vphase_theo = 				theo_freq[i] / k
	val_theo = 					vth_e / vphase_theo
	val_theo *= 				val_theo
	
	v_th__vphase__fit1[i] = val_fit1
	v_th__vphase__theo[i] = val_theo
fit1_gamma /= fit1_freq
theo_gamma /= theo_freq


clf()
plot( np.log10( v_th__vphase__fit1), np.log10( fit1_gamma) , label='fit 1' )
if len(fit2_gamma) > 0:
	plot( kvths, fit2_gamma, label='fit 2' )
plot( np.log10( v_th__vphase__theo), np.log10( theo_gamma), label='Theoretical (plasma dispersion)'  )
title('Dispersion Relation (Imag freq)')
xlabel('(vth / vphase)^2')
ylabel('Damping Rate/freq ( Imag freq / reall freq )')
legend( prop={'size':10} )
ylim( (-4, 0) )
xlim( ( -3, 0) )
savefig('dispersion__damping_experiment_all', dpi=200)


# print all fft of fixed ks... (for driver tunring)
clf()
for i in xrange(0, len(fixed_k_vs_time)):
	plot( fixed_k_vs_time[i], label='set %02d' % i)
legend(  prop={'size': 8} )
savefig('fft_of_fixed_k', dpi=200)

print fit1_freq




