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

for idx,dir in enumerate( valid_dirs ):
	
	try:
		os.chdir( dir )
		os.system("python /u/home/mori/muffinlo/utils/trunk1d/oshun1d_ld_e_vs_time__standing_wave.py")
	except:
		print "Error getting results for result in directory %s" % dir
		continue
	os.chdir("..")
	print "Success in running results from directory %s" % dir
	
	
print "Done."

