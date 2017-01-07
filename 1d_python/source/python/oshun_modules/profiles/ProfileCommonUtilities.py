# hi
#import oshun_modules.OshunGlobalState as GLOBAL
from ..OshunGlobalState import *

def createGlobalProfileArray( ):
	return np.zeros( STATE.spatial_axes__global[0].num_points )

def createLocalProfileArray( species ):
	return np.zeros( species.F.num_total_cells[0] )
