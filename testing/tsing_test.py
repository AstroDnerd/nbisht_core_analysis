

from starter2 import *
from collections import defaultdict
import scipy
import colors

from scipy.ndimage import gaussian_filter
import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 
from collections import defaultdict

import tsing
reload(tsing)
def cuml(ax,quan,color,label):
    the_x = nar(sorted(quan))
    the_y = np.arange( the_x.size)/the_x.size
    ax.plot( the_x, the_y, color=color,label=label)
sims=['u501', 'u502','u503']
#sims=[ 'u502','u503']
#sims=['u501']
sims=['u502']
if 1:
    for ns,sim in enumerate(sims):
        obj=tsing.te_tc(TL.loops[sim])
        obj.run(core_list=[74,368],make_plots=True)
