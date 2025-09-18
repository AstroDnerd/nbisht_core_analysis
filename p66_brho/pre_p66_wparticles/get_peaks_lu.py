
# Borrowed from Brho_particles66.py

from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
# - - - - - - - - - -

# NAZARE-an
directory_05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128'
directory_10 = '/archive2/luzlourdes/u10/DD0082/data0082'
directory_11 = '/archive2/luzlourdes/u11/DD0088/data0088'

# PEAKLIST files created: 
# u05 : 338  
# u10_082_peaklist.h5 : 256
# u11_088_peaklist.h5 : 365

ds = yt.load('%s'%(directory_11))

# CREATE .h5 peaklist for new data sets
if 0:
    new_indices = loop_tools.get_leaf_indices(ds,h5_name='u10_peaklist.h5')
    print("new .h5 file should have been created")


# CREATE .txt file of cores and # of particles 
if 1:
    coros = np.zeros(365,dtype=int)
    cores = loop_tools.get_leaf_indices(ds,h5_name='datasets_small/u11_088_peaklist.h5')

    for i in range(365): 
        coros[i] = len(cores.get(i,0)) 
        if i == 0:
            text = open("u11c_particles.txt", "w")
            text.write("\n%d %d"%(i,coros[i]))
            text.close
        if i > 0:
            text = open("u11c_particles.txt", "a")
            text.write("\n%d %d"%(i,coros[i]))
            text.close

