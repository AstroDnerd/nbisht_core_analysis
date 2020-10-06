
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
directory_10 = '/archive2/luzlourdes/u10/DD0082/data0082'
directory_11 = '/archive2/luzlourdes/u11/DD0088/data0088'

ds = yt.load('%s'%(directory_10))

new_indices = loop_tools.get_leaf_indices(ds,h5_name='u10_peaklist.h5')

print("new .h5 file should be created")

# continue smooching from my get_peaks.py on Stampede.
# but first make sure that the data_puller.py type of process
# is indeed using the target_frame to get the peaks.....




