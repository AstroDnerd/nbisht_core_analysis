from starter2 import *


import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array

from importlib import reload

import looper
reload(looper)
import loop_tools
reload(loop_tools)

if 'this_simname' not in dir():
    this_simname='u101'

#directory = 'u05-r4-l4-128-Beta0.2  u10_r4_l4_128-Beta2  u11_r4_l4_128-Beta20'
frame = dl.target_frames[this_simname]
ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
h5name = "%s_%04d_peaklist.h5"%(this_simname,frame)
master_clump = loop_tools.get_leaf_clumps(ds,small_test=False, h5_name = h5name) #,h5_name="NEW_PEAKS.h5" )
leaves = loop_tools.get_peak_indices(master_clump, ds, h5name) #,h5_name="NEW_PEAKS.h5" )
