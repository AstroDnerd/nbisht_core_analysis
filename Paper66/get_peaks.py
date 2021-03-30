"""
Getting peaks and particles on new sims.
All of this process should get streamlines.
1.)  In data_locations, fill out
    a.) target_frames
    b.) sims
    c.) peak_list (this will get made by get_peaks.py)
2.) Run get_peaks
3.) Move peaklist.h5 to its home.
4.) Run count_particles.py
5.) Make sure ??_n_particles.txt gets put in datasets_small.py
6.) Put ??_n_particless.txt in data_locations
6.) Now you should be mostly set.
7.) Run data_puller.
8.) Move the h5 files to /data/cb1/Projects/P19_CoreSimulations/CoreSets/
"""


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
    this_simname='u202'

#directory = 'u05-r4-l4-128-Beta0.2  u10_r4_l4_128-Beta2  u11_r4_l4_128-Beta20'
frame = dl.target_frames[this_simname]
ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
h5name = "%s_%04d_peaklist.h5"%(this_simname,frame)
master_clump = loop_tools.get_leaf_clumps(ds,small_test=False, h5_name = h5name) #,h5_name="NEW_PEAKS.h5" )
leaves = loop_tools.get_peak_indices(master_clump, ds, h5name) #,h5_name="NEW_PEAKS.h5" )
