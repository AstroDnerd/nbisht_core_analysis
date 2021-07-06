from starter2 import *
from yt.visualization.plot_modifications import *

#from yt.data_objects.level_sets.clump_handling import \
#            Clump, \
#            find_clumps, \
#            get_lowest_clumps
import yt.data_objects.level_sets.clump_handling as clump_handling
reload(clump_handling)

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
frame = 0 #dl.target_frames[this_simname]
ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
h5name = "%s_%04d_peaklist.h5"%(this_simname,frame)
if 'master_clump' not in dir():
    def get_potential_contours(ds,c_min=None,c_max=None,step=100,h5_name="NEW_PEAK_FILE.h5",pickle_name=None, 
                         subset=None, peak_radius=1.5,bad_particle_list=None, small_test=False):
        """get all the leaf indices for peaks in *ds*.
        If *pickle_name* is supplied, load from that, or if it doesn't exist, save to that.
        *subset*, if supplied, will restrict the indices returned.
        """
        if small_test:
            #center = ds.arr([0.07104492, 0.05688477, 0.1862793 ],'code_length')
            peak,center=ds.find_max('PotentialField')
            ad = ds.sphere(center,0.1)
        else:
            ad = ds.all_data()
        #ad  = ds.sphere([0.52075195, 0.74682617, 0.01196289], 0.1)
        master_clump = clump_handling.Clump(ad,'PotentialField')
        master_clump.add_validator("min_cells", 8)
        c_min = ad['PotentialField'].min()
        #c_max = 534069645. # ad["gas", "density"].max()
        c_max = ad['PotentialField'].max()

        step = (c_max-c_min)/20
        print(step)
        clump_handling.find_clumps_linear(master_clump, c_min, c_max, step)
# Write a text file of only the leaf nodes.
        #write_clumps(master_clump,0, "%s_clumps.txt" % ds)

        return master_clump
    master=get_potential_contours(ds)
proj = yt.ProjectionPlot(ds,0,'PotentialField')
proj.annotate_clumps(master.leaves)
proj.save('plots_to_sort/%s_clumps'%this_simname)
#master_clump = loop_tools.get_leaf_clumps(ds,small_test=False, h5_name = h5name) #,h5_name="NEW_PEAKS.h5" )
#leaves = loop_tools.get_peak_indices(master_clump, ds, h5name) #,h5_name="NEW_PEAKS.h5" )
