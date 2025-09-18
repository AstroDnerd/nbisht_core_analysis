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
    this_simname='u302'

#directory = 'u05-r4-l4-128-Beta0.2  u10_r4_l4_128-Beta2  u11_r4_l4_128-Beta20'

center=np.array([0.89331055, 0.1159668 , 0.4440918 ])#u202 core0258 center
rad=5e-2
if 'master_clump' not in dir():
    frame = dl.target_frames[this_simname]
    ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
    test_main = ds.sphere(center, rad)

    h5name = "%s_%04d_peaklist.h5"%(this_simname,frame)
    master_clump = loop_tools.get_leaf_clumps(ds,small_test=False, h5_name = h5name, master_region=test_main) #,h5_name="NEW_PEAKS.h5" )
    leaves = loop_tools.get_peak_indices(master_clump, ds, h5name) #,h5_name="NEW_PEAKS.h5" )

if 1:
    proj = ds.proj('density',1,center=center,data_source=test_main)
    pw = proj.to_pw()
    pw.set_cmap('density','Greys')

    #pw.annotate_clumps([master_clump]+master_clump.leaves)
    pw.annotate_these_particles2(1.0,col='r',positions= test_main['particle_position'])
    pw.zoom(0.5/rad)
    pw.save('plots_to_sort/%s_clump_test_4'%this_simname)
    #            proj = ds.proj(field,ax,center=center, data_source = sph) 

if 0:
    fig,ax=plt.subplots(1,1)
    ax.scatter(test_main['radius'], test_main['density'])
    axbonk(ax,xscale='log',yscale='log',ylabel='density',xlabel='radius', xlim=[1/2048,test_main['radius'].max()])
    fig.savefig('plots_to_sort/test_density.png')

if 1:
    from yt.data_objects.level_sets.clump_handling import \
                Clump, \
                find_clumps, \
                get_lowest_clumps
    mountain_top = Clump(test_main,('gas','density'))
    peak_density = test_main['density'].max()
    min_density = 1
    step = 1.2
    nzones=(mountain_top['density']>min_density).sum()
    #find_clumps(mountain_top, min_density,peak_density,step)
    mountain_top.find_children( 0.05*peak_density, peak_density)
    pw.annotate_clumps(mountain_top.leaves)
    pw.save('plots_to_sort/%s_clump_test_5'%this_simname)

