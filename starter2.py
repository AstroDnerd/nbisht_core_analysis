
#
# ALl the other things that aren't standard python packages.
#


import os
import sys
if not os.path.exists('dtools/starter1.py'):
    print("Need to get the submodule")
    print("git submodule init")
    print("git submodule update")
    sys.exit(-1)
from dtools.starter1 import *

#so we can import things from sub directories
path_list= ["./tools_data", "./tools_pdf", "./testing",\
            "./tools", "./tools_tracks", "./trash","./new_targets", "./p56_plots",
            "./other_particles", "./browser", "./sink_analysis"]
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]

import yt
from yt_names import *

import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
import tools_data.particle_ops as particle_ops
import tools_data.particle_grid_mask as particle_grid_mask
import looper2  

#import looper
#reload(looper)
import trackage
import tracks_read_write
import tracks_read_write as trw
from davetools import *
#import loop_apps

import data_locations as dl
import track_info
import track_loader

from collections import defaultdict
import colors
plot_dir = "%s/cdbreak_desktop/nikhilb_home/results"%os.environ['HOME']

def make_dir(dir_path):
        if not os.path.exists(dir_path):
                print("making directory:",dir_path)
                os.makedirs(dir_path)