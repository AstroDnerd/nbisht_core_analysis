
#
# ALl the other things that aren't standard python packages.
#

from starter1 import *
#so we can import things from sub directories
path_list= ["./tools_data", "./tools_pdf", "./testing",\
            "./tools", "./tools_tracks", "./trash","./new_targets", "./p56_plots",
            "./other_particles", "./browser"]
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]

import yt
from yt_names import *


import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
import particle_ops as particle_ops
import particle_grid_mask as particle_grid_mask

#import looper
#reload(looper)
import looper2
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
