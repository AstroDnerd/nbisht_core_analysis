from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

import three_loopers_mountain_top as TLM

#
# Project tracks
#
if 1:
    this_looper = TLM.loops['u301']

    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( this_looper.core_list)
    loop_apps.core_proj_multiple(this_looper,core_list=[10],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)
