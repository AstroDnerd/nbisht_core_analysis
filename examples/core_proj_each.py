from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
import projections
reload(projections)


#
# get mountain tops
#

import three_loopers_mountain_top as TLM

#
# Project tracks
#
if 1:
    if 'this_simname' not in dir():
        this_simname = 'u301'
    this_looper = TLM.loops[this_simname]

if 1:
    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( this_looper.core_list)
    projections.proj_cores_annotate_zoom(this_looper,
                                         #core_list=[10], 
                                         #frame_list=[0], 
                                         axis_list=[0],
                                         annotate_particles=False,
        field='density',color='r', zoom_level=4,plot_dir="./plots_to_sort", 
                                        )

