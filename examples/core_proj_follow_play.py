from starter2 import *
import colors
reload(colors)
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)

import three_loopers_1tff as tl

this_looper = tl.looper2
this_looper.plot_directory = "/home/dccollins/PigPen"

color_dict = colors.make_core_cmap(this_looper.core_list)
#loop_apps.core_proj_follow(this_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#import shift_snaps
#shift_snaps.shift_snaps(this_looper)
loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', 
                             core_list=[258], 
                             #frame_list=[0],
                             slab=False, 
                             zoom=True, 
                             only_sphere=True, center_on_sphere=True, annotate=False, p_size=7)
                           #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
