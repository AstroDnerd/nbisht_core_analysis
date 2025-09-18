from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)



this_simname = 'u301'

if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()

set_looper.plot_directory = "./plots_to_sort"

core_list = nar(ht1.cores_used)[ nar(ht1.overlaps[31])>0]
popper  = set([124,108, 66,65,62,61])
core_set = set(core_list)
core_list = list(core_set.difference(popper))
reload(colors)
color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
#loop_apps.core_proj_follow(set_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#import shift_snaps
#shift_snaps.shift_snaps(set_looper)
loop_apps.core_proj_multiple(set_looper,axis_list=[0], 
                             color_dict=color_dict,
                             field='density', 
                             core_list=core_list,
                             slab=False, zoom=True, only_sphere=False, center_on_sphere=True, 
                             annotate=True, frame_list = list(range(121)) + [
                            grids=False, float_positions=True, monotonic=True)
                           #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
