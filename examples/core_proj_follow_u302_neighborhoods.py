from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)

if 0:
    this_simname = 'u301'
    sim_list = ['u301']
    if 'set_looper' not in dir():
        savefile='u301_long_pos_only.h5'
        set_looper=looper.core_looper(directory= dl.sims[this_simname],savefile_only_trackage=savefile)
        thtr = set_looper.tr
        set_looper.out_prefix='XXXX'
        thtr.sort_time()
    set_looper.plot_directory = "./plots_to_sort"
if 1:
    import three_loopers_tenfour as TL4
    set_looper = TL4.loops['u401']
sim_list = ['u401']
import convex_hull_tools as CHT
reload(CHT)
if 'ht_set' not in dir() :
    ht_set = {}
    for this_simname in sim_list:
        ht_set[this_simname] = CHT.hull_tool(set_looper)
        ht_set[this_simname].make_hulls()
        ht_set[this_simname].make_overlaps()

import supersets
reload(supersets)
if 'st_set' not in dir():
    st_set={}
    for this_simname in sim_list:
        st_set[this_simname] = supersets.superset( set_looper, ht_set[this_simname])
        st_set[this_simname].find()

import colors
import core_proj
reload(core_proj)
for nset,superset in enumerate(st_set[this_simname].supersets):
    if nset != 1:
        continue

    core_list = list(superset)
    color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
#loop_apps.core_proj_follow(set_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#import shift_snaps
#shift_snaps.shift_snaps(set_looper)

    frame_list = set_looper.frame_list[1:-1]
    set_looper.out_prefix = '%s_S%02d'%(this_simname,nset)
    mono=core_proj.core_proj_multiple(set_looper,axis_list=[0], 
                                 color_dict=color_dict,
                                 frame_list = frame_list,
                                 field='density', 
                                 core_list=core_list,
                                 slab=False, zoom=True, only_sphere=True, center_on_sphere=True, 
                                 annotate=True,  plot_particles=True,
                                 grids=False, float_positions=True, monotonic=False, verbose=True,
                                      path_only=False, 
                                     plot_y_tracks=True, plot_points=False)
                               #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
