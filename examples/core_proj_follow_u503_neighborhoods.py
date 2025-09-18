from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)
this_simname = 'u503'
import three_loopers_u500 as TL5
set_looper = TL5.loops[this_simname]
sim_list = ['u503']

if 0:
    import convex_hull_tools as CHT
    reload(CHT)
    if 'ht' not in dir() :
        ht = {}
        for this_simname in sim_list:
            ht[this_simname] = CHT.hull_tool(set_looper)
            ht[this_simname].make_hulls()
            ht[this_simname].make_overlaps()

    import supersets
    reload(supersets)
    if 'st' not in dir():
        st={}
        for this_simname in sim_list:
            st[this_simname] = supersets.superset( set_looper, ht[this_simname])
            st[this_simname].find()
    do_supersets=True
else:
    do_supersets=False


import colors
import core_proj
reload(core_proj)
if do_supersets:
    stuff=st[this_simname].supersets
    for nset,superset in enumerate(stuff):
        if nset != 1:
            continue

        core_list = list(superset)
        color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
        frame_list = set_looper.tr.frames[::10]
        set_looper.out_prefix = '%s_full_S%02d'%(this_simname,nset)
else:
    core_list = np.unique(set_looper.core_ids)[:2]
    frame_list = set_looper.tr.frames[::10]
#frame_list=set_looper.tr.frames #[::10]
if 1:
    import movie_frames
    reload(movie_frames)
    movie_mask = movie_frames.quantized_mask(set_looper)
    frame_list=set_looper.tr.frames[movie_mask]
core_list = [52,53,55,56,57,58,74,76,77,78]
color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
color_dict = dict([(d,'r') for d in core_list])
import camera_path
reload(camera_path)

camera = camera_path.camera_1( set_looper, 'tight_float')
#camera = camera_path.camera_1( set_looper, 'domain')
#camera = camera_path.camera_1( set_looper, 'smooth_zoom')
frame_list=frame_list[-1:]
print('hey')
if 1:
    print('wut')
    mono=core_proj.core_proj_multiple(set_looper,camera,axis_list=[2], 
                                 color_dict=color_dict,
                                 frame_list = frame_list,
                                 field='density', 
                                 core_list=core_list,
                                 slab=False, zoom=True, only_sphere=False, center_on_sphere=True, 
                                 annotate=True,  plot_particles=False,
                                 grids=False, float_positions=True, monotonic=True, 
                                      verbose=True,
                                      path_only=False, 
                                     plot_y_tracks=True, plot_points=False)
                               #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])
