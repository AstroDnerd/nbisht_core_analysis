from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

this_simname = 'u302'
core_id = 258
mountain_top_name = "%s_mountain_top_core_c%04d.h5"%(this_simname,core_id)
do_mountain_projections=True
if 1:
    def verify_cores_very_dumb(top1):
        if top1.rhomax < 1100:
            return False
        else:
            return True
    outname = "u302_tracks_core_258.h5"
    #kludge={'peak_id':[0,1,364, 113,u302_was258]})
    #kludge={'peak_id':[258]}
    kludge={'peak_id':[core_id]}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, verify=verify_cores_very_dumb,kludge=kludge)
#
# Get tracks.  
#

if 1:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list = [target_frame] #list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    derived=[]

if 1:
    L258 = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname ,
                                     out_prefix = this_simname + "_c%04d"%core_id,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False,
                                     plot_directory = "./plots_to_sort"
                                  )
    print("READ TARGETS", mountain_top_name)
    L258.read_targets(mountain_top_name)
    L258.get_tracks()


#
# Project tracks
#
if 1:

    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( L258.core_list)
    loop_apps.core_proj_multiple(L258,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[0], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)


if 0:
    #read_write tests
    read_write_fname = "L258_trackage_only.h5"
    L258.save_trackage_only(read_write_fname)

if 0:
    L258_read = looper.core_looper(directory=dl.sims[this_simname],savefile=  read_write_fname)
    L258_read.out_prefix+='_read'

if 0:

    import colors
    reload(colors)
    reload(loop_apps)
    core_cmap = colors.make_core_cmap( L258.core_list)
    loop_apps.core_proj_multiple(L258_read,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)
