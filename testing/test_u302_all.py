from starter2 import *

import mountain_top
reload(mountain_top)
reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#
mountain_top_fname = "u302_all_test1_mountain_tops.h5"
read_write_fname = "L258_all_test1_trackage.h5"
this_simname = 'u302'

if 0:
    def verify_cores_very_dumb(leaf, peak_density):
        if peak_density < 1100:
            return False
        else:
            return True
    for this_simname in ['u302']:#['u301','u303']:# ['u302']:
        kludge={}
        #kludge={'peak_id':[0,1,364, 113,u302_was258]})
        #kludge={'peak_id':[258]}
        mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_fname, do_projections=False, verify=verify_cores_very_dumb,kludge=kludge)
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
    Lfull = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False,
                                     plot_directory = "./plots_to_sort"
                                  )
    Lfull.read_targets(mountain_top_fname)
    Lfull.get_tracks()


#
# Project tracks
#
if 1:

    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( Lfull.core_list)
    loop_apps.core_proj_multiple(Lfull,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)


if 0:
    #read_write tests
    Lfull.save_trackage_only(read_write_fname)

if 0:
    Lfull_read = looper.core_looper(directory=dl.sims[this_simname],savefile=  read_write_fname)
    Lfull_read.out_prefix+='_read'

if 0:

    import colors
    reload(colors)
    reload(loop_apps)
    core_cmap = colors.make_core_cmap( Lfull_read.core_list)
    loop_apps.core_proj_multiple(Lfull_read,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)
