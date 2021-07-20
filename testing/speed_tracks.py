#
# run get_mountain_tops first
#
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)

if 'this_simname' not in dir():
    this_simname = 'u302'

if 1:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list = [target_frame] #list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    derived=[]

disk_or_new = 'new'
if 'new_looper' not in dir() and disk_or_new:
    reload(looper)
    #output_name = "cores_%s_mountain34_take1.h5"%this_simname
    #target_fname = "%s_new_target_34.h5"%this_simname

    #output_name = "cores_%s_core14.h5"%this_simname
    #target_fname = "%s_core_14.h5"%this_simname
    #target_fname = "%s_core_0-20.h5"%this_simname
    #outname = 'u302_tracks_only_c0-20.h5')
    #target_fname = "%s_zero_test_mountain_tops.h5"%this_simname
    #outname = 'u302_zero_test_tracks_only.h5'
    target_fname = "u302_mountain_top_speed.h5"
    outname = "u302_tracks_speed.h5"


    new_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                    do_shift=False
                                  )
    new_looper.plot_directory = "./plots_to_sort"
    new_looper.read_targets(target_fname)
    import cProfile
    cProfile.run("new_looper.get_tracks()", "TrackProfile.txt") #tools/gprof2dot.py -f pstats TrackProfile.txt | dot -Tpng -o plots_to_sort/track_profile.png

    #loop_tools.re_shift_snaps(new_looper)
    #new_looper.save(outname)
    import tracks_read_write
    tracks_read_write.save_loop_only_trackage( new_looper, outname)
if 'new_looper' not in dir() and disk_or_new == 'disk':
    print("reload")
    file_list=['TEMP.h5']
    new_looper=looper.core_looper(directory=dl.sims[this_simname])
    new_looper.plot_directory = "/home/dccollins/PigPen"
    for nfile,fname in enumerate(file_list):
        new_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    new_looper.tr.sort_time()
    #import loop_tools
    #loop_tools.re_shift_snaps(new_looper)


import colors
reload(colors)
core_cmap = colors.make_core_cmap( new_looper.core_list)
reload(loop_apps)
loop_apps.core_proj_multiple(new_looper,#core_list=[258],
                             #frame_list=[0],
                             axis_list=[1], color_dict=core_cmap, particles=True,
                             only_sphere=False,zoom=False,
                             center_on_sphere=False, annotate=True,tracker_positions=True, shifted_tracker=False)
