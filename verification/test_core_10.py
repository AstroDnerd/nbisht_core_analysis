from starter2 import *

reload(looper)
import tracks_read_write
reload(tracks_read_write)
#
# get mountain tops
#

this_simname = 'u05'
core_id = 10
mountain_top_name = "%s_mountain_top_core_c%04d.h5"%(this_simname,core_id)
do_mountain_projections=True
if 0:
    import mountain_top
    reload(mountain_top)
    def verify_cores_very_dumb(top1):
        if top1.rhomax < 1100:
            return False
        else:
            return True
    outname = "u301_tracks_core_10_temp.h5"
    #kludge={'peak_id':[0,1,364, 113,u302_was258]})
    #kludge={'peak_id':[258]}
    kludge={'peak_id':[core_id]}
    MT=mountain_top.cut_mountain_top( this_simname, target_fname = mountain_top_name, do_projections=do_mountain_projections, verify=verify_cores_very_dumb,kludge=kludge)
#
# Get tracks.  
#

if 0:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list = list(range(0,target_frame,10))+[target_frame]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField']
    derived=[]


if 0:
    L258 = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname ,
                                     out_prefix = this_simname + "_c%04d"%core_id,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  [core_id],
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False,
                                     #plot_directory = "./plots_to_sort"
                                  )
    if 1:
        print("READ TARGETS", mountain_top_name)
        L258.read_targets(mountain_top_name)
    if 0:
        L258.get_target_indices(h5_name=dl.peak_list[this_simname],
                                         bad_particle_list=dl.bad_particles.get(this_simname,None))
    L258.plot_directory = "/home/dccollins/PigPen"
    L258.get_tracks()

#
# Read tracks
#

if 1:
    import three_loopers_mountain_top as TLM
    L258 = TLM.loops['u301']

if 0:
    #read_write tests
    read_write_fname = "L258_trackage_only.h5"
    L258.save_trackage_only(read_write_fname)

if 0:
    L258_read = looper.core_looper(directory=dl.sims[this_simname],savefile=  read_write_fname)
    L258_read.out_prefix+='_read'

if 1:

    import colors
    reload(colors)
    reload(loop_apps)
    core_cmap = colors.make_core_cmap( L258.core_list)
    loop_apps.core_proj_multiple(L258,core_list=[10],
                                 #frame_list=[0],
                                 axis_list=[0], color_dict=core_cmap, particles=True,
                                 only_sphere=True,zoom=True,
                                 center_on_sphere=True, annotate=False,tracker_positions=True, shifted_tracker=True)
