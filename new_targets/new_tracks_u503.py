#
# run get_mountain_tops first
#
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
reload(looper)

this_simname  = 'u503'
other_simname = 'u303'

mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%other_simname
outname = '%s_long_prim.h5'%this_simname
import xtra_energy
if 1:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list = list(range(0,target_frame+1,1))
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField','grav_x','grav_y','grav_z' ]
    fields += ['particle_pos_x', 'particle_pos_y', 'particle_pos_z', 'particle_index']
    derived=[xtra_energy.add_force_terms]

if 1:
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
    new_looper.read_targets(mountain_top_fname)

if 1:
    peaks = nar([new_looper.targets[core_id].peak_density for core_id in new_looper.core_list])
    remove = peaks < 1e4
    keep = peaks >= 1e4
    cores_to_cut = new_looper.core_list[ remove]
    for core_id in cores_to_cut:
        del new_looper.target_indices[core_id]
    new_looper.core_list = new_looper.core_list[ keep ]

if 0:
    #verify and remove bad particles.
    #Don't need to do this every time, it takes a minute.
    #need to have the bad_particles file 
    for frame in new_looper.frame_list:
        print('check particles, frame',frame)
        new_looper.verify_all_particles(frame)
    new_looper.save_bad_particles('datasets_small/%s_bad_particles_full.h5'%this_simname)

if 1:
    new_looper.read_bad_particles('datasets_small/%s_bad_particles_full.h5'%this_simname)
    new_looper.remove_bad_particles()
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


