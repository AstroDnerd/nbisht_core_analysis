#
# run get_mountain_tops first
#
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
reload(looper)
import looper2
reload(looper2)


LOOPER2 = True
looper_main = looper2.core_looper2
this_simname  = 'u502'
other_simname = 'u302'
mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%other_simname
outname = '%s_all_frame_all_prim.h5'%this_simname
bad_particle_fname_read='datasets_small/%s_bad_particles.h5'%this_simname
bad_particle_fname_save='datasets_small/%s_bad_particles.h5'%this_simname
#bad_particle_fname_read="datasets_small/u601_bad_particles_srsly.h5"

#bad_particle_fname_save='datasets_small/%s_bad_particles_take3.h5'%this_simname
#bad_particle_fname='datasets_small/%s_bad_particles_TEST.h5'%this_simname
import xtra_energy
target_frame = dl.target_frames[other_simname]
if 0:
    """Just One"""
    frame_list = list(range(0,target_frame,10)) + [target_frame]
    frame_list = [5]
if 1:
    """all frames"""
    target_frame = dl.target_frames[other_simname]
    frame_list = list(range(0,target_frame+1,1)) 

if 0:
    """first 4"""
    target_frame = dl.target_frames[other_simname]
    frame_list = list(range(0,3,1)) 

if 0:
    """Every 10"""
    target_frame = dl.target_frames[other_simname]
    frame_list = list(range(0,target_frame,10)) + [target_frame]


if 0:
    fields = ['x','y','z','density', 'cell_volume']
    derived=[]
if 1:
    fields = ['x','y','z','density', 'cell_volume']
    fields += ['velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    fields += ['PotentialField','grav_x','grav_y','grav_z' ]
    fields += ['particle_pos_x', 'particle_pos_y', 'particle_pos_z', 'particle_index']
    derived=[xtra_energy.add_force_terms]

if target_frame not in frame_list:
    print("YOU MUST HAVE THE LAST FRAME or the periodic unwrap fails")
    frame_list += [target_frame]

if 1:
    new_looper = looper_main(directory= dl.sims[other_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = target_frame,
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                    do_shift=False
                                  )
    new_looper.plot_directory = "./plots_to_sort"
    new_looper.read_targets(mountain_top_fname)

if 1:
    """get the right peaks"""
    """the right thing to do"""
    if LOOPER2:
        peaks = nar([new_looper.targets[core_id].peak_density for core_id in new_looper.core_list])
        remove = peaks < 1e4
        keep = peaks >= 1e4
        cores_to_cut = new_looper.core_list[ remove]
        mask = np.ones( new_looper.core_ids.size, dtype='bool')
        for core_id in cores_to_cut:
            mask = mask * (new_looper.core_ids != core_id)
        new_looper.core_ids = new_looper.core_ids[mask]
        new_looper.target_indices = new_looper.target_indices[mask]
        new_looper.core_list = new_looper.core_list[ keep ]
    else:
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

if 0:
    if os.path.exists(bad_particle_fname_save):
        print("File exists.  Overwrite? %s"%bad_particle_fname_save) 
        pdb.set_trace()
        print("File exists.  Overwrite? %s"%bad_particle_fname_save) 
    new_looper.save_bad_particles(bad_particle_fname_save)

#if 0:
#    new_looper.save_bad_particles('datasets_small/%s_bad_particles_full.h5'%this_simname)

if 1:
    new_looper.read_bad_particles(bad_particle_fname_read)
if 1:
    new_looper.remove_bad_particles()

if 1:
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


