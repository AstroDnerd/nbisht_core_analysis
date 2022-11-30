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
this_simname = 'u902'
other_simname = 'u302'
mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%other_simname
outname = '%s_vel_grad.h5'%this_simname
bad_particle_fname_read='datasets_small/%s_bad_particles.h5'%'u502'
bad_particle_fname_save='datasets_small/%s_bad_particles_NO_DONT_OVERWRITE.h5'%'u501'
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
    fields = [YT_x,YT_y,YT_z,YT_density, YT_cell_volume]
    fields += [YT_velocity_magnitude,YT_magnetic_field_strength, YT_velocity_divergence]
    fields += [YT_velocity_x,YT_velocity_y,YT_velocity_z]
    fields += [YT_dxvx,YT_dxvy,YT_dxvz]
    fields += [YT_dyvx,YT_dyvy,YT_dyvz]
    fields += [YT_dzvx,YT_dzvy,YT_dzvz]
    fields += [('gas','magnetic_field_%s'%s) for s in 'xyz']
    fields += [YT_potential_field]
    #fields += [YT_potential_field,YT_grav_x,YT_grav_y,YT_grav_z ]
    fields += [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z]
    fields += [YT_particle_pos_x, YT_particle_pos_y, YT_particle_pos_z, YT_particle_index]
    derived=[xtra_energy.add_force_terms, xtra_energy.add_v_grad]

if target_frame not in frame_list:
    print("YOU MUST HAVE THE LAST FRAME or the periodic unwrap fails")
    frame_list += [target_frame]

if 1:
    new_looper = looper_main(directory= dl.sims[other_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[other_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                    do_shift=False
                                  )
    new_looper.plot_directory = "./plots_to_sort"
    new_looper.read_targets(mountain_top_fname)

if 1:
    #Pick only these peaks above 1e4
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


if 1:
    new_looper.read_bad_particles(bad_particle_fname_read, core_hijack=0)
    new_looper.remove_bad_particles()


if 1:
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


