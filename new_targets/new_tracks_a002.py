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
this_simname  = 'b002'
other_simname = 'u302'
mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%other_simname
outname = '%s_all_particles.h5'%this_simname
bad_particle_fname_read='datasets_small/%s_bad_particles.h5'%'u501'
bad_particle_fname_save='datasets_small/%s_bad_particles_save.h5'%'u501'
#bad_particle_fname_read="datasets_small/u601_bad_particles_srsly.h5"

#bad_particle_fname_save='datasets_small/%s_bad_particles_take3.h5'%this_simname
#bad_particle_fname='datasets_small/%s_bad_particles_TEST.h5'%this_simname
import xtra_energy
target_frame = dl.target_frames[other_simname]
if 0:
    """Just One"""
    frame_list = list(range(0,target_frame,10)) + [target_frame]
    frame_list = [5]
if 0:
    """all frames"""
    target_frame = dl.target_frames[other_simname]
    frame_list = list(range(0,target_frame+1,1)) 

if 0:
    """first 4"""
    target_frame = dl.target_frames[other_simname]
    frame_list = list(range(0,3,1)) 

if 1:
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
    fields += [('gas','magnetic_field_%s'%s) for s in 'xyz']
    fields += [YT_potential_field,YT_acceleration_x, YT_acceleration_y, YT_acceleration_z]
    fields += [YT_particle_pos_x, YT_particle_pos_y, YT_particle_pos_z]
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

if 1:
    core_id=0
    new_looper.core_list=[core_id]
    i,j,k=np.mgrid[0:128:1,0:128:1,0:128:1]
    SL = tuple([slice(32,40)]*3)
    SL = tuple([slice(None)]*3)
    i_keep=i[SL]
    j_keep=j[SL]
    k_keep=k[SL]
    index = i_keep+128*(j_keep+128*k_keep)
    new_looper.target_indices=np.sort(index.flatten())
    new_looper.core_ids = np.zeros_like(new_looper.target_indices)



if 1:
    new_looper.read_bad_particles(bad_particle_fname_read)

if 1:
    new_looper.remove_bad_particles()

if 1:
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


