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
this_simname = 'u601'
other_simname = 'u301'
mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%other_simname
outname = '%s_every_10_all_prim.h5'%this_simname
bad_particle_fname_read='datasets_small/%s_bad_particles.h5'%'u501'
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

if 1:
    for missing_frame in [122, 123, 124]:
        if missing_frame in frame_list:
            frame_list.remove(missing_frame)

if 0:
    fields = ['x','y','z','density', 'cell_volume']
    derived=[]
if 1:
    fields = [YT_x,YT_y,YT_z,YT_density, YT_cell_volume]
    fields += [YT_velocity_magnitude,YT_magnetic_field_strength, YT_velocity_divergence]
    fields += [YT_velocity_x,YT_velocity_y,YT_velocity_z]
    fields += [('gas','magnetic_field_%s'%s) for s in 'xyz']
    fields += [YT_potential_field,YT_grav_x,YT_grav_y,YT_grav_z ]
    fields += [YT_particle_pos_x, YT_particle_pos_y, YT_particle_pos_z, YT_particle_index]
    derived=[xtra_energy.add_force_terms]

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

import bad_particle_hunt
if 1:

    print("Look for bad particles again.  Somehow we can't get ahead of this.")
    aaa = set(np.arange(128**3))
    badones=set()
    for frame in new_looper.frame_list:
        ds=new_looper.load(frame)
        ad=ds.all_data()
        pi = set(ad['all','particle_index'].v)
        badones.update(aaa-pi)
        for grid in ds.index.grids:
            these = set(grid['all','particle_index'].v)
            pi.difference_update( these)
        badones.update(pi)

        also_bad = bad_particle_hunt.check_particles(ds)
        badones.update(set(also_bad))
        print(frame, len(badones))

    new_looper.read_bad_particles(bad_particle_fname_read, core_hijack=0)

    bad_particle_id = [ 724134,  635702,  661226,  743270,  751995,  718196, 1354060,
                                1362500,  610123,  610189, 1930558, 1046537, 1841352, 1844125,
                                1845574, 1849410, 1853445, 1300291]
    badones.update(set(bad_particle_id))
    bad_particle_id = list(badones) #this is some confusing variable naming.
    bad_core_id = [0]*len(bad_particle_id)
    for bad_core, bad_part in zip(bad_core_id, bad_particle_id):
        new_looper.bad_particles[bad_core]=np.append(
                new_looper.bad_particles[bad_core], bad_part)
if 1:
    new_looper.remove_bad_particles()


if 1:
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


