#
# run get_mountain_tops first
#
from starter2 import *
import tracks_read_write
reload(dl)
looper_main = looper2.core_looper2

def get_tracks(trackname):
    track=track_info.tracks[trackname]
    target_frame=track.target_frame
    mountain_top_fname = track.mountain_top
    outname  = track.track_file

    if os.path.exists(outname):
        print("File exists, will not make tracks", outname)
        return 0


    if target_frame not in track.frame_list:
        print("YOU MUST HAVE THE LAST FRAME or the periodic unwrap fails")
        frame_list += [target_frame]

    new_looper = looper_main(directory= track.sim_directory,
                                     sim_name = trackname,
                                     out_prefix = trackname,
                                     target_frame = track.target_frame,
                                     frame_list = track.frame_list,
                                     core_list =  None,
                                     fields_from_grid=track.field_list,
                                     derived = track.derived_fields,
                                    do_shift=False
                                  )
    new_looper.plot_directory = track.plot_directory
    new_looper.read_targets(mountain_top_fname)

    if 0:
        #
        #  I think that this machinery is now in the mountain topper.
        #
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


    if track.bad_particles is not None:
        new_looper.read_bad_particles(bad_particle_fname_read, core_hijack=0)
        new_looper.remove_bad_particles()


    #
    # The actual work
    #
    new_looper.get_tracks()

    tracks_read_write.save_loop_trackage_only( new_looper, outname)


