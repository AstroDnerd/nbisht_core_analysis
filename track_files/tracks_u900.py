from starter2 import *
reload(track_info)


#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#


track_info.track('u902',
                 sim_directory=dl.sims['u202'],
                 target_frame = 118,
                 mountain_top = "datasets_small/u302_mountain_tops_take_9.h5",
                 peak_fname = 'datasets_small/u302_0118_peaklist.h5',
                 track_file = "%s/u900/u902_vel_grad.h5"%dl.coresets['ourset'],
                 mode_fname = "browser_data/core_formation_mode_new_u602.h5",
                 frame_list = "all_frames",
                 field_list = track_info.field_lists['everything'])


