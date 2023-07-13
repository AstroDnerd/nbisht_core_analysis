
from starter2 import *
reload(track_info)
#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#

track_info.track('t002',
                 sim_directory=dl.sims['u202'],
                 target_frame = 60,
                 mountain_top = "%s/t002_mountain_tops.h5"%dl.new_tracks,
                 peak_fname   = '%s/t002_0060_peaklist.h5'%dl.new_tracks,
                 track_file   = "%s/t002_every_ten_primitive.h5"%dl.new_tracks,
                 mode_fname   = "%s/t002_modes.h5"%dl.new_tracks,
                 clump_parameters={'small_test':True},
                 frame_list   = "every_ten",
                 field_list   = track_info.field_lists['primitive'])
