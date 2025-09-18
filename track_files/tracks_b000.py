

from starter2 import *
reload(track_info)
#
# Define tracks, including track location and build information.
# The constructor also registers each one in track_info.tracks
#

NAME = 'b002'
track_info.track(NAME,
                 sim_directory=dl.sims['u202'],
                 target_frame = 60,
                 mountain_top = "%s/%s_mountain_tops.h5"%(dl.new_tracks,NAME),
                 peak_fname   = '%s/%s_0060_peaklist.h5'%(dl.new_tracks,NAME),
                 track_file   = "%s/%s_every_one_primitive.h5"%(dl.new_tracks,NAME),
                 mode_fname   = "%s/%s_modes.h5"%(dl.new_tracks,NAME),
                 frame_list   = "all_frames",
                 field_list   = track_info.field_lists['primitive'])
