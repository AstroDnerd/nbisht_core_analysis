from starter2 import *
reload(looper)
reload(dl)

loops={}
sim_u201 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128-Beta0.2'
sim_u202 = '/scratch1/dcollins/Paper19/u202-Beta2'
sim_u203 = '/scratch1/dcollins/Paper19/u203-Beta20'

mountain_top = {'u301':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u301_new_tracks_take_9.h5",
                "u302":"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u302_new_tracks_take_9.h5",
                "u303":"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u303_new_tracks_take_9.h5"}
import tracks_read_write
loop_301 = looper.core_looper(directory= sim_u201)
tracks_read_write.load_trackage_only(loop_301,mountain_top['u301'])
loop_302 = looper.core_looper(directory= sim_u202)
tracks_read_write.load_trackage_only(loop_302,mountain_top['u302'])
#loop_302 = looper.core_looper(directory= sim_u202,
#                          savefile_only_trackage=mountain_top['u302'])
print(loop_301.frame_list)
#    loops[simname] = loop

#if 'u301' not in loops:
#    print('Load u301')
#    load_looper('u301')
#
#if 'u302' not in loops:
#    print('Load u302')
#    load_looper('u302')
#
#if 'u303' not in loops:
#    print('Load u303')
#    load_looper('u303')


