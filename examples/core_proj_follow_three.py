from starter2 import *

reload(looper)
import three_loopers as TL
reload(loop_apps)
import loop_tools
reload(loop_tools)
kill = []

TL.looper1.core_list= looper.get_all_nonzero(dl.n_particles['u05'])
TL.looper2.core_list= looper.get_all_nonzero(dl.n_particles['u10'])
TL.looper3.core_list= looper.get_all_nonzero(dl.n_particles['u11'])

TL.looper1.frame_list =list(range(0,130,10))+[125]
TL.looper2.frame_list =list(range(0,82,10))+[82]
TL.looper3.frame_list =list(range(0,88,10))+[88]

#for looper666 in [TL.looper1, TL.looper2, TL.looper3]:
loop_apps.core_proj_follow(TL.looper3,axis_list=[0])

