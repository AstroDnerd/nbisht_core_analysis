from starter2 import *

reload(looper)
import three_loopers as TL
reload(loop_apps)
import loop_tools
reload(loop_tools)

# LOAD LOOPER.PY????
kill = []


if TL.looper1.core_list not in dir():
    #TL.looper1.core_list= looper.get_all_nonzero(dl.n_particles['u201'])
    #TL.looper2.core_list= looper.get_all_nonzero(dl.n_particles['u202']) 
    #TL.looper3.core_list= looper.get_all_nonzero(dl.n_particles['u203'])
    
    TL.looper1.frame_list =dl.frame_list['u05']
    #TL.looper2.frame_list =dl.frame_list['u202']
    #TL.looper3.frame_list =dl.frame_list['u203']


print("BEFORE ALL LOOPERS loop")
#all_loopers = [TL.looper1, TL.looper2, TL.looper3]
#for i in range(len(all_loopers)):
#    loop_apps.core_proj_follow(all_loopers[i])

loop_apps.core_proj_follow(TL.looper1, core_list =[70])





