from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_3d_tool
reload(hair_3d_tool)


sim_list=['u301','u302','u303']
#sim_list=['u302','u303'] if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()
if 'lylist' not in dir() or True:
    import three_loopers_tenfour as TL4 
    tool = lyapunov_tool( TL4.['u401'])


    frame_list=tool.this_looper.frame_list
    frame_list=[0]
    rotato = True

    for frame in frame_list:
        tool.run( core_list = core_list, newplots = False, colors = color_dict,rotato=rotato, frame=frame)





if 0:
    import three_loopers_mountain_top as TLM
    lylist={}
    #for this_simname in sim_list:
    #    lylist[this_simname]= lyapunov_tool( TLM.loops[this_simname])
    lylist['u301'] = lyapunov_tool( set_looper)

#    for this_simname in  sim_list:
#        lylist[this_simname].run( )#core_list=[10,11])#core_list=[10])

    #lylist['u301'].run( core_list = [0,8,27,32,37,44,84,275])
    #lylist['u301'].run( core_list = [323])
    #lylist['u301'].run( core_list = [24])
    #lylist['u302'].run( core_list = [30,32])
    #lylist['u303'].run( core_list = [233,235])
    #lylist['u303'].run( core_list = [184])
    #lylist['u303'].run( core_list = [186])

    #core_list=[12,13,14,15,16,17,18,19,21] 
    target_core=31
    core_list = nar(ht1.cores_used)[ nar(ht1.overlaps[target_core   ])>0] 
    core_list = np.append(core_list,[target_core])
    popper  = set([124,108, 66,65,62,61])
    core_set = set(core_list)
    core_list = list(core_set.difference(popper))
    reload(colors)
    color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
    #lylist['u301'].run( core_list = core_list, newplots = False, colors = color_dict,rotato=False, frame=30)
    frame_list=lylist['u301'].this_looper.frame_list
    frame_list=[0]
    rotato = True
    for frame in frame_list:
        lylist['u301'].run( core_list = core_list, newplots = False, colors = color_dict,rotato=rotato, frame=frame)


