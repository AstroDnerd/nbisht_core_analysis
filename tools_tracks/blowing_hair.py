from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

sim_list=['u401']#,'u302','u303']
#sim_list=['u302','u303']

import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4
MOD = TL4
if 1:
    import convex_hull_tools as CHT
    reload(CHT)
    if 'ht_set' not in dir() :
        ht_set = {}
        for this_simname in sim_list:
            ht_set[this_simname] = CHT.hull_tool(MOD.loops[this_simname])
            ht_set[this_simname].make_hulls()
            ht_set[this_simname].make_overlaps()
    import supersets
    reload(supersets)
    if 'st_set' not in dir():
        st_set={}
        for this_simname in sim_list:
            st_set[this_simname] = supersets.superset(  MOD.loops[this_simname], ht_set[this_simname])
            st_set[this_simname].find(thresh=0.0)
reload(CHT)

if 1:
    for this_simname in sim_list:
        supersets = st_set[this_simname]
        for nset,this_superset in enumerate(supersets.supersets):
            if nset != 1:
                continue
            color_dict = colors.make_core_cmap(this_superset, cmap = 'tab20c', seed = -1)
            frame_list = MOD.loops[this_simname].frame_list
            CHT.plot_2d(ht_set[this_simname],core_list=this_superset,
                        frames=frame_list,accumulate=True,label_cores=[0], prefix = "S%02d"%nset,
                       color_dict=color_dict)

if 'lylist' not in dir() and False:
    lylist={}
    for this_simname in sim_list:
        lylist[this_simname]= hair_dryer.hair_tool( MOD.loops[this_simname])

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
    core_list = MOD.loops[this_simname].core_list
    for nset, superset in enumerate(st_set[this_simname].supersets):
        if nset != 1:
            continue
        core_list = list(superset)
        color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)

        lylist[this_simname]= lyapunov_tool( MOD.loops[this_simname])
        MOD.loops[this_simname].out_prefix = '%s_S%02d'%(this_simname, nset)
        lylist[this_simname].run( core_list = sorted(superset), newplots = False, colors = color_dict)


