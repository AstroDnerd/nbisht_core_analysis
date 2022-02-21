
from starter2 import *
import colors

import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']

if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL4.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = stay_close.close_tool( TL4.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL4.loops[this_simname], ht[this_simname])
        st[this_simname].find()


fig,ax=plt.subplots(1,3)
reload(CHT)
if 1:
    for ns,this_simname in enumerate(sim_list):
        tool=ht[this_simname]
        color_dict = colors.make_core_cmap(tool.cores_used)
        for nc,core_id in enumerate(ht[this_simname].cores_used):
            o = nar(tool.overlaps[core_id])
            od = overlap_numb[this_simname]
            in_me = od[nc,:]
            me_in = od[:,nc]
            yes_have_some = (in_me == 0 )*(me_in > 0 )
            butts = np.where(yes_have_some)[0]
            if yes_have_some.sum()==0:
                continue

            parts=nar(particles[this_simname])[butts]
            print(core_id, parts)
            core_ids = nar(tool.cores_used)[butts]
            other_core = core_ids[ np.argmax(parts)]

            core_list=[core_id,other_core]
            CHT.plot_2d(tool,core_list=core_list, color_dict=color_dict, prefix='BigButt_c%04d'%core_id,accumulate=True,axis_to_plot=[-1])



if 0:
    #kind of sucked
    for ns,this_simname in enumerate(sim_list):
        tool=ht[this_simname]
        color_dict = colors.make_core_cmap(tool.cores_used, cmap = 'tab20c', seed = -1)
        for nc,core_id in enumerate(ht[this_simname].cores_used):
            o = nar(tool.overlaps[core_id])

            if o[o>0].sum() == 0:
                core_list=[core_id]
            else:
                other_core = tool.cores_used[ np.argmax(o)]
                core_list=[core_id,other_core]
            CHT.plot_2d(tool,core_list=core_list, color_dict=color_dict, prefix='BFF_c%04d'%core_id,accumulate=True)


