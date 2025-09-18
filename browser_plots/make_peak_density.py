from starter2 import *


import three_loopers_six as TL
import mountain_top

if 1:
    for nsim,sim in enumerate(TL.loops):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)

        core_id_list=[]
        value_list=[]
        for core_id in sorted(core_list):
            core_id_list.append(core_id)
            value_list.append( this_looper.targets[core_id].peak_density)
        fptr = h5py.File('browser_data/peak_density_%s.h5'%sim,'w')
        fptr['core_ids']=nar(core_id_list)
        fptr['peak_density']=nar(value_list)
        fptr.close()




