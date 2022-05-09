from starter2 import *

import three_loopers_u500 as TL

for sim in TL.loops:
    this_looper=TL.loops[sim]
    core_list=sorted(np.unique(this_looper.tr.core_ids))


    values = []
    thtr=this_looper.tr
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        values.append( ms.nparticles)
    fptr=h5py.File('browser_data/%s_nparticles.h5'%sim,'w')
    fptr['core_ids']=core_list
    fptr['n_particles']=values
    fptr.close()
    print('wrote %s'%sim)

