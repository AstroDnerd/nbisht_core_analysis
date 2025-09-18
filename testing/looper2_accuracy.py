from starter2 import *

import convex_hull_tools as CHT 

if 1:
    fname1 = 'u701_long_prim.h5'
    directory=dl.sims['u301']
    L1 = looper2.core_looper2(directory=directory,savefile_only_trackage=fname1)

    HT1 = CHT.hull_tool(loop1)
    CHT.plot_2d(HT1, accumulate=True,core_list=[31,32])


if 0:
    porder1 = np.argsort(t1.particle_ids)
    porder2 = np.argsort(t2.particle_ids)

    d1 = t1.track_dict['density']
    d2 = t2.track_dict['density']
    print('density sorted', np.abs(d1[porder,:]-d2).sum())
    cores1_sorted=t1.core_ids[porder1]
    print( 'cores_sorted',cores1_sorted)
    cores2_sorted=t2.core_ids[porder2]
    print( 'cores_sorted',cores2_sorted)
    print( np.abs(cores1_sorted-cores2_sorted).sum())
