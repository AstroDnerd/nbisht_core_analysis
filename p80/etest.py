

from starter2 import *


from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')
import pcolormesh_helper as pch
reload(pch)
import progressbar
import time
import eigen
reload(eigen)

if 'ddd' not in dir() or clobber:
    ddd={}
if 1:
    for sim in sim_list:
        if sim not in ddd:
            ddd[sim]={}

        core_list=None
        #core_list=[7]
        core_list=[74]#, 112]
        #core_list=[112]
        #core_list=[74, 112]
        #core_list=[112]
        #core_list=[74]
        core_list=TL.loops[sim].core_by_mode['Alone']
        core_list=core_list[:15]

        #core_list=core_list[8:9]
        for core_id in core_list:
            if core_id in ddd[sim]:
                continue
            dddd = eigen.dq_dt2(TL.loops[sim])
            dddd.run(core_list=[core_id], OOM=False)
            #ddd[sim][core_id] = dddd
