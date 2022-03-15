from starter2 import *
from collections import defaultdict
import scipy
import colors

import convex_hull_tools as CHT
import hair_dryer
import otherones

import three_loopers_otherones as TLO
import three_loopers_six as TL6

import otherones_hair
reload(otherones_hair)
         
if 1:

    this_simname = 'x001'
    core_loop =  TL6.loops['u601']
    new_looper = TLO.loops['x001']
    #core_list=[323]
    core_list=[37]
    #core_list =np.unique(ht[this_simname].this_looper.tr.core_ids)[:1]
    IM = otherones_hair.images(new_looper, core_loop)
    IM.run(frames=[0, core_loop.target_frame], core_list=core_list, output_prefix=this_simname)#,core_list=[3])
    #IM.run(frames=[0, core_loop.target_frame])#,core_list=[3])
    #IM.run(frames=[0,118])
#hd = hair_dryer.hair_tool( this_looper )
#hd.run()
