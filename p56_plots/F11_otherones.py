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

    this_simname = 'x003'
    core_loop =  TL6.loops['u603']
    new_looper = TLO.loops['x003']
    #core_list=[323]
    core_list=[107]
    #core_list =np.unique(ht[this_simname].this_looper.tr.core_ids)[:1]
    fig,ax = plt.subplots(2,2, figsize=(6,6))
    fig.subplots_adjust(wspace=0, hspace=0)
    IM = otherones_hair.images(new_looper, core_loop)
    IM.run(frames=[0, core_loop.target_frame], core_list=core_list, output_prefix=this_simname, external_ax=ax)#,core_list=[3])
    ax[0][0].xaxis.tick_top()
    ax[0][1].xaxis.tick_top()
    ax[0][1].xaxis.set_label_position('top')
    ax[0][0].xaxis.set_label_position('top')
    ax[1][1].yaxis.tick_right()
    ax[0][1].yaxis.tick_right()
    ax[1][1].yaxis.set_label_position('right')
    ax[0][1].yaxis.set_label_position('right')
    outname = 'plots_to_sort/F11_otherones_u603_c0107.png'
    fig.savefig(outname)
    #IM.run(frames=[0, core_loop.target_frame])#,core_list=[3])
    #IM.run(frames=[0,118])
#hd = hair_dryer.hair_tool( this_looper )
#hd.run()
