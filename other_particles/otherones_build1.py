from starter2 import *
from collections import defaultdict
import scipy
import colors

import convex_hull_tools as CHT
import hair_dryer
import otherones
reload(hair_dryer)
import looper2
import three_loopers_tenfour as TL4
sim_list=['u402']
if 'this_simname' not in dir():
    this_simname = 'a002'
    other_simname = 'u302'
    save='a002_all_particles.h5'
    this_looper = looper2.core_looper2( directory = dl.sims[other_simname], savefile_only_trackage=save)
    print('make ms, takes about 90 sec')
    ms = trackage.mini_scrubber(this_looper.tr, 0, do_velocity=False)
    this_looper.ms = ms

    TL4.loops['u402'].big_loop=this_looper

if 'ht' not in dir():
    ht={}
for this_simname in sim_list:
    if this_simname not in ht:
        ht[this_simname] = CHT.hull_tool(TL4.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

if 'st' not in dir():
    import supersets
    reload(supersets)
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL4.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 'new_looper' not in dir():
    import otherones
    reload(otherones)
    core_list=[3,5]
    print('make new one')
    new_looper=otherones.find_other_ones('a002',ht['u402'],core_list=core_list)#, superset=st['u402'])
    outname = "otherones_b002_temp.h5"
if 0:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)
         
if 1:
    import otherones_hair
    reload(otherones_hair)
    IM = otherones_hair.images(new_looper, TL4.loops['u402'])
    IM.run(frames=[0],core_list=[3])
    #IM.run(frames=[0,118])
#hd = hair_dryer.hair_tool( this_looper )
#hd.run()
