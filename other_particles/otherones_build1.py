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
if 'this_simname' not in dir():
    new_simname = 'a001'
    other_simname = 'u301'
    this_simname = 'u401'
    save='a001_all_particles.h5'
    this_looper = looper2.core_looper2( directory = dl.sims[other_simname], savefile_only_trackage=save)
    print('make ms, takes about 90 sec')
    ms = trackage.mini_scrubber(this_looper.tr, 0, do_velocity=False)
    this_looper.ms = ms

    TL4.loops[this_simname].big_loop=this_looper
sim_list=[this_simname]

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

if 'new_looper' not in dir() :
    import otherones
    reload(otherones)
    print('make new one')
    core_list=None
    #core_list=np.unique(ht[this_simname].this_looper.tr.core_ids)

    for simname in sim_list:
        superset = None; name='has_neighbor'
        superset = st[simname]; name='no_neighbor'
        new_looper=otherones.find_other_ones(new_simname,ht[this_simname],core_list=core_list,superset=superset)
        outname = "otherones_%s_%s.h5"%(this_simname,name)
    if not os.path.exists(outname):
        print("SAVING FILE", outname)
        tracks_read_write.save_loop_trackage_only( new_looper, outname)

if 1:
    fig,axes=plt.subplots(2,2)
    ax=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0]

    n_particles = [(ht[this_simname].this_looper.tr.core_ids == core_id).sum() for core_id in new_looper.box.cores_used]

    Mtotal, Mother =   new_looper.box.mass_total_hull, new_looper.box.mass_other_ones
    ax2.scatter(n_particles, Mother/Mtotal)
    axbonk(ax2,xscale='log')
    ax.scatter(Mtotal, Mother)
    ax.plot([Mtotal.min(),Mtotal.max()],[Mtotal.min(),Mtotal.max()],c=[0.5]*3)
    axbonk(ax,xlabel=r'$M_{\rm{hull}}$', ylabel=r'$M_{\rm{otherones}}$', xscale='log',yscale='log')
    ax1.scatter(Mtotal,Mother/Mtotal)
    axbonk(ax1,xlabel=r'$M_{\rm{hull}}$', ylabel=r'$M_{\rm{otherones}}/M_{\rm{total}}$', xscale='log',yscale='linear')

    fig.savefig('plots_to_sort/otherones_mass_%s.png'%new_looper.sim_name)

         
if 1:
    import otherones_hair
    reload(otherones_hair)
    core_loop = TL4.loops[this_simname]
    core_list=None
    #core_list =np.unique(ht[this_simname].this_looper.tr.core_ids)[:1]
    IM = otherones_hair.images(new_looper, core_loop)
    IM.run(frames=[0, core_loop.target_frame], core_list=core_list, output_prefix=this_simname)#,core_list=[3])
    #IM.run(frames=[0, core_loop.target_frame])#,core_list=[3])
    #IM.run(frames=[0,118])
#hd = hair_dryer.hair_tool( this_looper )
#hd.run()
