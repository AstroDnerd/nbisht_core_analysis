from starter2 import *

import otherones
reload(otherones)

simlist = ['a001','a002','a003']
boxes={}
for sim in simlist:
    boxes[sim] = otherones.box_of_rain()
    fname = "box_of_masses_%s.h5"%sim
    print(fname)
    boxes[sim].read(fname)

if 1:
    #multi panel for exploring.
    fig,axes=plt.subplots(1,len(simlist), figsize=(12,8))
    axlist=axes
    for ns,sim in enumerate(simlist):
        box=boxes[sim]
        ax=axlist[ns]

        Mtotal, Mother =   box.mass_total_hull, box.mass_other_ones
        #ax.scatter(Mtotal,Mother/Mtotal)
        LLL = nar([Mtotal.min(), Mtotal.max()])
        ax.plot( LLL, LLL, c=[0.5]*3)
        ax.plot( LLL, 0.5*LLL, c=[0.5]*3)
        ax.plot( LLL, 0.25*LLL, c=[0.5]*3)

        ax.scatter(Mtotal,Mother, marker='+',c='k')
        axbonk(ax,xlabel=r'$M_{\rm{hull}}$', ylabel=r'$M_{\rm{otherones}}$', xscale='log',yscale='log')
        #axbonk(ax,xlabel=r'$M_{\rm{hull}}$', ylabel=r'$M_{\rm{otherones}}/M_{\rm{total}}$', xscale='log',yscale='linear')

    fig.savefig('plots_to_sort/otherones_mass.png')

if 0:
    #multi panel for exploring.
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
