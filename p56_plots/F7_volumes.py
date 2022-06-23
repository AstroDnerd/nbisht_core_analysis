from starter2 import *
import convex_hull_tools as CHT

import three_loopers_six as TL
sim_list=['u601','u602','u603']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls() 
        ht[this_simname].make_overlaps()

if 1:
    fig,ax=plt.subplots(1,3, figsize=(12,4))
    ax0=ax[0];ax1=atx[1];ax2=ax[2]
    for aaa in ax:
        aaa.set_aspect('equal')
    for ns,sim in enumerate(ht):
        the_y,the_x =  nar(ht[sim].cell_volumes), nar(ht[sim].hull_volumes)
        ax[ns].scatter( the_x,the_y)
        minmin = min([the_x.min(),the_y.min()])
        maxmax = max([the_x.max(),the_y.max()])
        ax[ns].plot([minmin,maxmax],[minmin,maxmax],c='k')
        axbonk(ax[ns],xscale='log',yscale='log',ylabel='Cell',xlabel='Hull',xlim=[minmin,maxmax],ylim=[minmin,maxmax])
    fig.savefig('plots_to_sort/volumes.png')

if 1:
    fig,ax=plt.subplots(1,3, figsize=(12,4))
    ax0=ax[0];ax1=ax[1];ax2=ax[2]
    for aaa in ax:
        aaa.set_aspect('equal')
    for ns,sim in enumerate(ht):
        the_y,the_x =  nar(ht[sim].cell_volumes), nar(ht[sim].hull_volumes)
        the_y = the_y/the_x
        ax[ns].scatter( the_x,the_y)
        minmin = min([the_x.min(),the_y.min()])
        maxmax = max([the_x.max(),the_y.max()])
        ax[ns].plot([minmin,maxmax],[minmin,maxmax],c='k')
        axbonk(ax[ns],xscale='log',yscale='log',ylabel='Cell/Hul',xlabel='Hull',xlim=[minmin,maxmax],ylim=[minmin,maxmax])
    fig.savefig('plots_to_sort/ratio.png')

