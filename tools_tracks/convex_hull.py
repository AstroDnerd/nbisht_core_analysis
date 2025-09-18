
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
import colors
plt.close('all')

import convex_hull_tools as CHT
reload(CHT)

#import three_loopers_1tff as tl
import three_loopers_mountain_top as TLM

if 'clobber' not in dir():
    clobber=False

if 'ht1' not in dir() or clobber: 
    ht1 = CHT.hull_tool(TLM.loops['u301'])
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(TLM.loops['u302'])
    #ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = CHT.hull_tool(TLM.loops['u303'])

#
# BIG PLOT
#
if 0:
    CHT.plot_2d(ht1,frames=[0], accumulate=True, label_cores=[323])
if 0:
    CHT.plot_2d(ht2,frames=[0], accumulate=True, label_cores=[])
if 0:
    CHT.plot_2d(ht3,frames=[0], accumulate=True, label_cores=[])

if 0:
    fractions,cores=get_overlapping_cores(ht3,185)
    catman = np.concatenate
    cores = catman([cores,[185]])[::-1]
    ht3b = hull_tool(TLM.loopes['u303'])
    ht3b.plot_2d(core_list = cores, accumulate=True)

if 'did_hulls' not in dir():
    did_hulls = False
if not did_hulls:
    did_hulls = True
    #
    # Compute overlaps
    #
    frame=0
    for htool in [ht1, ht2, ht3]:
        htool.overlaps=defaultdict(list)
        htool.make_hulls(frames=[frame])
        for core_1 in htool.cores_used:
            print("overlap li,", core_1)
            for core_2 in htool.cores_used:
                result = htool.check_hull_overlap(core_1,core_2)
                if core_1 == core_2:
                    result = -result
                htool.overlaps[core_1].append(result)


if 1:
    #
    # Distribution of Next Overlap
    #
    #fig3,ax3=plt.subplots(2,1, sharex=True)
    fig3,ax3=plt.subplots(1,1, sharex=True)
    ax3a=ax3
    #fig3.subplots_adjust(wspace=0, hspace=0)
    #ax3a=ax3[0]
    #ax3b=ax3[1]
    bins = np.linspace(0,1,21) #every 0.05
    for htool in [ht1, ht2, ht3]:
        c=colors.color[ htool.this_looper.out_prefix]
        next_fraction=[]
        no_overlap=0
        all_overlap = 0
        for core_1 in htool.cores_used:
            this_over =  nar(htool.overlaps[core_1])
            next_fraction.append(this_over.max())
            if this_over[this_over>=0].sum() < 1e-32:
                no_overlap  += 1
            all_overlap += (this_over > 0.99).any()
        print( "all ", all_overlap)
        print("no overlap",no_overlap)

        #ax3a.plot( sorted(next_fraction), np.arange(len(next_fraction))+1)
        ax3a.hist( next_fraction, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=bins)
        #ax3b.hist( next_fraction, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=16, cumulative=True, density=True)
        ax3a.scatter([0],[no_overlap],c=c,marker="*")#,s=1.)
        ax3a.scatter([1],[all_overlap],c=c,marker="*")#,s=1.)
        ax3a.legend(loc=2)
        #axbonk(ax3b,xlabel=r'$\rm{Overlap\ fraction\ with\ nearest\ neighbor}$',ylabel=r'$f_{\rm{cores}}$')
        axbonk(ax3a,xlabel=r'$\rm{Overlap\ fraction\ with\ nearest\ neighbor}$',ylabel=r'$N_{\rm{cores}}$',ylim=[0,100])
        #fig.savefig('plots_to_sort/%s_overlaps.png'%htool.this_looper.out_prefix)
    fig3.savefig('plots_to_sort/next_overlap_dist_n%04d.png'%frame)
                
if 1:
    #
    # Number of neghbors with more than X
    #
    frame=0

    if 1:
        fig3,ax3=plt.subplots(1,3, figsize=(12,4))
        for htool in [ht1, ht2, ht3]:
            c=colors.color[ htool.this_looper.out_prefix]
            for nf, frac in enumerate([.2, .5, .9]):
                N_over_M = []
                for core_1 in htool.cores_used:
                    this_over =  nar(htool.overlaps[core_1])
                    #next_fraction.append(this_over.max())
                    N_over_M.append( (this_over > frac).sum())
                N_bins = max(N_over_M)+1
                ax3[nf].hist( N_over_M, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=N_bins)
                #ax3[nf].legend(loc=1)
                axbonk(ax3[nf],ylabel=r'$N_{cores}$',xlabel=r'$F_{> %d %s }$'%(int(100*frac),"\\%"),xlim=[0,22])#,ylim=[0,70])
                #fig.savefig('plots_to_sort/%s_overlaps.png'%htool.this_looper.out_prefix)
        fig3.savefig('plots_to_sort/next_neighbor_over_M_%04d.pdf'%frame)

