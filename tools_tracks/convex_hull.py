
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})

import convex_hull_tools as CHT
reload(CHT)

import three_loopers_1tff as tl

if 'clobber' not in dir():
    clobber=False

if 'ht1' not in dir() or clobber: 
    ht1 = CHT.hull_tool(tl.looper1)
if 'ht2' not in dir() or clobber:
    ht2 = CHT.hull_tool(tl.looper2)
    #ht2.plot_2d(frames=[0])
if 'ht3' not in dir() or clobber:
    ht3 = CHT.hull_tool(tl.looper3)

if 0:
    CHT.plot_2d(ht1,frames=[0],core_list=[85,86, 306, 307, 308], accumulate=True)
if 0:
    CHT.plot_2d(ht1,frames=[0], accumulate=True, label_cores=[323])
if 1:
    CHT.plot_2d(ht2,frames=[0], accumulate=True, label_cores=[-1])
if 0:
    CHT.plot_2d(ht3,frames=[0], accumulate=True, label_cores=[])

if 0:
    fractions,cores=get_overlapping_cores(ht3,185)
    catman = np.concatenate
    cores = catman([cores,[185]])[::-1]
    ht3b = hull_tool(tl.looper3)
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


if 0:
    #
    # Distribution of Next Overlap
    #
    #fig3,ax3=plt.subplots(2,1, sharex=True)
    fig3,ax3=plt.subplots(1,1, sharex=True)
    ax3a=ax3
    #fig3.subplots_adjust(wspace=0, hspace=0)
    #ax3a=ax3[0]
    #ax3b=ax3[1]
    for htool in [ht1, ht2, ht3]:
        c=color[ htool.this_looper.out_prefix]
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
        ax3a.hist( next_fraction, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=16)
        #ax3b.hist( next_fraction, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=16, cumulative=True, density=True)
        ax3a.scatter([0],[no_overlap],c=c,marker="*")#,s=1.)
        #ax3a.scatter([1],[all_overlap],c=c,marker="*")#,s=1.)
        ax3a.legend(loc=2)
        #axbonk(ax3b,xlabel=r'$\rm{Overlap\ fraction\ with\ nearest\ neighbor}$',ylabel=r'$f_{\rm{cores}}$')
        axbonk(ax3a,xlabel=r'$\rm{Overlap\ fraction\ with\ nearest\ neighbor}$',ylabel=r'$N_{\rm{cores}}$',ylim=[0,70])
        #fig.savefig('plots_to_sort/%s_overlaps.png'%htool.this_looper.out_prefix)
    fig3.savefig('plots_to_sort/next_overlap_dist_n%04d.png'%frame)
                
if 0:
    #
    # Number of neghbors with more than X
    #
    frame=0
    for htool in [ht1, ht2, ht3]:
        htool.make_hulls(frames=[frame])

    if 1:
        fig3,ax3=plt.subplots(3,1, figsize=(4,12))
        for htool in [ht1, ht2, ht3]:
            c=color[ htool.this_looper.out_prefix]
            for nf, frac in enumerate([.2, .5, .9]):
                N_over_M = []
                for core_1 in htool.cores_used:
                    this_over =  nar(htool.overlaps[core_1])
                    #next_fraction.append(this_over.max())
                    N_over_M.append( (this_over > frac).sum())
                N_bins = max(N_over_M)+1
                ax3[nf].hist( N_over_M, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix,bins=N_bins)
                #ax3[nf].legend(loc=1)
                axbonk(ax3[nf],ylabel=r'$N_{cores}$',xlabel=r'$N_{> %d %s }$'%(int(100*frac),"\\%"),xlim=[0,22])#,ylim=[0,70])
                #fig.savefig('plots_to_sort/%s_overlaps.png'%htool.this_looper.out_prefix)
        fig3.savefig('plots_to_sort/next_neighbor_over_M_%04d.pdf'%frame)


if 0:
    import p56_plots.density_AC as AC
    reload(AC)
    if 'a1' not in dir():
        a1 = AC.ac_thing('u05'); a1.plot()
        a2 = AC.ac_thing('u10'); a2.plot()
        a3 = AC.ac_thing('u11'); a3.plot()
        acs={'u05':a1,'u10':a2,'u11':a3}

        
if 0:
    fig,ax=plt.subplots(1,1)
    #ax.hist(hull_lengths,histtype='step',color='k',normed=True)
    vals, bins = np.histogram(hull_lengths)
    bc = 0.5*(bins[1:]+bins[:-1])
    db = (bins[1:]-bins[:-1])
    ax.plot(bc,vals/vals.sum())
    ax.plot(AC.binned[1],AC.binned[2]/AC.binned[2][0],c='r')
    rect=patches.Rectangle((0,0),AC.L,AC.ACb[0],facecolor=[0.8]*3)
    ax.add_patch(rect)
    fig.savefig('plots_to_sort/%s_sizes.pdf'%this_simname)


if 0:
    #
    # hull lengths with density ac
    #
    fig,ax=plt.subplots(1,1)
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        c=color[ ht.this_looper.out_prefix]
        hull_lengths = nar(ht.hull_volumes)**(1./3)
        vals, bins = np.histogram(hull_lengths)
        bc = 0.5*(bins[1:]+bins[:-1])
        ax.plot(bc,vals/vals.sum(),color=c,label=ht.this_looper.out_prefix,linestyle="--")
        axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$',ylim=[0,0.6])

        ac = acs[ht.this_looper.out_prefix]
        ax.plot(ac.binned[1],ac.binned[2], c=c)
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/hull_lengths.pdf')
    plt.close('fig')

if 0:
    fig,ax = plt.subplots(1,1)
    def n_neighbors_above_f(self,fraction):
        N_neighbors=[]
        for core_id in self.cores_used:
            N_neighbors.append( (self.overlaps[core_id]>=fraction).sum())
        return N_neighbors

    Nn = 100
    Nf = 10
    N_nei_f = np.zeros([Nf, Nn])
    for htool in [ht1, ht2, ht3]:
        for iFr,fr in enumerate(np.arange(.1,1,.1)):
            print("do ", htool.this_looper.out_prefix, " fr",fr)
            N_nei = nar(n_neighbors_above_f(htool,fr))
            print( "Max neighbors", max(N_nei))
            for iN in range(max(N_nei)+1):
                N_nei_f[iFr,iN] = (N_nei >= iN).sum()

if 0:
    #
    # Hull volume vs total cell volume
    #
    frame = 10
    fig,axess=plt.subplots(1,3)
    if 'ext_hull' not in dir():
        ext_hull=extents()
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        if nrun != 2:
            continue
        ax=axess[nrun]
        ax.clear()
        name = ht.this_looper.out_prefix
        ht.make_hulls(frames=[frame])
        odd = nar(ht.hull_volumes) < nar(ht.cell_volumes)
        not_odd = nar(ht.hull_volumes) >= nar(ht.cell_volumes)
        ax.scatter(ht.hull_volumes,ht.cell_volumes,c='k')
        ax.scatter(nar(ht.hull_volumes)[odd],nar(ht.cell_volumes)[odd],c='r')
        ax.set_aspect('equal')
        #axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=[0,0.07],ylim=[0,0.07])
        ext_hull(nar(ht.hull_volumes))
        ext_hull(nar(ht.cell_volumes))
        ax.plot(ext_hull.minmax, ext_hull.minmax,c='g')
        axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',
               xlim=ext_hull.minmax,ylim=ext_hull.minmax, xscale='log',yscale='log')
    fig.savefig("plots_to_sort/hull_volume_n%04d.png"%(frame))

if 0:
    #
    # Hull volume vs total cell volume
    #
    hull_by_frame = {}
    hullvol = defaultdict(list)
    cellvol = defaultdict(list)
    looper_list=[tl.looper1,tl.looper2,tl.looper3]
    loopers = dict(zip([ looper.out_prefix for looper in looper_list], looper_list))

    for loop in looper_list:
        name = loop.out_prefix
        hull_by_frame[name]={}
        if name != 'u201':
            continue
        hvol = []
        cvol = []
        for nframe, frame in enumerate(loop.tr.frames):
            if name != 'u201':
                continue
            hull_by_frame[name][frame]=CHT.hull_tool(loop)
            hull_by_frame[name][frame].make_hulls(frames=[nframe])
            hvol.append( hull_by_frame[name][frame].hull_volumes)
            cvol.append( hull_by_frame[name][frame].cell_volumes)
        hullvol[name]=nar(hvol).transpose()
        cellvol[name]=nar(cvol).transpose()

if 0:
    #
    # Volumes and ratios by time.  
    #
    for name in hull_by_frame:
        fig,ax = plt.subplots(1,3)
        if name != 'u201':
            continue
        for nparticle, vols in enumerate(zip(hullvol[name],cellvol[name])):
            hvol,cvol = vols
            ax[0].plot( loopers[name].tr.times, hvol/hvol[0])
            ax[1].plot( loopers[name].tr.times, cvol/cvol[0])
            rat = hvol/cvol
            rat /= rat[:-2].mean()

            ax[2].plot( loopers[name].tr.times, rat)
        axbonk(ax[0],yscale='log')#, ylim=[1e-9,1e-1])
        axbonk(ax[1],yscale='log')#, ylim=[1e-9,1e-1])
        axbonk(ax[2],yscale='log')
        fig.savefig('plots_to_sort/ratio_time.png')





if 0:
    fig,axess=plt.subplots(1,3)
    if 'ext_hull' not in dir():
        ext_hull=extents()
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        if nrun != 2:
            continue
        ax=axess[nrun]
        ax.clear()
        name = ht.this_looper.out_prefix
        ht.make_hulls(frames=[frame])
        odd = nar(ht.hull_volumes) < nar(ht.cell_volumes)
        not_odd = nar(ht.hull_volumes) >= nar(ht.cell_volumes)
        ax.scatter(ht.hull_volumes,ht.cell_volumes,c='k')
        ax.scatter(nar(ht.hull_volumes)[odd],nar(ht.cell_volumes)[odd],c='r')
        ax.set_aspect('equal')
        #axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=[0,0.07],ylim=[0,0.07])
        ext_hull(nar(ht.hull_volumes))
        ext_hull(nar(ht.cell_volumes))
        ax.plot(ext_hull.minmax, ext_hull.minmax,c='g')
        axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',
               xlim=ext_hull.minmax,ylim=ext_hull.minmax, xscale='log',yscale='log')
    fig.savefig("plots_to_sort/hull_volume_n%04d.png"%(frame))







