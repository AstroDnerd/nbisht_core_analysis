

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
        c=colors.color[ ht.this_looper.out_prefix]
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
    hull_by_frame = {}
    hullvol = defaultdict(list)
    cellvol = defaultdict(list)
    looper_list=[TLH.loops['u301'],TLH.loops['u302'],TLH.loops['u303']]
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







