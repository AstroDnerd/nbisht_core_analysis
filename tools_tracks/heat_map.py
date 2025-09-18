from starter2 import *
import colors

def heat_map(quantity, these_times, ax=None,bins=None, hist_density=False,hist_norm=False, zlim=None):
    if bins is None:
        bins = np.linspace( quantity.min(), quantity.max(), 64)
    xbins = these_times
    ybins = 0.5*(bins[1:]+bins[:-1])
    dv = bins[1:]-bins[:-1]
    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    hist = np.zeros( [xbins.size,ybins.size])
    for ntime, time in enumerate(these_times):
        thishist,bins = np.histogram(quantity[:,ntime],bins=bins,density=hist_density)
        if hist_norm:
            thishist= thishist/thishist.sum()
        hist[ntime,:]=thishist


    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = hist[hist>0].min()
    if zlim is None:
        norm = mpl.colors.LogNorm( vmin =minmin, vmax=hist.max())
    else:
        norm = mpl.colors.LogNorm( vmin=zlim[0], vmax=zlim[1])
    #norm = mpl.colors.LogNorm(vmin=1,vmax=33)
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    return TheX, TheY, hist, dv, ploot

def heat_for_quantity(this_looper, field='density',core_list=None, bins=None, external_ax=None):
    if core_list is None:
        core_list = np.unique( this_looper.tr.core_ids)

    times = this_looper.tr.times/colors.tff
    for core_id in core_list:
        Q = this_looper.tr.c([core_id],field)

        if bins is None:
            bins = np.linspace( Q.min(), Q.max(), 64)
        elif bins == 'plog':
            bins = np.geomspace( Q[Q>0].min(), Q.max(), 64)

        if external_ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            ax = external_ax
            ax.clear()
        heat_map( Q, times, bins=bins, ax=ax)
        if external_ax is None:
            fig.savefig('plots_to_sort/heat_%s_%s_c%04d.png'%(field,this_looper.sim_name,core_id))



def plot_heat(times=None, cores_used=None, quan_dict=None, ax=None, bins=None, norm=False, positive=True, fixer=None, cut_last=False, tool=None):
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    if tool is not None:
        times = tool.this_looper.tr.times
        if cores_used is None:
            cores_used = np.unique(tool.this_looper.tr.core_ids)
    ntimes = len(times)
    ncores = len(cores_used)
    quantity = np.zeros([ncores,ntimes])
    these_times = times/tff_global
    outname='plots_to_sort/mass_time_heatmap_unique_particles.pdf'

    for ncore,core_id in enumerate(cores_used):
        this_quan = nar(quan_dict[core_id])
        if norm:
            this_quan /= this_quan[:6].mean()

        quantity[ncore,:]=this_quan
    if fixer:
        quantity=fixer(quantity)
    nc = len(cores_used)
    number_of_lines = min([100, nc])
    if nc > 10:
        take_a_few = ((nc-1)*np.random.random(number_of_lines)).astype('int')
    else:
        take_a_few = slice(None)
    a_few = np.arange(nc)[take_a_few]
    if 0:
        for ncore in a_few:
            ax.plot(these_times, quantity[ncore,:],c=[0.5]*4)
    TheX, TheY, hist, dv, plot=heat_map(quantity,  these_times,bins=bins,ax=ax)
    #axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
    #axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
    #axbonk(ax,ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')
    stuff = {'quantity':quantity, 'TheX':TheX,'TheY':TheY,'hist':hist, 'dv':dv}
    return stuff

