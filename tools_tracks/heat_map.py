from starter2 import *

def heat_map(quantity, these_times, ax=None,bins=None):
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
        thishist,bins = np.histogram(quantity[:,ntime],bins=bins,density=False)
        hist[ntime,:]=thishist


    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = hist[hist>0].min()
    norm = mpl.colors.LogNorm( vmin =minmin, vmax=hist.max())
    #norm = mpl.colors.LogNorm(vmin=1,vmax=33)
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    return TheX, TheY, hist, dv

def plot_heat(times, cores_used, quan_dict=None, ax=None, bins=None, norm=False, positive=True, fixer=None, cut_last=False):
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
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
    TheX, TheY, hist, dv=heatmap(quantity,  these_times,bins=bins,ax=ax)
    #axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
    #axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
    #axbonk(ax,ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')
    stuff = {'quantity':quantity, 'TheX':TheX,'TheY':TheY,'hist':hist, 'dv':dv}
    return stuff

