from starter2 import *


def plot_heat(tool=None, quan_dict=None, ax=None, bins=None, norm=False, positive=True, fixer=None):
    loop=tool.this_looper
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    ntimes = len(tool.times)
    ncores = len(tool.cores_used)
    quantity = np.zeros([ncores,ntimes])
    these_times = tool.times/tff_global
    outname='plots_to_sort/mass_time_heatmap_unique_particles.pdf'

    for ncore,core_id in enumerate(tool.cores_used):
        this_quan = nar(quan_dict[core_id])
        if norm:
            this_quan /= this_quan[:6].mean()

        quantity[ncore,:]=this_quan
    if fixer:
        quantity=fixer(quantity)
    nc = len(tool.cores_used)
    number_of_lines = min([10, nc])
    take_a_few = ((nc-1)*np.random.random(number_of_lines)).astype('int')
    #a_few = nar(tool.cores_used)[take_a_few]
    a_few = np.arange(nc)[take_a_few]
    #for ncore in a_few:
    #    ax.plot(these_times, quantity[ncore,:],c=[0.5]*4)
        #axes[3].plot(these_times, tool.dof[core_id]/tool.dof[core_id][-1],c=[0.5]*4)
    if bins is None:
        bins = np.linspace( quantity.min(), quantity.max(), 64)
    xbins = these_times
    ybins = 0.5*(bins[1:]+bins[:-1])
    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    hist = np.zeros( [xbins.size,ybins.size])
    for ntime, time in enumerate(these_times):
        thishist,bins = np.histogram(quantity[:,ntime],bins=bins)
        hist[ntime,:]=thishist
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    #minmin = hist[hist>0].min()
    norm = mpl.colors.LogNorm(vmin=1,vmax=33)
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    #axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
    #axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
    #axbonk(ax,ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')
    return quantity

