from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class mass_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.unique_mass=defaultdict(list)
        self.total_mass=defaultdict(list)
        self.dof =defaultdict(list)
        self.cores_used=[]
    def run(self,core_list=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times


            for nf,frame in enumerate(thtr.frames):
                density = thtr.c([core_id],'density')[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                mask2 = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                unique_mass = (density[mask2]*cell_volume[mask2]).sum()
                total_mass =  (density*cell_volume).sum() 
                self.total_mass[core_id].append( total_mass)
                self.unique_mass[core_id].append(unique_mass)
                self.dof[core_id].append( mask2.sum())
            print("WWWWWW", self.unique_mass[core_id])


def plot_mass_tracks(tool, ax):
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    ntimes = tool.times.size
    ncores = len(tool.cores_used)
    masses = np.zeros([ntimes,ncores])
    these_times = tool.times/tff_global
    outname='plots_to_sort/mass_time_heatmap_unique_particles.pdf'
    for ncore,core_id in enumerate(tool.cores_used):
        this_mass = tool.unique_mass[core_id]
        if len(this_mass) == 0:
            pdb.set_trace()
        masses[:,ncore]=this_mass/nar(this_mass[:6]).mean()
        dof = nar(tool.dof[core_id])
        #masses[:,ncore] *= dof[0]/dof
        if np.isnan(masses).any():
            print("NANN")
            pdb.set_trace()
    nc = len(tool.cores_used)
    number_of_lines = min([10, len(tool.cores_used)])
    take_a_few = ((nc-1)*np.random.random(number_of_lines)).astype('int')
    print(take_a_few)
    for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]):
        rel_mass=masses[:,ncore]
        ax.plot(these_times, rel_mass,c=[0.5]*4)
        #axes[3].plot(these_times, tool.dof[core_id]/tool.dof[core_id][-1],c=[0.5]*4)
    mass_bins_edge = np.logspace(-3,4,101)
    mass_bins_edge = np.logspace(-2,1,101)
    xbins = these_times
    ybins = 0.5*(mass_bins_edge[1:]+mass_bins_edge[:-1])
    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    hist = np.zeros( [xbins.size,ybins.size])
    for ntime, time in enumerate(these_times):
        thishist,bins = np.histogram(masses[ntime,:],bins=mass_bins_edge)
        hist[ntime,:]=thishist
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = hist[hist>0].min()
    norm = mpl.colors.LogNorm(vmin=1,vmax=33)
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    #axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
    #axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
    axbonk(ax,ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')

