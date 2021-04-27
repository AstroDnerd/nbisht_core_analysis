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
        core_list=all_cores
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
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

import three_loopers_1tff as tl
if 'clobber' not in dir():
    clobber=True
if 'mass_tool1' not in dir() or clobber:
    mass_tool1=mass_tool(tl.looper1)
    mass_tool1.run()

if 'mass_tool2' not in dir() or clobber:
    mass_tool2=mass_tool(tl.looper2)
    mass_tool2.run()
if 'mass_tool3' not in dir() or clobber:
    mass_tool3=mass_tool(tl.looper3)
    mass_tool3.run()

fig,ax=plt.subplots(1,3, figsize=(12,4))
axes=ax.flatten()
if 1:
    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        ntimes = tool.times.size
        ncores = len(tool.cores_used)
        masses = np.zeros([ntimes,ncores])
        these_times = tool.times/tff_global
        outname='plots_to_sort/mass_time_heatmap_unique_particles.pdf'
        this_mass = tool.unique_mass[core_id]
        for ncore,core_id in enumerate(tool.cores_used):
            masses[:,ncore]=this_mass/nar(this_mass[:6]).mean()

        nc = len(tool.cores_used)
        take_a_few = ((nc-1)*np.random.random(10)).astype('int')
        print(take_a_few)
        for ncore,core_id in enumerate(nar(tool.cores_used)[take_a_few]):
            rel_mass=masses[:,ncore]

            axes[nt].plot(these_times, rel_mass,c=[0.5]*4)
            #axes[3].plot(these_times, tool.dof[core_id]/tool.dof[core_id][-1],c=[0.5]*4)

        mass_bins_edge = np.logspace(-3,4,101)
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
        ploot=axes[nt].pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
        axes[nt].plot(these_times, [2]*tool.times.size,c='r')#,lw=0.2)
        axes[nt].plot(these_times, [1./2]*tool.times.size,c='r')#,lw=0.2)
        axbonk(axes[nt],ylabel=None,xlabel=r'$t/t_{\rm{ff}}$',yscale='log')
    axes[0].set_ylabel(r'$M/\bar{M_{\rm{early}}}$')
    fig.colorbar(ploot)
    fig.savefig(outname)
    print(outname)

