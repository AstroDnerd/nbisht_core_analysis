from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class mass_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.dof =defaultdict(list)
        self.cores_used=[]
        self.mass=defaultdict(list)
    def run(self,core_list=None):
        #dx=1./2048
        dx=1/1024
        #dx=1./128
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles < 10:
                continue
            print('mass ',self.this_looper.sim_name, core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times

            nzones=nar([0]+list(np.cumsum((128*2**np.array([0,1,2,3,4]))**3)))
            for nf,frame in enumerate(thtr.frames):
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                dx = cell_volume**(1./3)
                nx = 1/dx
                dx.shape = dx.size,1
                ix =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
                iy =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
                iz =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
                density = thtr.c([core_id],'density')[:,nf]
                dx.shape=dx.size
                level = ((np.log(dx)+np.log(128))/np.log(1/2)).astype('int')
                #print(level)
                index = (ix + nx*(iy * nx*iz)+nzones[level]).astype('int')
                ar = np.argsort(index)
                rs = np.argsort(ar)
                isorted=index[ar]
                mask = np.ones_like(density,dtype='bool')
                mask[1:] = isorted[1:]-isorted[:-1] != 0
                mask2 = mask[ rs]
                #mass = (density[mask2]*cell_volume[mask2]).sum()
                mass = (density*cell_volume).sum()
                self.mass[core_id].append(mass)
                self.dof[core_id].append( mask2.sum())

import track_loader as TL
sims=['u501']#,'u502','u503']
TL.load_tracks(sims)
import tsing
#reload(tsing)
tsing_tool=tsing.get_tsing(TL.loops)
if 'clobber' not in dir():
    clobber=False
if 'MT' not in dir() or clobber:
    MT={}
    for sim in sims:
        temp=mass_tool(TL.loops[sim])
        core_list = TL.loops[sim].core_by_mode[mode]
        core_list=None
        temp.run(core_list=core_list)
        MT[sim]=temp

if 1:
    max_max=43494
    rm = rainbow_map(np.log10(max_max))
    fig,axes=plt.subplots(3,1)
    for nt,sim in enumerate(sims):
        tool=MT[sim]
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        these_times = tool.times
        for mode in ['Alone']:
            core_list = tool.cores_used
            for ncore,core_id in enumerate(core_list):
                #rel_mass = (tool.mass[core_id]-tool.mass[core_id][0])/tool.mass[core_id][0]
                #rel_mass=tool.mass[core_id]/nrm#/tool.mass[core_id][-1]
                ts=tsing_tool[sim].tsing_core[core_id]
                these_times = tool.times/colors.tff/ts
                dof = nar(tool.dof[core_id])
                rel_dof = (dof - dof.min())/(dof.max()-dof.min())
                max_max=max([max_max,dof.max()])
                axes[nt].plot(these_times, rel_dof,c=[0.5]*4, linewidth=0.1)
        axes[nt].set(xlabel=r'$t/t_{\rm{ff}}$', xlim=[0,2])
        axes[nt].axvline(1)
    axes[0].set_ylabel(r'$\rm{Relative\ Degrees\ Of\ Freedom}$')


    outname='plots_to_sort/degrees_of_freedom.pdf'
    fig.savefig(outname)
    print(outname)

