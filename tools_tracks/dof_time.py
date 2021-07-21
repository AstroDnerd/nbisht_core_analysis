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
                ix =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
                iy =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
                iz =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
                density = thtr.c([core_id],'density')[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                index = ix + nx*(iy * nx*iz)
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
    max_max=43494
    rm = rainbow_map(np.log10(max_max))
    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))
        these_times = tool.times/tff_global
        for ncore,core_id in enumerate(nar(tool.cores_used)):
            #rel_mass = (tool.mass[core_id]-tool.mass[core_id][0])/tool.mass[core_id][0]
            #rel_mass=tool.mass[core_id]/nrm#/tool.mass[core_id][-1]
            dof = nar(tool.dof[core_id])
            rel_dof = (dof - dof.min())/(dof.max()-dof.min())
            max_max=max([max_max,dof.max()])
            axes[nt].plot(these_times, rel_dof,c=rm(np.log10(dof[0])))
            axes[nt].set_xlabel(r'$t/t_{\rm{ff}}$')
    axes[0].set_ylabel(r'$\rm{Relative\ Degrees\ Of\ Freedom}$')


    outname='plots_to_sort/degrees_of_freedom.pdf'
    fig.savefig(outname)
    print(outname)

