from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class mass_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.mass=defaultdict(list)
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
            if ms.nparticles < 3:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            self.times = thtr.times

            for nf,frame in enumerate(thtr.frames):
                x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]#or whatever the number of zones is
                y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
                z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
                density = thtr.c([core_id],'density')[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                index = x + nx*(y * nx*z)
                ar = np.argsort(index)
                rs = np.argsort(ar)
                isorted=index[ar]
                mask = np.ones_like(density,dtype='bool')
                mask[1:] = isorted[1:]-isorted[:-1] != 0
                mask2 = mask[ rs]
                mass = (density[mask2]*cell_volume[mask2]).sum()
                self.mass[core_id].append(mass)

import three_loopers as tl
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

fig,ax=plt.subplots(2,2)
axes=ax.flatten()
if 1:
    for nt,tool in enumerate([mass_tool1,mass_tool2,mass_tool3]):



        for core_id in [9]: #tool.cores_used:
            rel_mass = (tool.mass[core_id]-tool.mass[core_id][0])/tool.mass[core_id][0]
            if rel_mass[8] > 0:
                print("got one", core_id)
                print(tool.this_looper.tr.c([core_id],'density').shape)


            axes[nt].plot(tool.times, rel_mass,c=[0.5]*4)
        axes[nt].set_yscale('symlog',linthresh=1)

    fig.savefig('plots_to_sort/mass_time.pdf')

if 0:
    for nt,tool in enumerate([mass_tool1]):
        rel_mass = (1-tool.mass[core_id]/tool.mass[core_id][0])

        for core_id in tool.cores_used[4::10]:
            axes[nt].plot(tool.times,x  ,c=[0.5]*4)
        axes[nt].set_yscale('symlog',linthresh=1)

    fig.savefig('plots_to_sort/mass_time.pdf')

