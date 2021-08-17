from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')


class density_mean_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.peak_density=defaultdict(list)
        self.mean_density=defaultdict(list)
        self.masses=defaultdict(list)
        self.volumes=defaultdict(list)

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
                mass = (density[mask2]*cell_volume[mask2]).sum()
                volume = (cell_volume[mask2]).sum()
                self.peak_density[core_id].append( density.max())
                self.masses[core_id].append(mass)
                self.volumes[core_id].append(volume)
                self.mean_density[core_id].append( mass/volume)



import three_loopers_mountain_top as TLM

if 'dmt_u301' not in dir():
    dmt_u301 = density_mean_tool(TLM.loops['u301'])
    dmt_u301.run()

for dmt in [dmt_u301]:
    name = dmt.this_looper.out_prefix
    fig,ax=plt.subplots(2,2)
    ax0=ax[0][0]; ax1=ax[0][1]
    ax2=ax[1][0]; ax2=ax[1][1]

    peak_density = [ dmt.peak_density[core_id][-1] for core_id in dmt.cores_used]
    mean_density = [ dmt.mean_density[core_id][-1] for core_id in dmt.cores_used]

    ax0.scatter( peak_density, mean_density)
    axbonk(ax0,xscale='log',yscale='log',xlabel='peak_density',ylabel='mean_density')
    print('woo')

    fig.savefig('plots_to_sort/%s_peak_mean_density.png'%name)
