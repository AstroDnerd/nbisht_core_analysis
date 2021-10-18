from starter2 import *

import three_loopers_mountain_top as TLM

simlist = ['u301','u302','u303']

fig,ax = plt.subplots(1,1)

for sim in simlist:
    loop = TLM.loops[sim]
    if loop.targets is None:
        mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%sim
        loop.read_targets(mountain_top_fname)
    peaks = nar([loop.targets[core_id].peak_density for core_id in loop.core_list])
    ax.hist(np.log10(peaks), histtype='step')
fig.savefig('plots_to_sort/peak_thing.png')



