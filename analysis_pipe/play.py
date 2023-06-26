
from starter2 import *

import track_loader as TL

track_list=['t002']
TL.load_tracks(track_list)
trackname='t002'
this_looper=TL.loops[trackname]
thtr=this_looper.tr
mountain_top_fname = track_info.tracks[trackname].mountain_top
if this_looper.targets is None:
    this_looper.read_targets_only(mountain_top_fname)


times=thtr.times
for core_id in this_looper.core_list:
    ms = trackage.mini_scrubber(thtr,core_id)
    peak = this_looper.targets[core_id]
    x,y,z=ms.raw_x,ms.raw_y,ms.raw_z
    if 0:
        print( np.abs((x[:,-1] - peak.peak_location[0])).max()*2048)
        print( np.abs(y[:,-1] - peak.peak_location[1]).max()*2048)
        print( np.abs(z[:,-1] - peak.peak_location[2]).max()*2048)
    print(ms.density[:,-1])
    fig,ax=plt.subplots(1,1)
    ax.plot(times, ms.density.transpose())
    ax.set(yscale='log')
    fig.savefig('plots_to_sort/fu_c%04d'%core_id)

    





