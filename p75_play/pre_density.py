

from starter2 import *
import track_loader as TL
import ucore
#reload(ucore)


import radial_sphere 
reload(radial_sphere)
import radial_plots
reload(radial_plots)

import convex_hull_tools as CHT

sim_list=['u502']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

mon = monster.closet['u502']

if 'ht' not in dir():
    ht={}
    for sim in ucore.etrack_list:
        ht[sim] = CHT.hull_tool(TL.loops[sim])
        ht[sim].make_hulls()

times=[]
for net,et in enumerate(ucore.etrack_list):
    times.append( float(TL.tracks[et].tr.times[-1]))
times=nar(times)
import hull_density as hd
reload(hd)
if 1:
    ucore_volumes={}
    ucore_mean={}
    ucore_std={}
    for uc in ucore.ucore_list:
        print('uc',uc.uid)
        ucore_mean[uc.uid]=np.zeros(len(ucore.etrack_list))
        ucore_std[uc.uid]=np.zeros(len(ucore.etrack_list))
        ucore_volumes[uc.uid]=np.zeros(len(ucore.etrack_list))
        for net,et in enumerate(ucore.etrack_list):
            if et in uc.core_id_by_et:
                core_id = uc.core_id_by_et[et]
            else:
                continue
            rho_mean,rho_var =  hd.density(mon,ht[et].points_3d[core_id])
            ucore_mean[uc.uid][net]=rho_mean
            ucore_std[uc.uid][net]=rho_var
            ind = np.where(ht[et].cores_used == core_id)[0][0]
            volume = ht[et].hull_volumes[ind]
            ucore_volumes[uc.uid][net]=volume
if 1:
    fig,axes=plt.subplots(1,3)
    ax0=axes[0];ax1=axes[1];ax2=axes[2]
    cum=0
    for uc in ucore.ucore_list:
        uid = uc.uid
        if uc.uid not in ucore_mean:
            continue
        #cum = ucore_volumes[uc.uid] + cum
        #ax.plot(times/colors.tff,cum)
        ok = ucore_mean[uc.uid]>0
        ax0.plot(times[ok]/colors.tff, ucore_mean[uc.uid][ok])
        ok = ucore_std[uc.uid]>0
        ax1.plot(times[ok]/colors.tff, ucore_std[uc.uid][ok])
        ok = ucore_volumes[uc.uid]>0
        ax2.plot(times[ok]/colors.tff, ucore_volumes[uc.uid][ok])
    ax2.set(yscale='log')
    fig.savefig('%s/rho'%(plot_dir))


