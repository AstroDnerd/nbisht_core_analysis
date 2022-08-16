
from starter2 import *
from collections import defaultdict
import scipy
import colors

from scipy.ndimage import gaussian_filter
import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def v_hair(this_looper,core_list=None, suffix=''):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    fig,axes=plt.subplots(3,1)
    fig.subplots_adjust(hspace=0)
    ax=axes[0]; ax1=axes[1]; ax2=axes[2]
    for core_id in core_list:
        print('do',core_id)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

        sl=slice(None)
        c=[0.5]*4
        c='k'

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]

        dv = ms.cell_volume[sl].transpose()[mask,:]
        vv = dv.sum(axis=1)
        #vx = ms.rel_vx[sl].transpose()[mask,:]
        #vy = ms.rel_vy[sl].transpose()[mask,:]
        #vz = ms.rel_vz[sl].transpose()[mask,:]
        v22 = ms.rel_vmag[sl].transpose()[mask,:]

        vr = ms.vr_rel[sl].transpose()[mask,:]
        vt = (ms.vt2_rel[sl].transpose()[mask,:])**0.5

        vrm=vr.mean(axis=1)
        v2 = v22.mean(axis=1)
        vtm=vt.mean(axis=1)

        UB = gaussian_filter(rho.max(axis=1),1)
        tf = times.flatten()
        dt = tf[1:]-tf[:-1]
        dU = UB[1:]-UB[:-1]
        tc = 0.5*(times[1:]+times[:-1])
        dU = dU[1:-1]
        dt = dt[1:-1]
        tc = tc[1:-1]
        dudt=dU/dt
        thresh = 1e5
        singularity = np.where( dudt >thresh)[0][0]
        tsing = times[singularity]

        #velocity plots
        if 1:
            ok = vrm>0
            ax1.plot(times[ok]/tsing, vrm[ok], 'r--', linewidth=0.1)
            ax1.plot(times[~ok]/tsing, np.abs(vrm[~ok]), c='r',label=r'$v_r$', linewidth=0.1)
        ax2.plot(times/tsing, vtm, c='c', label=r'$v_t$', linewidth=0.1)
        ax.plot(times/tsing, v2, c=c, label=r'$v$', linewidth=0.1)
    labels = [r'$v_{rms}$', r'$v_{radial}$', r'$v_{tangent}$']
    for na,aaa in enumerate(axes):
        aaa.set( ylim=[0,11], ylabel=labels[na])
        if na==2:
            aaa.set(xlabel=r'$t/t_{singularity}$')
        else:
            aaa.set(xlabel='', xticks=[])
        aaa.axhline(1,c=[0.5]*4)
        aaa.axvline(1, c=[0.5]*4)



    outname='plots_to_sort/%s_velocity_mean_hair_t%s.png'%(this_looper.sim_name,suffix)
    fig.savefig(outname)
    print(outname)




sims=['u501', 'u502','u503']
#sims=[ 'u502','u503']
for sim in sims:
    for mode in ['One','Binary','Cluster']:
        core_list = TL.loops[sim].core_by_mode[mode]
        v_hair(TL.loops[sim], core_list=core_list,suffix=mode)

