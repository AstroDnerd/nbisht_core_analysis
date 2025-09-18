
from starter2 import *
from collections import defaultdict
import scipy
import colors

from scipy.ndimage import gaussian_filter
import hair_dryer
reload(hair_dryer)

#import three_loopers_u500 as TL
import track_loader as TL
import movie_frames 

def v_hair(this_looper,core_list=None, suffix='', norm=False, tsing_in=None):

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
    fig,axes=plt.subplots(4,1, figsize=(6,6))
    fig.subplots_adjust(hspace=0)
    ax=axes[1]; ax1=axes[2]; ax2=axes[3]; ax3=axes[0]
    #ax=axes[4]; ax1=axes[2]; ax2=axes[3]; ax3=axes[0]; ax4=axes[1]
    rho_ext=extents()
    eng_ext=extents()
    fig2,axB=plt.subplots(1,1)
    mean_rho_mean=[]
    for core_id in core_list:
        print('do',core_id)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)
        ms.compute_ke_rel(core_id)
        ms.compute_ge(core_id)


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

        ke = ms.ke_rel[sl].transpose()[mask,:]
        ge = ms.ge[sl].transpose()[mask,:]
        ketotal = (ke*dv).sum(axis=1)
        getotal = np.abs((ge*dv).sum(axis=1))
        if 1:
            eng_ext(ketotal)
            eng_ext(getotal)
        if 0:
            eng_ext(getotal/ketotal)

        rhomean = rho.mean(axis=1)
        rhomax = rho.max(axis=1)
        rhomin = rho.min(axis=1)
        mean_rho_mean.append(rhomean[0])
        rho_ext(rhomean)
        rho_ext(rhomin)

        #tsing = tsing_in.tend_core[core_id]
        tsing = tsing_in.tsing_core[core_id]


        #velocity plots
        ax.plot(times/tsing, v2, c=c, label=r'$v$', linewidth=0.1)
        if 1:
            ok = vrm>0
            ax1.plot(times[ok]/tsing, vrm[ok], 'r--', linewidth=0.1)
            ax1.plot(times[~ok]/tsing, np.abs(vrm[~ok]), c='r',label=r'$v_r$', linewidth=0.1)
        ax2.plot(times/tsing, vtm, c='c', label=r'$v_t$', linewidth=0.1)
        #ax3.plot(times/tsing, rhomin, c='r', linewidth=0.1)
        ax3.plot(times/tsing, rhomean, c='g', linewidth=0.1)
        #ax3.plot(times/tsing, rhomax, c='b', linewidth=0.1)

        #ax4.plot(times/tsing, ketotal,c='r', linewidth=0.1)
        #ax4.plot(times/tsing, getotal,c='g', linewidth=0.1)
        #ax4.plot(times/tsing, getotal/ketotal,c='g', linewidth=0.1)
        axB.plot(getotal/getotal[0],ketotal/ketotal[0])

    test_t = np.arange(0,1,0.01)
    a=1.8614
    mean_rho_mean=nar(mean_rho_mean)
    if 0:
        for rho_mean in [mean_rho_mean.min(),mean_rho_mean.max()]:
            #rho_mean = np.mean(mean_rho_mean)
            rhot = (1-test_t**2)**(-a)
            rho_ff = rho_mean*rhot
            ax3.plot(test_t,rho_ff)

    axB.set(xscale='log',yscale='log', xlim=eng_ext.minmax,ylim=eng_ext.minmax)
    fig2.savefig('plots_to_sort/grasping')
    labels = [ r'$\overline{v_{\rm{rms}}}/c_s$', r'$|\overline{v_{\rm{R}}}|/c_s$', r'$|\overline{v_{\rm{T}}}|/c_s$']
    for na,aaa in enumerate(axes[1:]):
        ylim = [0,7.8]
        xlim = [0,1.2]
        aaa.set( ylim=ylim, ylabel=labels[na], xlim=xlim)
        if na==2:
            aaa.set(xlabel=r'$t/t_{\rm{sing}}$')
        else:
            aaa.set(xlabel='', xticks=[])
        aaa.axhline(1,c=[0.5]*4)
        aaa.axvline(1, c=[0.5]*4)
    ax3.axvline(1,c=[0.5]*4)
    ax3.set(yscale='log', ylim=rho_ext.minmax, xlim=ax.get_xlim(), ylabel=r'$\overline{n}~[\rm{cm^{-3}}]$', xticks=[])
    #ax4.set(yscale='log',ylim=eng_ext.minmax)



    outname='plots_to_sort/%s_velocity_mean_hair_t%s.pdf'%(this_looper.sim_name,suffix)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)

    fig.savefig(outname)
    print(outname)


if 1:
    #First Paper version.
    #don't touch.
    #probably replaced.
    sims=['u501']
    TL.load_tracks(sims)
    #sims=[ 'u502','u503']
    for sim in sims:
        for mode in ['A']:
            core_list = TL.loops[sim].core_by_mode[mode]
            #core_list = np.unique(TL.loops[sim].tr.core_ids)

            #core_list=core_list[:5]

            v_hair(TL.loops[sim], core_list=core_list,suffix=mode,norm=False, tsing_in=tsing_tool[sim])

if 0:
    #play with tsing
    #don't touch.
    sims=['u502']
    TL.load_tracks(sims)
    #sims=[ 'u502','u503']
    import tsing
    tsing_tool=tsing.get_tsing(TL.loops)
    for sim in sims:
        for mode in ['A']:
            core_list = TL.loops[sim].core_by_mode[mode]
            #core_list = np.unique(TL.loops[sim].tr.core_ids)

            #core_list=core_list[:5]

            v_hair(TL.loops[sim], core_list=core_list,suffix=mode,norm=False, tsing_in=tsing_tool[sim])


if 0:
    #loop over everything.  
    #don't play with this one.
    sims=['u501', 'u502','u503']
#sims=[ 'u502','u503']
    for sim in sims:
        for mode in ['Alone','Binary','Cluster']:
            core_list = TL.loops[sim].core_by_mode[mode]
            v_hair(TL.loops[sim], core_list=core_list,suffix=mode)

