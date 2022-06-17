
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def velocity_things(this_looper,core_list=None, do_plots=True):

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

    sigma_3d = np.zeros([len(times), len(core_list)])
    vr_all   = np.zeros_like(sigma_3d)
    vt_all   = np.zeros_like(sigma_3d)
    for nc,core_id in enumerate(core_list):
        print('V %s %d'%(this_looper.sim_name,core_id))
        fig,ax=plt.subplots(1,1)
        ax2=ax.twinx()

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]

        dv = ms.cell_volume[sl].transpose()[mask,:]
        vv = dv.sum(axis=1)
        vx = ms.rel_vx[sl].transpose()[mask,:]
        vy = ms.rel_vy[sl].transpose()[mask,:]
        vz = ms.rel_vz[sl].transpose()[mask,:]
        v22 = ms.rel_vmag[sl].transpose()[mask,:]

        vr = ms.vr_rel[sl].transpose()[mask,:]
        vt = (ms.vt2_rel[sl].transpose()[mask,:])**0.5

        vrm=vr.mean(axis=1)
        v2 = v22.mean(axis=1)
        vtm=vt.mean(axis=1)

        sigma_3d[:,nc]=v2
        vr_all[:,nc]=vrm
        vt_all[:,nc]=vtm
        if do_plots:
            if 0:
                #4th moment.  Doesn't do much.
                Q=4
                v4 = (((vx**Q+vy**Q+vz**Q)).mean(axis=1))**(1/Q)
                ax2.plot(times, v4, c='g')
            #v4 = ((vx**4+vy**4+vz**4)**0.25*dv).sum(axis=1)/vv
            ax2.plot(times, v22.mean(axis=1), c='k', label=r'$\sigma_v$')
            if 1:
                ok = vrm>0
                ax2.plot(times[ok], vrm[ok], 'r--')
                ax2.plot(times[~ok], np.abs(vrm[~ok]), c='r',label=r'$v_r$')
            ax2.plot(times, vtm, c='c', label=r'$v_t$')
            axbonk(ax2, ylim=[0,10], ylabel=r'$(V^q)^1/q$')

            ax.plot(times , rho, c=c, linewidth=0.1)
            axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

            outname='plots_to_sort/%s_rho_vnorm_t_c%04d.png'%(this_looper.sim_name,core_id)
            fig.savefig(outname)

    return {'sigma_3d':sigma_3d,'vr':vr_all,'vt':vt_all, 'times':times}




sims=['u501', 'u502','u503']
if 'vels' not in dir():
    vels={}
    for sim in sims:
        vels[sim]=velocity_things(TL.loops[sim], do_plots=True)#, core_list=[323])

for sim in vels:
    fig,ax=plt.subplots(1,3)
    fig.subplots_adjust(wspace=0)

    ax[0].plot( vels[sim]['times'], vels[sim]['sigma_3d'], c=[0.5]*4, linewidth=0.3)
    ax[1].plot( vels[sim]['times'], np.abs(vels[sim]['vr']), c=[0.5]*4, linewidth=0.3)
    ax[2].plot( vels[sim]['times'], vels[sim]['vt'], c=[0.5]*4, linewidth=0.3)
    ax[0].set_title(r'$\sigma_{3D}$')
    ax[1].set_title(r'$v_r$')
    ax[2].set_title(r'$v_t$')
    ax[0].set_ylim([0,11])
    ax[1].set_ylim(ax[0].get_ylim())
    ax[2].set_ylim(ax[0].get_ylim())
    ax[1].set_ylabel('')
    ax[2].set_ylabel('')
    ax[1].set_yticks([])
    ax[2].set_yticks([])
    ax[0].set_xlabel(r'$t/t_{\rm{ff}}$')
    ax[1].set_xlabel(r'$t/t_{\rm{ff}}$')
    ax[2].set_xlabel(r'$t/t_{\rm{ff}}$')
    for aaa in ax:
        aaa.plot( vels[sim]['times'], vels[sim]['times']*0+1,'k--')
    fig.savefig('plots_to_sort/vels_%s'%sim)



