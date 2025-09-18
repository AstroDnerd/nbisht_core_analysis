
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import track_loader as TL

def simple_hair(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    all_times=thtr.times
    all_frames=thtr.frames
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff

    for core_id in core_list:
        fig,axes=plt.subplots(3,1, figsize=(6,10))
        ax=axes[0];ax1=axes[2]; ax3=axes[1]
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)

            
        LOS = 0
        x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
        y = [2,2,0][LOS] # unfolds nicely.

        if False:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,30)
            #c=[0,0,0,0.1]
            c=[0.1]*4
            Linewidth1=Linewidth2=0.2

        print(sl)
        print(c)
        rho = ms.density[sl].transpose()#*colors.density_units
        rho = rho[mask,:]

        dv = ms.cell_volume[sl].transpose()[mask,:]
        vv = dv.sum(axis=1)
        vx = ms.rel_vx[sl].transpose()[mask,:]
        vy = ms.rel_vy[sl].transpose()[mask,:]
        vz = ms.rel_vz[sl].transpose()[mask,:]

        v22_all = ms.rel_vmag[:].transpose()[mask,:]
        vr_all =  ms.vr_rel[:].transpose()[mask,:]
        vt_all = (ms.vt2_rel[:].transpose()[mask,:])**0.5

        vrm=vr_all.mean(axis=1)
        v2 = v22_all.mean(axis=1)
        vtm=vt_all.mean(axis=1)

        
        rho_plot=ax1.twinx()
        print(rho.shape,c)
        rho_plot.plot(times, rho*colors.density_units, c=c, linewidth=Linewidth1)
        rho_plot.set(yscale='log',ylabel=r'$n~ [cm^{-3}]$')

        ax1.plot(times, v2, c='k')
        ax1.plot(times, vtm, c='c')
        ax1.plot(times, np.abs(vrm), c='r')
        ax1.set(ylabel=r'$v/c_s$', xlabel=r'$t/t_{\rm{ff}}$')



        p = [ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(),ms.particle_z[sl].transpose()]

        for aaa in [ax,ax3]:
            aaa.scatter( p[x][0,:].flatten(),p[y][0,:].flatten(),c='k',s=0.1)
            aaa.scatter( p[x][-1,:].flatten(),p[y][-1,:].flatten(),c='r',s=0.1)
            aaa.plot( p[x], p[y], c=c, linewidth=0.3)
            aaa.set(xlabel=r'$z/L$', ylabel=r'$y/L$')

        x0,x1=[0.090,0.175]
        y0,y1=[0.15,0.25]
        ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],c='r',linewidth=Linewidth1)
        ax3.set(xlim=[x0,x1],ylim=[y0,y1])
        outname='plots_to_sort/%s_mosh_%s_c%04d.png'%(this_looper.sim_name,'xyz'[LOS],core_id)
        fig.tight_layout()
        fig.savefig(outname)
        print(outname)




sims=[ 'u502']
for sim in sims:
    core_list=[9]
    print('word')
    simple_hair(TL.loops[sim],core_list=core_list)

