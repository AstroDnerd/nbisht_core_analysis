
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 
import multiplots
reload(multiplots)
def simple_rho(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    G = colors.G
    gx = thtr.track_dict['grav_x']
    gy = thtr.track_dict['grav_y']
    gz = thtr.track_dict['grav_z']
    vx = thtr.track_dict['velocity_x']
    vy = thtr.track_dict['velocity_y']
    vz = thtr.track_dict['velocity_z']
    GE2 = (gx*vx+gy*vy+gz*vz)
    v2 = vx*vx+vy*vy+vz*vz
    g2 = gx*gx+gy*gy+gz*gz
    GE2/=np.sqrt(v2*g2)
    ge_min=GE2.min()
    ge_max=GE2.max()
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)
        #fig,ax=multiplots.four_on_the_floor()
        #axlist = ax.flatten()

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        #ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]
        gx = thtr.c([core_id],'grav_x')[sl].transpose()[mask,:]
        gy = thtr.c([core_id],'grav_y')[sl].transpose()[mask,:]
        gz = thtr.c([core_id],'grav_z')[sl].transpose()[mask,:]
        vx = thtr.c([core_id],'velocity_x')[sl].transpose()[mask,:]
        vy = thtr.c([core_id],'velocity_y')[sl].transpose()[mask,:]
        vz = thtr.c([core_id],'velocity_z')[sl].transpose()[mask,:]
        v2 = vx*vx+vy*vy+vz*vz
        EG = (gx*gx+gy*gy+gz*gz)/(8*np.pi*colors.G)
        EK = 0.5*rho*v2



        names = ['density','EK','EG','EG/EK']
        arrz = [rho, EK, EG, EG/EK]
        scales=['log','log','log','log']
        for n , fname in enumerate(names):
            arr = arrz[n]
            yscale=scales[n]
            #axlist[n].plot(times,arr,c=c,linewidth=0.1)
            #axbonk(axlist[n], xlabel='t/tff', ylabel=fname, yscale=yscale)
            c = [0.5]*4
            c[n]=1.
            c0=[0.5]*4
            c1=[1.0,0.0,0.0,0.5]
            c2=[0.0,1.0,0.0,0.5]
            c3=[1.0,0.0,1.0,0.5]
            c=[c0,c1,c2,c3][n]
            ax.plot( times, arr, c=c, linewidth=0.1)

        axbonk(ax,xlabel='t/tff',ylabel='Q',yscale='log')

        #GE2 = (gx*vx+gy*vy+gz*vz)
        #GE2/=np.sqrt(v2*g2)

        #ax.plot(times , GE2, c=c, linewidth=0.1)
        #axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$(\nabla \phi)^2/8 pi G$',yscale='log', ylim=[ge_min,ge_max])
        #ax.set_yscale('symlog',linthresh=100)
        #ax2=ax.twinx()
        #c=[1.0,0.1,0.1,0.1]
        #ax2.plot(times , rho, c=c, linewidth=0.1)
        #axbonk(ax2,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log')
        #outname='plots_to_sort/%s_gdotrhov_t_c%04d.png'%(this_looper.sim_name,core_id)
        outname='plots_to_sort/%s_eng_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
        fig3,ax3=plt.subplots(1,1)
        ax3.boxplot( GE2)
        fig3.savefig('plots_to_sort/derp.png')





sims=['u501', 'u502','u503']
for sim in sims:
    core_list = np.unique(TL.loops[sim].tr.core_ids)
    #core_list=core_list[7:8] #[323]
    simple_rho(TL.loops[sim],core_list=core_list)


