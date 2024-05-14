
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import track_loader as TL
import movie_frames 

from scipy.optimize import root
def simple_rho(this_looper,core_list=None, tsing_tool=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    time_mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[time_mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    fig,axes=plt.subplots(2,1, figsize=(6,8))
    ax=axes[0];ax1=axes[1]#; ax2=axes[2]
    time=extents()
    rho_max_mean=0
    rho_avg_mean=0
    rho_min_mean=0
    tsing_mean=0
    ext = extents()
    for core_id in core_list:
        print('plot',core_id)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4


        rho = ms.density[sl].transpose()
        rho_max = rho[time_mask,:].max(axis=1)
        rho_max /= rho_max[0]
        rho_min = rho[time_mask,:].min(axis=1)
        rho_min /= rho_min[0]
        rho_avg = rho[time_mask,:].mean(axis=1)
        rho_avg /= rho_avg[0]
        rho_max_mean += rho_max[0]
        rho_min_mean += rho_min[0]
        rho_avg_mean += rho_avg[0]


        
        if tsing_tool is not None:
            tsing=tsing_tool.tsing_core[core_id]
            tsung=tsing_tool.tend_core[core_id]
            tsing_mean+=tsung
        else:
            tsung=1

        a=1.8614
        free_fall = (1-(times.flatten()/tsung)**2)**-a
        #ax.plot(times/tsung, free_fall)
        y1=rho_max#/free_fall
        ax.plot(times/tsung, y1, c=c)
        y2=rho_avg#/free_fall
        ax1.plot(times/tsung, y2, c=c)
        if 0:
            ax2.plot( times/tsung, free_fall)
        #ax2.plot(times, rho_min, c=c)
        ext(rho)
        #ext(nar([0.01,100]))

    time_lim=time.minmax
    time_lim=[0,1.2]
    ax.set(xlabel=r'$t/t_{\rm{SUNG}}$',    ylabel=r'$n_{\rm{max}}/n_{\rm{max}}(t=0)$',yscale='log', ylim=ext.minmax, xlim=time_lim)
    ax1.set(xlabel=r'$t/t_{\rm{SUNG}}$', ylabel=r'$n_{\rm{avg}}/n_{\rm{avg}}(t=0)$',yscale='log', ylim=ext.minmax, xlim=time_lim)
    if 0:
        ax2.set(xlabel=r'$t/t_{SUNG}$', ylabel=r'$\rm{min} \rho$',yscale='log', ylim=ext.minmax, xlim=time_lim)

    #pdb.set_trace()
    #rho_mean=rho_mean/len(core_list)
    #tsing_mean /= len(core_list)
    rho_max_mean /= len(core_list)
    rho_min_mean /= len(core_list)
    rho_avg_mean /= len(core_list)

    if 0:
        #check that the actual ODE is reproduced.
        b = 4*np.pi/3
        def f(u,x):
            return (u[1], -b/u[0]**2)
        tff = np.sqrt(3*np.pi/32)
        r0=1
        t=np.arange(0,tff,tff/1000)
        init=[1,0]
        from scipy.integrate import odeint
        sol=odeint(f,init,t)

        r_t = sol[:,0]
        rho_t = rho_mean/r_t**3
        ax.plot(t/tff, rho_t)
        ax1.plot(t/tff, rho_t)

    a=1.8614
    tff = np.sqrt(3*np.pi/32)
    t=np.arange(0,tff,tff/1000)

    r2 = (1-(t/tff)**2)**(a/3)
    rho_g_max = rho_max_mean/r2**3
    ax.plot(t/tff,rho_g_max,c='g')
    rho_g_avg = rho_avg_mean/r2**3
    ax1.plot(t/tff,rho_g_avg,c='g')
    rho_g_min = rho_min_mean/r2**3
    #ax2.plot(t/tff,rho_g_min,c='g')
    #pdb.set_trace()

    if 0:
        #play around to see how many orders you need to take.
        #turns out all of them.
        x=t/tff
        r3 = 1-a/3*(x)**2 + 1/18*(a-3)*a*(x)**4 - 1/162*((a-6)*(a-3)*a)*x**6 + (a-9)*(a-6)*(a-3)*a*x**3/1944
        rho3 = rho_mean/r3**3
        ax1.plot(t/tff,rho3,c='r')

    outname='plots_to_sort/%s_free_fall.pdf'%(this_looper.sim_name)
    fig.tight_layout()
    fig.savefig(outname)
    print(outname)




sims=['u501', 'u502','u503']
#sims=['u502']
import tsing
TL.load_tracks(sims)
tsing_tool = tsing.get_tsing(TL.loops)
for sim in sims:
    core_list = TL.loops[sim].core_by_mode['A']
    #core_list=[74]
    simple_rho(TL.loops[sim],core_list=core_list, tsing_tool=tsing_tool[sim])

