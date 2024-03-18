
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
    fig,axes=plt.subplots(1,2, figsize=(6,4))
    ax=axes[0];ax1=axes[1]
    time=extents()
    rho_mean=0
    tsing_mean=0
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
        rho = rho[time_mask,:].max(axis=1)


        
        if tsing_tool is not None:
            tsing=tsing_tool.tsing_core[core_id]
            tsung=tsing_tool.tend_core[core_id]
            tsing_mean+=tsung
        else:
            tsung=1
        ax.plot(times, rho, c=c)
        ax1.plot(times/tsung, rho, c=c)
        ax1.axvline(tsing/tsung)
        time(times/tsung)
        index = np.argmin( np.abs(times/tsung-1))
        rho_mean += rho[0]

    time_lim=time.minmax
    time_lim=[0,1.2]
    ax.set(xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max], xlim=time_lim)
    ax1.set(xlabel=r'$t/t_{SING}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max], xlim=time_lim)

    #pdb.set_trace()
    rho_mean=rho_mean/len(core_list)
    tsing_mean /= len(core_list)

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
    r2 = (1-(t/tff)**2)**(a/3)
    rho2 = rho_mean/r2**3
    ax1.plot(t/tff,rho2,c='g')

    x=t/tff
    r3 = 1-a/3*(x)**2 + 1/18*(a-3)*a*(x)**4 - 1/162*((a-6)*(a-3)*a)*x**6 + (a-9)*(a-6)*(a-3)*a*x**3/1944
    rho3 = rho_mean/r3**3
    ax1.plot(t/tff,rho3,c='r')

    outname='plots_to_sort/%s_free_fall.png'%(this_looper.sim_name)
    fig.savefig(outname)
    print(outname)




sims=['u501', 'u502','u503']
sims=['u502']
import tsing
TL.load_tracks(sims)
tsing_tool = tsing.get_tsing(TL.loops)
for sim in sims:
    core_list = TL.loops[sim].core_by_mode['Alone']
    simple_rho(TL.loops[sim],core_list=core_list, tsing_tool=tsing_tool[sim])

