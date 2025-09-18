
from starter2 import *
import track_loader as TL
import ucore
#reload(ucore)

sim_list=['u502']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

def trace_movie(mon,core_id,prefix='film_'):
    ms = mon.get_ms(core_id)
    ms.particle_pos(core_id)
    t = mon.all_times/colors.tff
    t.shape=t.size,1
    nskip=100
    rho = ms.density.transpose()[:,::nskip]
    px = ms.particle_x.transpose()[:,::nskip]
    py = ms.particle_y.transpose()[:,::nskip]

    cube_min_x = ms.this_x.min()
    cube_max_x = ms.this_x.max()
    cube_min_y = ms.this_y.min()
    cube_max_y = ms.this_y.max()
    cube_min_z = ms.this_z.min()
    cube_max_z = ms.this_z.max()
    left = np.array([cube_min_x,cube_min_y,cube_min_z])
    right= np.array([cube_max_x,cube_max_y,cube_max_z])



    timelim = [0,t.max()]
    rholim = [rho.min(),rho.max()]
    rho_norm = mpl.colors.LogNorm(vmin=rholim[0],vmax=rholim[1])
    for frame in mon.frames[::10]:
        nf = mon.get_frame_index(frame)
        fig,ax=plt.subplots(1,2)
        ax[0].plot(t[:nf,:], rho[:nf,:], c=[0.5]*4,zorder=2)
        ax[0].set(xlim=timelim,ylim=rholim,yscale='log')

        ax[1].plot(px[:nf,:],py[:nf,:],c=[0.5]*4,zorder=1)
        ax[1].scatter(px[nf,:],py[nf,:],c=rho[nf,:],norm=rho_norm,s=20,zorder=2,alpha=0.5, ec='None')



        fig.savefig('%s/%sf%04d'%(plot_dir,prefix,nf))
        plt.close(fig)




if 'things' not in dir():
    things={}

for sim in sim_list:
    last_etrack = ucore.etrack_list[-1]
    mon = monster.closet[sim]
    for uc in ucore.ucore_list[:1]:
        core_id = uc.core_id_by_et[last_etrack]
        trace_movie(mon,core_id,prefix='film_')
