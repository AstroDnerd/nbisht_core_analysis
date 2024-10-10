
from starter2 import *
import colors

import hair_dryer
reload(hair_dryer)

import track_loader as TL
import movie_frames 

def simple_rho(this_looper,core_list=None, tsing_dict = None, tend_dict = None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    mask = np.ones_like(thtr.times,dtype='bool')
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()
    for core_id in core_list:
        fig,ax=plt.subplots(1,1)

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        ms.particle_pos(core_id)

        #if ms.nparticles < 1000:
        sl=slice(None)
        c=[0.5]*4
        #else:
        #    sl = slice(None,None,10)
        #    #c=[0,0,0,0.1]
        #    c=[0.1]*4

        rho = ms.density[sl].transpose()
        rho = rho[mask,:]

        ax.plot(times , rho, c=c, linewidth=0.1)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log', ylim=[rho_min,rho_max])

        #add tsing and tend lines
        if tsing_dict != None:
            if core_id in list(tsing_dict.keys()):
                ax.axvline(x=tsing_dict[core_id], c = 'green')
        
        if tend_dict != None:
            if core_id in list(tend_dict.keys()):
                ax.axvline(x=tend_dict[core_id], c = 'red')
            

        outname='plots_to_sort/%s_rho_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)


def run(trackname, core_list=None, tsing_dict = None, tend_dict = None):
    TL.load_tracks([trackname])
    this_looper=TL.tracks[trackname]
    simple_rho(this_looper, core_list=core_list, tsing_dict = tsing_dict, tend_dict = tend_dict)

if 0:
#sims=['u501', 'u502','u503']
    sims=['t002']
    TL.load_tracks(sims)
    for sim in sims:
        #core_list = TL.loops[sim].core_list[0:1]
        core_list=None
        simple_rho(TL.loops[sim], core_list=core_list)

