
from starter2 import *
from collections import defaultdict
import scipy
import colors
reload(colors)

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def div_v_hair(this_looper,core_list=None):

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
    B_all = thtr.track_dict['magnetic_field_strength']
    B_min=B_all.min()
    B_max=B_all.max()
    fig,ax=plt.subplots(1,1)
    RM=np.zeros([len(times),len(core_list)])
    BM=np.zeros([len(times),len(core_list)])
    DivV=np.zeros([len(times),len(core_list)])
    AllDivV = thtr.track_dict['velocity_divergence']
    MinDivV = AllDivV.min()
    MaxDivV = AllDivV.max()
    for nc,core_id in enumerate(core_list):
        ax.clear()
        #print("%s %d"%(this_looper.sim_name,core_id))


        dv = thtr.c([core_id],'cell_volume').transpose()[mask,:]
        Div = thtr.c([core_id],'velocity_divergence').transpose()[mask,:]
        mean_div = (Div*dv).sum(axis=1)/dv.sum(axis=1)
        DivV[:,nc]=mean_div
        ax.plot(times, Div, c=[0.5]*4,linewidth=0.1)

        axbonk(ax,xlabel='t',ylabel='DivV')#,xscale='log',yscale='log')
        ax.set_yscale('symlog',linthresh=250)
        ax.plot(times,times*0+250, c=[0.5]*4)
        ax.plot(times,times*0-250, c=[0.5]*4)
        ax.set_ylim([MinDivV,MaxDivV])
        outname='plots_to_sort/%s_divv_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
    plt.close('all')




plt.close('all')
sims=['u501', 'u502','u503']
for sim in sims:
    div_v_hair(TL.loops[sim])#, core_list=[323])
    #break

