
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
    DivV=np.zeros([len(times),len(core_list)])
    for nc,core_id in enumerate(core_list):

        dv = thtr.c([core_id],'cell_volume').transpose()[mask,:]
        Div = thtr.c([core_id],'velocity_divergence').transpose()[mask,:]
        mean_div = (Div*dv).sum(axis=1)/dv.sum(axis=1)
        DivV[:,nc]=mean_div

    fig,ax=plt.subplots(1,1)
    DivVall=thtr.track_dict['velocity_divergence']
    DivVall_early=DivVall[:,:50].flatten()
    YYY = np.mgrid[0:1:DivVall_early.size*1j]
    print(YYY.size)
    ax.plot(sorted(DivVall_early),YYY)
    outname='plots_to_sort/%s_divv_early_hist.png'%(this_looper.sim_name)
    axbonk(ax,xlim=[-1000,1000])
    fig.savefig(outname)
    print(outname)

    if 0:
        fig,ax=plt.subplots(1,1)
        rho_min=RM.min()
        rho_max=RM.max()
        B_min=BM.min()
        B_max=BM.max()/2
        print(B_min,B_max)
        for iframe in range(times.size):
            ax.clear()
            ax.plot( RM[:iframe,:], BM[:iframe,:], c=[0.5]*4,linewidth=0.1)
            axbonk(ax,xlabel='rho',ylabel='B',xscale='log',yscale='log',xlim=[rho_min,rho_max],ylim=[B_min,B_max])
            outname='plots_to_sort/MeanBvsRho_%s_i%04d.png'%(this_looper.sim_name,iframe)
            fig.savefig(outname)
            print(outname)
    if 0:
        fig,ax=plt.subplots(1,1)
        ax.plot(times,RM,c=[1,0.5,0.5,0.5],linewidth=0.1)
        ax.plot(times,BM,c=[0.5,1.0,0.5,0.5],linewidth=0.1)
        ax.plot(times,np.abs(DivV),c=[0.5,0.5,1.0,0.5],linewidth=0.1)
        #ax.plot(times,BM/RM,c=[0.5,0.5,1.0,0.5],linewidth=0.1)

        X = nar([times.flatten()]*BM.shape[1]).transpose().flatten()
        Y = np.log((BM/RM).flatten())
        pfit = np.polyfit(X,np.exp(Y),1)
        ax.set_title('log10 B/Rho = B0 e^(%.1f t)'%(pfit[0]))
        #print(pfit)
        #ax.clear()
        #ax.scatter(X,10**Y)
        ax.plot(times, times*pfit[0]+pfit[1],c='k')


        #ax.plot(times,RM/BM**2,c=[0.5,0.5,1.0,0.5],linewidth=0.1)
        outname='plots_to_sort/b_and_rho_%s.png'%(this_looper.sim_name)
        axbonk(ax,xlabel='t/tff',ylabel='B,rho',yscale='log')
        fig.savefig(outname)
        print(outname)







plt.close('all')
sims=['u501', 'u502','u503']
for sim in sims:
    div_v_hair(TL.loops[sim], core_list=[323])

