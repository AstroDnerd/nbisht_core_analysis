
from starter2 import *
from collections import defaultdict
import scipy
import colors
reload(colors)

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def simple_rho(this_looper,core_list=None):

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
    for nc,core_id in enumerate(core_list):
        #print("%s %d"%(this_looper.sim_name,core_id))

        rho = thtr.c([core_id],'density').transpose()[mask,:]

        dv = thtr.c([core_id],'cell_volume').transpose()[mask,:]
        B=thtr.c([core_id],'magnetic_field_strength').transpose()[mask,:]/colors.mean_field[this_looper.sim_name]
        rho_mean = (rho*dv).sum(axis=1)/dv.sum(axis=1)
        RM[:,nc]=rho_mean
        B_mean = (B*dv).sum(axis=1)/dv.sum(axis=1)
        BM[:,nc]=B_mean

        Div = thtr.c([core_id],'velocity_divergence').transpose()[mask,:]
        mean_div = (Div*dv).sum(axis=1)/dv.sum(axis=1)
        DivV[:,nc]=mean_div

        ax.plot(rho_mean,B_mean,c=[0.5]*4,linewidth=0.1)
        axbonk(ax,xlabel='rho',ylabel='B',xscale='log',yscale='log')
    outname='plots_to_sort/%s_mean_b_and_rho.png'%(this_looper.sim_name)
    fig.savefig(outname)
    print(outname)
    plt.close('all')

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
    if 1:
        fig,ax=plt.subplots(1,1)
        ax.plot(times,RM,c=[1,0.5,0.5,0.5],linewidth=0.1)
        ax.plot(times,BM,c=[0.5,1.0,0.5,0.5],linewidth=0.1)
        #ax.plot(times,np.abs(DivV),c=[0.5,0.5,1.0,0.5],linewidth=0.1)
        ax.plot(times,BM/RM,c=[0.5,0.5,1.0,0.5],linewidth=0.1)

        if 1:
            X = nar([times.flatten()]*BM.shape[1]).transpose().flatten()
            Y = np.log((BM/RM).flatten())
            pfit = np.polyfit(X,np.exp(Y),1)
            ax.set_title('log10 B/Rho = B0 e^(%.1f t)'%(pfit[0]))
            #print(pfit)
            #ax.clear()
            #ax.scatter(X,10**Y)
            ax.plot(times, times*pfit[0]+pfit[1],c='k')
        if 0:
            X = nar([times.flatten()]*BM.shape[1]).transpose().flatten()
            Y = np.log((BM/RM).flatten())
            pfit = np.polyfit(X,Y,1)
            ax.set_title('log10 B/Rho = B0 e^(%.1f t)'%(pfit[0]))
            print(pfit)
            #ax.clear()
            #ax.scatter(X,10**Y)
            ax.plot(times, np.exp(times*pfit[0]+pfit[1]),c='k')


        #ax.plot(times,RM/BM**2,c=[0.5,0.5,1.0,0.5],linewidth=0.1)
        outname='plots_to_sort/b_and_rho_%s.png'%(this_looper.sim_name)
        axbonk(ax,xlabel='t/tff',ylabel='B,rho',yscale='log')
        fig.savefig(outname)
        print(outname)







plt.close('all')
sims=['u501', 'u502','u503']
for sim in sims:
    simple_rho(TL.loops[sim])#, core_list=[323])

