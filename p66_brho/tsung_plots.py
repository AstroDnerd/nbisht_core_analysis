from starter2 import *
import xtra_energy

from scipy.optimize import curve_fit
from scipy import stats
from scipy.ndimage import gaussian_filter
import core_proj_three
reload(core_proj_three)
import other_scrubber
reload(other_scrubber)
#import three_loopers_six as TL
import camera_path
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
#sim_list=['u502']
plt.close('all')

            
def replotter(obj,suffix='', redlines=False, subset=0):
    OneExt = False
    suffix=""

    if subset==1:
        #synthetic los obs
        row_dict={'b_los_cyl':0,'b_los_r1':1,'b_los_r2':2}
        profs=['b_los_cyl','b_los_r1','b_los_r2']
    if subset==2:
        #synthetic obs try one first
        row_dict={'B_los_r1':0, 'N_r1':0,'B_los_r1log':0, 'N_r1log':0}
        profs=['B_los_r1','N_r1'] 

    Nplots=len(profs) -1
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(16,5))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles,axes): #axes[0] 
        ax.set_title(title)

    ext = [extents() for n in range(len(profs))]
    args = {'linewidth':0.2, 'c':'green', 'alpha':0.4} #'c':[0.5]*3


    # FOR EACH PROFILE
    #for nprofile, profile in enumerate(profs):
    profiles = obj.profiles_gas
    core_list=list(profiles.keys())
    
    # CORES
    Nlog = defaultdict(list)
    Blog = defaultdict(list)
    kappas = []
    for core_id in core_list:
        frames = sorted(list(profiles[core_id].keys()))
        row = row_dict['N_r1']
        
        # FRAMES
        for nframe,frame in enumerate(frames):
            ax = axes[nframe]
            #args['c']=[0.5]*3

            X = nar(profiles[core_id][frame]['N_r1'])
            Y = nar(profiles[core_id][frame]['B_los_r1']) 

            ext[-1](X)
            ext[row](Y)
            ax.scatter(X, Y, **args)

            if 1:
                if nframe == 0:
                    Nlog[0].extend(nar(profiles[core_id][frame]['N_r1log']))
                    Blog[0].extend(nar(profiles[core_id][frame]['B_los_r1log']))
                if nframe == 1: 
                    Nlog[1].extend(nar(profiles[core_id][frame]['N_r1log']))
                    Blog[1].extend(nar(profiles[core_id][frame]['B_los_r1log']))
                if nframe == 2: 
                    Nlog[2].extend(nar(profiles[core_id][frame]['N_r1log']))
                    Blog[2].extend(nar(profiles[core_id][frame]['B_los_r1log']))
                if nframe == 3: 
                    Nlog[3].extend(nar(profiles[core_id][frame]['N_r1log']))
                    Blog[3].extend(nar(profiles[core_id][frame]['B_los_r1log']))
        '''            
        if core_id == core_list[-1]:
            for i in range(4):
                print('this should be printed 4 times')
                pfit = np.polyfit(Nlog[i],Blog[i],1)
                alpha = pfit[0]
                kappas.append(alpha)
                b_o = pfit[1]
                n_Rho = np.linspace(np.array(Nlog[i]).min(),np.array(Nlog[i]).max(),num=len(Nlog[i]))
                N_Rho = 10 ** n_Rho
                B_two = 10 ** (alpha*n_Rho + b_o)
                plt.plot(N_Rho,B_two,color='k',linestyle='dashed')
        '''

    if 'B_los_r1' in profs:
        row = row_dict['B_los_r1']
        for i, ax in enumerate(axes):
            print('this should be printed 4 times')
            pfit = np.polyfit(Nlog[i],Blog[i],1)
            alpha = pfit[0]
            kappas.append(alpha)
            b_o = pfit[1]
            n_Rho = np.linspace(np.array(Nlog[i]).min(),np.array(Nlog[i]).max(),num=len(Nlog[i]))
            N_Rho = 10 ** n_Rho
            B_two = 10 ** (alpha*n_Rho + b_o)
            ax.plot(N_Rho,B_two,color='k',linestyle='dashed')
            ax.set(xscale='log',yscale='log',ylabel=r'$B_{los}$',xlabel=r'$N$',xlim=ext[-1].minmax, ylim=ext[row].minmax)
    print('kappas',kappas)
    print('saving')
    timescale = ['0_tsing','tsing_tsung','4_panel'][obj.timescale]
    #fig.savefig('b_los_r1_%s_%s_exts_kaps'%(obj.this_looper.sim_name,timescale))
    fig.savefig('testing')


import tsing
reload(tsing)
sim_list=['u501','u502','u503']
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

import tsung_cores
reload(tsung_cores)
if 1:
    sim_list=['u501']
    if 'mp' not in dir():
        for sim in sim_list:
            all_cores=np.unique(TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            #core_list = TL.loops[sim].core_by_mode['Alone']
            #core_list.pop(2)
            #core_list=core_list[:2]
            #core_list=None

            mp=tsung_cores.multipro(TL.loops[sim])
            timescale = 2 #0= 0-tsing, 1=tsing-tsing 2=4 panel
            mp.run(core_list=core_list,tsing=tsing_tool[sim], timescale=timescale,get_particles=False,save_sorts=True )
            replotter(mp,subset=2)
