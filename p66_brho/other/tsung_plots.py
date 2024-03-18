
from starter2 import *
import tsing
reload(tsing)
import tsung_cores
reload(tsung_cores)
           

def replotter(obj, suffix='', redlines=False, subset=0):
    suffix=""

    if subset==1:
        #synthetic los obs
        row_dict={'b_los_cyl':0,'b_los_r1':1,'b_los_r2':2}
        profs=['b_los_cyl','b_los_r1','b_los_r2']
    if subset==2:  #NEED TO CHANGE THESE TO THE MORE COMMON NAMES
        #synthetic obs try one first
        row_dict={'B_los_r1':0, 'Rho':0,'B_los_r1log':0, 'Rholog':0}
        profs=['B_los_r1','Rho'] 
    if subset==3:
        row_dict={'B':0, 'Rho':0, 'Blog':0, 'Rholog':0}
        profs=['B','Rho'] 

    Nplots=len(profs) -1
    Ntimes=len(obj.titles)
    fig,axes = plt.subplots(Nplots,Ntimes, figsize=(16,5))
    fig.subplots_adjust(hspace=0,wspace=0)
    for title,ax in zip(obj.titles,axes): #axes[0] 
        ax.set_title(title)

    ext = [extents() for n in range(len(profs))]
    args = {'linewidth':0.2, 'c':'red', 'alpha':0.4} #'c':[0.5]*3


    # FOR EACH PROFILE
    if subset==2:
        profiles = obj.profiles_gas
        core_list=list(profiles.keys())
    
    # CORES
    Rholog = defaultdict(list)
    Bmaglog = defaultdict(list)
    kappas = []
    for core_id in core_list:
        frames = sorted(list(profiles[core_id].keys()))
        row = row_dict['Rho']
        
        # FRAMES
        for nframe,frame in enumerate(frames):
            ax = axes[nframe]

            X = nar(profiles[core_id][frame]['Rho_sph'])
            Y = nar(profiles[core_id][frame]['Bmag_sph']) 

            ext[-1](X)
            ext[row](Y)
            ax.scatter(X, Y, **args)

            if 1:
                if nframe == 0:
                    Rholog[0].extend(nar(profiles[core_id][frame]['Rho_sph_log']))
                    Bmaglog[0].extend(nar(profiles[core_id][frame]['Bmag_sph_log']))
                if nframe == 1: 
                    Nlog[1].extend(nar(profiles[core_id][frame]['Rholog']))
                    Blog[1].extend(nar(profiles[core_id][frame]['Blog']))
                if nframe == 2: 
                    Nlog[2].extend(nar(profiles[core_id][frame]['Rholog']))
                    Blog[2].extend(nar(profiles[core_id][frame]['Blog']))
                if nframe == 3: 
                    Nlog[3].extend(nar(profiles[core_id][frame]['Rholog']))
                    Blog[3].extend(nar(profiles[core_id][frame]['Blog']))

    if 'B' in profs:
        row = row_dict['B']
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
            #ax.set(xscale='log',yscale='log',ylabel=r'$B_{los}$',xlabel=r'$N$',xlim=ext[-1].minmax, ylim=ext[row].minmax)
            ax.set(xscale='log',yscale='log',ylabel=r'$B$',xlabel=r'$\rho$',xlim=ext[-1].minmax, ylim=ext[row].minmax)
    print('kappas',kappas)
    print('saving')
    timescale = ['0_tsing','tsing_tsung','4_panel'][obj.timescale]
    fig.savefig('b_particles_%s_%s_exts_kaps'%(obj.this_looper.sim_name,timescale))
    #fig.savefig('testing_particles')


plt.close('all')
import three_loopers_u500 as TL
sim_list=['u501']  #02, 03
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

if 1:
    if 'mp' not in dir():
        for sim in sim_list:
            all_cores=np.unique(TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            core_list=core_list[:1]  #debug

            mp=tsung_cores.multipro(TL.loops[sim])
            timescale = 2 
            mp.run(core_list=core_list, tsing=tsing_tool[sim], timescale=timescale, get_particles=True, save_sorts=True )

            #replotter(mp,subset=3)  #2:gas, 3:particles



