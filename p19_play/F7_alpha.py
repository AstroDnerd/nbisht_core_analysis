from starter2 import *
from collections import defaultdict
import scipy
import colors
reload(tsing)

class alpha_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.alpha_c=defaultdict(list)
        self.pearson_c=defaultdict(list)
        self.alpha_t=defaultdict(list)
        self.pearson_t=defaultdict(list)
        self.cores_used=[]
        self.name = self.this_looper.out_prefix

    def run(self,core_list=None,do_plots=True):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        tmap = rainbow_map( len(thtr.frames))
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue

            self.cores_used.append(core_id)

            print('go ', self.name,core_id)
            fig,ax = plt.subplots(1,3, figsize=(12,4))
            ax0=ax[0]; ax1=ax[1]; ax2=ax[2]
            for nf,frame in enumerate(thtr.frames):
                time = thtr.times[nf]
                density = thtr.c([core_id],'density')[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                this_r = ms.r[:,nf]
                mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                #HAVE to do fits and pearson on the log.
                r_to_fit = np.log10( this_r[mask])
                rho_to_fit = np.log10( density[mask])
                #r_to_fit =  this_r[mask]
                #rho_to_fit = density[mask]

                if mask.sum() > 3:
                    pfit = np.polyfit( r_to_fit, rho_to_fit, 1)
                    pearson = scipy.stats.pearsonr( r_to_fit, rho_to_fit)
                else:
                    pearson = [0.0,0]
                    pfit = [0.0,0]
                        
                self.alpha_c[core_id].append(pfit)
                self.pearson_c[core_id].append(pearson)

                self.alpha_t[nf].append(pfit[0])
                self.pearson_t[nf].append(pearson[0])



                if do_plots:
                    ax0.scatter( 10**r_to_fit, 10**rho_to_fit, c=[tmap(nf)]*mask.sum())
                    ax1.scatter( time, pfit[0], c=[tmap(nf)], marker="*")
                    ax2.scatter( time, pearson[0], c=[tmap(nf)], marker="*")
            if do_plots:
                axbonk(ax0,xlabel=r'$r$', ylabel=r'$\rho$',yscale='log',xscale='log')
                for LLL in [-0.5,-1,-1.5,-2]:
                    ax1.plot( thtr.times, [LLL]*len(thtr.times), c=[0.5]*4)

                fig.savefig('plots_to_sort/%s_alpha_time_c%04d'%(self.name, core_id))
                plt.close(fig)

#import three_loopers_mountain_top as TLM
#sim_list=['u301','u302','u303']
#import three_loopers_tenfour as TL4
import three_loopers_u500 as TL
sim_list=['u401','u402','u403']
sim_list=['u601','u602','u603']
sim_list=['u501','u502','u503']
sim_list=['u501']
loops = TL.loops
if 'atlist' not in dir():
    atlist={}
    for this_simname in sim_list:
        atlist[this_simname]= alpha_tool( loops[this_simname])

    for this_simname in  sim_list:
        atlist[this_simname].run(do_plots=False )

if 'tsing_tool' not in dir() or True:
    tsing_tool={}
    for ns,sim in enumerate(sims):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

if 1:
    fig,axes=plt.subplots(3,3, figsize=(12,12))
    for ns,sim in enumerate(sim_list):
        if ns>0:
            continue
        at = atlist[sim]
        times = at.this_looper.tr.times/colors.tff
        frames = at.this_looper.tr.frames
        alpha_square = np.zeros([len(at.alpha_t[0]),len(times)])
        tstool=tsing_tool[sim]
        for ntime, time in enumerate(times):
            alpha_square[:,ntime]=at.alpha_t[ntime]

        for nc, core_id in enumerate(at.core_list):
            modes = at.this_looper.mode_dict[core_id]
            for mode in ['One','Binary','Cluster']:
                if mode in modes:
                    mode_mask[mode][nc]=1

        ax=axes[ns]
        for nc, core_id in enumerate(at.core_list):
            if at.core_list[nc] != tstool.core_list[nc]:
                raise
            if core_id  in at.this_looper.core_by_mode['One']:
                nrow=0
            elif core_id in at.this_looper.core_by_mode['Binary']:
                nrow=1
            elif core_id in at.this_looper.core_by_mode['Cluster']:
                nrow=2
            else:
                raise
            ax = axes[nrow][ns]

            ax.plot( times/tstool.tsing_core[core_id], alpha_square[nc,:], c=[0.1]*4)
            #ax.set(xlim=[0,3])

    fig.savefig('plots_to_sort/alpha_stretch.pdf')

if 0:
    fig,axes=plt.subplots(3,3)
    for ns,sim in enumerate(sim_list):
        if ns>0:
            continue
        at = atlist[sim]
        times = at.this_looper.tr.times
        frames = at.this_looper.tr.frames
        alpha_square = np.zeros([len(at.alpha_t[0]),len(times)])
        for ntime, time in enumerate(times):
            alpha_square[:,ntime]=at.alpha_t[ntime]

        mode_mask={}
        for mode in ['One','Binary','Cluster']:
            #KLUDGE
            mode_mask[mode]=np.zeros_like( at.core_list)
        for nc, core_id in enumerate(at.core_list):
            modes = at.this_looper.mode_dict[core_id]
            for mode in ['One','Binary','Cluster']:
                if mode in modes:
                    mode_mask[mode][nc]=1
        for nm, mode in enumerate(['One','Binary','Cluster']):
            print('box')
            axes[ns][nm].boxplot( alpha_square[mode_mask[mode]==1, :])
            #axes[ns][nm].boxplot( alpha_square )
    fig.savefig('plots_to_sort/alpha_mode.pdf')



if 1:
    #violin plots

    for this_simname in  ['u501']:#,'u302','u303']:

        fig,ax=plt.subplots(1,2)
        ax0=ax[0]; ax1=ax[1]
        at = atlist[this_simname]
        times = at.this_looper.tr.times
        frames = at.this_looper.tr.frames
        alpha_square = np.zeros([len(at.alpha_t[0]),len(times)])
        for ntime, time in enumerate(times):
            alpha_square[:,ntime]=at.alpha_t[ntime]
        ax0.boxplot(alpha_square)

        pearson_square = np.zeros([len(at.pearson_t[0]),len(times)])
        for ntime, time in enumerate(times):
            pearson_square[:,ntime]=at.pearson_t[ntime]
        ax1.boxplot(pearson_square)

        #this is an annoying amount of code to set the x-tick-labels.
        G = 1620/(4*np.pi)
        tff  = np.sqrt(3*np.pi/(32*G*1))
        these_times = times[::3]/tff
        ticks = ax0.xaxis.get_ticklocs()
        new_ticks = ticks[::3]
        ax0.xaxis.set_ticks(new_ticks)
        ax0.xaxis.set_ticklabels(["%0.2f"%s for s in these_times])

        ax1.xaxis.set_ticks(new_ticks)
        ax1.xaxis.set_ticklabels(["%0.2f"%s for s in these_times])
        fig.savefig('plots_to_sort/%s_violins.pdf'%at.name)








