from starter2 import *
from collections import defaultdict
import scipy
import colors

class random_pearson():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.name = self.this_looper.out_prefix
        self.pearson_randos=[]
        self.alpha_randos=[]
        self.fits=[]
        self.radius=[]

    def run(self,core_list=None,do_plots=True, frame=0):
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
        ds = self.this_looper.load(frame)

        frame_ind = frame
        if frame > 0:
            print("Need to fix the frame index")
            return

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue

            self.cores_used.append(core_id)

            this_r = ms.r[:,frame_ind]
            density =thtr.c([core_id],'density')[:,frame_ind]
            cell_volume  =thtr.c([core_id],'cell_volume')[:,frame_ind]

            radius = np.sqrt((this_r**2 * density * cell_volume).sum()/(density*cell_volume).sum())

            center = np.random.random(3)
            print("Do ", self.name, core_id, radius)

            sphere = ds.sphere(center, radius)
            r_to_fit = np.log10( sphere['radius'])
            rho_to_fit = np.log10( sphere['density'])
            #r_to_fit =  sphere['radius']
            #rho_to_fit =  sphere['density']
            pfit = np.polyfit( r_to_fit, rho_to_fit,1)
            self.fits.append( pfit)
            self.alpha_randos.append( pfit[0])
            pearson = scipy.stats.pearsonr( r_to_fit, rho_to_fit)
            self.pearson_randos.append(pearson[0])







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
import three_loopers_six as TL
sim_list=['u401','u402','u403']
sim_list=['u601','u602','u603']
loops = TL.loops
if 'atlist' not in dir():
    atlist={}
    for this_simname in sim_list:
        atlist[this_simname]= alpha_tool( loops[this_simname])

    for this_simname in  sim_list:
        atlist[this_simname].run(do_plots=False )

if 'spheres' not in dir():
    spheres={}
    for this_simname in sim_list:
        spheres[this_simname]= random_pearson( loops[this_simname])
    for this_simname in  sim_list:
        spheres[this_simname].run()

import means_etc
reload(means_etc)
if 1:

    for ns,this_simname in  enumerate(sim_list):
        Pearson_rand=[]
        Pearson_sim=[]
        alpha_rand=[]
        alpha_sim=[]
        if ns == 0:
            lab = sim_list[ns]
        else:
            lab = None
        Pearson_rand += (spheres[this_simname].pearson_randos)
        Pearson_sim += (atlist[this_simname].pearson_t[0])
        alpha_rand += (spheres[this_simname].alpha_randos)
        alpha_sim += (atlist[this_simname].alpha_t[0])

        plt.clf()
        fig, ax, ax_alpha,ax_r =means_etc.three_way_bean()
        ax_r.hist( Pearson_sim, histtype='step', color='r',linestyle='-',density=True, orientation='horizontal')
        ax_r.hist( Pearson_rand, histtype='step', color='k', linestyle='-',density=True, label=lab, orientation='horizontal')
        ax.scatter( alpha_sim, Pearson_sim, color='r', s=0.1, label='core preimages')
        ax.scatter( alpha_rand, Pearson_rand, color='k', s=0.1, label='random spheres')
        ax_alpha.hist( alpha_sim, histtype='step', color='r',linestyle='-',density=True)
        ax_alpha.hist( alpha_rand, histtype='step', color='k', linestyle='-',density=True, label=lab)

        axbonk(ax,xlabel=r'$\alpha_1$',ylabel=r'Pearson $r_P$')
        fig.savefig('plots_to_sort/pearson_randos_%s.pdf'%this_simname)
        plt.close(fig)





if 0:
    #violin plots

    for this_simname in  ['u301']:#,'u302','u303']:
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








