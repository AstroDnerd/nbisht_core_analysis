
from starter2 import *
from collections import defaultdict
import scipy
import colors
import movie_frames 

from scipy.ndimage import gaussian_filter

tsing_tool={}
def get_tsing(loop_dict):
    for ns,sim in enumerate(loop_dict):
        if sim in tsing_tool:
            continue
        else:
            print("tsing, %s"%(sim))

        loop=loop_dict[sim]
        obj=te_tc(loop)
        sim=loop.sim_name
        tsing_tool[sim]=obj
        tsing_tool[sim].run()
    return tsing_tool


class te_tc:
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.tsing = defaultdict(list)
        self.tend = defaultdict(list)
        self.tsing_core={}
        self.cores_used=[]
        self.tend_core={}
        self.fsung=[]
        self.mode=[]
        self.tsing_frame={}
        self.tend_frame={}


        
    def run(self,core_list=None, suffix='',make_plots=False):
        this_looper=self.this_looper
        print('Do',this_looper.sim_name)

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        times=thtr.times[mask]+0 #the zero makes a copy
        frames=thtr.frames[mask]+0
        times.shape=times.size,1
        times=times/colors.tff
        self.core_list=core_list
        for core_id in core_list:

                
            #ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
            #ms.particle_pos(core_id)

            sl=slice(None)
            c=[0.5]*4
            c='k'

            #rho = ms.density[sl].transpose()
            rho = thtr.c([core_id], 'density').transpose()[mask,:]+0
            rho_sorted=rho+0
            rho_sorted.sort(axis=1)
            rho_max = rho_sorted[:,-5]
            #rho = rho[mask,:]
            #rho_max=rho.max(axis=1)

            UB = gaussian_filter(rho_max,1)
            tf = times.flatten()
            dt = tf[1:]-tf[:-1]
            dU = UB[1:]-UB[:-1]
            tc = 0.5*(times[1:]+times[:-1])
            dU = dU[1:-1]
            dt = dt[1:-1]
            tc = tc[1:-1]
            dudt=dU/dt
            thresh = 1e5
            singularity = np.where( dudt >thresh)[0][0]
            if (dudt[singularity:]<=0).any():
                collapse_done = np.where( dudt[singularity:] < 0)[0][0]
                collapse_done += singularity
            else:
                collapse_done=-1
            self.tsing_frame[core_id]=frames[singularity]
            self.tend_frame[core_id]=frames[collapse_done]
            this_tsing = times[singularity][0]  #the extra [0] is for the plottable shape of times.
            this_tend = times[collapse_done][0]
            self.tsing_core[core_id]=this_tsing
            self.tend_core[core_id]=this_tend
            '''
            if core_id in this_looper.core_by_mode['Alone']:
                self.mode.append(1)
                self.tsing['Alone'].append( this_tsing)
                self.tend['Alone'].append( this_tend)
            elif core_id in this_looper.core_by_mode['Binary']: 
                self.mode.append(2)
                self.tsing['Binary'].append( this_tsing)
                self.tend['Binary'].append( this_tend)

            elif core_id in this_looper.core_by_mode['Cluster']: 
                self.mode.append(3)
                self.tsing['Cluster'].append( this_tsing)
                self.tend['Cluster'].append( this_tend)
            else:
                print("You broke something.")
                raise
            '''
            self.fsung.append( (rho[collapse_done,:] > 5e3).sum()/rho[collapse_done,:].size)


            if make_plots:
                print('save',core_id)

                plt.close('all')
                fig,ax=plt.subplots(2,1)
                ax0=ax[0];ax1=ax[1]
                ax[0].plot( times, rho, color=[0.5]*4, linewidth=0.1)
                ax[0].set(yscale='log',xlabel='t/tff',ylabel='rho')
                ax[0].axvline(self.tsing_core[core_id])
                ax[0].axvline(self.tend_core[core_id])
                ax[0].plot(times,rho_max,c='r')
                #ax[1].plot(
                fig.savefig('plots_to_sort/tsing_%s_c%04d'%(this_looper.sim_name,core_id))
