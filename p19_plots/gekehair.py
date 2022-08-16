
from starter2 import *
from collections import defaultdict
import scipy
import colors
import xtra_energy
import camera_path
reload(camera_path)
import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

from scipy.ndimage import gaussian_filter

class gekeratio:
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.avg_ge=[]
        self.avg_ke=[]
        self.rho_max=[]
        self.rho_avg=[]
        self.times=None
        self.tsing=[]
    def run(self,core_list=None, do_plots=True):
        this_looper=self.this_looper

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        all_times=thtr.times
        all_frames=thtr.frames
        times=thtr.times[mask]+0 #the zero makes a copy
        times.shape=times.size,1
        times=times/colors.tff
        frames=all_frames[mask]
        rho_all = thtr.track_dict['density']
        rho_min=rho_all.min()
        rho_max=rho_all.max()
        self.times=times.flatten()

        mini_scrubbers={}
        for nc,core_id in enumerate(core_list):
            print('V %s %d'%(this_looper.sim_name,core_id))
                
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
            ms.particle_pos(core_id)
            ms.compute_ge(core_id)
            ms.compute_ke(core_id)
            ms.compute_ke_rel(core_id)
            mini_scrubbers[core_id]=ms

            if ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4

            rho = ms.density[:,mask].transpose()
            boo= np.abs(ms.ge.mean(axis=0)[mask])
            self.avg_ge.append(boo)
            moo = ms.ke_rel.mean(axis=0)[mask]
            self.avg_ke.append(moo)

            self.rho_max.append( rho.max(axis=1))
            self.rho_avg.append( rho.mean(axis=1))


            #
            # Pick a bunch of times to analyze.
            # Especially the start and end of singularity
            #
            UB = gaussian_filter(rho.max(axis=1),1)
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
            self.tsing.append( times[singularity])

            if (dudt[singularity:]<=0).any():
                collapse_done = np.where( dudt[singularity:] < 0)[0][0]
                collapse_done += singularity
            else:
                collapse_done=-1


sims=['u501', 'u502','u503']
import three_loopers_u500 as TL
import mass_tools
if 0:
    if 'mt' not in dir():
        mt={}
    for sim in sims:
        if sim not in mt:
            mt[sim]=mass_tools.mass_tool(TL.loops[sim])
            mt[sim].run()

sims=['u501', 'u502','u503']
#sims=['u502']#, 'u501']
for sim in sims:
    #core_list=[381]
    #core_list={'u501':[323], 'u502':[381]}[sim]
    #core_list={'u501':[323], 'u502':[112]}[sim]
    #core_list=[31,32]
    #core_list=[114]

    core_list=None
    annotate_phases=False
    core_list = [112]
    #annotate_phases=True
    for mode in ['One','Binary','Cluster']:
        core_list = TL.loops[sim].core_by_mode[mode]
        obj=gekeratio(TL.loops[sim])
        obj.run(do_plots=False, core_list=core_list)#, mass=mt[sim].unique_mass, dof=mt[sim].dof, volume=mt[sim].volume)

        fig,ax=plt.subplots(1,1)
        for ng,ge in enumerate(obj.avg_ge):
            ax.plot( obj.times/obj.tsing[ng], ge/obj.avg_ke[ng], c='k', linewidth=0.2)
            ax.axhline(1,c=[0.5]*4, linewidth=0.2)
            ax.axhline(0.5,c=[0.5]*4, linewidth=0.2)
            #ax.plot( obj.times/obj.tsing[ng], ge, c='k', linewidth=0.2)
            #ax.plot( obj.times/obj.tsing[ng], obj.rho_max[ng], c='r', linewidth=0.2)
            #ax.plot( obj.times/obj.tsing[ng], obj.rho_avg[ng], c='g', linewidth=0.2)
        ax.set(yscale='log',xscale='linear',xlabel='t/tsing', ylim=[1e-2,100], ylabel=r'$E_G/E_K$')
        fig.savefig('plots_to_sort/gekehair_%s_%s.pdf'%(sim,mode))

