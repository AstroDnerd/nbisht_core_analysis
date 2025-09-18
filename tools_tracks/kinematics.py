

"""

Basic kinematics 

"""




from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
reload(trackage)
plt.close('all')


class kinematics():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = self.this_looper.sim_name
        self.cores_used=[]
        self.mass_flux=defaultdict(list)
        self.mean_density=defaultdict(list)
        self.mean_radial=defaultdict(list)

    def go(self,core_list=None, do_plots=False):
        self.r=[]
        self.vr=[]
        self.v2=[]

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        rmin, rmax = 1./2048, 0.4
        self.times = tsorted[:-1]

        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 3:
                continue
            print('kinematics ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            for nt, time in enumerate(self.times):
                density = thtr.c([core_id],'density')[:,nt]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
                this_r = ms.r[:,nt]
                this_r[ this_r < rmin] = rmin
                asort = np.argsort(this_r)
                unsort = np.argsort(asort)
                rsort = this_r[asort]
                self.r.append(rsort)
                dv = cell_volume[asort]

                ms.get_central_velocity(core_id,nt)
                #ms.get_central_velocity2(core_id,nt)
                if 0:
                    self.lab='rel'
                    #radial relative to raw-mean position and raw-mean velocity
                    #vr_rel = ( v - <v>)dot (r - <r>)/||r-<r>|| (\hat{r_rel})
                    #<v> = raw_vx.mean() no weighteing
                    #<r> = this_x.mean()^2 + this_y.mean()^2
                    vr = ms.vr_rel[:,nt]  #the radial one, works ok
                if 0:
                    #total magnitude || v - v0||
                    #rel_vmag = || v - <v>||
                    #<v> as above
                    self.lab='mag' 
                    vr = ms.rel_vmag[:,nt]  #testing
                if 0:
                    #cen_vmag = ||raw_v - vcentral||
                    self.lab='cen'
                    vr = ms.cen_vmag[:,nt]  #testing
                if 1:
                    #the one that make S2.
                    #rc_vmag = ( v - v_central)dot \hat{r_rel}
                    self.lab='RC'
                    vr = ms.rc_vmag[:,nt]  #testing
                if 0:
                    # mass-weighted centroid
                    self.lab='RCd'
                    vr = ms.rcd_vmag[:,nt]  #testing

                #vrs = vr[asort]
                #sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
                #self.vr.append(sigma_vr2)

                this_density = (density*dv).sum()/dv.sum()
                this_flux = (vr*density*dv).sum()/dv.sum()
                this_vr   = (vr*dv).sum()/dv.sum()
                if np.isnan(this_flux).any():
                    pdb.set_trace()
                self.mean_density[core_id].append( this_density)
                self.mass_flux[core_id].append(this_flux)
                self.mean_radial[core_id].append( this_vr)




if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4
import sf2
frame=0
sim_list = ['u401']#,'u402','u403']
if 'Dtool' not in dir():
    Dtool={}
    for this_simname in sim_list:
        Dtool[this_simname] = kinematics( TL4.loops[this_simname])
        Dtool[this_simname].go()
        #Ltool[this_simname].v_hist(frame, core_list = [13,14,15,16,17,18,19,21])
#   run2 = sub_trial(TLM.loops['u302'])
#   run2.v_hist(frame)#, core_list = [10,32,323])
#   run3 = sub_trial(TLM.loops['u303'])
#   run3.v_hist(frame)#, core_list = [10,32,323])
import heat_map
reload(heat_map)

if 1:
    for this_simname in sim_list:
        fig,ax=plt.subplots(1,3, figsize=(12,8))
        if 0:
            thresh = 1e-2
            top = 0.3
            linbins=16
            logbins=16
            bins=np.unique(np.concatenate([ np.linspace(0,thresh,linbins), np.geomspace(thresh, top, logbins)]))
        #bins = np.geomspace(1e-4,0.3,32)
        bins=np.geomspace( 1,1e8,64)
        qmatrix=heat_map.plot_heat( tool=Dtool[this_simname], quan_dict = Dtool[this_simname].mean_density, ax=ax[0],bins=bins)
        axbonk(ax[0],xlabel=r'$t/t_{\rm{ff}}$', yscale='log', ylabel='<rho>')
        ##bins=np.linspace(-100,100,64)
        #thresh=1
        #end=1e8
        #N1 = 64
        #N2 = 16
        #bins=np.unique(np.concatenate([-np.geomspace(thresh,end,N1)[::-1],
        #                               np.linspace(-thresh,thresh,N2),
        #                               np.geomspace(thresh,end,N1)]))
        bins = -np.geomspace(0.4,1e7,64)[::-1]
        #bins = np.geomspace(0.4,1e7,64)
        #def fixer(matt):
        #    return np.abs(matt)
        fixer=None
        fmatrix=heat_map.plot_heat( tool=Dtool[this_simname], quan_dict = Dtool[this_simname].mass_flux, ax=ax[1],bins=bins, fixer=fixer)
        ax[1].set_yscale('symlog',linthresh=0.4)
        ax[1].set_ylabel(' rho v_r')
        vmatrix=heat_map.plot_heat( tool=Dtool[this_simname], quan_dict = Dtool[this_simname].mean_radial, ax=ax[2])
        outname = "plots_to_sort/%s_mean_density.png"%this_simname
        fig.savefig(outname)
        plt.close(fig)

