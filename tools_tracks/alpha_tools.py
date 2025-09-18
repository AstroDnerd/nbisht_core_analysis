


from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict

G = 1620/(np.pi*4)
class sub_trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = self.this_looper.sim_name
        self.cores_used=[]
        self.subsonic_length=defaultdict(list)
        self.subvirial_length=defaultdict(list)
        self.max_length=defaultdict(list)
        self.rms_length=defaultdict(list)

        self.mass_fraction_interior=defaultdict(list)
        self.EG_interior=defaultdict(list)
        self.EK_interior=defaultdict(list)
        self.virial_parameter=defaultdict(list)


    def run(self,velocity_method='vct',core_list=None, do_plots=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        rmin, rmax = 1./2048, 0.4
        nx=tsorted.size; ny=64
        self.rbins = np.logspace(np.log10(rmin), np.log10(rmax),ny+1)
        #self.vbins = np.linspace(vmin, vmax, ny+1)
        self.times = tsorted
        self.hist = np.zeros([nx,ny])
        def cen(arr):
            return 0.5*(arr[1:]+arr[:-1])
        #self.TheX = np.r_[(ny)*[cen(self.times)]].transpose()
        #self.TheY = np.r_[(nx)*[cen(self.rbins)]]
        #self.x_del = (self.rbins[1:]-self.rbins[:-1])
        #self.y_del = (self.times[1:]-self.times[:-1])
        #self.x_del.shape = (self.x_del.size,1)
        #self.dv = 1./(self.x_del*self.y_del)

        if do_plots is not None:
            fig3,ax3=plt.subplots(1,1)
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            ms.compute_ge()
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 10:
                continue
            print('SubLength ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            fig,ax=plt.subplots(1,2, figsize=(8,4))
            rmap = rainbow_map(self.times.size)
            if 'virial_radius' in do_plots:
                fig2,ax2=plt.subplots(2,2)
                ax2list=ax2.flatten()

            for nt, time in enumerate(self.times):
                mask = ms.compute_unique_mask(core_id,1/2048,nt)
                ms.get_central_velocity(core_id,nt)

                density = thtr.c([core_id],'density')[:,nt]
                phi = thtr.c([core_id],'PotentialField')[:,nt]
                ge = ms.ge[:,nt]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
                this_r = ms.r[:,nt]
                this_r[ this_r < rmin] = rmin
                asort = np.argsort(this_r)
                unsort = np.argsort(asort)
                rsort = this_r[asort]
                dv = cell_volume[asort]

                #ms.get_central_velocity2(core_id,nt)
                if velocity_method == 'vrm':
                    self.lab='vrm'
                    #radial relative to raw-mean position and raw-mean velocity
                    #vr_rel = ( v - <v>)dot (r - <r>)/||r-<r>|| (\hat{r_rel})
                    #<v> = raw_vx.mean() no weighteing
                    #<r> = this_x.mean()^2 + this_y.mean()^2
                    vr = ms.vr_rel[:,nt]  #the radial one, works ok
                if velocity_method == 'vtm':
                    #total magnitude || v - v0||
                    #rel_vmag = || v - <v>||
                    #<v> as above
                    self.lab='vtm'
                    vr = ms.rel_vmag[:,nt]  #testing
                if velocity_method == 'vtc':
                    #cen_vmag = ||raw_v - vcentral||
                    self.lab='vtc'
                    vr = ms.cen_vmag[:,nt]  #testing
                if velocity_method == 'vrc':
                    #THIS ONE.
                    #the one that make S2.
                    #rc_vmag = ( v - v_central)dot \hat{r_rel}
                    self.lab='vrc'
                    vr = ms.rc_vmag[:,nt]  #testing
                if velocity_method == 'vtcd':
                    #the one that make S2.
                    #r_rel now relative to density-weighted center
                    #rc_vmag = ( v - v_central)dot \hat{r_rel}
                    self.lab='vtcd'
                    vr = ms.rcd_vmag[:,nt]  #testing
                vrs = vr[asort]
                sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)

                v2_sorted=sigma_vr2

                subsonic_length=0
                if (v2_sorted <1).any():
                    subsonic_index = max(np.where( v2_sorted < 1)[0])
                    subsonic_length = rsort[subsonic_index]
                self.subsonic_length[core_id].append( subsonic_length)

                #
                # mass fraction interior
                #
                mass_total = (density[mask]*cell_volume[mask]).sum()
                sorted_mask = mask[asort]

                mass_interior = np.cumsum( density[asort][sorted_mask]*cell_volume[asort][sorted_mask])
                self.mass_fraction_interior.append(mass_interior/mass_total)
                #
                # virial
                #

                EG_interior = np.cumsum(G*ge[asort]*cell_volume[asort])
                EK_interior = np.cumsum(G*density[asort]*vrs**2*cell_volume[asort])
                virial = EK_interior/EG_interior
                self.virial_parameter[core_id].append( virial)
                subvirial_length=0
                if (virial<=1).any():
                    subvirial_index = max(np.where( virial <= 1)[0])
                    subvirial_length = rsort[subvirial_index]
                self.subvirial_length[core_id].append(subvirial_length)

                self.max_length[core_id].append(rsort.max())
                rms_length = np.sqrt(( density*this_r**2*cell_volume).sum()/(density*cell_volume).sum())
                self.rms_length[core_id].append(rms_length)

                if 'virial_radius' in do_plots:
                    ccc=[rmap(nt)]*vrs.size
                    ax2list[0].scatter(vrs, EK_interior, color=ccc)
                    ax2list[1].scatter(vrs, EG_interior, color=ccc)
                    ax2list[2].scatter(vrs, virial, color=ccc)
                    frame = thtr.frames[nt]

                if 'vec' in do_plots:
                    fig2,ax2=plt.subplots(1,1)
                    #the_y=ms.this_y[:,nt]
                    #the_z=ms.this_z[:,nt]
                    the_y=ms.ry_rel[:,nt]
                    the_z=ms.rz_rel[:,nt]
                    the_vy = ms.rel_vy[:,nt]
                    the_vz = ms.rel_vz[:,nt]
                    #ax2.scatter(the_y,the_z, c=[rmap(nt)]*the_y.size)
                    ax2.quiver(the_y,the_z,the_vy,the_vz, color=[rmap(nt)]*the_y.size)
                    frame = thtr.frames[nt]
                    fig2.savefig('plots_to_sort/%s_vmag_1_c%04d_n%04d.png'%(self.name,core_id,frame))
                    plt.close(fig2)



                #this_hist, edges= np.histogram(rsort, weights=v2_sorted, density=False, bins=self.rbins)
                #rcen = 0.5*(self.rbins[1:]+self.rbins[:-1])
                #ok = this_hist > 0
                #ax.plot( rcen[ok], this_hist[ok], c=rmap(nt))
                #self.hist[nt,:] = this_hist.astype('float')
                if 'vr' in do_plots:
                    ax[0].plot( rsort, v2_sorted, c=rmap(nt))

            if 'virial_radius' in do_plots:
                fig2.savefig('plots_to_sort/%s_virial_radius_c%04d_n%04d.png'%(self.name,core_id,frame))
                plt.close(fig2)

            if 'Lt' in do_plots:
                ax3.plot( self.times, self.subsonic_length[core_id],c=[0.5]*3)
                ax[0].plot( self.rbins, self.rbins*0+1, c=[0.5]*3)
                c=[rmap(nt) for nt, time in enumerate(self.times)] 
                ax[1].scatter(self.times, self.subsonic_length[core_id], c=c)
                axbonk(ax[0],xscale='log',yscale='linear',xlabel='r',ylabel=r'$\sigma_v^2(r)$')
                ax[0].set_yscale('symlog',linthresh=1)
                fig.savefig('plots_to_sort/%s_vhist_time_c%04d.png'%(self.name,core_id))
        if 'vr' in do_plots:
            fig3.savefig('plots_to_sort/%s_subsonic_length.png'%(self.name))
