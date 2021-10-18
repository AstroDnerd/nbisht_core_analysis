
"""

THIS ALMOST WORKS.
I need to tweak the structure function in tools_tracks/sf2 in some
odd ways to make it work, though.  There's a wierd extra 5.

"""




from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
reload(trackage)
plt.close('all')


class sub_trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = self.this_looper.sim_name
        self.cores_used=[]
        self.subsonic_length=defaultdict(list)

    def v_hist(self,core_list=None, do_plots=False):
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

        fig3,ax3=plt.subplots(1,1)
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 10:
                continue
            print('SubLength ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            fig,ax=plt.subplots(1,2, figsize=(8,4))
            rmap = rainbow_map(self.times.size)
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

                #ms.get_central_velocity(core_id,nt)
                ms.get_central_velocity2(core_id,nt)
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
                if 0:
                    #the one that make S2.
                    #rc_vmag = ( v - v_central)dot \hat{r_rel}
                    self.lab='RC'
                    vr = ms.rc_vmag[:,nt]  #testing
                if 1:
                    #the one that make S2.
                    #r_rel now relative to density-weighted center
                    #rc_vmag = ( v - v_central)dot \hat{r_rel}
                    self.lab='RC'
                    vr = ms.rcd_vmag[:,nt]  #testing
                vrs = vr[asort]
                sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
                self.vr.append(sigma_vr2)

                v2_sorted=self.vr[-1]

                subsonic_length=0
                if (v2_sorted <1).any():
                    subsonic_index = max(np.where( v2_sorted < 1)[0])
                    subsonic_length = rsort[subsonic_index]
                self.subsonic_length[core_id].append( subsonic_length)

                if do_plots:
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
                if do_plots:
                    ax[0].plot( rsort, v2_sorted, c=rmap(nt))
            if do_plots:
                ax3.plot( self.times, self.subsonic_length[core_id],c=[0.5]*3)
                ax[0].plot( self.rbins, self.rbins*0+1, c=[0.5]*3)
                c=[rmap(nt) for nt, time in enumerate(self.times)] 
                ax[1].scatter(self.times, self.subsonic_length[core_id], c=c)
                axbonk(ax[0],xscale='log',yscale='linear',xlabel='r',ylabel=r'$\sigma_v^2(r)$')
                ax[0].set_yscale('symlog',linthresh=1)
                fig.savefig('plots_to_sort/%s_vhist_time_c%04d.png'%(self.name,core_id))
        if do_plots:
            fig3.savefig('plots_to_sort/%s_subsonic_length.png'%(self.name))



if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4
import sf2
frame=0
sim_list = ['u401']
if 'Ltool' not in dir():
    Ltool={}
    for this_simname in sim_list:
        Ltool[this_simname] = sub_trial( TL4.loops[this_simname])
        Ltool[this_simname].v_hist(do_plots=True)
        #Ltool[this_simname].v_hist(frame, core_list = [13,14,15,16,17,18,19,21])
#   run2 = sub_trial(TLM.loops['u302'])
#   run2.v_hist(frame)#, core_list = [10,32,323])
#   run3 = sub_trial(TLM.loops['u303'])
#   run3.v_hist(frame)#, core_list = [10,32,323])
import heat_map
reload(heat_map)

if 1:
    for this_simname in sim_list:
        fig,ax=plt.subplots(1,1)
        if 0:
            thresh = 1e-2
            top = 0.3
            linbins=16
            logbins=16
            bins=np.unique(np.concatenate([ np.linspace(0,thresh,linbins), np.geomspace(thresh, top, logbins)]))
        bins = np.geomspace(1e-4,0.3,32)
        qmatrix=heat_map.plot_heat( tool=Ltool[this_simname], quan_dict = Ltool[this_simname].subsonic_length, ax=ax,bins=bins)
        outname = "plots_to_sort/%s_subsonic_heatmap.png"%this_simname
        fig.savefig(outname)
        plt.close(fig)


if 0:
    fig,ax=plt.subplots(1,1)
    nsub = []
    for rrr in [run1]:
        for core_id in rrr.cores_used:
            ax.plot( rrr.times, rrr.subsonic_length[core_id],  c=[0.5]*3, linewidth=0.1)
            nsub.append( (nar(rrr.subsonic_length[core_id])>0).sum())
        axbonk(ax, xlabel='t',ylabel=r'$L_{sub}$')
        ax.set_yscale('symlog',linthresh=1./2048)
        fig.savefig('plots_to_sort/%s_subsonic_length.png'%(rrr.name))
        ax.clear()
        ax.hist(nsub)
        fig.savefig('plots_to_sort/%s_n_subsonic.png'%(rrr.name))
