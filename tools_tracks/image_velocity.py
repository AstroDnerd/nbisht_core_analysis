
"""

THIS ALMOST WORKS.
I need to tweak the structure function in tools_tracks/sf2 in some
odd ways to make it work, though.  There's a wierd extra 5.

"""




from starter2 import *
import matplotlib.image as mpimg
reload(trackage)
from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')


class trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]

    def by_frame(self,frame,core_list=None, do_plot=False):
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
        if do_plot:
            fig,ax=plt.subplots(1,1)
        rmin, rmax = 1./2048, 0.4
        vmin, vmax = 0, 8
        nx=ny=32
        self.rbins = np.logspace(np.log10(rmin), np.log10(vmax),nx+1)
        self.vbins = np.linspace(0, vmax, ny+1)
        self.hist = np.zeros([nx,ny])
        def cen(arr):
            return 0.5*(arr[1:]+arr[:-1])
        self.TheX = np.r_[(ny)*[cen(self.rbins)]].transpose()
        self.TheY = np.r_[(nx)*[cen(self.vbins)]]
        self.x_del = (self.rbins[1:]-self.rbins[:-1])
        self.y_del = (self.vbins[1:]-self.vbins[:-1])
        self.x_del.shape = (self.x_del.size,1)
        self.dv = 1./(self.x_del*self.y_del)

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms=ms
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 15:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            nt = np.where(thtr.frames == frame)[0][0]

            ms.get_central_velocity(core_id,nt)
            density = thtr.c([core_id],'density')[:,nt]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
            this_r = ms.r[:,nt]
            this_r[ this_r < rmin] = rmin
            asort = np.argsort(this_r)
            unsort = np.argsort(asort)
            rsort = this_r[asort]
            self.r.append(rsort)
            dv = cell_volume[asort]

            if 0:
                vr = ms.vr_rel[:,nt]  #the radial one, works ok
            if 1:
                vr = ms.rel_vmag[:,nt]  #testing
            if 1:
                vr = ms.cen_vmag[:,nt]
            vrs = vr[asort]
            sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
            self.vr.append(sigma_vr2)
            fig,ax=plt.subplots(1,2)
            ax0=ax[0];ax1=ax[1]
            ax0.scatter(this_r, ms.cen_vmag[:,nt])
            ax0.plot(rsort, np.sqrt(sigma_vr2), c='k')
            fig.savefig("%s/velocity_radius_%s_c%04d_n%04d.pdf"%(dl.output_directory, self.this_looper.out_prefix,core_id, frame))
            plt.close(fig)

            fig, ax=plt.subplots(1,1)
            ax.set_aspect('equal')
            #ax.plot([0,1,1,0,0],[0,0,1,1,0])
            #delta=0.1   
            #ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)

            the_x  = ms.this_x[:,nt]
            the_y  = ms.this_y[:,nt]
            the_z  = ms.this_z[:,nt]
            the_vx = ms.cen_vx[:,nt]
            the_vy = ms.cen_vy[:,nt]
            the_vz = ms.cen_vz[:,nt]
            ax.quiver(the_y,the_z,the_vy,the_vz)
            ax.scatter(ms.mean_y[nt], ms.mean_z[nt],marker='*',c='r',s=4)

            ax.set_title('test 1')
            #pdb.set_trace()
            fig.savefig("%s/velocity_image_%s_c%04d_n%04d.pdf"%(dl.output_directory, self.this_looper.out_prefix,core_id, frame))
            plt.close(fig)
            fig, ax=plt.subplots(1,1)
            ax.hist(the_vx,bins=100, histtype='step',label='vx')
            ax.hist(the_vy,bins=100, histtype='step',label='vy')
            ax.hist(the_vz,bins=100, histtype='step',label='vz')
            ax.hist(vrs,bins=100, histtype='step',label='vr')
            ax.legend(loc=0)
            fig.savefig("%s/velocity_hist_%s_c%04d_n%04d.pdf"%(dl.output_directory, self.this_looper.out_prefix,core_id, frame))                             

            plt.close('all')


if 'do_all_plots' not in dir():
    do_all_plots = False


#import three_loopers as TL
import three_loopers_1tffPlus as TL
import sf2
frame=0
for frame in [0]:#,10,20,30,40,50,60]:
    run1 = trial(TL.looper1)
    run1.by_frame(frame,core_list=[308])
if 0:
    run1 = trial(TL.looper1)
    run1.by_frame(frame)
    run2 = trial(TL.looper2)
    run2.by_frame(frame)
    run3 = trial(TL.looper3)
    run3.by_frame(frame)
