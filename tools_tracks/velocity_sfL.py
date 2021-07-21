
"""

THIS ALMOST WORKS.
I need to tweak the structure function in tools_tracks/sf2 in some
odd ways to make it work, though.  There's a wierd extra 5.

"""




from starter2 import *
import matplotlib.image as mpimg

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
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 15:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            nt = np.where(thtr.frames == frame)[0][0]

            density = thtr.c([core_id],'density')[:,nt]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
            this_r = ms.r[:,nt]
            this_r[ this_r < rmin] = rmin
            asort = np.argsort(this_r)
            unsort = np.argsort(asort)
            rsort = this_r[asort]
            self.r.append(rsort)
            dv = cell_volume[asort]

            if 1:
                vr = ms.vr_rel[:,nt]  #the radial one, works ok
                vrs = vr[asort]
                sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
            if 0:
                vr = ms.rel_vmag[:,nt]  #testing
                vrs = vr[asort]
                sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
            self.vr.append(sigma_vr2)

            v2_sorted=self.vr[-1]



            if do_plot:
                ax.plot(rsort,sigma_vr2, c=tmap(nt))

            this_hist, xedge, yedge= np.histogram2d(rsort, v2_sorted, bins=[self.rbins,self.vbins])
            self.hist+= this_hist.astype('float')
            #rquan = np.digitize(rsort, rbins)
            #vquan = np.digitize(sigma_vr, vbins)
            #self.hist[rquan,vquan] += 1
        

        if do_plot:
            axbonk(ax,xlabel='r',ylabel=r'$\sigma_v(r)$',xscale='log',xlim=[rmin,0.4], ylim=[0,10])
            #ax.set_yscale('symlog',linthresh=10)

            outname = "%s/%s_sigma_vr%04d.png"%(dl.output_directory, self.this_looper.out_prefix, frame)
            print(outname)
            fig.savefig(outname)
            fig,ax=plt.subplots(1,1)

            ax.pcolormesh(self.rbins, self.vbins, self.hist)
            fig.savefig("%s/%s_hist_cXXXX_n%04d.png"%(dl.output_directory, self.this_looper.out_prefix, frame))

            plt.close('all')


if 'do_all_plots' not in dir():
    do_all_plots = False


#import three_loopers as TL
import three_loopers_1tffPlus as TL
import sf2
frame=0
if 0:
    run1 = trial(TL.looper1)
    run1.by_frame(frame)
    run2 = trial(TL.looper2)
    run2.by_frame(frame)
    run3 = trial(TL.looper3)
    run3.by_frame(frame)
#for frame in [0]: #range(0,110,10):
#    run1.plot(frame)


def plot(self,frame, my_sf2=None):
    pdf = self.hist/(self.dv)
    pdf /= (pdf*self.dv).sum()
    fig,ax=plt.subplots(1,1)
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    for r,v in zip(self.r,self.vr):
        ax.plot(r,v,c=[0.5,0.5,0.5,0.3],lw=0.1)
    for r,v in zip(self.r,self.vr):
        ax.scatter(r[0],v[0],c='k')
    #ax.plot( self.rbins, [1.0]*self.rbins.size,c=[0.5]*3)
    ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.5)
    axbonk(ax,yscale='linear',xscale='log', xlim=[1./2048,0.4], ylim=[0,9], ylabel=r'$\sigma_{v,total}^2$',xlabel=r'$r$')
    fig.colorbar(ploot,ax=ax)
    if my_sf2 is not None:
        ax.plot(my_sf2[0],my_sf2[1],c='k')
    ax.set_title('TAKE 3')
    fig.savefig("%s/long_velocity_sf_%s_hist_cXXXX_n%04d.pdf"%(dl.output_directory, self.this_looper.out_prefix, frame))
    plt.close(fig)

import sf2
reload( sf2)
for this_run in [run1, run2, run3]:
    #if 'msf' not in dir() or True:
        #msf = sf2.make_sf(this_run.this_looper,0)
        #rbins,SS = msf.bin_take3(); SS/=2*np.pi
    rbins, SS = dpy('RadialSF.h5',['R','SF2L'])
    SS /= 128**3
    SS /= 4
        #rbins,SS = msf.bin_kludged()


    plot(this_run,0, my_sf2=[rbins,SS])
#plot(run1,0, my_sf2=[rbins,SS])
