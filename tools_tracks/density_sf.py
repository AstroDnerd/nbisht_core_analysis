
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


class dtrial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.dsf_list=[]

    def by_frame(self,frame,core_list=None, do_plot=False):
        self.r=[]
        self.rho=[]

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
        #dmin, dmax = -1000,1000
        self.dmin, self.dmax = 0, 1000
        nx=ny=64
        self.rbins = np.logspace(np.log10(rmin), np.log10(self.dmax),nx+1)
        self.dbins = np.linspace(self.dmin, self.dmax, ny+1)
        self.hist = np.zeros([nx,ny])
        def cen(arr):
            return 0.5*(arr[1:]+arr[:-1])
        self.TheX = np.r_[(ny)*[cen(self.rbins)]].transpose()
        self.TheY = np.r_[(nx)*[cen(self.dbins)]]
        self.x_del = (self.rbins[1:]-self.rbins[:-1])
        self.y_del = (self.dbins[1:]-self.dbins[:-1])
        self.x_del.shape = (self.x_del.size,1)
        self.dv = 1./(self.x_del*self.y_del)

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 100:
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

            dbar = mean_density = 1.0
            delta = (density - dbar)/dbar

            delta2 = (delta*delta*cell_volume).sum()/cell_volume.sum()
            delta_rms = np.sqrt(delta2)
            delta_prime = (delta-delta_rms)**2
            dsf_sorted = delta_prime[asort]
            dsf = np.cumsum(dsf_sorted**2*dv)/np.cumsum(dv)
            self.dsf_list.append(dsf)

            d2_sorted=self.dsf_list[-1]


            this_hist, xedge, yedge= np.histogram2d(rsort, d2_sorted, bins=[self.rbins,self.dbins])
            self.hist+= this_hist.astype('float')
            #rquan = np.digitize(rsort, rbins)
            #vquan = np.digitize(sigma_vr, vbins)
            #self.hist[rquan,vquan] += 1



def plot(self,frame, my_sf2=None):
    pdf = self.hist/(self.dv)
    pdf /= (pdf*self.dv).sum()
    fig,ax=plt.subplots(1,1)
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    for r,v in zip(self.r,self.dsf_list):
        ax.plot(r,v,c=[0.5,0.5,0.5,0.6],lw=0.1)
    for nc,rv in enumerate(zip(self.r,self.dsf_list)):
        r,v=rv
        core_id = self.cores_used[nc]
        ax.scatter(r[0],v[0],c='k')
        ax.text(r[0],v[0],"%d"%core_id)
    #ax.plot( self.rbins, [1.0]*self.rbins.size,c=[0.5]*3)
    ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.2)
    axbonk(ax,yscale='linear',xscale='log', xlim=[1./2048,0.4], ylim=[self.dmin,self.dmax], ylabel=r'$\sigma_{v,total}^2$',xlabel=r'$r$')
    fig.colorbar(ploot,ax=ax)
    print("============")
    print(my_sf2)
    print("============")

    if my_sf2 is not None:
        ax.plot(my_sf2[0],my_sf2[1],c='k')
    ax.set_title('TAKE 1')
    self.lab=''
    fig.savefig("%s/density_sf_%s_%s_hist_cXXXX_n%04d.pdf"%(dl.output_directory,self.lab, self.this_looper.out_prefix, frame))
    plt.close(fig)

def plot_each_line(self,frame, my_sf2=None):
    pdf = self.hist/(self.dv)
    ok = pdf>0
    pdf[ok] /= (pdf[ok]*self.dv[ok]).sum()
    fig,ax=plt.subplots(1,1)
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.2)
    axbonk(ax,yscale='linear',xscale='log', xlim=[1./2048,0.4], ylim=[0,10], ylabel=r'$\sigma_{v,total}^2$',xlabel=r'$r$')
    self.vr=nar(self.vr)
    vr_max = max([v.max() for v in self.vr])
    if my_sf2 is not None:
        ax.plot(my_sf2[0],my_sf2[1],c='k')
    for nc,rv in enumerate(zip(self.r,self.vr)):
        r,v=rv
        core_id=self.cores_used[nc]
        ax.clear()
        ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.2,shading='nearest')
        axbonk(ax,yscale='linear',xscale='log', xlim=[1./2048,0.4], ylim=[0,vr_max], ylabel=r'$\sigma_{v,total}^2$',xlabel=r'$r$')
        ax.plot(my_sf2[0],my_sf2[1],c='k')
        #ax.plot(r,v,c=[0.5,0.5,0.5,0.6],lw=0.1)
        ax.plot(r,v,c='k')
        ax.scatter(r[0],v[0],c='k')
        outname = "%s/velocity_sf_%s_hist_c%04d_n%04d.png"%(dl.output_directory, self.this_looper.out_prefix, core_id,frame)
        print(outname)
        fig.savefig(outname)
    #ax.plot( self.rbins, [1.0]*self.rbins.size,c=[0.5]*3)
    fig.colorbar(ploot,ax=ax)
    ax.set_title('TAKE 3')
    outname="%s/velocity_sf_%s_hist_cXXXX_n%04d.pdf"%(dl.output_directory, self.this_looper.out_prefix, frame)
    print(outname)
    fig.savefig(outname)
    plt.close(fig)
#import sf2
#reload( sf2)
if 'do_all_plots' not in dir():
    do_all_plots = False


#import three_loopers as TL
import three_loopers_1tff as TL
import sf2
reload(sf2)
frame=0
if 'drun1' not in dir():
    drun1 = dtrial(TL.looper1)
    drun1.by_frame(frame)
    #run2 = trial(TL.looper2)
    #run2.by_frame(frame)
    #run3 = trial(TL.looper3)
    #run3.by_frame(frame)
#for frame in [0]: #range(0,110,10):
#    run1.plot(frame)

if 0:
    for this_run in [drun1]:#, run2, run3]:
        if 'msf' not in dir() or True:
            msf = sf2.make_sf(this_looper=this_run.this_looper,frame=0, field='density')
            rbins,SS = msf.bin_take3(); SS/=2*np.pi
        ##plot_each_line(this_run,0, my_sf2=[rbins,SS])
        #plot(this_run,0, my_sf2=[rbins,SS])
        plot(this_run,0, my_sf2=[rbins,SS])
    #plot(run1,0, my_sf2=[rbins,SS])
