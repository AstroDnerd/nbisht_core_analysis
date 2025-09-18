
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
    def __init__(self,this_looper, my_sf):
        self.this_looper=this_looper
        self.my_sf = my_sf
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
        vmin, vmax = 0.1, 100
        nx=ny=64
        self.rbins = np.logspace(np.log10(rmin), np.log10(vmax),nx+1)
        #self.vbins = np.linspace(vmin, vmax, ny+1)
        self.vbins = np.logspace(np.log10(vmin), np.log10(vmax), ny+1)
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

            ms.get_central_velocity(core_id,nt)
            if 0:
                self.lab='rel'
                vr = ms.vr_rel[:,nt]  #the radial one, works ok
            if 0:
                self.lab='mag' 
                vr = ms.rel_vmag[:,nt]  #testing
            if 0:
                self.lab='cen'
                vr = ms.cen_vmag[:,nt]  #testing
            if 1:
                self.lab='RC'
                vr = ms.rc_vmag[:,nt]  #testing
            vrs = vr[asort]
            sigma_vr2 = np.cumsum(vrs**2*dv)/np.cumsum(dv)
            self.vr.append(sigma_vr2)

            v2_sorted=self.vr[-1]


            this_hist, xedge, yedge= np.histogram2d(rsort, v2_sorted, bins=[self.rbins,self.vbins])
            self.hist+= this_hist.astype('float')
            #rquan = np.digitize(rsort, rbins)
            #vquan = np.digitize(sigma_vr, vbins)
            #self.hist[rquan,vquan] += 1


if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
import three_loopers_mountain_top as TLM
import sf2
frame=0
if 0:
    run1 = trial(TLM.loops['u301'])
    run1.by_frame(frame)
    run2 = trial(TLM.loops['u302']) 
    run2.by_frame(frame)
    run3 = trial(TLM.loops['u303'])
    run3.by_frame(frame)
#for frame in [0]: #range(0,110,10):
#    run1.plot(frame)


def plot(self,frame, my_sf2=None,longorno=''):
    pdf = self.hist/(self.dv)
    pdf /= (pdf*self.dv).sum()
    fig,ax=plt.subplots(1,1)
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    vr_max = max([v.max() for v in self.vr])
    vr_min = min([ 0.1, my_sf2[1].min()])
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    for r,v in zip(self.r,self.vr):
        #the index at r=.1
        i_r01 = np.argmin( np.abs( r-0.1))
        sf_r01 = np.argmin( np.abs(my_sf2[0] -0.1))
        scale = my_sf2[1][sf_r01]/v[i_r01] 
        ax.plot(r,v*scale,c=[0.5,0.5,0.5,0.6],lw=0.1)
    for nc,rv in enumerate(zip(self.r,self.vr)):
        r,v=rv
        core_id = self.cores_used[nc]
        #ax.scatter(r[0],v[0],c='k')
        #ax.text(r[0],v[0],"%d"%core_id)
    #ax.plot( self.rbins, [1.0]*self.rbins.size,c=[0.5]*3)
    #ax.scatter(nar([1,2,3,4])/128,[10]*4)
    ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.2)
    axbonk(ax,yscale='log',xscale='log', xlim=[1./128,0.4], ylim=[vr_min,vr_max], ylabel=r'$\sigma_{v,total}^2$',xlabel=r'$r$')
    ax.set_yscale('symlog',linthresh=1)
    from matplotlib.ticker import MultipleLocator
    ml = MultipleLocator(10)
    #ax.tick_params(axis='y',which='minor')
    #ax.yaxis.set_minor_locator(ml)
    ax.set_yticks( np.concatenate([ np.arange(0,1,0.1), np.arange(1,10), np.arange(10,100,10)]))
    fig.colorbar(ploot,ax=ax)
    if my_sf2 is not None:
        ax.plot(my_sf2[0],my_sf2[1],c='k')
    outname = "%s/%svelocity_sf_%s_%s_hist_cXXXX_n%04d.pdf"%(dl.output_directory,longorno,self.lab, self.this_looper.out_prefix, frame)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)

def clean(X,Y):
    ok = (np.isnan(Y)==False)*(Y>0)
    return X[ok], Y[ok]

if 'R' not in dir():
    R = [[],[],[]]
    SF2L = [[],[],[]]
    Things = []
    R[0], SF2L[0] = dpy(dl.sf_path + '/SF2_L_subset_2_u05.h5',['Rmag','SF2_L'])
    R[1], SF2L[1] = dpy(dl.sf_path + '/SF2_L_subset_2_u10.h5',['Rmag','SF2_L'])
    R[2], SF2L[2] = dpy(dl.sf_path + '/SF2_L_subset_2_u11.h5',['Rmag','SF2_L'])
    SCALE = 128**6
    SCALE = 1/128**3
    import radial_binner as rb
    reload(rb)
    for n in range(3):
        #bins = np.linspace(0,128,129)/128
        bins = np.linspace(0,64,65)/64
        bin_edge, bin_center, values =  rb.rb2( R[n]/128, SF2L[n], bins=bins)
        bin_center,values = clean(bin_center,values)
        Things.append([bin_center,values])

SCALE = 128**6
SCALE = 1/128**3


if 1:
    import sf2
    reload( sf2)
    for nrun,this_run in enumerate( [run1, run2, run3]):
        if 'msf' not in dir() or True:
            msf = sf2.make_sf(this_run.this_looper,0)
            rbins,SS = msf.bin_sf2(); SS/=2*np.pi
        #plot(this_run,0, my_sf2=[rbins,SS])
        plot(this_run,0, my_sf2=[Things[n][0],Things[n][1]*SCALE],longorno='long')
