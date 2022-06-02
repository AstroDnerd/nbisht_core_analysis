
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
        self.rmin, self.rmax = 1./2048, 0.4
        self.vmin, self.vmax = 0.1, 100
        nx=67; ny=64
        if 0:
            self.xscale='linear'
            self.yscale='linear'
            self.rbins = np.linspace(self.rmin, self.rmax,nx+1)
            self.vbins = np.linspace(0, self.vmax, ny+1)
        if 1:
            self.xscale='log'
            self.yscale='log'
            self.rbins = np.logspace(np.log10(self.rmin), np.log10(self.rmax),nx+1)
            self.vbins = np.logspace(np.log10(self.vmin), np.log10(self.vmax), ny+1)
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
            #self.rmin = 1./2048
            nt = np.where(thtr.frames == frame)[0][0]

            density = thtr.c([core_id],'density')[:,nt]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
            this_r = ms.r[:,nt]
            this_r[ this_r < self.rmin] = self.rmin
            asort = np.argsort(this_r)
            unsort = np.argsort(asort)
            rsort = this_r[asort]
            if 1:
                rsort = rsort[1:]
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
            if 1:
                sigma_vr2=sigma_vr2[1:]
            self.vr.append(sigma_vr2)

            v2_sorted=self.vr[-1]


            this_hist, xedge, yedge= np.histogram2d(rsort, v2_sorted, bins=[self.rbins,self.vbins])
            self.hist+= this_hist.astype('float')
            #rquan = np.digitize(rsort, rbins)
            #vquan = np.digitize(sigma_vr, vbins)
            #self.hist[rquan,vquan] += 1


def plot(self,frame, my_sf2=None,longorno='', external_ax=None, external_norm=None):
    pdf = self.hist/(self.dv)
    pdf /= (pdf*self.dv).sum()
    if external_ax is None:
        fig,ax=plt.subplots(1,1)
    else:
        ax=external_ax
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = pdf[pdf>0].min()
    vr_max = max([v.max() for v in self.vr])
    vr_min = min([ 0.1, my_sf2[1].min()])
    print("MIN MIN", minmin, pdf.max())
    print("MIN MIN", self.hist[ self.hist>0].min(), self.hist.max())
    if external_norm is None:
        norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
    else:
        norm = external_norm

    #plot individuals
    for r,v in zip(self.r,self.vr):
        #the index at r=.1
        #does not work.  Makes a bowtie.
        #i_r01 = np.argmin( np.abs( r-0.1))
        #sf_r01 = np.argmin( np.abs(my_sf2[0] -0.1))
        #scale = my_sf2[1][sf_r01]/v[i_r01] 
        #ax.plot(r,v*scale,c=[0.5,0.5,0.5,0.6],lw=0.1)
        ax.plot(r,v,c=[0.5,0.5,0.5,0.1],lw=0.3)


    mean_of_sigma = (self.TheY*pdf).sum(axis=1)/pdf.sum(axis=1)
    mean_of_sigma = (self.TheY*pdf*self.x_del).sum(axis=1)/(pdf*self.x_del).sum(axis=1)
    ax.plot(self.TheX[:,0], mean_of_sigma, 'k--')
    ploot=ax.pcolormesh(self.TheX, self.TheY, self.hist,cmap=cmap,norm=norm,alpha=0.2, shading='nearest')

    xlim=self.rmin,self.rmax
    xlim=1./128,self.rmax
    ylim=self.vmin,self.vmax
    axbonk(ax,yscale=self.yscale,xscale=self.xscale, xlim=xlim, ylim=ylim, ylabel=r'$\sigma_{v,L}^2, S_{2,L}$',xlabel=r'$r$')
    if external_ax is None:
        fig.colorbar(ploot,ax=ax)
    if my_sf2 is not None:
        print(my_sf2)
        the_x = my_sf2[0][1:]
        the_y = my_sf2[1][1:]
        ax.plot(the_x,the_y,c='k')
        pfit = np.polyfit(np.log10(the_x),np.log10(the_y), 1)
        print("================",pfit)
        ax.plot( the_x, 10**(np.log10(the_x)*pfit[0]+pfit[1]),c='r')

    outname = "%s/%svelocity_sf_%s_%s_hist_cXXXX_n%04d.pdf"%(dl.output_directory,longorno,self.lab, self.this_looper.out_prefix, frame)
    if external_ax is None:
        fig.savefig(outname)
        print(outname)
        plt.close(fig)
    return ploot


if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
#import three_loopers_mountain_top as TLM
import three_loopers_six as TL6
MOD = TL6
import sf2
frame=0
if 'run1' not in dir():
    run1 = trial(MOD.loops['u601'])
    run1.by_frame(frame)
    run2 = trial(MOD.loops['u602']) 
    run2.by_frame(frame)
    run3 = trial(MOD.loops['u603'])
    run3.by_frame(frame)
#for frame in [0]: #range(0,110,10):
#    run1.plot(frame)

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
    fig,ax=plt.subplots(1,1)
    for n in range(3):
        the_x = Things[n][0][1:]
        the_y = Things[n][1][1:]*SCALE
        ax.plot(the_x,the_y,c='k')
        pfit = np.polyfit(np.log10(the_x),np.log10(the_y), 1)
        print("================",pfit)
        ax.plot( the_x, 10**(np.log10(the_x)*pfit[0]+pfit[1]),c='r')
        ax.set_yscale('log')
        ax.set_xscale('log')
    fig.savefig('plots_to_sort/messup.png')
if 1:
    import sf2
    reload( sf2)
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    delta=0.07
    fig.subplots_adjust(hspace=0, wspace=0, right=1-delta/2, left=delta)

    norm = mpl.colors.LogNorm(vmin=2.4e-11,vmax=0.004)
    norm = mpl.colors.LogNorm(vmin=1,vmax=4000)
    for nrun,this_run in enumerate( [run1, run2, run3]):
        #if 'msf' not in dir() or True:
        #    msf = sf2.make_sf(this_run.this_looper,0)
        #    rbins,SS = msf.bin_sf2(); SS/=2*np.pi
        #plot(this_run,0, my_sf2=[rbins,SS])
        plot_obj=plot(this_run,0, my_sf2=[Things[nrun][0],Things[nrun][1]*SCALE],longorno='long',
             external_ax=ax[nrun], 
             external_norm = norm)
        if nrun > 0:
            ax[nrun].set_ylabel('')
            ax[nrun].set_yticks([])
    fig.colorbar(plot_obj)
    fig.savefig('plots_to_sort/SF2.pdf')
