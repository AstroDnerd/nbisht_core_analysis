
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

form='pdf'
if 'looper1' not in dir():
    from three_loopers import *


tm = rainbow_map(15)
#all_cores = all_cores[:10]
def make_prof(ds,fields,weight_field=None,accumulation=False,fractional=True,n_bins=64,extrema=None):
    reg = ds.all_data()
    prof = yt.create_profile(reg,fields[0],fields[1] ,weight_field=weight_field,accumulation=accumulation,
                            fractional=fractional, n_bins=n_bins, extrema=extrema)
    the_x = 0.5*(prof.x_bins[1:]+prof.x_bins[0:-1])
    the_y = prof[fields[1]]
    #if units[0] is not None:
    #    the_x = the_x.in_units(units[0])
    #if units[1] is not None and fractional is not True:
    #    the_y = the_y.in_units(units[1])
    output={}
    output['prof']=prof
    output['the_x']=the_x
    output['the_y']=the_y
    return output

plt.clf()
odir=os.environ['HOME']+'/PigPen/'
odir = "./plots_to_sort/"


def return_core_count(self,n_min=1):
    output=[]
    for core_id in self.target_indices:
        if self.target_indices[core_id].size >= n_min:
            output.append(core_id)
    return output


class means_etc():
    def __init__(self,this_looper,n_min=3):
        thtr=this_looper.tr
        core_list = return_core_count(this_looper,n_min)
        self.dmeans = np.zeros_like(core_list,dtype='float')
        self.dstds = np.zeros_like(core_list,dtype='float')
        self.d_logmeans = np.zeros_like(core_list,dtype='float')
        self.d_logstds  = np.zeros_like(core_list,dtype='float')
        self.v_logmeans = np.zeros_like(core_list,dtype='float')
        self.v_logstds  = np.zeros_like(core_list,dtype='float')
        self.variance = np.zeros_like(core_list,dtype='float')
        self.vmeans    = np.zeros_like(core_list,dtype='float')
        self.vstds = np.zeros_like(core_list,dtype='float')
        self.npart = np.zeros_like(core_list,dtype='float')
        self.vrel  =  np.zeros_like(core_list,dtype='float')
        self.volume =  np.zeros_like(core_list,dtype='float')
        self.temp_all_rho=[]

        for i,nc in enumerate(core_list):
            this_density = thtr.c(int(nc),'density')[:,0]
            self.temp_all_rho.append(this_density)
            this_vel = thtr.c(int(nc),'velocity_magnitude')[:,0]
            this_volume = thtr.c(int(nc),'cell_volume')[:,0]
            self.npart[i] = this_density.size
            self.volume[i] = this_volume.sum()
            self.dmeans[i]=this_density.mean()
        
            self.dstds[i] = this_density.std()
            self.d_logmeans[i]=np.exp(np.log(this_density).mean())
            self.d_logstds[i] =np.exp(np.log(this_density).std())
            self.v_logmeans[i]=np.exp(np.log(this_vel).mean())
            self.v_logstds[i] =np.exp(np.log(this_vel).std())
            self.vmeans[i]= this_vel.mean()
            self.vstds[i] = this_vel.std()
            ms = trackage.mini_scrubber(thtr,[nc],do_velocity=True)
            self.variance[i] =  np.mean(ms.rel_vx**2)+np.mean(ms.rel_vy**2)+np.mean(ms.rel_vz**2)

    def make_profiles(self,looper):
        #note to future self, I put this here but didn't finish.
        if 'prof0' not in dir() and False:
            frame = dl.target_frames[this_simname]
            ds0 = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],0,0))
            ds1 = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
            prof0=make_prof(ds0,['density','cell_volume'])
            prof1=make_prof(ds1,['density','cell_volume'])

        if 'prof_vel' not in dir():
            prof_vel=make_prof(ds1,['velocity_magnitude','cell_volume'])



m1 = means_etc( looper1 )
m2 = means_etc( looper2 )
m3 = means_etc( looper3 )

looper1.c='r'
looper2.c='g'
looper3.c='b'

def three_way_bean():
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    return axScatter,axHistx, axHisty

#    # the scatter plot
#    axScatter.scatter(x, y)
#
#    # now determine nice limits by hand
#    binwidth = 0.25
#    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
#    lim = (int(xymax/binwidth) + 1) * binwidth
#
#    axScatter.set_xlim((-lim, lim))
#    axScatter.set_ylim((-lim, lim))
#
#    bins = np.arange(-lim, lim + binwidth, binwidth)
#    axHistx.hist(x, bins=bins)
#    axHisty.hist(y, bins=bins, orientation='horizontal')
#
#    axHistx.set_xlim(axScatter.get_xlim())
#    axHisty.set_ylim(axScatter.get_ylim())
##
#    plt.savefig(outname)

   


if 1:
    if 1:
        ax, ax_den_hist,ax_vel_hist=three_way_bean()
        ok = slice(None)
        ax.scatter(np.log10(m1.dmeans[ok]),np.log10(m1.variance[ok])/2,c=looper1.c,label=looper1.out_prefix, s=0.1)
        ax.scatter(np.log10(m2.dmeans[ok]),np.log10(m2.variance[ok])/2,c=looper2.c,label=looper2.out_prefix, s=0.1)
        ax.scatter(np.log10(m3.dmeans[ok]),np.log10(m3.variance[ok])/2,c=looper3.c,label=looper3.out_prefix, s=0.1)


        ax_den_hist.hist(np.log10(m1.dmeans[ok]), histtype='step',color=looper1.c,label=looper1.out_prefix)
        ax_den_hist.hist(np.log10(m2.dmeans[ok]), histtype='step',color=looper2.c,label=looper2.out_prefix)
        ax_den_hist.hist(np.log10(m3.dmeans[ok]), histtype='step',color=looper3.c,label=looper3.out_prefix)

        ax_vel_hist.hist(np.log10(m1.variance[ok])/2, histtype='step', orientation='horizontal',color=looper1.c)
        ax_vel_hist.hist(np.log10(m2.variance[ok])/2, histtype='step', orientation='horizontal',color=looper2.c)
        ax_vel_hist.hist(np.log10(m3.variance[ok])/2, histtype='step', orientation='horizontal',color=looper3.c)
        axbonk(ax,yscale='linear', xscale='linear',  ylabel=r'$\log_{10} \sigma_v$', xlabel=r'$\log_{10} \langle\rho\rangle$')
        axbonk(ax_vel_hist,yscale='linear', xscale='linear',  ylabel=None, xlabel=r'$N$')
        axbonk(ax_den_hist,yscale='linear', xscale='linear',  ylabel=r'$N$', xlabel=None)
        ax_den_hist.legend(loc=1)
        ax_den_hist.set_xticks([])
        ax_vel_hist.set_yticks([])
        ax.set_xlim(-1,2)
        ax_den_hist.set_xlim(ax.get_xlim())
        ax_vel_hist.set_ylim(ax.get_ylim())

        plt.savefig(odir+'/pre_rho_rms_v_rms.%s'%form)

    if 0:
        ax.clear()
        ax.scatter(dmeans[ok],vstds[ok])
        axbonk(ax,yscale='log', xscale='log', ylabel=r'$\sigma_v$', xlabel=r'$\langle\rho\rangle$')
        fig.savefig(odir+'/pre_rho_mean_v_rms.%s'%form)
    if 0:
        ax.clear()
        ax.scatter(dmeans,vmeans)
        axbonk(ax,yscale='log', xscale='log', ylabel=r'$\langle v \rangle$', xlabel=r'$\langle \rho \rangle$')
        fig.savefig(odir+'/pre_rho_v_mean.%s'%form)

    if 0:
        ax.clear()
        ax.errorbar(d_logmeans,v_logmeans, xerr=d_logstds, yerr=v_logstds)
        axbonk(ax,yscale='log', xscale='log', ylabel=r'$\langle v \rangle$', xlabel=r'$\langle \rho \rangle$')
        fig.savefig(odir+'/%s_pre_rho_v_log_errb.%s'%(this_simname,form))

if 0:
    """density: all the PDFs in one plot"""
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    ax1.plot( prof0['the_x'],prof0['the_y']*128**3,c='k')
    #ax1.plot( prof1['the_x'],prof1['the_y'],c='k')
    for i,nc in enumerate(all_cores):
        this_density = thtr.c([int(nc)],'density')[:,0]
        if this_density.size > 1:
            hist_raw, edges =np.histogram(np.log10(this_density),bins=np.linspace(-2,2,16))
            hhh = hist_raw#/hist_raw.sum()
            centers = 0.5*(edges[1:]+edges[:-1])
            ax1.plot( 10**centers, hhh)
    axbonk(ax1,yscale='log',xscale='log', xlabel=r'$\log_{10}\rho$',ylabel=r'$N$')
    fig.savefig(odir+"density_pdf.%s"%form)
    plt.close(fig)

if 0:
    """velocity: all the PDFs in one plot"""
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    ax1.plot( prof_vel['the_x'],prof_vel['the_y']*128**3,c='k')
    for i,nc in enumerate(all_cores):
        this_vel = thtr.c([int(nc)],'velocity_magnitude')[:,0]
        if this_vel.size > 1:
            hist_raw, edges =np.histogram(this_vel)
            hhh = hist_raw#/hist_raw.sum()
            centers = 0.5*(edges[1:]+edges[:-1])
            ok = hhh>0
            ax1.plot( centers[ok], hhh[ok])
    form = 'pdf'
    axbonk(ax1,yscale='log',xscale='log', xlabel=r'$v$',ylabel=r'$N$')
    fig.savefig(odir+"velocity_pdf.%s"%form)

if 0:
    #relative and tangential histograms
    fig, axes=plt.subplots(2,2)
    axv0 = axes[0][0]
    axv1 = axes[0][1]
    axv2 = axes[1][0]
    axv3 = axes[1][1]
    #for ni,frame in enumerate(thtr.frames):
    for ni,frame in enumerate(thtr.frames):
        for x in axes.flatten():
            x.clear()
        for nc,core_id in enumerate(core_list):
            density = thtr.c([core_id],'density')
            if density.size < 3:
                continue


            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            relative_v = np.sqrt(ms.rel_v2[:,ni])
            #ax3.clear()
            ok = (relative_v >0) * ( relative_v > 1e-10)
            vel_bins = np.logspace(-3,3)
            nums, bin_edges, rect=axv0.hist(relative_v[ok], histtype='step', color=[0.5]*3, bins=vel_bins)
            vr = ms.vr_rel[:,ni]
            vt = np.sqrt( ms.vt2_rel[:,ni])
            nums, bin_edges, rect=axv1.hist(vr[ok], histtype='step', color=[0.5]*3)
            nums, bin_edges, rect=axv2.hist(vt[ok], histtype='step', color=[0.5]*3,bins=vel_bins)

            axv3.scatter(vr[ok],vt[ok],c=[ [0.5]*3]*ok.sum(),s=0.1)
            #ext_x(bin_edges)
            #ext_y(nums)
            #oname = 'plots_to_sort/vels_ind_c%04d_n%04d.png'%(core_id, frame)
        axbonk(axv0,xlabel=r'$\log_{10} v_{rel}$', ylabel=r'$N$',yscale='log',xscale='log')
        axbonk(axv1,xlabel=r'$v_{r}$', ylabel=r'$N$',xscale='linear',yscale='log', xlim=[-15,15])
        axbonk(axv2,xlabel=r'$v_{t}$', ylabel=r'$N$',xscale='log',yscale='log')
        outname = 'plots_to_sort/vel_histograms_n%04d'%frame
        fig.savefig(outname)
        print(outname)

if 0:
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    hist_raw, edges =np.histogram(np.log10(volume**(1./3)))
#,bins=np.linspace(-2,2,16))
    hhh = hist_raw#/hist_raw.sum()
    centers = 0.5*(edges[1:]+edges[:-1])
    ax1.plot( centers, hhh,c=tm(ni),label=ni, marker='*')
    ax1.legend(loc=1)
    ax1.set_title('test')
    fig.savefig('plots_to_sort/sizes.png')
    plt.close(fig)


if 0:
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    tm = rainbow_map(len(thtr.frames))
    for ni,frame in enumerate(thtr.frames):
        hist_raw, edges =np.histogram(thtr.track_dict['velocity_magnitude'][:,ni])
#,bins=np.linspace(-2,2,16))
        hhh = hist_raw#/hist_raw.sum()
        centers = 0.5*(edges[1:]+edges[:-1])
        ax1.plot( centers, hhh,c=tm(ni),label=ni)
    ax1.legend(loc=1)
    ax1.set_title('test')
    fig.savefig('plots_to_sort/vels.png')

#core_list = this_looper.core_list


if 0:
    #velocity PDFs
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    fig2, axes2 = plt.subplots(1,1)
    ax2 = axes2 
    tm = rainbow_map(len(thtr.frames))

    fig3, ax3 = plt.subplots(2,2)
    ax3.clear()
    ext_x=extents()
    ext_x(nar([-2.56e+00, 3.07e+00]))

    ext_y=extents()
    ext_y(nar([0,4.6e3]))
    for ni,frame in enumerate(thtr.frames):
        print(frame)
        relative_vels=[]
        ax3.clear()
        for nc,core_id in enumerate(core_list):
            density = thtr.c([core_id],'density')
            if density.size < 3:
                continue


            #miniscrubber computes distance, r^2, several other quantities
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            relative_v = np.sqrt(ms.rel_v2[:,ni])
            #ax3.clear()
            ok = (relative_v >0) * ( relative_v > 1e-10)
            nums, bin_edges, rect=ax3.hist(np.log(relative_v[ok]), histtype='step', color=[0.5]*3)
            #ext_x(bin_edges)
            #ext_y(nums)
            #oname = 'plots_to_sort/vels_ind_c%04d_n%04d.png'%(core_id, frame)
            #fig3.savefig(oname)
            #print(oname)
            relative_vels += list(  relative_v)
            print( ni, nc)


        bins = np.linspace(0,30,16)
        hist_raw, edges =np.histogram(thtr.track_dict['velocity_magnitude'][:,ni],bins=bins)
        hhh = hist_raw#/hist_raw.sum()
        centers = 0.5*(edges[1:]+edges[:-1])
        ax1.plot( centers, hhh/hhh.sum(),c=tm(ni),label=ni)

        hist_rel, edges2 =np.histogram(relative_vels,bins=bins)
        centers2 = 0.5*(edges2[1:]+edges2[:-1])
        ax2.plot( centers2, hist_rel/hist_rel.sum(),c=tm(ni),label=ni)
        #axbonk(ax3, xlabel=r'$\log_{10}v_{rel}$', ylabel=r'$N$', xlim=ext_x.minmax, ylim=ext_y.minmax, yscale='log')
        axbonk(ax3, xlabel=r'$\log_{10}v_{rel}$', ylabel=r'$N$', xlim=[-3,3], ylim=ext_y.minmax, yscale='log')
        outname='plots_to_sort/vels_ind_n%04d.png'%ni; print(outname)
        fig3.savefig(outname)
    axbonk(ax1, xlabel=r'$v$',ylabel=r'$V(v_{total}|*)$')
    axbonk(ax2, xlabel=r'$v$',ylabel=r'$V(v_{rel}|*)$')

    ax1.legend(loc=1)
    fig.savefig('plots_to_sort/vels.png')
    fig2.savefig('plots_to_sort/vels_boosted.png')


#ok=npart>1
#fig,ax=plt.subplots()

if 0:
    """Several different PDFs"""
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    #ax1.plot( prof1['the_x'],prof1['the_y'],c='k')
    for i,nc in enumerate(all_cores):
        this_density = thtr.c(int(nc),'density')[:,0]
        if this_density.size > 1:
            ax1.clear()
            ax1.plot( prof0['the_x'],prof0['the_y']*128**3,c='k')
            hist_raw, edges =np.histogram(np.log10(this_density),bins=np.linspace(-2,2,16))
            hhh = hist_raw#/hist_raw.sum()
            centers = 0.5*(edges[1:]+edges[:-1])
            ax1.plot( 10**centers, hhh)
            form = 'pdf'
            axbonk(ax1,yscale='log',xscale='log', xlabel=r'$\log_{10}\rho$',ylabel=r'$N$', ylim=[0.5,1e5])
            oname = odir+"density_pdf_c%04d.%s"%(nc,form)
            fig.savefig(oname)
            print(oname)
    plt.close(fig)

if 0:
    if 1:
        ax.clear()
        ax.scatter(vmeans[ok],vstds[ok])
        axbonk(ax,xscale='linear',yscale='linear',xlim=[0,20],ylim=[0,20],xlabel=r'$\langle v\rangle$', ylabel=r'$\sigma_v$')
        fig.savefig(odir+'/pre_meanv_rms_v.%s'%form)
    if 1:
        ax.clear()
        ax.scatter(logdstds[ok],vstds[ok],c='k')
        #ax.yscale('log')
        #ax.xscale('log')
        #ax.ylabel(r'$\log_{10} ||v||$')
        x = extents()
        x(logdstds[ok])
        y = extents()
        y(vstds[ok])
        axbonk(ax,xscale='linear',yscale='linear',xlim=x.minmax,ylim=y.minmax, xlabel=r'$\sigma_{\ln \rho}$', ylabel=r'$\sigma_v$')
        fig.savefig(odir+'/pre_logrho_rms_v_rms.%s'%form)


if 0:
    ax.clear()
    for i,nc in enumerate(all_cores):
        this_density = thtr.c(int(nc),'velocity_magnitude')[:,0]
        #ax.hist(np.log10(this_density),histtype='step',color='k')
        ax.hist(this_density,histtype='step',color='k')
    ax.yscale('log')
    #ax.xlabel(r'$\log_{10} ||v||$')
    ax.xlabel(r'$||v||$')
    ax.ylabel(r'$N(v)$')
    fig.savefig(odir+'/p56_allpdf_vel_notlog.%s'%form)

if 0:
    ax.clear()
    for i,nc in enumerate(all_cores):
        this_density = thtr.c(int(nc),'density')[:,0]
        ax.hist(np.log10(this_density),histtype='step',color='k')
    ax.yscale('log')
    fig.savefig(odir+'/p47_allhist.png')


if 0:
    means = np.zeros_like(all_cores,dtype='float')
    stds = np.zeros_like(all_cores,dtype='float')
    npart = np.zeros_like(all_cores,dtype='float')
    for i,nc in enumerate(all_cores):
        this_density = np.log(thtr.c(int(nc),'density')[:,0])
        means[i]=this_density.mean()
        npart[i] = this_density.size
        stds[i]=np.sqrt(np.sum((this_density-means[i])**2)/npart[i])
    ax.clear()
    ax.scatter(means,stds)
    ax.ylabel(r'$\sqrt{\langle \ln \rho^2 \rangle}$')
    ax.xlabel(r'$\langle \ln \rho \rangle$')
    fig.savefig(odir+'/p47_log_mean_std.%s'%form)

