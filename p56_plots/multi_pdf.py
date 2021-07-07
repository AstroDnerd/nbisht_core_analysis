
from starter2 import *

import data_locations as dl
reload(dl)
reload(trackage)
plt.close('all')

if 'sim_list' not in dir():
    this_simname = 'u11'
    out_prefix='u11'
    form = 'pdf'
    files_u05=  glob.glob("/scratch1/dcollins/Paper19/Datasets/all_primitives/*h5")
    files_u11 = ["../Datasets/u11_primitives_cXXXX_n0000.h5"]; out_prefix='u11'
    files_u10 = ["../Datasets/u10_primitives_cXXXX_n0000.h5"]; out_prefix='u11'


    class sim_details():
        def __init__(self,name,file_list,beta="", loop=None):

            self.name=name
            self.file_list=file_list
            self.my_looper=None
            self.tracker = None
            self.density_profiles=[]
            self.velocity_profiles=[]
            self.means=None
            self.loop=loop
            self.beta=beta


    import three_loopers_1tff as tl
    u05s = sim_details('u05',files_u05,beta=0.2, loop=tl.looper1)
    u10s = sim_details('u10',files_u10,beta=2.0, loop=tl.looper2)
    u11s = sim_details('u11',files_u11,beta=20, loop=tl.looper3)


        
    sim_list=[u05s,u10s,u11s]

     
    if 0:
        for sim in sim_list:
            if sim.my_looper is None:
                sim.my_looper=looper.core_looper(directory=dl.sims[this_simname])
                for nfile,fname in enumerate(sim.file_list):
                    sim.my_looper.load_loop(fname)
                    print( "File %d of %d"%(nfile,len(sim.file_list)))
                sim.tracker = sim.my_looper.tr
                sim.tracker.sort_time()

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
#rm = rainbow_map(len(all_cores))
#tm = rainbow_map(15)


plt.clf()
odir=os.environ['HOME']+'/PigPen/'
odir = "./plots_to_sort/"

class means_etc():
    def __init__(self,thtr,core_list=None):
        if core_list is None:
            core_list = np.unique(thtr.core_ids)
        self.core_list=core_list
        self.thtr = thtr
        self.dmeans = np.zeros_like(core_list,dtype='float')
        self.dstds = np.zeros_like(core_list,dtype='float')
        self.d_logmeans = np.zeros_like(core_list,dtype='float')
        self.d_logstds  = np.zeros_like(core_list,dtype='float')
        self.v_logmeans = np.zeros_like(core_list,dtype='float')
        self.v_logstds  = np.zeros_like(core_list,dtype='float')
        self.vmeans    = np.zeros_like(core_list,dtype='float')
        self.vstds = np.zeros_like(core_list,dtype='float')
        self.vstds_xyz = np.zeros_like(core_list,dtype='float')
        self.npart = np.zeros_like(core_list,dtype='float')
        self.vrel  =  np.zeros_like(core_list,dtype='float')
        self.volume =  np.zeros_like(core_list,dtype='float')

        self.dmean_volume = np.zeros_like(core_list,dtype='float')
        self.dstd_volume_1 = np.zeros_like(core_list,dtype='float')
        self.dstd_volume_2 = np.zeros_like(core_list,dtype='float')
        self.mass = np.zeros_like(core_list,dtype='float')
        self.Bmean = np.zeros_like(core_list,dtype='float')
        self.Brms = np.zeros_like(core_list,dtype='float')

        self.densities=[]
        self.volumes=[]
        self.bfields=[]



        for i,nc in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,int(nc),do_velocity=True)
            frame = 1
            this_density = thtr.c(int(nc),'density')[:,frame]
            this_vel = thtr.c(int(nc),'velocity_magnitude')[:,frame]
            this_volume = thtr.c(int(nc),'cell_volume')[:,frame]
            this_magnetic = thtr.c(int(nc),'magnetic_field_strength')[:,frame]
            self.densities.append(this_density)
            self.volumes.append(this_volume)
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

            rvx = ms.rel_vx #vx - <vx>
            rvy = ms.rel_vy
            rvz = ms.rel_vz
            sigma_2 = (rvx**2).sum()/rvx.size+ (rvy**2).sum()/rvx.size+ (rvz**2).sum()/rvx.size
            self.vstds_xyz[i] =np.sqrt(sigma_2)

            self.dmean_volume[i] = (this_density*this_volume).sum()/this_volume.sum()
            self.dstd_volume_1[i] =np.sqrt( ( (this_density-1)**2*this_volume).sum()/this_volume.sum())
            self.dstd_volume_2[i] =np.sqrt( ( (this_density-self.dmean_volume[i])**2*this_volume).sum()/this_volume.sum())

            self.mass[i] = (this_density*this_volume).sum()
            self.bfields.append(this_magnetic)
            self.Bmean[i] = (this_magnetic*this_volume).sum()/this_volume.sum()
            self.Brms[i] = np.sqrt(((this_magnetic)**2*this_volume).sum()/this_volume.sum())


for sim in sim_list:
    if sim.means is None:
        sim.means = means_etc( sim.loop.tr )

if 0:
    fig, ax = plt.subplots(1,3, figsize=(12,4))
    ax0=ax[0]
    ax1=ax[1]
    ax2=ax[2]
    frame=0
    rm = rainbow_map(10)
    for sim in sim_list[:1]:
        order = np.argsort( sim.means.dmean_volume )
        for nc, core_order in enumerate(order):
            core_id = sim.means.core_list[core_order]
            this_rho = sim.means.thtr.c(int(core_id), 'density')[:,frame]
            this_log_rho = np.log10(this_rho)
            if this_rho.size < 10:
                continue

            ax0.hist( this_log_rho, histtype='step')
            ax0.set_yscale('log')
            ax1.scatter( sim.means.dmean_volume[core_order],
                        sim.means.dstd_volume_2[core_order])
            axbonk(ax1, xlim=[1e-1, 10],ylim=[1e-1,10], yscale='log',xscale='log')

        ok = sim.means.npart 
        ratio = sim.means.dstd_volume_2/sim.means.dmean_volume 
        ax2.scatter( sim.means.npart, sim.means.dmean_volume)#, histtype='step')
        axbonk(ax2,xscale='log', yscale='log', xlabel=r'$N$', ylabel=r'$\rho$')
        ax1.plot( sim.means.dmean_volume, sim.means.dmean_volume,c='k')
        outname=odir+'/multi_hist_X%03d.png'%(nc)
        fig.savefig(outname)
        print(outname)


outname='x'



fig,ax=plt.subplots(1,1)

if 1:
    rho_ext = extents()
    v_ext = extents()
    only_once=True
    for sim in sim_list:
        ok = sim.means.npart > 10
        kwargs={}
        kwargs['label']=label=r'$\beta=%0.1f$'%sim.beta
        kwargs['marker']= {'u05':'*','u10':'v','u11':'^'}[sim.name]
        cset={'u05':0.0,'u10':0.5,'u11':0.7}
        kwargs['color'] = [cset[sim.name]]*3
        if 0:
            # rho sigma-v
            if only_once:
                ax.plot( [1.0,1.0],[0,100],c=[0.5]*4)
                ax.plot( [0.1,100],[1,1],c=[0.5]*4)
                only_once=False
            ax.scatter(sim.means.dmeans[ok],sim.means.vstds_xyz[ok],**kwargs)
            rho_ext(sim.means.dmeans[ok])
            v_ext(sim.means.vstds_xyz[ok])
            axbonk(ax,yscale='log', xscale='log', ylabel=r'$\sigma_v$', xlabel=r'$\langle\rho\rangle$',
                   xlim=rho_ext.minmax,ylim=v_ext.minmax)
            ax.legend(loc=3)
            outname=odir+'/pre_rho_mean_v_rms.%s'%(form)
        if 0:
            # rho sigma-rho, dumb way
            the_x = sim.means.dmeans[ok]
            the_y = sim.means.dstds[ok]
            pfit = np.polyfit( the_x, the_y,1)
            new_y = the_y
            ax.scatter(the_x,new_y)
            #ax.plot( the_x, pfit[0]*the_x+pfit[1],c='k')
            axbonk(ax,yscale='log', xscale='log',  xlabel=r'$\langle\rho\rangle$', ylabel=r'$\sigma_{\rho}$')
            outname=odir+'/pre_rho_mean_rho_rms.%s'%(form)
        if 0:
            #this is the good one.
            # rho sigma-rho 
            the_x = sim.means.dmean_volume[ok]
            the_y = sim.means.dstd_volume_2[ok]
            pfit = np.polyfit( the_x, the_y,1)
            new_y = the_y
            ax.scatter(the_x,new_y)
            ax.plot( the_x, the_x,c='k')
            axbonk(ax,yscale='log', xscale='log',  xlabel=r'$\langle\rho\rangle$', ylabel=r'$\sigma_{\rho}$')
            outname = odir+'/pre_rho_mean_rho_rms_avg2.%s'%(form)
        if 1:
            # brms vs v
            Bmean = sim.means.Bmean
            Brms  = sim.means.Brms
            the_y = Brms
            the_x = sim.means.vstds_xyz**2#/Bmean
            the_x =the_x[ok]; the_y=the_y[ok]
            #pfit = np.polyfit( the_x, the_y,1)
            ax.scatter(the_x,the_y,**kwargs)
            ax.set_yscale('log')
            ax.set_xscale('log')
            #xlim = [the_x.min(),the_x.max()]
            #ylim = [the_y.min(),the_y.max()]
            #ax.set_xlim(xlim); ax.set_ylim(ylim)
            #axbonk(ax,yscale='log', xscale='log',  ylabel=r'$\sigma_B$', xlabel=r'$\sigma_{v}^2/B$')
            outname = odir+'/pre_v_rms_brms.%s'%(form)
        if 0:
            #
            # sigma-v sigma-rho
            the_x = sim.means.vstds_xyz[ok]
            the_y = sim.means.dstd_volume_2[ok]
            pfit = np.polyfit( the_x, the_y,1)
            ax.scatter(the_x,the_y)
            axbonk(ax,yscale='log', xscale='log',  ylabel=r'$\sigma_\rho$', xlabel=r'$\sigma_{v}$')
            outname = odir+'/pre_v_rms_rho_rms_avg2.%s'%(form)

        if 0:
            the_x = sim.means.mass[ok]
            the_y = sim.means.dstd_volume_2[ok]
            pfit = np.polyfit( the_x, the_y,1)
            ax.scatter(the_x,the_y)
            axbonk(ax,yscale='log', xscale='log',  xlabel=r'$M$', ylabel=r'$\sigma_{\rho}$')
            outname = odir+'/pre_mass_rho_rms_avg2.%s'%(form)

        if 0:
            the_x = sim.means.volume[ok]
            the_y = sim.means.dstd_volume_2[ok]
            pfit = np.polyfit( the_x, the_y,1)
            ax.scatter(the_x,the_y)
            axbonk(ax,yscale='log', xscale='log',  xlabel=r'$Volume$', ylabel=r'$\sigma_{\rho}$')
            outname=odir+'/pre_volume_rho_rms_avg2.%s'%(form)


        if 0:
            ax.scatter(sim.means.dmeans,sim.means.vmeans)
            axbonk(ax,yscale='log', xscale='log', ylabel=r'$\langle v \rangle$', xlabel=r'$\langle \rho \rangle$')
            outname = odir+'/pre_rho_v_mean.%s'%(form)

        if 0:
            ax.errorbar(sim.means.d_logmeans,sim.means.v_logmeans, 
                        xerr=sim.means.d_logstds, yerr=sim.means.v_logstds)
            axbonk(ax,yscale='log', xscale='log', ylabel=r'$\langle v \rangle$', xlabel=r'$\langle \rho \rangle$')
            outname = odir+'/pre_rho_v_log_errb.%s'%(form)
    ax.legend(loc=0)
    fig.savefig(outname)
    print(outname)

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

fig,ax=plt.subplots()

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

