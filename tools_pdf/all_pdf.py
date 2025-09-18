
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
file_list = glob.glob("/scratch1/dcollins/Paper19/Datasets/all_primitives/*h5")[:10]
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='test'
#file_list = file_list[:2]


if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    del dmeans, prof0, prof1, prof_vel
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))

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

if 'prof0' not in dir():
    ds0 = yt.load("%s/DD%04d/data%04d"%(dl.enzo_directory,0,0))
    ds1 = yt.load("%s/DD%04d/data%04d"%(dl.enzo_directory,125,125))
    prof0=make_prof(ds0,['density','cell_volume'])
    prof1=make_prof(ds1,['density','cell_volume'])

if 'prof_vel' not in dir():
    ds = yt.load("%s/DD%04d/data%04d"%(dl.enzo_directory,0,0))
    prof_vel=make_prof(ds1,['velocity_magnitude','cell_volume'])

#all_cores=all_cores[:5]
all_densities={}
if 'dmeans' not in dir():
    dmeans = np.zeros_like(all_cores,dtype='float')
    dstds = np.zeros_like(all_cores,dtype='float')
    logdmeans = np.zeros_like(all_cores,dtype='float')
    logdstds = np.zeros_like(all_cores,dtype='float')
    vmeans = np.zeros_like(all_cores,dtype='float')
    vstds = np.zeros_like(all_cores,dtype='float')
    npart = np.zeros_like(all_cores,dtype='float')
    for i,nc in enumerate(all_cores):
        this_density = thtr.c(int(nc),'density')[:,0]
        npart[i] = this_density.size
        this_vel = thtr.c(int(nc),'velocity_magnitude')[:,0]
        dmeans[i]=this_density.mean()
        dstds[i] = this_density.std()
        logdmeans[i]= np.log(this_density).mean()
        logdstds[i] =  np.log(this_density).std()
        vmeans[i]= this_vel.mean()
        vstds[i] = this_vel.std()

if 1:
    """all the PDFs in one plot"""
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    ax1.plot( prof0['the_x'],prof0['the_y']*128**3,c='k')
    #ax1.plot( prof1['the_x'],prof1['the_y'],c='k')
    for i,nc in enumerate(all_cores):
        this_density = thtr.c(int(nc),'density')[:,0]
        if this_density.size > 1:
            hist_raw, edges =np.histogram(np.log10(this_density),bins=np.linspace(-2,2,16))
            hhh = hist_raw#/hist_raw.sum()
            centers = 0.5*(edges[1:]+edges[:-1])
            ax1.plot( 10**centers, hhh)
    form = 'pdf'
    axbonk(ax1,yscale='log',xscale='log', xlabel=r'$\log_{10}\rho$',ylabel=r'$N$')
    fig.savefig(odir+"density_pdf.%s"%form)
    plt.close(fig)

if 0:
    """all the PDFs in one plot"""
    fig, axes = plt.subplots(1,1)
    ax1 = axes
    ax1.plot( prof_vel['the_x'],prof_vel['the_y']*128**3,c='k')
    for i,nc in enumerate(all_cores):
        this_vel = thtr.c(int(nc),'velocity_magnitude')[:,0]
        if this_vel.size > 1:
            hist_raw, edges =np.histogram(this_vel)
            hhh = hist_raw#/hist_raw.sum()
            centers = 0.5*(edges[1:]+edges[:-1])
            ok = hhh>0
            ax1.plot( centers[ok], hhh[ok])
    form = 'pdf'
    axbonk(ax1,yscale='log',xscale='log', xlabel=r'$v$',ylabel=r'$N$')
    fig.savefig(odir+"velocity_pdf.%s"%form)
ok=npart>1
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
        ax.scatter(logdstds[ok],vstds[ok])
        #ax.yscale('log')
        #ax.xscale('log')
        #ax.ylabel(r'$\log_{10} ||v||$')
        axbonk(ax,xscale='linear',yscale='linear',xlim=[0,20],ylim=[0,20], xlabel=r'$\sigma_{\ln \rho}$', ylabel=r'$\sigma_v$')
        fig.savefig(odir+'/pre_logrho_rms_v_rms.%s'%form)
if 0:
    if 1:
        ax.clear()
        ax.scatter(dstds[ok],vstds[ok])
        ax.yscale('log')
        ax.xscale('log')
        #ax.ylabel(r'$\log_{10} ||v||$')
        ax.ylabel(r'$\sigma_v$')
        ax.xlabel(r'$\sigma_{\rho}$')
        fig.savefig(odir+'/pre_rho_rms_v_rms.%s'%form)
    if 1:
        ax.clear()
        ax.scatter(dmeans[ok],vstds[ok])
        ax.yscale('log')
        ax.xscale('log')
        #ax.ylabel(r'$\log_{10} ||v||$')
        ax.ylabel(r'$\sigma_v$')
        ax.xlabel(r'$\langle\rho\rangle$')
        fig.savefig(odir+'/pre_rho_mean_v_rms.%s'%form)
    if 1:
        ax.clear()
        ax.scatter(dmeans,vmeans)
        ax.yscale('log')
        ax.xscale('log')
        #ax.ylabel(r'$\log_{10} ||v||$')
        ax.ylabel(r'$\langle v \rangle$')
        ax.xlabel(r'$\langle \rho \rangle$')
        fig.savefig(odir+'/pre_rho_v_mean.%s'%form)
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

