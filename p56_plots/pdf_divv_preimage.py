
"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import data_locations as dl
import xtra_energy
reload(loop_apps)
from scipy.optimize import curve_fit


#
#
#



#import three_loopers_1tff as tl
#this_looper = tl.looper1
#this_simname = this_looper.out_prefix
this_simname = 'u201'


#import testing.cic_test as cic
import testing.early_mask as em
reload(em)


def toplot(prof,quan = 'cell_volume'):
    xbins = prof.x_bins
    bin_center = 0.5*(xbins[1:]+xbins[:-1])
    bin_widths = xbins[1:]-xbins[:-1]
    pdf = prof[quan]
    pdf = pdf/bin_widths
    return xbins, bin_center,pdf,bin_widths
def gaussian(the_x,norm,x0,sigma):
    #return norm*np.exp( -(the_x-x0)**2/(2*sigma**2))
    return norm/np.sqrt(2*np.pi*sigma**2)*np.exp( -(the_x-x0)**2/(2*sigma**2))

if 'ds' not in dir():
    frame=100
    #ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
    ds = yt.load("%s/DD%04d/data%04d"%(dl.sims['u201'],frame,frame))
    em.add_tracer_density(ds)
    ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    #all_target_indices = np.concatenate( list(this_looper.target_indices.values()))
    all_target_indices = nar(dpy('target_ind.h5','target_ind')).flatten()
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    deposit_tuple = ("deposit","target_particle_volume")
    #ad[deposit_tuple]
        
if 'prof_mask' not in dir():
    bins1 = np.linspace(-1000,1000,100)
    bins2 = np.logspace(np.log10(1000), np.log10(1e5),100)
    bins = np.concatenate([-bins2[::-1], bins1, bins2])
    bins = ds.arr(np.unique(bins), ad['velocity_divergence'].units)
    override_bins = {'velocity_divergence':bins}

    prof_all  = yt.create_profile(ad,bin_fields=['velocity_divergence'],fields=['cell_volume'],weight_field=None, override_bins=override_bins)
    prof_mask = yt.create_profile(ad,bin_fields=['velocity_divergence'],fields=[deposit_tuple],weight_field=None, override_bins=override_bins)
    prof_all.save_as_dataset('%s_divv_pdf_all_n%04d.h5'%(this_simname,frame))
    prof_mask.save_as_dataset('%s_divv_pdf_mask_n%04d.h5'%(this_simname,frame))


if 1:
    #get data
    bbb1, bcen1, vals1, db= toplot(prof_all)
    bbb2, bcen2, vals2, db = toplot(prof_mask,quan=deposit_tuple[1])
    db1 = bbb1[1:]-bbb1[:-1]
    db2 = bbb2[1:]-bbb2[:-1]

if 1:
    fig,ax=plt.subplots(1,1)
    ax.plot( bcen1,vals1,'k',linewidth=2, label=r'$V(\nabla \cdot v)$')
    ax.plot( bcen2,vals2,'k--',linewidth=2, label=r'$V(\nabla \cdot v|*)$')
    axbonk(ax,yscale='log',ylabel='V',xlabel=r'$\nabla \cdot v$')
    ax.set_xscale('symlog',linthresh=1000)
    outname = 'plots_to_sort/%s_hist_divv_n%04d.pdf'%(this_simname,frame)
    fig.savefig(outname)
    ax.clear()
    ax.plot( bcen1,vals1/db1,'k',linewidth=2, label=r'$V(\nabla \cdot v)$')
    ax.plot( bcen2,vals2/db2,'k--',linewidth=2, label=r'$V(\nabla \cdot v|*)$')
    axbonk(ax,yscale='log',ylabel='V',xlabel=r'$\nabla \cdot v$')
    ax.set_xscale('symlog',linthresh=1000)
    outname = 'plots_to_sort/%s_pdf_divv_n%04d.pdf'%(this_simname,frame)
    fig.savefig(outname)

if 0:
    fig2, ax2=plt.subplots(1,1)
    cuml_all  = np.cumsum(vals1)
    cuml_mask = np.cumsum(vals2)
    ax2.plot( bcen1, cuml_all/cuml_all[-1], c='k')
    ax2.plot( bcen2, cuml_mask/cuml_mask[-1], 'k--')
    #ax2.plot( bcen1,vals1,c='k')
    #ax2.plot( bcen2,vals2,'k--')
    axbonk(ax2,xlabel=r'$\nabla \cdot v$',ylabel=r'$\int V(nabla \cdot v)$',xscale='log',yscale='log')
    outname = 'plots_to_sort/%s_cuml_divv_n%04d.pdf'%(this_simname,frame)
    fig2.savefig(outname)
    print(outname)


if 0:
    #fit for gaussians
    fits1, cov1 = curve_fit(gaussian,np.log10(bcen1),vals1, p0=[1,1,1])
    a1, mu1, sig1 = fits1
    fits2, cov2 = curve_fit(gaussian,np.log10(bcen2),vals2, p0=[1,1,1])
    a2, mu2, sig2 = fits2

if 0:
    #overplot gaussians
    ax.plot( bcen2, gaussian(np.log10(bcen2), *fits2),'k--',linewidth=0.5)
    ax.plot( bcen1, gaussian(np.log10(bcen1), *fits1), 'k--', linewidth=2)

if 0:
    #comput and plot the ratio
    ratio = vals2/vals1
    ok = bcen1>0.1
    ok = np.logical_and(ratio>0, ok)
    ax.plot( bcen1[ratio>0], ratio[ratio>0],label=r"$V(*|\rho)$",c=[0.5]*4)

if 0:
    #fit powerlaw for ratio
    from scipy.optimize import curve_fit
    def powerlaw(r,rho0, r0, alpha):
        return alpha*np.log10(r/r0) + np.log10(rho0)
    popt, pcov = curve_fit(powerlaw, bcen1[ok], np.log10(ratio[ok]), p0=[1,1,-2])

if 0:
    #plot powerlaws
    ax.plot( bcen1[ok], 10**powerlaw(bcen1[ok], popt[0], popt[1],popt[2]),label=r'$\rho^{%0.2f}$'%(popt[2]))
    #ax.plot( bcen1[ok], 10**powerlaw(bcen1[ok], popt[0], popt[1],0.5))
    #ax.plot( bcen1, bcen1**popt[2]*gaussian(np.log10(bcen1), *fits1)*vals2.max(), 'k--', linewidth=2)

if 0:
    #plot rho^2 P(rho)
    n=0
    for core in this_looper.target_indices:
        n+=this_looper.target_indices[core].size
    eta1 = n/128**3
    p_star_given_rho = vals1*bcen1**popt[2]
    eta = vals2.max()/p_star_given_rho.max()
    p_star_given_rho *= eta1
    ax.plot( bcen1,p_star_given_rho,linewidth=2, label=r'$\eta_1\rho^{%0.2f}P(\rho)$'%popt[2], linestyle='--',c=[0.5]*4)




if 0:
    #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
    outname = "plots_to_sort/%s_pdf_density_preimage_fits.pdf"%this_simname
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='linear',yscale='linear')
    axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
    ax.legend(loc=3)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)
    

if 0:
    ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
    ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    fp={}
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    fp['target_indices']=all_target_indices
    fp['mask_to_get']=np.zeros_like(all_target_indices,dtype='int32')
    for ns, xs in enumerate([0.,0.5,0.75]):
        sp=yt.SlicePlot(ds,fields=['deposit_target_particles_boolean'],axis=0,center=[0.5]*3,field_parameters=fp)
        sp.save("plots_to_sort/u05_n%04d_sl%04d"%(frame,ns))
        sp=yt.SlicePlot(ds,fields=['PotentialField'],axis=0,center=[0.5]*3,field_parameters=fp)
        sp.save("plots_to_sort/u05_n%04d_sl%04d"%(frame,ns))

