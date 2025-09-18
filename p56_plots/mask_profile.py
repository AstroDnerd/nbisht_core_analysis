
"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import xtra_energy
reload(loop_apps)
from scipy.optimize import curve_fit
core_list=looper.get_all_nonzero()
frame_list=[0]# range(0,130,10)
fields=['density']  
#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
save_field = '../Datasets/all_cores_n0000.h5'
if 'this_looper' not in dir() and os.path.exists(save_field) and True:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,savefile=save_field)
    this_looper.derived = [xtra_energy.add_force_terms] #for some reason derived quantities aren't saved.

if 'this_looper' not in dir() and False:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='datasets_small/bad_particles.h5')
    this_looper.get_tracks()
    this_looper.save('all_cores_n%04d.h5'%0)


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

if 1:
    frame=0
    ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
    em.add_tracer_density(ds)
    ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
    deposit_tuple = ("deposit","target_particle_volume")
    #ad[deposit_tuple]
        
if 'prof_mask_density' not in dir():
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    #bins={'velocity_x':np.linspace(-25,25,64)}
    #bins['PotentialField']= np.linspace(-50,50,64)
    bins=None
    print('PROFILE density')
    prof_all_density  = yt.create_profile(ad,bin_fields=['density'],fields=['cell_volume'],weight_field=None, override_bins=bins)
    print('PROFILE deposit target')
    prof_mask_density = yt.create_profile(ad,bin_fields=['density'],fields=[deposit_tuple],weight_field=None, override_bins=bins)
    print('PROFILE DONE')

if 1:
    fig,ax=plt.subplots(1,1)
    ax.plot( bcen1,vals1,'k',linewidth=2, label=r'$V(\rho)$')
    ax.plot( bcen2,vals2,'k--',linewidth=2, label=r'$V(\rho|*)$')


if 0:
    #fit for gaussians
    bbb1, bcen1, vals1, db= toplot(prof_all_density)
    fits1, cov1 = curve_fit(gaussian,np.log10(bcen1),vals1, p0=[1,1,1])
    a1, mu1, sig1 = fits1
    bbb2, bcen2, vals2, db = toplot(prof_mask_density,quan=deposit_tuple[1])
    fits2, cov2 = curve_fit(gaussian,np.log10(bcen2),vals2, p0=[1,1,1])
    a2, mu2, sig2 = fits2

if 0:
    #overplot gaussians
    ax.plot( bcen2, gaussian(np.log10(bcen2), *fits2),'k--',linewidth=0.5)
    ax.plot( bcen1, gaussian(np.log10(bcen1), *fits1), 'k--', linewidth=2)

if 1:
    #comput and plot the ratio
    ratio = vals2/vals1
    Pstar = 1./128**3*211  
    ok = bcen1>0.1
    ok = np.logical_and(ratio>0, ok)
    ax.plot( bcen1[ratio>0], ratio[ratio>0],label=r"$V(*|\rho)$",c=[0.5]*4)

if 1:
    #plot rho^2 P(rho)
    p_star_given_rho = vals1*bcen1**0.5
    p_star_given_rho *= vals2.max()/p_star_given_rho.max()
    ax.plot( bcen1,p_star_given_rho,linewidth=2, label=r'$\rho^{1/2}P(\rho)$', linestyle='--',c=[0.5]*4)

if 0:
    #fit powerlaw for ratio
    from scipy.optimize import curve_fit
    def powerlaw(r,rho0, r0, alpha):
        return alpha*np.log10(r/r0) + np.log10(rho0)
    popt, pcov = curve_fit(powerlaw, bcen1[ok], np.log10(ratio[ok]), p0=[1,1,-2])

if 0:
    #plot powerlaws
    ax.plot( bcen1[ok], 10**powerlaw(bcen1[ok], popt[0], popt[1],popt[2]))
    ax.plot( bcen1[ok], 10**powerlaw(bcen1[ok], popt[0], popt[1],0.5))
    ax.plot( bcen1, bcen1**popt[2]*gaussian(np.log10(bcen1), *fits1)*vals2.max(), 'k--', linewidth=2)



if 1:
    #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
    outname = "plots_to_sort/pdf_density_preimage_fits.pdf"
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

