
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

#this_simname = 'u11'


frame_list=[0]# range(0,130,10)
core_list=looper.get_all_nonzero(dl.n_particles[this_simname])
fields=['density']  
#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
#save_field = '../Datasets/all_cores_n0000.h5'
#this_simname = 'u10'
#save_field = '../Datasets/u10_primitives_cXXXX_n0000.h5'
#save_field = '../Datasets/u11_primitives_cXXXX_n0000.h5'
save_field = '../Datasets/%s_prim_phi_cXXXX_n0000.h5'%this_simname
if 'this_looper' not in dir() and os.path.exists(save_field):
    directory = dl.sims[this_simname]
    this_looper = looper.core_looper(directory= directory,savefile=save_field)
    this_looper.derived = [xtra_energy.add_force_terms] #for some reason derived quantities aren't saved.

if 'this_looper' not in dir() and False:
    directory = dl.sims[this_simname]
    this_looper = looper.core_looper(directory= directory,
                                     derived=[],
                                     sim_name = this_simname,
                                     out_prefix = 'plots_to_sort/%s_test'%this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name=dl.peaks[this_simname],
                                     bad_particle_list=dl.bad_particles[this_simname])
    this_looper.get_tracks()


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
        
if 'prof_all_pot' not in dir():
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    bins={'PotentialField':np.linspace(-32,32,64)}
    #bins['PotentialField']= np.linspace(-50,50,64)
    #bins=None
    prof_all_pot  = yt.create_profile(ad,bin_fields=['PotentialField'],fields=['cell_volume'],weight_field=None, override_bins=bins)
    prof_mask_pot = yt.create_profile(ad,bin_fields=['PotentialField'],fields=[deposit_tuple],weight_field=None, override_bins=bins)
    prof_all_pot.save_as_dataset('%s_potential_pdf_all_n0000.h5'%this_simname)
    prof_mask_pot.save_as_dataset('%s_potential_pdf_mask_n0000.h5'%this_simname)

if 1:
    plt.close('all')
    bbb1, bcen1, vals1, db= toplot(prof_all_pot)
    bbb2, bcen2, vals2, db = toplot(prof_mask_pot,quan=deposit_tuple[1])
    fig,ax=plt.subplots(1,1)

    ok1 = vals1>0
    ok2 = vals2>0
    ax.plot( bcen1[ok1],vals1[ok1],'k',linewidth=2, label=r'$V(v)$')


    mass_fraction = 1./128**3*211  
    particle_fraction_scale= 1# all_target_indices.size/128.**3
    ax.plot( bcen2[ok2],vals2[ok2]/particle_fraction_scale,'k--',linewidth=2, label=r'$V(v|*)$')

if 1:
    #comput and plot the ratio
    ratio = vals2/vals1
    ok = bcen1>0.1
    ok = np.logical_and(ratio>0, ok)
    ax.plot( bcen1[ratio>0], ratio[ratio>0],label=r"$V(*|v)$",c=[0.5]*4)

if 0:
    #fit powerlaw for ratio
    #note this doesn't work.
    from scipy.optimize import curve_fit
    def powerlaw(r,rho0, r0, alpha):
        return alpha*np.log10(r/r0) + np.log10(rho0)
    popt, pcov = curve_fit(powerlaw, bcen1[ok], np.log10(ratio[ok]), p0=[1,1,-2])
    ax.plot( bcen1[ok], 10**powerlaw(bcen1[ok], popt[0], popt[1],popt[2]),label=r'play %0.2f'%(popt[2]))


if 1:
    #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
    outname = "plots_to_sort/%s_pdf_potential_preimage_fits.pdf"%this_simname
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='linear',yscale='linear')
    axbonk(ax,xlabel=r'$\Phi$',ylabel=r'$V(\Phi)$',xscale='linear',yscale='log')
    ax.legend(loc=3)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)
    
if 1:
    fig2, ax2=plt.subplots(1,1)
    cuml_phi_all  = np.cumsum(vals1)
    cuml_phi_mask = np.cumsum(vals2)
    ax2.plot( bcen1, cuml_phi_all, c='k')
    ax2.plot( bcen2, cuml_phi_mask/cuml_phi_mask[-1], 'k--')
    #ax2.plot( bcen1, vals1, c='k')
    #ax2.plot( bcen2, vals2, 'k--')
    axbonk(ax2,xlabel=r'$\Phi$',ylabel=r'$\int P(\Phi)$',xscale='linear',yscale='log')
    fig2.savefig('plots_to_sort/%s_cuml_phi_n%04d.pdf'%(this_simname,frame))
