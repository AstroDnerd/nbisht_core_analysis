
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

if 'this_simname' not in dir():
    this_simname = 'u11'


frame_list=[0]# range(0,130,10)
core_list=looper.get_all_nonzero(dl.n_particles[this_simname])
fields=['density']  
#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
save_field = '../Datasets/all_cores_n0000.h5'
save_field = '../Datasets/u10_primitives_cXXXX_n0000.h5'
save_field = '../Datasets/u11_primitives_cXXXX_n0000.h5'
if 'this_looper' not in dir() and os.path.exists(save_field) and True:
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
        
if 'prof_mask_vel' not in dir():
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    #bins={'velocity_x':np.linspace(-25,25,64)}
    #bins['PotentialField']= np.linspace(-50,50,64)
    bins=None
    prof_all_vel  = yt.create_profile(ad,bin_fields=['velocity_magnitude'],fields=['cell_volume'],weight_field=None, override_bins=bins)
    prof_mask_vel = yt.create_profile(ad,bin_fields=['velocity_magnitude'],fields=[deposit_tuple],weight_field=None, override_bins=bins)

if 'vrms' not in dir():
    ad = ds.all_data()
    vx = ad['velocity_x']
    vy = ad['velocity_y']
    vz = ad['velocity_z']
    vrms = np.sqrt( (vx**2+vy**2+vz**2).mean() )/np.sqrt(3)

if 1:
    plt.close('all')
    bbb1, bcen1, vals1, db= toplot(prof_all_vel)
    bbb2, bcen2, vals2, db = toplot(prof_mask_vel,quan=deposit_tuple[1])

if 1:
    fig2, ax2=plt.subplots(1,1)
    cuml_all  = np.cumsum(vals1)
    cuml_mask = np.cumsum(vals2)
    ax2.plot( bcen1, cuml_all/cuml_all[-1], c='k')
    ax2.plot( bcen2, cuml_mask/cuml_mask[-1], 'k--')
    #ax2.plot( bcen1,vals1,c='k')
    #ax2.plot( bcen2,vals2,'k--')
    axbonk(ax2,xlabel=r'$\rho$',ylabel=r'$\int V(rho)$',xscale='log',yscale='log')
    outname = 'plots_to_sort/%s_cuml_v_n%04d.pdf'%(this_simname,frame)
    fig2.savefig(outname)
    print(outname)

if 1:
    fig,ax=plt.subplots(1,1)

    ok1 = vals1>0
    ok2 = vals2>0
    ax.plot( bcen1[ok1],vals1[ok1],'k',linewidth=2, label=r'$V(v)$')


    ax.plot( bcen2[ok2],vals2[ok2],'k--',linewidth=2, label=r'$V(v|*)$')

if 1:
    #maxwellian
    particle_fraction_scale= all_target_indices.size/128.**3
    v = bcen1[ok2]
    sigmav = vrms
    maxwell = 1./(2*np.pi*sigmav**2)*v**2*np.exp(-(v-0)**2/(2*sigmav**2))
    maxwell *= particle_fraction_scale
    ok3 = maxwell > 1e-7
    ax.plot( v[ok3], maxwell[ok3] , label=r'$\sigma_{v,1d}=%0.1f$'%vrms, linestyle='--',c=[0.5]*4)


if 1:
    #comput and plot the ratio
    ratio = vals2/vals1
    ok = bcen1>0.1
    ok = np.logical_and(ratio>0, ok)
    ax.plot( bcen1[ratio>0], ratio[ratio>0],label=r"$V(*|v)$",c=[0.5]*4)

if 1:

    if 'prof_mask_veln' not in dir():
        prof_mask_veln = {}
        adn={}
if 0:   
    for n in [10,20,30,40,50,60,70,80,90,100]:
        if n not in prof_mask_veln:
            if n not in adn:
                ds = this_looper.load(frame=n,derived=[em.add_tracer_density])
                em.add_tracer_density(ds)
            if n not in adn:
                adn[n] = ds.all_data()
                adn[n].set_field_parameter('target_indices',all_target_indices)
                adn[n].set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
            this_ad = adn.get(n,ds.all_data())

            prof_mask_veln[n] = yt.create_profile(this_ad,bin_fields=['velocity_magnitude'],fields=[deposit_tuple],\
                                                  weight_field=None, override_bins=bins)
    for n in prof_mask_veln:
        bbbx, bcenx, valsx, dbx= toplot(prof_mask_veln[n],quan=deposit_tuple[1])
        dx = bbbx[1:]-bbbx[:-1]
        total = (dx*valsx).sum()
        print("TOTAL",total)
        ax.plot( bcenx[ok1],valsx[ok1],linewidth=2, label=r'$V(v)$')


if 1:
    #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
    outname = "plots_to_sort/%s_pdf_velocity_preimage_fits.pdf"%this_simname
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
    #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='linear',yscale='linear')
    axbonk(ax,xlabel=r'$||v||$',ylabel=r'$V(||v||)$',xscale='log',yscale='log')
    ax.legend(loc=3)
    fig.savefig(outname)
    print(outname)
    plt.close(fig)
    
