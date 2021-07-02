
"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import data_locations as dl
import xtra_energy
reload(loop_apps)
from scipy.optimize import curve_fit

import testing.early_mask as em
reload(em)
#fit powerlaw for ratio
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)


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

if 'loop_dict' not in dir():
    if 0:
        import three_loopers_1tff as tl
        loop_dict = {'u201':tl.looper1, 'u202':tl.looper2, 'u203':tl.looper3}
    else:
        loop_dict = {}
        core_set={'u201':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/all_cores_n0000.h5',
                  'u202':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u10_primitives_cXXXX_n0000.h5',
                  'u203':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u11_primitives_cXXXX_n0000.h5'}
    for this_simname in  ['u201','u202','u203']:
        directory = dl.sims[this_simname]
        save_field = core_set[this_simname]
        print("LOAD", this_simname)
        loop_dict[this_simname] = looper.core_looper(directory= directory,savefile=save_field)




vrms = {'u201':5.2, 'u202':5.1, 'u203':5.4}

for this_simname in ['u201','u202','u203']:
    #if this_simname != 'u201' and skipper==True:
    #    continue

    this_looper = loop_dict[this_simname]


    #FIELD = 'density'
    #FIELD = 'velocity_magnitude'
    FIELD = 'magnetic_field_strength'
    #FIELD = 'PotentialField'
    eta1 = 0.040253639221191406 # = number of core particles / total
        #for core in this_looper.target_indices:
        #    n+=this_looper.target_indices[core].size
        #eta1 = n/128**3
    frame=0
    deposit_tuple = ("deposit","target_particle_volume")
    xscale = {'PotentialField':'linear'} 

    prof_dir = "."
    version = 't1_'
    version = ''
    prof_full_fname = "%s/%spreimage_pdf_full_%s_%s_n%04d.h5"%(prof_dir, version,this_simname, FIELD, frame)
    prof_part_fname = "%s/%spreimage_pdf_part_%s_%s_n%04d.h5"%(prof_dir, version,this_simname, FIELD, frame)
    print(prof_full_fname)

    if os.path.exists(prof_full_fname):
        bbb1, bcen1, vals1, db1 = dpy(prof_full_fname, ['bin_edge','bin_center','vals','db'])
        bbb2, bcen2, vals2, db2 = dpy(prof_part_fname, ['bin_edge','bin_center','vals','db'])
    else:
        if 1:
            ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
            em.add_tracer_density(ds)
            ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
        #ad[deposit_tuple]
        all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in this_looper.target_indices])
        ad.set_field_parameter('target_indices',all_target_indices)
        ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
        #bins={'velocity_x':np.linspace(-25,25,64)}
        #bins['PotentialField']= np.linspace(-50,50,64)
        bins = {'PotentialField':np.linspace(-32,32,64),'density':None,'magnetic_field_strength':None,'velocity_magnitude':None}
        prof_all_density  = yt.create_profile(ad,bin_fields=[FIELD],fields=['cell_volume'],weight_field=None, override_bins=bins)
        prof_mask_density = yt.create_profile(ad,bin_fields=[FIELD],fields=[deposit_tuple],weight_field=None, override_bins=bins)

        bbb1, bcen1, vals1, db1= toplot(prof_all_density)
        bbb2, bcen2, vals2, db2 = toplot(prof_mask_density,quan=deposit_tuple[1])

        fptr1 = h5py.File(prof_full_fname,'w')
        fptr1.create_dataset( 'bin_edge', data=bbb1)
        fptr1.create_dataset( 'bin_center', data=bcen1)
        fptr1.create_dataset( 'vals', data=vals1)
        fptr1.create_dataset( 'db', data=db1)
        fptr1.close()

        fptr2 = h5py.File(prof_part_fname,'w')
        fptr2.create_dataset( 'bin_edge', data=bbb2)
        fptr2.create_dataset( 'bin_center', data=bcen2)
        fptr2.create_dataset( 'vals', data=vals2)
        fptr2.create_dataset( 'db', data=db2)
        fptr2.close()




    #
    # Get rid of zeros in PDF.
    #
    ok1 = vals1 > 0
    ok2 = vals2 > 0
    ok_both=ok1*ok2

    field_latex = {'density':r"\rho", 'velocity_magnitude':'v','magnetic_field_strength':'B','PotentialField':r'\Phi'}[FIELD]
    PDF_LABEL = r'$V(%s)$'%(field_latex)
    PDF_LABEL_C = r'$V(%s|*)$'%(field_latex)
    PDF_LABEL_R = r'$V(*|%s)$'%(field_latex)


    if 1:
        fig,ax=plt.subplots(1,1)
        ax.plot( bcen1[ok1],vals1[ok1],'k',linewidth=2, label=PDF_LABEL)
        ax.plot( bcen2[ok2],vals2[ok2],'k--',linewidth=2, label=PDF_LABEL_C)

    if 1:
        ratio = vals2/vals1
        ok = bcen1>0.1
        ok = np.logical_and(ratio>0, ok)
        ax.plot( bcen1[ratio>0], ratio[ratio>0],label=PDF_LABEL_R,c=[0.5]*4)


    ks_string = ''
    def do_ks(test_dist, original_dist):
        KS_output = ks_2samp(test_dist, original_dist)
        crit_stat = 1.36*np.sqrt( (original_dist.size + test_dist.size)/(original_dist.size*test_dist.size))
        ks_string = r'$D=%0.2f (D_c=%0.2f, p=%0.2f)$'%(KS_output.statistic,crit_stat, KS_output.pvalue) 
        return ks_string

    if FIELD == 'density':
        popt, pcov = curve_fit(powerlaw, bcen1[ok], np.log10(ratio[ok]), p0=[1,1,-2])

        p_star_given_rho = vals1*bcen1**popt[2]
        p_star_given_rho *= eta1
        ok_p = ok_both
        ks_string = do_ks( vals2[ok_p], p_star_given_rho[ok_p])
        lab = r'$\eta_1\rho^{%0.2f}V(\rho)$: %s'%(popt[2], ks_string)
        ax.plot( bcen1,p_star_given_rho,linewidth=2, label=lab, linestyle='--',c=[0.5]*4)

    if FIELD == 'velocity_magnitude':
        #maxwellian
        v = bcen1
        sigmav = vrms[this_simname]
        maxwell = 1./(2*np.pi*sigmav**2)*v**2*np.exp(-(v-0)**2/(2*sigmav**2))
        maxwell *= eta1
        ok_p = (maxwell > 1e-7)*ok2

        ks_string = do_ks( vals2[ok_p], maxwell[ok_p])

        lab = r'$\eta_1 \mathcal{M}(\sigma_{1d}=%0.1f)$: %s'%(sigmav, ks_string)


        ax.plot( v[ok_p], maxwell[ok_p] , label=lab, linestyle='--',c=[0.5]*4)
        p_star_given_rho = maxwell



    if FIELD == 'magnetic_field_strength':
        p_star_given_rho = vals1#*bcen1**0.5
        p_star_given_rho *= eta1

        ks_string = do_ks( vals2[ok_both], p_star_given_rho[ok_both])
        ax.plot( bcen1,p_star_given_rho,linewidth=2, label=r'$\eta_1 V(B)$ %s'%ks_string, linestyle='--',c=[0.5]*4)
        ok_p = ok_both




    if 1:
        #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
        outname = "plots_to_sort/%s_pdf_%s_preimage_fits.pdf"%(this_simname,FIELD)
        #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
        #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='linear',yscale='linear')
        axbonk(ax,xlabel=r'$%s$'%field_latex,ylabel=r'$V(%s)$'%field_latex,xscale=xscale.get(FIELD,'log'),yscale='log')
        ax.legend(loc=3)
        fig.savefig(outname)
        print(outname)
        plt.close(fig)


        
    if 0:
        # Do the KS test
        # https://abhyankar-ameya.medium.com/kolmogorov-smirnov-two-sample-test-with-python-70c309107c78
        #
        fig2, ax2=plt.subplots(1,1)
        print(KS_output)
        cuml_all  = np.cumsum(original_dist)
        cuml_mask = np.cumsum(test_dist)
        ax2.plot( bcen1[ok_p], cuml_all/cuml_all[-1], c='k')
        ax2.plot( bcen2[ok_p], cuml_mask/cuml_mask[-1], 'k--')
        #ax2.plot( bcen1,vals1,c='k')
        #ax2.plot( bcen2,vals2,'k--')
        axbonk(ax2,xlabel=r'$\rho$',ylabel=r'$\int V(rho)$',xscale='log',yscale='log')
        outname = 'plots_to_sort/%s_cuml_%s_n%04d.pdf'%(this_simname,FIELD,frame)
        fig2.savefig(outname)
        print(outname)
