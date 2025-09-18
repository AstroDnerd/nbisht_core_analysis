
from starter2 import *
import data_locations as dl
import xtra_energy
from scipy.optimize import curve_fit

import testing.early_mask as em
reload(em)
#fit powerlaw for ratio
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp

import mean_field
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)


def toplot(prof,quan = YT_cell_volume):
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
        #take 1, not enough particles.
        import three_loopers_1tff as tl
        loop_dict = {'u201':tl.looper1, 'u202':tl.looper2, 'u203':tl.looper3}
    elif 0:
        #take 1.5, a slightly different core set
        loop_dict = {}
        core_list={}
        core_set={'u05':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/all_cores_n0000.h5',
                  'u10':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u10_primitives_cXXXX_n0000.h5',
                  'u11':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/u11_primitives_cXXXX_n0000.h5'}
        for this_simname in ['u05','u10','u11']:
            directory = dl.sims[this_simname]
            save_field = core_set[this_simname]
            print("LOAD", this_simname)
            loop_dict[this_simname] = looper.core_looper(directory= directory,savefile=save_field)
            core_list[this_simname] =looper.get_all_nonzero(dl.n_particles[this_simname])
    elif 0:
        import three_loopers_mountain_top as TLM
        reload(TLM)
        loop_dict = TLM.loops
    elif 0:
        import three_loopers_tenfour as TL4
        loop_dict = TL4.loops
    elif 0:
        import three_loopers_six as TL6
        loop_dict = TL6.loops
    else:
        import track_loader as TL
        loop_dict = TL.loops




vrms = {'u201':5.2, 'u202':5.1, 'u203':5.4}
vrms.update( {'u301':5.2, 'u302':5.1, 'u303':5.4} )
vrms.update( {'u401':5.2, 'u402':5.1, 'u403':5.4} )
vrms.update( {'u601':5.2, 'u602':5.1, 'u603':5.4} )
vrms.update( {'u501':5.2, 'u502':5.1, 'u503':5.4} )
vrms.update( {'u05':5.2, 'u10':5.1, 'u11':5.4})
sims_to_use = ['u301', 'u302','u303']
sims_to_use = ['u401', 'u402','u403']
sims_to_use = ['u601', 'u602','u603']
sims_to_use = ['u501', 'u502','u503']
#sims_to_use = ['u05','u10','u11']
#sims_to_use = ['u201', 'u202','u203']

#vrms['b002']=vrms['u202']
#sims_to_use=['b002']

fig, axes = plt.subplots(4,3, figsize=(10,12))
for nsim,this_simname in enumerate(sims_to_use):
    #if this_simname != 'u11':# and skipper==True:
    #    continue

    TL.load_tracks([this_simname])

    this_looper = loop_dict[this_simname]
    for fnum in range(4):
        FIELD = [YT_density, YT_velocity_magnitude, YT_magnetic_field_strength, YT_potential_field][fnum]
        labelset  = ['abc','def','ghi','jkl'][fnum]
        label=labelset[nsim]
        ax = axes[fnum][nsim]


        #FIELD = YT_density; label = 'abc'[nsim]
        #FIELD = YT_velocity_magnitude; label='def'[nsim]
        #FIELD = YT_magnetic_field_strength; label='ghi'[nsim]
        #FIELD = YT_potential_field; label='jkl'[nsim]
        #FIELD = YT_grav_energy_2; label='jkl'[nsim]

        #core_list = [10,32,84]
        core_list = np.unique(this_looper.tr.core_ids)
        if 0:
            #looper 1 version
            #if the next line fails, use this one.
            all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
        if 1:
            #the looper 2 version.  If this fails, you have looper 1.
            all_target_indices = this_looper.target_indices
        all_target_indices = all_target_indices.astype('int64')
        eta1 = len(all_target_indices)/128**3
        print("ETA 1 %s %0.4f"%(this_simname, eta1))
            #for core in this_looper.target_indices:
            #    n+=this_looper.target_indices[core].size
            #eta1 = n/128**3
        frame=0
        deposit_tuple = ("deposit","target_particle_volume")
        xscale = {'PotentialField':'linear'} 

        prof_dir = "."
        version = 't2_'
        version = ''
        prof_full_fname = "%s/%spreimage_pdf_full_%s_%s_n%04d.h5"%(prof_dir, version,this_simname, FIELD[1], frame)
        prof_part_fname = "%s/%spreimage_pdf_part_%s_%s_n%04d.h5"%(prof_dir, version,this_simname, FIELD[1], frame)

        if os.path.exists(prof_full_fname):
            print("READ FROM DISK")
            bbb1, bcen1, vals1, db1 = dpy(prof_full_fname, ['bin_edge','bin_center','vals','db'])
            bbb2, bcen2, vals2, db2 = dpy(prof_part_fname, ['bin_edge','bin_center','vals','db'])
        else:
            print("MAKE NEW PROFILES")
            if 1:
                ds = this_looper.load(frame=frame,derived=[em.add_tracer_density])
                em.add_tracer_density(ds)
                xtra_energy.add_energies(ds)
                ad = ds.all_data() #ds.region([0.5]*3,[0.4]*3,[0.6]*3)
            #all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list[this_simname]])
            ad.set_field_parameter('target_indices',all_target_indices)
            ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
            #dep = ad[deposit_tuple] make sure you uncomment the mask_to_get flag, too.
            #ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
            #bins={'velocity_x':np.linspace(-25,25,64)}
            #bins['PotentialField']= np.linspace(-50,50,64)
            bins = {'PotentialField':np.linspace(-32,32,65),
                    'density':np.logspace(np.log10(5e-3),2,65),
                    'magnetic_field_strength':np.logspace( np.log10( 0.4), np.log10(110), 65),
                    'velocity_magnitude':np.logspace(np.log10(0.03356773506543828), 
                                                     np.log10(30.661519625727557),65),
                    'grav_energy_2':np.geomspace(0.01,50,64)}


            #density1 = np.abs(ad[FIELD[1]])/ad[YT_kinetic_energy]
            density1 = np.abs(ad[FIELD[1]])#/ad[YT_kinetic_energy]
            #bins['grav_energy_2']=np.geomspace(density1.min(),density1.max(),64)
            bbb1 = bins[FIELD[1]]
            cell_volume1 = ad[YT_cell_volume]
            vals1, bbb1 = np.histogram(density1, weights=cell_volume1, bins=bbb1)
            bcen1=0.5*(bbb1[1:]+bbb1[:-1])
            db1 = bbb1[1:]-bbb1[:-1]
            vals1/=db1

            frame_ind = np.where(this_looper.tr.frames == frame)[0][0]
            #density2 = this_looper.tr.track_dict[FIELD[1]][:,frame_ind]
            #cell_volume2 = this_looper.tr.track_dict[YT_cell_volume][:,frame_ind]
            if FIELD[1] == 'grav_energy_2':
                gx     = this_looper.tr.c(core_list, 'grav_x')[:,frame_ind]
                gy     = this_looper.tr.c(core_list, 'grav_y')[:,frame_ind]
                gz     = this_looper.tr.c(core_list, 'grav_z')[:,frame_ind]
                #vx     = this_looper.tr.c(core_list, 'velocity_x')[:,frame_ind]
                #vy     = this_looper.tr.c(core_list, 'velocity_y')[:,frame_ind]
                #vz     = this_looper.tr.c(core_list, 'velocity_z')[:,frame_ind]
                #de     = this_looper.tr.c(core_list, 'density')[:,frame_ind]
                density2 = 1/(np.pi*8*colors.G)*(gx*gx+gy*gy+gz*gz)
                #density2/=0.5*de*(vx*vx+vy*vy+vz*vz)

            else:
                density2     = this_looper.tr.c(core_list, FIELD[1])[:,frame_ind]


            cell_volume2 = this_looper.tr.c(core_list, 'cell_volume')[:,frame_ind]
            vals2, bbb2 = np.histogram(density2, weights=cell_volume2, bins=bbb1)
            bcen2=0.5*(bbb2[1:]+bbb2[:-1])
            db2 = bbb2[1:]-bbb2[:-1]
            vals2 /= db2


        if FIELD[1] == 'magnetic_field_strength':
            bcen1 /= mean_field.B[this_simname]
            bcen2 /= mean_field.B[this_simname]
        #
        # Get rid of zeros in PDF.
        #
        ok1 = vals1 > 0
        ok2 = vals2 > 0
        ok_both=ok1*ok2

        field_latex = {'density':r"\rho/\rho_0", 'velocity_magnitude':'v/c_s','magnetic_field_strength':'B/B0','PotentialField':r'\Phi', 'grav_energy_2':"EG"}[FIELD[1]]
        PDF_LABEL = r'$V(%s)$'%(field_latex)
        PDF_LABEL_C = r'$V(%s|*)V(*)$'%(field_latex)
        PDF_LABEL_R = r'$V(*|%s)$'%(field_latex)


        if 1:
            #fig,ax=plt.subplots(1,1)
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

        if FIELD[1]== 'density':

            popt, pcov = curve_fit(powerlaw, bcen1[ok], np.log10(ratio[ok]), p0=[1,1,-2])

            #p_star_given_rho = vals1*bcen1**popt[2]
            #p_star_given_rho *= eta1
            rhomax = bcen1.max()
            p_star_given_rho = (bcen1/rhomax)**(popt[2])*vals1
            ok_p = ok_both
            ks_string = do_ks( np.log(vals2[ok_p])[::-1], np.log(p_star_given_rho[ok_p])[::-1])
            lab = r'$V(*) \rho^{%0.2f}V(\rho)$: %s'%(popt[2], ks_string)
            print('=========')
            print("DENSITY %s %s"%(this_simname,ks_string))
            ax.plot( bcen1,p_star_given_rho,linewidth=2, label=lab, linestyle='--',c=[0.5]*4)

        if FIELD[1] == 'velocity_magnitude':
            #maxwellian
            v = bcen1
            sigmav = vrms[this_simname]
            maxwell = 1./(2*np.pi*sigmav**2)*v**2*np.exp(-(v-0)**2/(2*sigmav**2))
            maxwell *= eta1
            ok_p = (maxwell > 1e-7)*ok2

            ks_string = do_ks( vals2[ok_p], maxwell[ok_p])
            print('=========')
            print("SPEED %s %s"%(this_simname,ks_string))

            lab = r'$V(*) \mathcal{M}(v;\sigma_{1d}=%0.1f)$: %s'%(sigmav, ks_string)


            ax.plot( v[ok_p], maxwell[ok_p] , label=lab, linestyle='--',c=[0.5]*4)
            p_star_given_rho = maxwell



        if FIELD[1] == 'magnetic_field_strength':
            p_star_given_rho = vals1#*bcen1**0.5
            p_star_given_rho *= eta1

            ks_string = do_ks( vals2[ok_both], p_star_given_rho[ok_both])
            ax.plot( bcen1,p_star_given_rho,linewidth=2, label=r'$V(*) V(B)$ %s'%ks_string, linestyle='--',c=[0.5]*4)
            ok_p = ok_both


            print('=========')
            print("FIELD %s %s"%(this_simname,ks_string))


        if 1:


            #ax.plot( bcen2,vals2*vals1.max()/vals2.max(),'r:')
            outname = "plots_to_sort/%s_pdf_%s_preimage_fits.pdf"%(this_simname,FIELD[1])
            #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='log',yscale='log')
            #axbonk(ax,xlabel=r'$\rho$',ylabel='V(rho)',xscale='linear',yscale='linear')
            
            axbonk(ax,xlabel=r'$%s$'%field_latex,ylabel=r'$V(%s)$'%field_latex,xscale=xscale.get(FIELD[1],'log'),yscale='log')

            if 1:
                #figure sublabel
                lab = r'$3%s$'%label
                #xloc = ax.get_xlim()[1]*0.5
                #yloc = ax.get_ylim()[0]*1.2
                xloc = 0.5
                yloc = 0.05
                ax.text(xloc,yloc,lab, transform=ax.transAxes)
            #ax.legend(loc=3)
            #fig.savefig(outname)
            print(outname)
            plt.close(fig)

for fnum in range(4):
    ext = extents()
    for nsim,this_simname in enumerate(sims_to_use):
        ax = axes[fnum][nsim]
        ext(nar(ax.get_ylim()))
    for nsim,this_simname in enumerate(sims_to_use):
        ax = axes[fnum][nsim]
        ax.set( ylim=ext.minmax)
        if nsim:
            ax.set(ylabel='',yticks=[])
fig.subplots_adjust(wspace=0, hspace=0)
fig.tight_layout()
fig.savefig('plots_to_sort/Figure3.pdf')


            
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
    outname = 'plots_to_sort/%s_cuml_%s_n%04d.pdf'%(this_simname,FIELD[1],frame)
    fig2.savefig(outname)
    print(outname)
