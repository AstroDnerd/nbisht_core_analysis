
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
np.set_printoptions(threshold=sys.maxsize)
# --- --- --- --- --- --- ---

class kapparadial(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.bad_cores = []

    def qtyRun(self, sim, core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # FIELDS STORED
        self.columnrho_sph = [np.zeros(len(core_list)) for x in range(6)]
        self.blos_sph = [np.zeros(len(core_list)) for x in range(6)]

        # THE FRAMES
        the_frame = thtr.frames[-1:] 
        ds = self.this_looper.load(the_frame[0]) 
        radius=1e-2
        radius = ds.arr(radius,'code_length')

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)

            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)

            p49_fields.add_QU(ds)
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_radius = (1/128) * 3.5 
            the_area = np.pi * (the_radius**2) 
     
            xyz = [0,1,2]
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']

            # MAKE ONE SPHERE ONLY
            the_sphere = ds.sphere(the_center, the_radius) 
            j=0

            # HM...ideas from alpha_r.py
            rrr = the_sphere['radius']
            dv = the_sphere['cell_volume']
            B  = the_sphere[B[j]]

            args=np.argsort(rrr)
            rho = sph['density'][args]
            rrrs = rrr[args]
            dvs  = dv[args]
            Bs   = B[args]

            # means
            column_rho = np.cumsum(rho*dvs)/np.cumsum(dvs)
            B_los = np.cumsum(Bs * rho * dvs)/np.cumsum(dvs)

            one = 1/128
            rbins = [int(args.size/(0.5*one)), int(args.size/one), int(args.size/(1.5*one)), int(args.size/(2*one)), int(args.size/(2.5*one)), int(args.size/(3*one))]
            for here in [rbins]:
                self.columnrho_sph[nc][rbins] = column_rho[:here] 
                self.blos_sph[nc][rbins] = B_los[:here]

            '''
            for num in range(len(rbins)):  #one direction for now 
                self.columnrho_sph[nc][rbins] = (the_sphere['density'] * the_sphere['cell_volume']).sum()/the_area
                self.blos_sph[nc][rbins] = (the_sphere['density'] * the_sphere[B[j]] * the_sphere['cell_volume']).sum()/the_sphere['gas','cell_mass'].sum()
            '''

        # THEN
        kappas = []
        for num in range(6)
            # get kappa per radius size
            fit = np.polyfit(np.log(self.columnrho_sph, np.log(self.blos_sph),1)  
            kappas.append(fit[0])

        fig,ax=plt.subplots(1,1)
        ax.scatter(rbins,kappas, c='r')
        fig.savefig('kappas_r_%s'%(sim))



        print("check-in")
        # CACHE TO DISK 
        if 0: 
            hfivename = 'p66_brho/kappa_r_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')

            Fptr['core_id']=self.cores_used 
            Fptr['B_los_sph']=self.blos_sph
            Fptr['N_sph']=self.columnrho_sph
            Fptr.close()


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=kapparadial(TL6.loops['u601'])
    #core_list1 = TL6.loops['u601'].core_by_mode['Alone']
if 'scope2' not in dir() or clobber:
    scope2=kapparadial(TL6.loops['u602'])
    #core_list2 = TL6.loops['u602'].core_by_mode['Alone']
if 'scope3' not in dir() or clobber:
    scope3=kapparadial(TL6.loops['u603'])
    #core_list3 = TL6.loops['u603'].core_by_mode['Alone']

simnames = ['u601','u602', 'u603']


for nt,tool in enumerate([scope1,scope2,scope3]):

    # WHICH CORES & RUN IF NOT CACHED
    if 0:
        all_cores = np.unique(tool.this_looper.tr.core_ids)
        core_list = all_cores
        #core_list = all_cores[:1]  #DEBUG
        #core_list = core_list1[:3] 
        tool.qtyRun(nt,core_list=core_list) 

    # ONCE CACHED, HISTOGRAMS SCATTERS & REL ERROR REGARDING POLARIZATION 
    if 1: 
        hfivename = 'p66_brho/h5files/b_sphere_%s.h5'%(nt)  #EDIT
        Fptr = h5py.File(hfivename,'r')

        cores_used = Fptr['core_id']

        b_los_sph = Fptr['B_los_sph']
        b_LOS_sph = np.concatenate((b_los_sph[0],b_los_sph[1],b_los_sph[2]))
        b_lossphlog = np.log10(abs(b_LOS_sph))

        n_sph = Fptr['N_sph']
        N_sph = np.concatenate((n_sph[0],n_sph[1],n_sph[2]))
        n_sphlog = np.log10(N_sph)

        b_tot_sph = Fptr['B_tot_sph']
        b_TOT_sph = np.concatenate((b_tot_sph[0],b_tot_sph[1],b_tot_sph[2]))
        b_totsphlog = np.log10(abs(b_TOT_sph))
        
        rho_sph = Fptr['Rho_avg_sph']
        Rho_sph = np.concatenate((rho_sph[0],rho_sph[1],rho_sph[2]))
        rho_sphlog = np.log10(Rho_sph)
        


        the_bins = 128  #default
        color=['b','g','orange']
        # ONE DIRECTION AT A TIME 
        for i in range(3):
            # 3MT PLOTS 
            if 0: 
                the_bins = 16 
                #the_weights = np.ones(len(np.delete(relerr_UQfrac[i],index)))/len(np.delete(relerr_UQfrac[i],index))
                pdf, bins = np.histogram(np.delete(relerr_UQfrac[i],index), bins=the_bins, density=True)#, weights =the_weights)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
                #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                plt.xlabel(r'$P_{error: x, y, z}$')
                plt.ylabel('PDF')
                outname = 'UQfrac_relerr_PDF_%d_%s'%(i,simnames[nt])
                #outname = 'UQfrac_relerr_%s'%(simnames[nt])

            if 0: 
                plt.scatter(n_rho_mcyl_one[i],abs(b_dcf_mcyl_one[i]),c='orange',alpha=0.4)
                plt.scatter(n_rho_mcyl_one[i],abs(b_pos_mcyl_one[i]),c='k',alpha=0.4)
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel(r'$N$ $(cm^{-2})$')
                plt.ylabel(r'$B_{dcf}$, $B_{pos}$ $(\mu G)$')
                outname ='B_dcfpos_vsN_%d_%s'%(i,simnames[nt])

            # SAVE ONE DIRECTION AT A TIME
            if 0:
                plt.savefig(outname)
                print('plotted_%d'%i)
                plt.close('all')


        # ALL DIRECTIONS AT ONCE  
        if 1: 
            if 1:
                pfit = np.polyfit(rho_sphlog,b_totsphlog, 1)
                alpha = pfit[0]
                b_poso = pfit[1]
                n_Rho = np.linspace(rho_sphlog.min(),rho_sphlog.max(),num=len(rho_sphlog))
                N_Rho = 10 ** n_Rho
                B_two = 10 ** (alpha*n_Rho + b_poso)
                plt.plot(N_Rho,B_two,color='k',linestyle='dashed')
            if 1:
                pearX,pearY = scipy.stats.pearsonr(Rho_sph,abs(b_TOT_sph))
            if 0:
                pfit2 = np.polyfit(n_sphlog,b_lossphlog, 1)
                alpha2 = pfit2[0]
                b_poso2 = pfit2[1]
                n_Rho2 = np.linspace(n_sphlog.min(),n_sphlog.max(),num=len(n_sphlog))
                N_Rho2 = 10 ** n_Rho2
                B_two2 = 10 ** (alpha2*n_Rho2 + b_poso2)
                plt.plot(N_Rho2,B_two2,color='g',linestyle='dashed')
            if 0:
                pearX2,pearY2 = scipy.stats.pearsonr(N_sph,abs(b_LOS_sph))

            #plt.scatter(n_RHO_mcyl_one, abs(b_LOS_mcyl_one),c = 'k', alpha=0.4)
            plt.scatter(Rho_sph, abs(b_TOT_sph),c = 'k', alpha=0.4)
            #plt.scatter(rho_AVE, abs(b_PAR),c = 'k', alpha=0.4)

            #plt.scatter(abs(b_DCF_cyl), abs(b_DCF_midcyl), c = 'orange', alpha=0.4)
            #plt.scatter(abs(b_TOT_midcyl)[ok_los], abs(b_DCF_midcyl)[ok_los], c = 'orange', alpha=0.4)
            #plt.axline((0, 0), slope=1, c='k', linewidth=0.5)

            plt.xscale('log')
            plt.yscale('log')
            #plt.xlabel(r'$N_{sph,los}$ $(cm^{-2})$')
            plt.xlabel(r'$rho_{ave}$ $(cm^{-3})$')
            plt.ylabel(r'$B_{sph,tot}$ $(\mu G)$')
            #plt.ylabel(r'$B_{particles}$ $(\mu G)$')
            #plt.title(r'$\kappa_{1,2} = %f,%f$ $R_{1,2}=%f,%f$'%(alpha,alpha2,pearX,pearX2))
            plt.title(r'$\kappa= %f$ $R = %f$'%(alpha,pearX))
            outname ='b_sphtot_nook_%s'%(simnames[nt])

        #color=['g','orange']  #b_tor blue, b_pol orange
        #blosm1_frac = abs(b_LOS_mcyl_one)/abs(b_TOT_mcyl_one) 
        #blosm2_frac = abs(b_LOS_cyl)/abs(b_TOT_cyl)
        if 0: 
            the_bins = 64 #16 
            pdf1, bins1 = np.histogram(blosm1_frac, bins=the_bins, density=True) 
            pdf2, bins2 = np.histogram(blosm2_frac, bins=the_bins, density=True)
            #pdf3, bins3 = np.histogram(dTHETA_cyl, bins=the_bins, density=True)
            cdf1= np.cumsum(pdf1)
            cdf2= np.cumsum(pdf2)
            #cdf3= np.cumsum(pdf3)
            bin_centers1 = 0.5*(bins1[1:]+bins1[:-1]) 
            bin_centers2 = 0.5*(bins2[1:]+bins2[:-1]) 
            #bin_centers3 = 0.5*(bins3[1:]+bins3[:-1]) 
            plt.plot(bin_centers1,pdf1,alpha=0.7,color='k')
            plt.plot(bin_centers2,pdf2,alpha=0.7,color='g')
            #plt.plot(bin_centers3,pdf3,alpha=0.7,color='k')
            plt.xlabel(r'$B_{los}/B{tot}$ mcyl1(black), cyl(green)')
            plt.ylabel(r'PDF')
            outname ='B_losVStotfrac_mcyl1_cyl_%s'%(simnames[nt])

        # SAVING ALL AT ONCE
        if 1:
            plt.savefig(outname)
            print('plotted_%d'%i)
            plt.close('all')

        # AND CLOSE H5 FILE
        Fptr.close()


