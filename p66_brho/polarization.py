
'''
polarization into Enzo
'''

from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
# --- --- --- --- --- --- ---

class polarization(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.bad_cores = []

    def qtyRun(self,sim,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # FIELDS STORED
        self.Q_long = [np.zeros(len(core_list)) for x in range(3)]
        self.Q_short = [np.zeros(len(core_list)) for x in range(3)]
        self.U_long = [np.zeros(len(core_list)) for x in range(3)]
        self.U_short = [np.zeros(len(core_list)) for x in range(3)]
        self.Q_long = [np.zeros(len(core_list)) for x in range(3)]

        self.theta_long = [np.zeros(len(core_list)) for x in range(3)]
        self.UQfrac_long = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho = [np.zeros(len(core_list)) for x in range(3)]

        self.theta_short = [np.zeros(len(core_list)) for x in range(3)]
        self.UQfrac_short = [np.zeros(len(core_list)) for x in range(3)]
        self.columnrho_short = [np.zeros(len(core_list)) for x in range(3)]


        # THE FRAMES
        the_frame = thtr.frames[-1:] 

        # CORE-LOOP
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)

            ds = self.this_looper.load(the_frame[0]) 
            #pdb.set_trace()  #TAKE OFF
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)

            p49_fields.add_QU(ds)
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]
            the_radius = 1/128
            the_area= np.pi * (the_radius**2) 
            
            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            the_cyl = {}
            the_mid_cyl = {}
            for i in range(3):
                the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))
                the_mid_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=the_radius) 
            
            Q = ['Qx','Qy','Qz']
            U = ['Ux','Uy','Uz']
            for j in range(3):  #REVISE!!
                self.Q_long[j][nc] = (the_cyl[j][Q[j]] * the_cyl[j]['cell_volume']).sum()
                self.Q_short[j][nc] = (the_mid_cyl[j][Q[j]] * the_mid_cyl[j]['cell_volume']).sum()

                self.U_long[j][nc] = (the_cyl[j][U[j]] * the_cyl[j]['cell_volume']).sum()
                self.U_short[j][nc] = (the_mid_cyl[j][U[j]] * the_mid_cyl[j]['cell_volume']).sum()

                Q_long =self.Q_long[j][nc] 
                U_long =self.U_short[j][nc]
                Q_short =self.Q_short[j][nc] 
                U_short =self.U_short[j][nc]

                self.theta_long[j][nc] = np.arctan2(U_long,Q_long) * 180 / np.pi  
                self.theta_short[j][nc] = np.arctan2(U_short,Q_short) * 180 / np.pi  

                self.columnrho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()
                self.columnrho_short[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()

                N_long = self.columnrho[j][nc]  
                N_short = self.columnrho_short[j][nc]  
                self.UQfrac_long[j][nc] = np.sqrt(Q_long**2 + U_long**2)/N_long 
                self.UQfrac_short[j][nc] = np.sqrt(Q_short**2 + U_short**2)/N_short 

                if math.isnan(self.UQfrac_long[j][nc]) == True:
                    print('a nan value')
                    self.bad_cores.append(core_id)
                    self.UQfrac_long[j][nc] = 0.0
                    self.UQfrac_short[j][nc] = 0.0



        # CACHE TO DISK 
        if 0: 
            hfivename = 'polarization_%s_ii.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['core_id']=self.cores_used
            Fptr['theta_long']=self.theta_long
            Fptr['theta_short']=self.theta_short
            Fptr['pol_frac_long']=self.UQfrac_long 
            Fptr['pol_frac_short']=self.UQfrac_short
            Fptr.close()



# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=polarization(TL6.loops['u601'])
if 'scope2' not in dir() or clobber:
    scope2=polarization(TL6.loops['u602'])
if 'scope3' not in dir() or clobber:
    scope3=polarization(TL6.loops['u603'])

simnames = ['u601','u602', 'u603']


for nt,tool in enumerate([scope1,scope2,scope3]):

    # WHICH CORES & RUN IF NOT CACHED
    if 0:
        all_cores = np.unique(tool.this_looper.tr.core_ids)
        #core_list = all_cores[:1]  #DEBUG
        core_list = all_cores
        tool.qtyRun(nt,core_list=core_list) 

    # ONCE CACHED, HISTOGRAMS SCATTERS & REL ERROR REGARDING POLARIZATION 
    if 1: 
        hfivename = 'p66_brho/polarization_%s.h5'%(nt)
        Fptr = h5py.File(hfivename,'r')
        cores_used = Fptr['core_id'][()]
        
        theta_long = Fptr['theta_long'][()]
        theta_short = Fptr['theta_short'][()]
        diff_theta = abs(theta_short - theta_long)  # EDIT FOR > 90 BELOW 
        relerr_theta = abs((theta_short - theta_long)/theta_short)  

        UQfrac_long = Fptr['pol_frac_long'][()]
        UQfrac_short = Fptr['pol_frac_short'][()] 
        relerr_UQfrac = abs((UQfrac_short - UQfrac_long)/UQfrac_short) 

        ''' 
        # WITHOUT CACHING
        theta_long = tool.theta_long
        theta_short = tool.theta_short
        UQfrac_long = tool.UQfrac_long
        UQfrac_short = tool.UQfrac_short
        '''
  
        the_bins = 128  #default
        dxyz = ['x','y','z']
        color=['b','g','orange']
        relerr_UQfrac_all = []
        diff_theta_all = []
        for i in range(3):

            index = []
            # CORRECTED FOR THETA DIFF > 90
            for j,theta in enumerate(diff_theta[i]):
                if theta > 90:
                    #print('theta',theta)
                    minusoneeighty = 180 - theta
                    #print('corrected theta',minusoneeighty)
                    diff_theta[i][j] = minusoneeighty

            # CORRECTED FOR POL FRAC ERROR > 5
            for k,frac in enumerate(relerr_UQfrac[i]):
                if frac > 1:
                    print('polfrac',frac)
                    print('which core, which direction',cores_used[k],i)
                    index.append(k)

            # HISTOGRAMS
            if 0:  # THETA PDFS & CDFS 
                pdf_long, bins_long = np.histogram(theta_long[i], bins=the_bins, density=True)
                cdf_long = np.cumsum(pdf_long)
                bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
                plt.plot(bin_centers_long,pdf_long,c=color[i])

                pdf_short, bins_short = np.histogram(theta_short[i], bins=the_bins, density=True)
                cdf_short = np.cumsum(pdf_short)
                bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
                plt.plot(bin_centers_short,pdf_short,c=color[i],linestyle='dashed')            
            if 0:  # POL FRACS LONG & SHORT PDFS & CDFS 
                pdf_long, bins_long = np.histogram(UQfrac_long[i], bins=the_bins, density=True)
                cdf_long = np.cumsum(pdf_long)
                bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
                plt.plot(bin_centers_long,cdf_long,c=color[i])#alpha=0.5)

                pdf_short, bins_short = np.histogram(UQfrac_short[i], bins=the_bins, density=True)
                cdf_short = np.cumsum(pdf_short)
                bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
                plt.plot(bin_centers_short,cdf_short,c=color[i],linestyle='dashed')#,alpha=0.5)

            # SCATTERS
            if 0:  # THETA 
                plt.scatter(theta_long[i], theta_short[i], alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=1.0)
                outname = 'polarization_fraction_relerr_%d_%s'%(i,simnames[nt])
            if 0:  # POLARIZATION FRACS 
                plt.scatter(UQfrac_long[i], UQfrac_short[i], alpha=0.7)
                plt.axline((0, 0), slope=1, c='k', linewidth=1.0)

            # RELATIVE ERRORS
            if 0:  # THETA 
                pdf, bins = np.histogram(relerr_theta[i], bins=the_bins, density=True)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
            if 0:  # SAVE, LABEL
                plt.xlabel(r'$\theta_%d$ rel error'%i)
                plt.ylabel('PDF')
                outname = 'polarization_theta_relerr_%d_%s'%(i,simnames[nt])

            if 0:  # POLARIZATION FRACS 
                #the_bins = 64 
                pdf, bins = np.histogram(relerr_UQfrac[i], bins=the_bins, density=True)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,cdf,c=color[i],alpha=0.7)
            if 0:  # SAVE, LABEL
                plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ rel error %d'%i)
                plt.ylabel('PDF')
                outname = 'polarization_fraction_relerr_%d_%s'%(i,simnames[nt])

            # 3MT PLOTS 
            if 0: 
                the_bins = 16 
                #the_weights = np.ones(len(np.delete(diff_theta[i],index)))/len(np.delete(diff_theta[i],index))
                pdf, bins = np.histogram(np.delete(diff_theta[i],index), bins=the_bins, density=True)#, weights=the_weights)
                cdf= np.cumsum(pdf)
                bin_centers = 0.5*(bins[1:]+bins[:-1]) 
                plt.plot(bin_centers,pdf,c=color[i],alpha=0.7)
                #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                plt.xlabel(r'$\Delta\theta_{x,y,z}$')
                plt.ylabel('PDF')
                outname ='theta_diff_%d_%s'%(i,simnames[nt])
                #outname ='theta_diff_%s'%(simnames[nt])
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
                plt.scatter(relerr_UQfrac[i],diff_theta[i],c=color[i],alpha=0.7)
                plt.xlabel(r'$(P_{short}-P_{long})/P_{short}$')
                plt.ylabel(r'$\theta_{short:x,y,z} - \theta_{long:x,y,z}$')
                plt.xlim(-0.05,1.05)
                plt.ylim(-2,92)
                outname ='thetadiff_vs_relerrUQfrac_%d_%s'%(i,simnames[nt])
                #outname ='thetadiff_vs_relerrUQfrac_%s'%(simnames[nt])
            if 0:  # OR HEAT MAP 
                the_bins = 32 
                fig, ax = plt.subplots(1,1)
                output = pch.simple_phase(np.delete(relerr_UQfrac[i],index), np.delete(diff_theta[i],index), log=False, ax=ax, nBins=the_bins)
                fig.colorbar(output['plot'])
                ax.set_xlim(-0.05,1.05)
                ax.set_ylim(-2,92)
                outname ='thetadiff_vs_relerrUQfrac_heatmap_%d_%s'%(i,simnames[nt])
                fig.savefig(outname)
                plt.close('all')
            if 0:  # TO PLOT ALL AT ONCE
                relerr_UQfrac_all = np.concatenate((relerr_UQfrac_all,np.delete(relerr_UQfrac[i],index)))
                diff_theta_all = np.concatenate((diff_theta_all,np.delete(diff_theta[i],index)))


            # SAVE ONE DIRECTION AT A TIME
            if 0:
                plt.savefig(outname)
                print('plotted_%d'%i)
                plt.close('all')


        # ALL DIRECTIONS AT ONCE  
        if 0:
            plt.xlabel(r'$\theta$')
            plt.ylabel('PDF')
            plt.ylim(0,0.40)
            outname = 'polarization_theta_shortlong_xyz_%s'%(simnames[nt])
        if 0:
            plt.xlabel(r'$\theta$ rel error')
            plt.ylabel('PDF')
            plt.xlim(0.0,1.2)
            plt.ylim(-0.2,7.0)
            outname = 'polarization_theta_relerr_zoomed_xyz_%s'%(simnames[nt])

        if 0:
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$')
            plt.ylabel('PDF')
            outname = 'polarization_fraction_CDFshortlong_xyz_%s'%(simnames[nt])
        if 0: 
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ rel error')
            plt.ylabel('PDF')
            #plt.xlim(0.0,1.2)
            #plt.ylim(-0.2,9.0)
            outname = 'polarization_fraction_relerr_cdf_bins64_xyz_%s'%(simnames[nt])

        if 0: 
            plt.xlabel(r'$\theta$ LONG')
            plt.ylabel(r'$\theta$ SHORT')
            outname = 'pol_theta_short_vs_long_%s'%(simnames[nt])
        if 0: 
            plt.xlabel(r'$(Q^2 + U^2)^{1/2}/N$ LONG')
            plt.ylabel(r'$(Q^2 + U^2)^{1/2}/N$ SHORT')
            plt.ylim(-0.02,0.85)
            plt.xlim(-0.02,0.6)
            outname = 'pol_fraction_short_vs_long_%s'%(simnames[nt])

        # SAVING ALL AT ONCE
        if 0:
            plt.savefig(outname)
            print('plotted_%d'%i)
            plt.close('all')

        if 0:  # HEAT MAP 
            the_bins = 32 
            fig, ax = plt.subplots(1,1)
            output = pch.simple_phase(relerr_UQfrac_all, diff_theta_all, log=False, ax=ax, nBins=the_bins)
            fig.colorbar(output['plot'])
            ax.set_xlim(-0.05,1.05)
            ax.set_ylim(-2,92)
            plt.xlabel(r'$P_{error:x,y,z}$')
            plt.ylabel(r'$\Delta\theta_{x,y,z}$')
            outname ='thetadiff_vs_relerrUQfrac_heatmap_%s'%(simnames[nt])
            fig.savefig(outname)
            plt.close('all')

        # AND CLOSE H5 FILE
        Fptr.close()


