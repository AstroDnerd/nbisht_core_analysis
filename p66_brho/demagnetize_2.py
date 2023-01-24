
'''
synthetic observations, version 2
'''

from cmath import log
from re import X
from turtle import width
from starter2 import *
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)
import multiplots 
reload(multiplots)
import davetools
reload(davetools)
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from matplotlib.ticker import PercentFormatter

np.set_printoptions(threshold=sys.maxsize)
from icecream import ic

class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.demagnetized_long = []
        self.demagnetized_short = []
        self.demagnetized_long_ratioXYZ = []
        self.demagnetized_short_ratioXYZ = []

    def qtyRun(self,sim,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        # CORE_LIST
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        self.Rho = [np.zeros(len(core_list)) for x in range(3)]
        self.Field = [np.zeros(len(core_list)) for x in range(3)]

        self.synthRho = [np.zeros(len(core_list)) for x in range(3)]
        self.synthRho_mid = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_mid = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_abs = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_mid_abs = [np.zeros(len(core_list)) for x in range(3)]
        
        self.radii = [np.zeros(len(core_list)) for x in range(3)]

        self.demagnetized_long_ratio = [np.zeros(len(core_list)) for x in range(3)]
        self.demagnetized_short_ratio = [np.zeros(len(core_list)) for x in range(3)]

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:] 


        # CORES
        for nc,core_id in enumerate(core_list):

            # OPEN THE FIGURES
            if 0:
                fig, ax = plt.subplots(1,1)  #see this option below
                #m_fig, m_ax = plt.subplots(1,2)

            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # PIECES FOR THE OBJECTS
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            #self.radii[sim][nc] = rinf[sim][nc] # OR a fixed 1/128; what is most realistic for a telescope
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]

            twoeight = 1/128
            fivesix = 1/256
            radiuses = [twoeight]#, fivesix] # 0, or 1
            if 0:  #if looping the two radiuses
                fig, ax = plt.subplots(1,2)

            #cylinders = ['short','long']
            cylinders = ['short']  #doesn't choose anything, simply iterates once
            for val in range(len(cylinders)):
                the_radius = twoeight 
                the_area= np.pi * (the_radius**2) 

                # MAKE THE OBJECTS:
                xyz = [0,1,2]
                the_cyl = {}
                the_mid_cyl = {}
                for i in range(3):
                    the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))
                    the_mid_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=the_radius) 


                # THE DENSITY & FIELD: 3D 
                if 0:
                    all_frames = thtr.frames
                    for nf,frame in enumerate(all_frames): 
                        if frame == all_frames[-1]:
                            mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf) 
                            cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                            density = thtr.c([core_id],'density')[mask,nf]
                            bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                            by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                            bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                            b = [bx, by, bz]
                            for j in range(3):
                                self.Field[j][nc] = (b[j] * cell_volume).sum()/cell_volume.sum()  #try weight by density?
                                self.Rho[j][nc] = (density * cell_volume).sum()/(cell_volume.sum())  


                # THE DENSITY & FIELD, ZONE METHOD:  2D 
                B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
                for j in range(3):
                    self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                    #self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['cell_volume'].sum()
                    #self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['density']).sum()/the_cyl[j]['density'].sum()  #TRY..
                    self.synthField[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()
                    the_cyl_Babs = np.abs(the_cyl[j][B[j]])
                    self.synthField_abs[j][nc] = (the_cyl[j]['density'] * the_cyl_Babs * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()

                    self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['cell_volume'].sum()
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['density']).sum()/the_mid_cyl[j]['density'].sum()
                    self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()
                    the_mid_cyl_Babs = np.abs(the_mid_cyl[j][B[j]])
                    self.synthField_mid_abs[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl_Babs * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()


                    # DEMAGNETIZATION
                    # HOW MUCH DEMAGNETIZATION, TAKE II
                    if 1:  # accounts all cores
                        notabs_long = np.abs(self.synthField[j][nc])
                        notabs_short = np.abs(self.synthField_mid[j][nc])

                        abs_long = self.synthField_abs[j][nc]
                        abs_short = self.synthField_mid_abs[j][nc]
                        
                        ratio_long = notabs_long/abs_long 
                        ratio_short = notabs_short/abs_short 

                        self.demagnetized_long_ratio[j][nc] = ratio_long
                        self.demagnetized_short_ratio[j][nc] = ratio_short
                        self.demagnetized_short_ratioXYZ.append(ratio_long)
                        self.demagnetized_short_ratioXYZ.append(ratio_short)


                    # HOW MUCH DEMAGNETIZATION, TAKE I
                    if 0:
                        fig, ax = plt.subplots(1,1)  #just want B field for now
                    dxyz = ['x','y','z']
                    if 0:
                        m_order = np.argsort(the_mid_cyl[j][dxyz[j]]) 
                        zm = the_mid_cyl[j][dxyz[j]][m_order]
                        sortedmBrho = (the_mid_cyl[j]['density']*the_mid_cyl[j][B[j]])[m_order]
                        sortedmrho = (the_mid_cyl[j]['density'])[m_order]
                        cumsum_bm = np.cumsum(sortedmBrho)
                        mfieldsum = cumsum_bm.v
                        cumsum_rhom = np.cumsum(sortedmrho) 
                        mrhosum = cumsum_rhom.v

                        f_order = np.argsort(the_cyl[j][dxyz[j]]) 
                        zf = the_cyl[j][dxyz[j]][f_order]
                        sortedBrho = (the_cyl[j]['density']*the_cyl[j][B[j]])[f_order]
                        sortedBrho_abs = (the_cyl[j]['density']*the_cyl_Babs)[f_order]
                        sortedrho = (the_cyl[j]['density'])[f_order]
                        cumsum_bf = np.cumsum(sortedBrho)
                        cumsum_bf_abs = np.cumsum(sortedBrho_abs)
                        fieldsum = cumsum_bf.v
                        fieldsum_gauss = gaussian_filter(fieldsum,5)
                        fieldsum_abs = cumsum_bf_abs.v
                        cumsum_rhof = np.cumsum(sortedrho)
                        rhosum = cumsum_rhof.v

                    # revise, this may not be needed, can do both at once w/o loop
                    if cylinders[val] == 'long':
                        logged = np.log10(abs(self.synthField[j][nc]))
                    if cylinders[val] == 'short':
                        logged = np.log10(abs(self.synthField_mid[j][nc]))

                    if logged <= 1: 
                        print("an outlier core! ",j,core_id)

                        if 0:  #HEADS UP, cylinders[val] doesn't mean anything right now
                            if cylinders[val] == 'long':
                                bf = fieldsum
                                bmax = max(fieldsum) 
                                bmin = min(fieldsum)
                                i_max = 0
                                i_min = 0
                                i_zero = 0
                                for i, b in enumerate(fieldsum):
                                    if b == bmax:
                                        i_max = i
                                        print('i_max found')
                                    if b == bmin:
                                        i_min = i
                                        print('i_min found')

                                # part i, ii: 
                                deb_one = 0
                                deb_two = 0
                                demagnet = 0
                                if i_max < i_min:
                                    deb_one = abs(bf[i_max] - bf[0]) 
                                    deb_two = abs(bf[i_min] - bf[-1])
                                if i_min < i_max:
                                    deb_one = abs(bf[i_min] - bf[0])
                                    deb_two = abs(bf[i_max] - bf[-1])
                                demagnet = deb_one + deb_two
                                if demagnet == 0:
                                    print('a zero!')
                                if demagnet != 0.0:
                                    self.demagnetized_long.append(demagnet)

                            if cylinders[val] == 'short':
                                bf = mfieldsum
                                bmax = max(mfieldsum) 
                                bmin = min(mfieldsum)
                                i_max = 0
                                i_min = 0
                                i_zero = 0
                                for i, b in enumerate(mfieldsum):
                                    if b == bmax:
                                        i_max = i
                                        print('i_max found')
                                    if b == bmin:
                                        i_min = i
                                        print('i_min found')

                                # part i, ii: 
                                deb_one = 0
                                deb_two = 0
                                demagnet = 0
                                if i_max < i_min:
                                    deb_one = abs(bf[i_max] - bf[0]) 
                                    deb_two = abs(bf[i_min] - bf[-1])
                                if i_min < i_max:
                                    deb_one = abs(bf[i_min] - bf[0])
                                    deb_two = abs(bf[i_max] - bf[-1])
                                demagnet = deb_one + deb_two
                                if demagnet == 0:
                                    print('a zero!')
                                if demagnet != 0.0:
                                    self.demagnetized_short.append(demagnet)


                    # TRY SMOOTHING THE CURVE
                    def smooth(y, box_pts):
                        box = np.ones(box_pts)/box_pts
                        y_smooth = np.convolve(y, box, mode='same')
                        return y_smooth
                    # TRY SMOOTHING THE CURVE ANOTHER TAKE
                    #ysg = savgol_filter(fieldsum, 51, 10)

                    # LET'S PLOT
                    save_full = 'magnetization_long_%d_%d_%s_plot'%(j,core_id,sim)
                    save_mid = 'magnetization_mid_%d_%d_%s'%(j,core_id,sim)
                    if 0:
                        #m_ax[0].scatter(zm,cumsum_rhom,color='k',alpha=0.5)
                        #m_ax[1].scatter(zm,cumsum_bm,color='b',alpha=0.5)

                        ax.plot(zf,fieldsum,color='b',alpha=0.7)
                        #ax.plot(zf,smooth(fieldsum,1000),'r',alpha=0.5)
                        #ax.plot(zf,ysg,'g')
                        if 1:  #plot abs vs no abs
                            ax.plot(zf,fieldsum_abs,color='g',alpha=0.7)

                        ax.set_xlabel('i')
                        ax.set_ylabel('int rho * Bi di')

                        if core_id in [121,122,0]:
                            ax.set_xlim(0.0,0.1)
                        if core_id == 121 and j==1:
                            ax.set_xlim(0.450,0.525)
                        if core_id == 121 and j==2:
                            ax.set_xlim(0.45,0.5)
                        if core_id == 1:
                            ax.set_xlim(0.06,0.1)
                        if core_id in [21,127]:
                            ax.set_xlim(0.0,0.2)
                        if core_id in [171]:
                            ax.set_xlim(0.0,0.3)
                        if core_id == 32:
                            ax.set_xlim(0.1,0.25)
                        if core_id in [44,64]:
                            ax.set_xlim(0.1,0.3)
                        if core_id == 44 and j==2:
                            ax.set_xlim(0.125,0.2)
                        if core_id in [180]:
                            ax.set_xlim(0.2,0.6)
                        if core_id in [75,74]:
                            ax.set_xlim(0.3,0.5)
                        if core_id == 135:
                            ax.set_xlim(0.350,0.425)
                        if core_id == 85:
                            ax.set_xlim(0.3,0.6)
                        if core_id in [148]:
                            ax.set_xlim(0.4,0.6)
                        if core_id in [67,83,292]:
                            ax.set_xlim(0.5,0.7)
                        if core_id in [308]:
                            ax.set_xlim(0.6,0.7)
                        if core_id == 286 and j==0:
                            ax.set_xlim(0.64,0.67)
                        if core_id == 286 and j==1:
                            ax.set_xlim(0.5,0.6)
                        if core_id == 93:
                            ax.set_xlim(0.6,0.8)
                        if core_id == 123:
                            ax.set_xlim(0.52,0.56)
                        if core_id == 248 and j==0:
                            ax.set_xlim(0.86,0.9)
                        if core_id in [231]:
                            ax.set_xlim(0.8,0.9)
                        if core_id == 248 and j == 2:
                            ax.set_xlim(0.0,0.6)
                        if core_id == 294:
                            ax.set_xlim(0.73,0.76)
                        if core_id == 297:
                            ax.set_xlim(0.65,0.7)
                        if core_id == 276:
                            ax.set_xlim(0.350,0.4)
                        if core_id == 201:
                            ax.set_xlim(0.450,0.525)
                        if core_id == 126:
                            ax.set_xlim(0.54,0.58)


                    if 0:  #GO AHEAD, SAVE 
                        fig.savefig(save_full)
                        #m_fig.savefig(save_mid)
                        plt.close(fig)
                        #plt.close(m_fig)




# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
if 'scope1' not in dir() or clobber:
    scope1=telescope(TL6.loops['u601'])
if 'scope2' not in dir() or clobber:
    scope2=telescope(TL6.loops['u602'])
if 'scope3' not in dir() or clobber:
    scope3=telescope(TL6.loops['u603'])

simnames = ['u601','u602', 'u603']

for nt,tool in enumerate([scope1,scope2,scope3]): 
    # WHICH CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    #core_list = all_cores[:5]  #DEBUG
    core_list = all_cores

    # RUN
    tool.qtyRun(nt,core_list=core_list) 
    if 0:  #HISTOGRAM THE DEMAGNETIZATION TAKE I
        demagnet_long= tool.demagnetized_long
        demagnet_long_log = np.log10(demagnet_long)
        the_bins_long = math.isqrt(len(demagnet_long)) 
        plt.hist(demagnet_long_log, bins=the_bins_long, density=True, histtype='step')

        demagnet_short= tool.demagnetized_short
        demagnet_short_log = np.log10(demagnet_short)
        the_bins_short = math.isqrt(len(demagnet_short)) 
        plt.hist(demagnet_short_log, bins=the_bins_short, density=True, histtype='step',linestyle=('dashed'))

        plt.xlabel('log(B) G*g/cm**3')
        plt.ylabel('PDF')
        outname = 'demagnetized_%s_v2'%simnames[nt]
        
    if 1:  #HISTOGRAM THE DEMAGNETIZATION TAKE II; xyz
        the_bins = 128
        
        demagnet_long = tool.demagnetized_long_ratio
        demagnet_short = tool.demagnetized_short_ratio
  
        dxyz = ['x','y','z']
        color=['g','b','r']
        for i in range(3):
            pdf_long, bins_long = np.histogram(demagnet_long[i], bins=the_bins, density=True)
            cdf_long = np.cumsum(pdf_long)
            bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
            #plt.plot(bin_centers_long,pdf_long,c=color[i],alpha=0.5)
            plt.plot(bin_centers_long,cdf_long,c=color[i])
        
            pdf_short, bins_short = np.histogram(demagnet_short[i], bins=the_bins, density=True)
            cdf_short = np.cumsum(pdf_short)
            bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
            #plt.plot(bin_centers_short,pdf_short,c=color[i],linestyle='dashed',alpha=0.5)
            plt.plot(bin_centers_short,cdf_short,c=color[i],linestyle='dashed')
            
            if 0:
                plt.xlabel('(int B%s * rho * d%s)/(int |B%s| * rho * d%s)'%(dxyz[i],dxyz[i],dxyz[i],dxyz[i]))
                plt.ylabel('PDF')

        if 1:
            plt.xlabel('(int Bi * rho * di)/(int |Bi| * rho * di)')
            plt.ylabel('CDF')

        outname = 'demagnetized_inverse_cdfxyz_%s'%(simnames[nt])
        plt.savefig(outname)
        print('plotted_%d'%i)
        plt.close('all')

    if 0:  #HISTOGRAM THE DEMAGNETIZATION TAKE II
        the_bins = 128

        demagnet_long= tool.demagnetized_long_ratio
        the_bins_long = math.isqrt(len(demagnet_long))  #18
        pdf_long, bins_long = np.histogram(demagnet_long, bins=the_bins, density=True)#, histtype='step')#, weights=np.ones(len(demagnet_long))/len(demagnet_long))
        cdf_long = np.cumsum(pdf_long)
        bin_centers_long = 0.5*(bins_long[1:]+bins_long[:-1]) 
        plt.plot(bin_centers_long,pdf_long,c='r')
        #plt.plot(bin_centers_long,cdf_long,c='b')
        
        demagnet_short= tool.demagnetized_short_ratio
        the_bins_short = math.isqrt(len(demagnet_short)) 
        pdf_short, bins_short = np.histogram(demagnet_short, bins=the_bins, density=True)#, histtype='step',linestyle='dashed')#, weights=np.ones(len(demagnet_short))/len(demagnet_short))
        cdf_short = np.cumsum(pdf_short)
        bin_centers_short = 0.5*(bins_short[1:]+bins_short[:-1]) 
        plt.plot(bin_centers_short,pdf_short,c='r',linestyle='dashed')
        #plt.plot(bin_centers_short,cdf_short,c='b',linestyle='dashed')
  
        print('simulation', nt)
        plt.xlabel('(int Bi * rho * di)/(int |Bi| * rho * di)')
        #plt.xlim(-0.25,1.25)
        plt.ylabel('PDF')

        #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        outname = 'demagnetized_inverse_pdf_%s'%simnames[nt]

    if 0:  #HISTOGRAM THE <B_particles>/<B_los>  
        pField3D = tool.Field
        Field3D = np.concatenate((pField3D[0],pField3D[1],pField3D[2]))
        Field_3D = abs(Field3D)

        pField = tool.synthField
        Field = np.concatenate((pField[0],pField[1],pField[2]))
        Field_long = abs(Field)
        pField_mid = tool.synthField_mid
        Field_mid = np.concatenate((pField_mid[0],pField_mid[1],pField_mid[2]))
        Field_short = abs(Field_mid)

        ratio_long = Field_3D/Field_long 
        ratio_short = Field_3D/Field_short 

        the_bins_long = math.isqrt(len(ratio_long)) 
        the_bins_short = math.isqrt(len(ratio_short)) 
        plt.hist(ratio_long, bins=64, density=True, histtype='step')
        plt.hist(ratio_short, bins=64, density=True, histtype='step',linestyle=('dashed'))

        plt.xlabel(r'$ \left\langle B_{particles} \right\rangle / \left\langle B_{los} \right\rangle $')
        plt.xlim(-1,15)
        plt.ylabel('PDF')
        outname = 'ratio_p_los_xaxiscut2%s'%simnames[nt]

    if 0: 
        plt.savefig(outname)
        print('plotted')
        plt.close('all')
        #pdb.set_trace()



