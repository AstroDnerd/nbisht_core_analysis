
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
np.set_printoptions(threshold=sys.maxsize)

from icecream import ic

class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []
        self.demagnetized_long = []
        self.demagnetized_short = []

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
        
        self.radii = [np.zeros(len(core_list)) for x in range(3)]

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:] 


        # CORES
        for nc,core_id in enumerate(core_list):

            # OPEN THE FIGURES
            if 0:
                fig, ax = plt.subplots(1,1)  #just want B field for now
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
            cylinders = ['short']
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
                if 1:
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
                    self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['cell_volume'].sum()
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['density']).sum()/the_mid_cyl[j]['density'].sum()
                    self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()


                    # HOW MUCH DEMAGNETIZATION
                    if 0:
                        dxyz = ['x','y','z']
                        if cylinders[val] == 'long':
                            logged = np.log10(abs(self.synthField[j][nc]))
                        if cylinders[val] == 'short':
                            logged = np.log10(abs(self.synthField[j][nc]))
                        if logged <= 1: 
                            print("an outlier core! ",j,core_id)

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
                            sortedrho = (the_cyl[j]['density'])[f_order]
                            cumsum_bf = np.cumsum(sortedBrho)
                            fieldsum = cumsum_bf.v
                            cumsum_rhof = np.cumsum(sortedrho)
                            rhosum = cumsum_rhof.v

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

                        # LET'S PLOT
                        save_full = 'magnetization_full_%d_%d_%s_zoom2'%(j,core_id,sim)
                        save_mid = 'magnetization_mid_%d_%d_%s'%(j,core_id,sim)
                        if 0:
                            #m_ax[0].scatter(zm,cumsum_rhom,color='k',alpha=0.5)
                            #m_ax[1].scatter(zm,cumsum_bm,color='b',alpha=0.5)

                            #ax[0].scatter(zf,cumsum_rhof,color='k',alpha=0.5)
                            ax.scatter(zf,cumsum_bf,color='b',alpha=0.5)
                            ax.set_xlabel('i')
                            ax.set_xlim(0.0750,0.08)
                            ax.set_ylabel('int rho * Bi di')
                        
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
    #core_list = all_cores[:3]  #DEBUG
    core_list = all_cores

    # RUN
    tool.qtyRun(nt,core_list=core_list) 
    if 0:  #HISTOGRAM THE DEMAGNETIZATION
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
        
    if 1:  #HISTOGRAM THE <B_particles>/<B_los>  
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

    plt.savefig(outname)
    print('plotted the histogram')
    plt.close('all')
    #pdb.set_trace()



