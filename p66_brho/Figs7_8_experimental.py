
'''
synthetic observations, version 4
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

from icecream import ic

class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

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
        
        self.synthRho_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.synthRho_frbmid = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_frbmid = [np.zeros(len(core_list)) for x in range(3)]

        self.radii = [np.zeros(len(core_list)) for x in range(3)]

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:] 

        # CORES
        fig, ax = plt.subplots(1,2)
        #m_fig, m_ax = plt.subplots(1,2)
        for nc,core_id in enumerate(core_list):
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

            for val in range(len(radiuses)):
                the_radius = radiuses[val] 
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

                    self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['cell_volume'].sum()
                    #self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['density']).sum()/the_mid_cyl[j]['density'].sum()

                    self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()


                    # MORE INVESTIGATIONS - for outlier "high" cores: magfield and density
                    dxyz = ['x','y','z']
                    if 0:  #furthermore, create x, y, z lists 
                        if sim == 0: 
                            thirtylist = [67, 180, 286, 126] 
                            fmlist = [44, 64, 67, 73, 74, 83, 85, 121, 122, 171, 180, 201, 231, 248, 286, 294, \
                                    0, 1, 21, 93, 135, 276, 308, 32, 123, 126, 127, 148, 292, 297, \
                                    70, 80, 81, 176, 201, 286, 0, 1, 17, 121, 135, 231, 276, 286, 18, 32, 44, 92, 121, 123, 126, 148, 248, 292, 297]
                            if j == 0:
                                flist =[44, 64, 67, 73, 74, 83, 85, 121, 122, 171, 180, 201, 231, 248, 286, 294]
                                mlist =[]
                            if j ==1: 
                                flist =[]
                                mlist =[]
                            if j == 2:
                                flist = []
                                mlist =[]
                            corelist = np.unique(fmlist)
                        if sim == 1:
                            thirtylist = [193, 278]
                            fmlist = [7, 16, 37, 105, 191, 193, 270, 300, 301, 315, 329, 378, 381, 9, 84, 268, 373, \
                                    9, 49, 81, 83, 113, 211, 271, 278, 356, 377, \
                                    37, 63, 69, 105, 191, 270, 300, 347, 9, 79, 193, 301, 370, 373, 9,37, 51, 81, 113, 191, 203, 211, 271, 278, 283, 377]
                            corelist = np.unique(fmlist)
                        if sim == 2:
                            thirtylist = [3, 91, 165]
                            fmlist = [3, 16, 24, 74, 84, 88, 90, 91, 130, 136, 138, 185, 222, 258, 270, 271, 274, 275, \
                                    23, 27, 57, 89, 93, 103, 119, 124, 128, 182, 208, 219, 223, 230, 233, 237, 248, \
                                    21, 58, 83, 85, 86, 104, 140, 191, 220, 228, 229, 231, 246, 269, 291, 295, \
                                    74, 84, 136, 185, 216, 258, 270, 271, 274, 275, 27, 57, 84, 85, 93, 103, 124, 128, \
                                    136, 138, 182, 208, 237, 248,3, 8, 27, 83, 84, 85, 86, 103, 104, 140, 186, 191, 223, 231, 246, 269, 291]   
                            corelist = np.unique(fmlist)
                        if j == 0:  #change to doing full or mid one at a time and overlaying all full & mid x, y, z respectively 
                            if core_id in flist:
                                print("an outlier core! ",j)

                                m_order = np.argsort(the_mid_cyl[j][dxyz[j]]) 
                                zm = the_mid_cyl[j][dxyz[j]][m_order]
                                sortedmBrho = (the_mid_cyl[j]['density']*the_mid_cyl[j][B[j]])[m_order]
                                sortedmrho = (the_mid_cyl[j]['density'])[m_order]
                                cumsum_mb = np.cumsum(sortedmBrho)
                                cumsum_mrho = np.cumsum(sortedmrho) 

                                #the_By = cumsum_mb/abs(self.Field[j][nc]) 
                                #the_rhoy = cumsum_mrho/self.Rho[j][nc]
                                #m_ax[0].scatter(zm,the_rhoy,color='k',alpha=0.5)
                                #m_ax[1].scatter(zm,the_By,color='b',alpha=0.5)
                                #m_ax[0].scatter(zm,cumsum_mrho,color='k',alpha=0.5)
                                #m_ax[1].scatter(zm,cumsum_mb,color='b',alpha=0.5)

                                order = np.argsort(the_cyl[j][dxyz[j]]) 
                                z = the_cyl[j][dxyz[j]][order]
                                sortedBrho = (the_cyl[j]['density']*the_cyl[j][B[j]])[order]
                                sortedrho = (the_cyl[j]['density'])[order]
                                cumsum_b = np.cumsum(sortedBrho)
                                cumsum_rho = np.cumsum(sortedrho)
                                ax[0].scatter(z,cumsum_rho,color='k',alpha=0.5)
                                ax[1].scatter(z,cumsum_b,color='b',alpha=0.5)
                                
                                #save_cs3D = 'cs3Dratio_cvmid_rho3DB_abs_%d_%d_%s'%(j,core_id,sim)
                                #save_mid = 'magnetization_mid_%d_%d_%s'%(j,core_id,sim)
                                #m_fig.tight_layout() 
                                #m_fig.savefig(save_mid)
                                #fig.tight_layout() 
                                
                                #plt.close(m_fig)
                                #pdb.set_trace()
                
                                # CURRENTLY SET UP SO THAT IT IS OVERLAYED! it saves at the end of the core loop.


                # PROJECTIONS & THE DENSITY & FIELD, FRB: 2D
                if 0:
                    # PROJECT DENSITY 
                    projs_cyl_rho = []
                    projs_mid_cyl_rho = []
                    frbs_cyl_rho = []
                    frbs_mid_cyl_rho = []
                    for k in range(3):
                        projs_cyl_rho.append(yt.ProjectionPlot(ds,k,('gas','density'),center=the_center,\
                                                width=the_radius, data_source = the_cyl[k])) 
                        projs_mid_cyl_rho.append(yt.ProjectionPlot(ds,k,('gas','density'),center = the_center,\
                                                width=the_radius, data_source = the_mid_cyl[k])) 
                        frbs_cyl_rho.append(projs_cyl_rho[k].frb)
                        frbs_mid_cyl_rho.append(projs_mid_cyl_rho[k].frb)
                        
                        # TO SEE THE ANNOTATED MAGNETIC FIELDS
                        if 0:  #ONE WAY - independent
                            projs_cyl_rho[k].annotate_magnetic_field()
                            projs_cyl_rho[k].save('core%d_%s_rad%s'%(core_id,sim,val))
                            print('magfields annotated')

                    if 0:
                        # PROJECT FIELD 
                        projs_cyl_B = []
                        projs_mid_cyl_B = []
                        frbs_cyl_B = []
                        frbs_mid_cyl_B = []
                        for m in range(3):
                            projs_cyl_B.append(yt.ProjectionPlot(ds,m,B[m],weight_field =('gas','density'),center=the_center,\
                                                  width=the_radius, data_source = the_cyl[m]))
                            projs_mid_cyl_B.append(yt.ProjectionPlot(ds,m,B[m],weight_field =('gas','density'),center=the_center,\
                                                      width=the_radius, data_source = the_mid_cyl[m]))
                            frbs_cyl_B.append(projs_cyl_B[m].frb)
                            frbs_mid_cyl_B.append(projs_mid_cyl_B[m].frb)


                        # GET DENSITY AND FIELD
                        length = len(frbs_mid_cyl_rho[0]['gas','density'])  
                        for n in range(3):
                            self.synthRho_frb[n][nc] = (frbs_cyl_rho[n]['gas','density']).sum()/length
                            self.synthRho_frbmid[n][nc] = (frbs_mid_cyl_rho[n]['gas','density']).sum()/length

                            self.synthField_frb[n][nc] = (frbs_cyl_B[n][B[n]]).sum()/length
                            self.synthField_frbmid[n][nc] = (frbs_mid_cyl_B[n][B[n]]).sum()/length

                # IMAGES: PICKED Z PROJ ONLY FOR NOW
                    if 0:  
                        theArray_rho = frbs_cyl_rho[2]['gas','density'].v
                        theArray_midrho = frbs_mid_cyl_rho[2]['gas','density'].v 

                        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                        cmap.set_under('w')
                        minmin = theArray_rho[theArray_rho>0].min()
                        norm = mpl.colors.LogNorm(vmin=minmin,vmax=theArray_rho.max())
                        
                        # TO OVERLAY B FIELD
                        X = np.unique(frbs_cyl_rho[2]['x'].v)
                        Y = np.unique(frbs_cyl_rho[2]['y'].v)
                        xx,yy = np.meshgrid(X,Y)
                        #pdb.set_trace()

                        if the_radius == 1/128: 
                            ax[0].imshow(theArray_rho, cmap=cmap,norm=norm)#,shading='nearest')  #try pcolormesh
                            ax[0].streamplot(xx,yy,frbs_cyl_B[2][B[0]].v,frbs_cyl_B[2][B[1]].v)
                            print('inside 128')
                        if the_radius == 1/256: 
                            ax[1].imshow(theArray_rho, cmap=cmap,norm=norm)#,shading='nearest')
                            ax[1].streamplot(xx,yy,frbs_cyl_B[2][B[0]].v,frbs_cyl_B[2][B[1]].v)
                            print('inside 256')

                            plt.savefig("frbrho_128_256_z_core%d_%s.png"%(core_id,sim))
                            plt.close('all')

        #save_full = 'magnetization_full_%d_%s'%(j,sim)
        #fig.savefig(save_full)
        #plt.close(fig)


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

if 0:  # TO OBTAIN INFLECTION RADIUS
    if 'rinf1' not in dir():
        rinf1 = r_inflection.R_INFLECTION(TL6.loops['u601'])
        rinf_1 = rinf1.run()
    if 'rinf2' not in dir():
        rinf2 = r_inflection.R_INFLECTION(TL6.loops['u602'])
        rinf_2 = rinf2.run()
    if 'rinf3' not in dir():
        rinf3 = r_inflection.R_INFLECTION(TL6.loops['u603'])
        rinf_3 = rinf3.run()
    rinfs = [rinf_1, rinf_2, rinf_3]


simnames = ['u601','u602', 'u603']
DB = {}
# IF ALL SIMS AT ONCE
if 0: 
    fig,ax = plt.subplots(2,3)#, figsize=(6,2))
    fig.subplots_adjust(wspace=0, hspace=0)

for nt,tool in enumerate([scope1]):#,scope2,scope3]): 
    # WHICH CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    core_list = all_cores[2:4]  #DEBUG
    #core_list = all_cores

    # RUN
    tool.qtyRun(nt,core_list=core_list) 

    # WHAT TO PLOT
    pRho3D = tool.Rho
    Rho3D = np.concatenate((pRho3D[0],pRho3D[1],pRho3D[2]))

    pRho = tool.synthRho
    Rho = np.concatenate((pRho[0],pRho[1],pRho[2]))
    pRho_mid = tool.synthRho_mid
    Rho_mid = np.concatenate((pRho_mid[0],pRho_mid[1],pRho_mid[2]))

    pRho_frb = tool.synthRho_frb
    Rho_frb = np.concatenate((pRho_frb[0],pRho_frb[1],pRho_frb[2]))
    pRho_midfrb = tool.synthRho_frbmid
    Rho_midfrb = np.concatenate((pRho_midfrb[0],pRho_midfrb[1],pRho_midfrb[2]))
    
    pField3D = tool.Field
    Field3D = np.concatenate((pField3D[0],pField3D[1],pField3D[2]))

    pField = tool.synthField
    Field = np.concatenate((pField[0],pField[1],pField[2]))
    pField_frb = tool.synthField_frb
    Field_frb = np.concatenate((pField_frb[0],pField_frb[1],pField_frb[2]))
    pField_mid = tool.synthField_mid
    Field_mid = np.concatenate((pField_mid[0],pField_mid[1],pField_mid[2]))
    pField_midfrb = tool.synthField_frbmid
    Field_midfrb = np.concatenate((pField_midfrb[0],pField_midfrb[1],pField_midfrb[2]))    

    # SET THE DESIRED X & Y
    the_x = Rho
    the_y = abs(Field)
    the_xx = Rho_mid 
    the_yy = abs(Field_mid)

    rho2D_ratio = the_y/the_x 
    field2D_ratio = the_yy/the_xx
    high = field2D_ratio > 30

    RHO = np.log10(the_x) 
    BLOS = np.log10(abs(the_y))  
    RHOO = np.log10(the_xx) 
    BLOSS = np.log10(abs(the_yy))  
    ok = BLOS > 1
    okk = BLOSS > 1

    if 0:
        # WHICH CORES ARE LOW
        core_low = []
        three_core_list = np.concatenate((core_list,core_list,core_list))
        for j,k in enumerate(three_core_list):
            if ok[j] == False:
                core_low.append(k)
        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d BLow_cores 256: %s \n"%(nt,core_low))
        alphaFile.close()
        print('written')
    if 0:
        # WHICH CORES ARE HIGH
        core_high = []
        direction = []
        three_core_list = np.concatenate((core_list,core_list,core_list))
        print('length 3 core lists',len(three_core_list))
        for j,k in enumerate(three_core_list):
            if j == 112:
                print('not in x')
            if j == 225:
                print('not in y')
            if high[j] == True:
                core_high.append(k)
                direction.append(j)
                print('core_id',k)
        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d Bhigh_ratio_cores 128: %s %s \n"%(nt,core_high,direction))
        alphaFile.close()
        print('written')

    # LABELS
    label_rho3D =r'$\rho_{3D}$' 
    label_rho2D =r'$\rho_{2D_full}$' 
    label_rho2Dmid =r'$\rho_{2D,mid}$' 
    label_rho2Dmf = r'$\rho_{2D,mid}, \rho_{2D,full}$'
    label_b3D = r'$B_{3D}$'
    label_b2D = r'$B_{2D,full}$'
    label_b2Dmid = r'$B_{2D,mid}$'
    label_b2Dmf = r'$B_{2D,mid}, B_{2D,full}$'
    label_blos = r'$\left\langle\mid B_i \mid\right\rangle (\mu G)$'
    label_Nlos = r'$\left\langle N_i \right\rangle$'
    label_rratio =r'$\rho_{2Dmid}/\rho_{2D}$' 
    label_fratio =r'$B_{2Dmid}/B_{2D}$' 
    # MINS MAXS AND BINS
    the_rhomin = the_xx.min()
    the_rhomax = the_x.max()
    the_rhbins = math.isqrt(len(the_x))
    the_rhobins = np.linspace(the_rhomin,the_rhomax,num=the_rhbins)
    the_fieldmin = the_y.min()
    the_fieldmax = the_yy.max()
    the_fibins = math.isqrt(len(the_y))
    the_fieldbins = np.linspace(the_fieldmin,the_fieldmax,num=the_fibins) 

    the_rhofracmin = rho2D_ratio.min()
    the_rhofracmax = rho2D_ratio.max()
    the_rbins = math.isqrt(len(rho2D_ratio))
    the_rhofracbins = np.linspace(the_rhofracmin,the_rhofracmax,num=the_rbins)
    the_fieldfracmin = field2D_ratio.min()
    the_fieldfracmax = field2D_ratio.max()
    the_fbins = math.isqrt(len(field2D_ratio))
    the_fieldfracbins = np.linspace(the_fieldfracmin,the_fieldfracmax,num=the_fbins) 

    # WHAT IS THE ALPHA VALUE
    if 0:
        pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        alpha = pfit[0]
        BLOS_o = pfit[1] 

        mpfit = np.polyfit(RHOO[okk],BLOSS[okk],1) 
        malpha = mpfit[0]
        BLOSS_o = mpfit[1] 

        DB[nt] = []
        DB[nt].append(alpha)

    # WRITE it to the records        
    if 0:
        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,DB[nt][0]))
        alphaFile.close()
        print('sims alphas ',DB[nt])

    # SINGLE PANELS
    if 0:
        fig,ax = plt.subplots(1,1)
        ax.scatter(the_x,the_y,color='purple',alpha=0.4,marker='*')
        ax.scatter(the_x[ok],the_y[ok],color='k',marker='*')
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)
        ax.plot(RHO_X,BLOS_Y,color='grey',linestyle='dashed') 

        if 1: #OVERLAYED
            ax.scatter(the_xx,the_yy,color='olive',alpha=0.4,marker='*')
            ax.scatter(the_xx[okk],the_yy[okk],color='g',marker='*')  # could try red to match final frame of 3D
            RHO_xx = np.linspace(RHOO[okk].min(),RHOO[okk].max(),num=len(RHOO[okk]))
            RHO_XX = 10 ** RHO_xx
            BLOS_YY = 10 ** (alpha*RHO_xx + BLOSS_o)
            ax.plot(RHO_XX,BLOS_YY,color='grey',linestyle='dotted')  
            
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(40e-2,20e1)
        ax.set_ylim(1e-1,1e4)
        ax.set_xlabel(label_rho2Dmf)
        ax.set_ylabel(label_b2Dmf)

        savename = 'BxyzVsNxyz_fm_cyl128_synth_%s.png'%nt
        fig.savefig(savename)
        plt.close(fig)
        
    # PLOTTINGS - 3 WAY:  short and long; maybe a way to merge the current Fig 7 & 8.  
    if 1:
        plt.clf()
        figa, axa, axtop, axright = multiplots.three_way_bean()  #figsize=(6,6)
        '''
        axa.scatter(the_x,the_y,color='red',alpha=0.4,marker='*')
        axa.scatter(the_x[ok],the_y[ok],color='k',marker='*')
        axa.scatter(the_xx,the_yy,color='blue',alpha=0.4,marker='*')
        axa.scatter(the_xx[okk],the_yy[okk],color='g',marker='*')

        axtop.hist(the_x,bins=the_rhobins,histtype='step',color='k')
        axtop.hist(the_xx,bins=the_rhobins,histtype='step',color='g')
        axright.hist(the_y,bins=the_fieldbins,histtype='step',color='k',orientation='horizontal')
        axright.hist(the_yy,bins=the_fieldbins,histtype='step',color='g',orientation='horizontal')
        
        axa.set_xlabel(label_rho2Dmf)
        axa.set_ylabel(label_b2Dmf)
        axa.set_xscale('log')
        axa.set_yscale('log')
        axtop.set_ylabel(r'$N$') 
        axright.set_xlabel(r'$N$') 

        #axa.set_xlim([-0.1,nmax+0.1])
        #axa.set_xlim(-5,50)
        axright.set_ylim(axa.get_ylim())
        axtop.set_xlim(axa.get_xlim())
        axtop.set_xticks([])
        axright.set_yticks([])

        '''
        savename3 = 'xxx_long_short_threeway128_%s'%nt
        #figa.tight_layout()
        figa.savefig(savename3)
        
    # PLOTTINGS - 3 WAY: ratios   #maybe a way to merge the current Fig 7 & 8.  
    if 0:
        plt.clf()
        figa, axa, axtop, axright = multiplots.three_way_bean()

        axa.scatter(rho2D_ratio, field2D_ratio,c='g')
        axtop.hist(rho2D_ratio,bins=the_rhobins,histtype='step',color='k')
        axright.hist(field2D_ratio,bins=the_fieldbins,histtype='step',color='b',orientation='horizontal')
        
        axa.set_xlabel(label_rratio)
        axa.set_ylabel(label_fratio)
        axtop.set_ylabel(r'$N$') 
        axright.set_xlabel(r'$N$') 

        #axa.set_xlim([-0.1,nmax+0.1])
        #axa.set_xlim(-5,50)
        axright.set_ylim(axa.get_ylim())
        axtop.set_xlim(axa.get_xlim())
        axtop.set_xticks([])
        axright.set_yticks([])

        savename3 = 'Rhomidfull_vs_Bmidfull_threeway128_%s'%nt
        figa.savefig(savename3)


    rho_extents=davetools.extents()
    rho_extents(the_x)
    mrho_extents=davetools.extents()
    mrho_extents(the_xx)
    magfield_extents = davetools.extents()
    magfield_extents(the_y)
    mmagfield_extents = davetools.extents()
    mmagfield_extents(the_yy)
    # PLOTTINGS - 6 WAY - let's see
    if 0: 
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)
        RHOO_x = np.linspace(RHOO[okk].min(),RHOO[okk].max(),num=len(RHOO[okk]))
        RHOO_X = 10 ** RHOO_x
        BLOSS_Y = 10 ** (malpha*RHOO_x + BLOSS_o)
        if nt == 0:
            ax[0][0].scatter(the_x,the_y,color='r',alpha=0.4)
            ax[0][0].scatter(the_x[ok],the_y[ok],color='b',alpha=0.4)
            ax[0][0].plot(RHO_X,BLOS_Y,color='k') 
            ax[0][0].set_ylabel(label_b2D)
            ax[1][0].scatter(the_xx,the_yy,color='r',alpha=0.4)
            ax[1][0].scatter(the_xx[okk],the_yy[okk],color='b',alpha=0.4)
            ax[1][0].plot(RHO_X,BLOS_Y,color='k') 
            ax[1][0].set_ylabel(label_b2Dmid)
        if nt == 1: 
            ax[0][1].scatter(the_x,the_y,color='r',alpha=0.4)
            ax[0][1].scatter(the_x[ok],the_y[ok],color='b',alpha=0.4)
            ax[0][1].plot(RHO_X,BLOS_Y,color='k') 
            ax[0][1].tick_params(axis='y',labelleft=False)
            ax[1][1].scatter(the_xx,the_yy,color='r',alpha=0.4)
            ax[1][1].scatter(the_xx[okk],the_yy[okk],color='b',alpha=0.4)
            ax[1][1].plot(RHO_X,BLOS_Y,color='k') 
            ax[1][1].set_xlabel(label_rho2D)
            ax[1][1].tick_params(axis='y',labelleft=False)
                    
        if nt == 2: 
            ax[0][2].scatter(the_x,the_y,color='r',alpha=0.4)
            ax[0][2].scatter(the_x[ok],the_y[ok],color='b',alpha=0.4)
            ax[0][2].plot(RHO_X,BLOS_Y,color='k') 
            ax[0][2].tick_params(axis='y',labelleft=False)
            ax[1][2].scatter(the_xx,the_yy,color='r',alpha=0.4)
            ax[1][2].scatter(the_xx[okk],the_yy[okk],color='b',alpha=0.4)
            ax[1][2].plot(RHO_X,BLOS_Y,color='k') 
            ax[1][2].tick_params(axis='y',labelleft=False)

            for i in range(3):
                ax[0][i].set_xscale('log')
                ax[0][i].set_yscale('log')
                ax[0][i].set_aspect('equal')
                ax[0][i].set_xlim(rho_extents.minmax)
                ax[0][i].set_ylim(magfield_extents.minmax)

                ax[1][i].set_xscale('log')
                ax[1][i].set_yscale('log')
                ax[1][i].set_aspect('equal')
                ax[1][i].set_xlim(mrho_extents.minmax)
                ax[1][i].set_ylim(mmagfield_extents.minmax)
                saving = 'full_mid_128_sixway'
                fig.savefig(saving)
                plt.close(fig)

    # PLOTTINGS - TWO PANELS
    if 0: 
        fig, ax = plt.subplots(1,2)
        #fig.subplots_adjust(wspace=0, hspace=0)
        ax[0].scatter(the_x,the_xx,color='k',alpha=0.5)
        ax[0].axline((0, 0), slope=1, c='k')
        ax[1].scatter(the_y,the_yy,color='b',alpha=0.5)
        ax[1].axline((0, 0), slope=1, c='k')

        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')

        #ax[0].set_xlim(1e-1,1e1)
        #ax[0].set_ylim(1e-1,1e1)
        #ax[1].set_xlim(1e-1,1e1)
        #ax[1].set_ylim(1e-1,1e1)

        ax[0].set_xlabel(label_rho2D)
        ax[0].set_ylabel(label_rho2Dmid)
        ax[1].set_xlabel(label_b2D)
        ax[1].set_ylabel(label_b2Dmid)
        ax[0].set_xlim(1e0,30e2)
        ax[0].set_ylim(1e0,30e2)
        ax[1].set_xlim(1e0,30e2)
        ax[1].set_ylim(1e0,30e2)

        fig.tight_layout() 
        saving = '2Dmidvsfull_128_%s'%nt
        fig.savefig(saving)
        plt.close(fig)

    # PLOTTINGS - 4 PANELS
    if 0: 
        fig,ax = plt.subplots(2,2, figsize=(6,6))
        fig.subplots_adjust(wspace=0, hspace=0)

        ax[0][0].scatter(the_x,the_y,color='k',alpha=0.4)
        ax[0][0].axline((0, 0), slope=1, c='k')
        ax[0][1].scatter(the_xx,the_yy,color='b',alpha=0.4)
        ax[0][1].axline((0, 0), slope=1, c='k')
        ax[1][0].hist(rho2D_ratio, bins=the_rhobins, density=True, histtype='step', color='k')
        ax[1][1].hist(field2D_ratio, bins=the_fieldbins, density=True, histtype='step', color='b')
        ax[1][1].set_xlim(-1,50)

        y_rvals = ax[0][0].get_yticks()
        y_fvals = ax[0][1].get_yticks()
        ax[0][0].set_yticklabels(['{:.3f}'.format(x/len(rho2D_ratio)) for x in y_rvals])
        ax[0][1].set_yticklabels(['{:.3f}'.format(x/len(field2D_ratio)) for x in y_fvals])
         
        ax[0][0].set_xscale('log')
        ax[0][0].set_yscale('log')
        ax[0][0].set_xlabel(label_rho2D)
        ax[0][0].set_ylabel(label_rho2Dmid)

        ax[0][1].set_xscale('log')
        ax[0][1].set_yscale('log')
        ax[0][1].set_xlabel(label_b2D)
        ax[0][1].set_ylabel(label_b2Dmid)

        ax[1][0].set_xlabel(label_rratio)
        ax[1][1].set_xlabel(label_fratio)

        ax[0][0].xaxis.tick_top()
        ax[0][1].xaxis.tick_top()
        ax[0][1].xaxis.set_label_position('top')
        ax[0][0].xaxis.set_label_position('top')
        ax[1][1].yaxis.tick_right()
        ax[0][1].yaxis.tick_right()
        ax[1][1].yaxis.set_label_position('right')
        ax[0][1].yaxis.set_label_position('right')

        savename4 = 'RB2Dmidvsful_ratiohistos_abs_fourway128_%s'%nt
        fig.savefig(savename4)
        plt.close(fig)



# EDIT
# call synth_obs_plot.py
# pdb.set_trace()
