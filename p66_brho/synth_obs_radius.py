
'''
synthetic observations, version 4
'''

from cmath import log
from re import X
from turtle import width
from starter2 import *
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)


class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

    def qtyRun(self,sim,rinf,core_list=None): 
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
        for nc,core_id in enumerate(core_list):
            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # PIECES FOR THE OBJECTS
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            self.radii[sim][nc] = rinf[sim][nc] # OR a fixed 1/128; what is most realistic for a telescope
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]

            twoeight = 1/128
            fivesix = 1/256
            radiuses = [twoeight, fivesix] # 0, or 1
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


                # THE DENSITY & FIELD, ZONE METHOD:  2D 
                B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
                for j in range(3):
                    self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                    self.synthField[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()

                    self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area
                    self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()
                #print('length of appended array', len(the_cyl[0]))
                print('length of x appended array', len(the_cyl[0]['density']))

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
                                self.Field[j][nc] = (b[j] * cell_volume).sum()/cell_volume.sum()
                                self.Rho[j][nc] = (density * cell_volume).sum()/(cell_volume.sum())  

                # PROJECTIONS & THE DENSITY & FIELD, FRB: 2D
                if 1:
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

                    if 1:
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
                        print('length frb appended array',len(frbs_cyl_B))  #this should be 3, yep


                        # GET DENSITY AND FIELD
                        length = len(frbs_mid_cyl_rho[0]['gas','density'])  
                        print('length of x of frb appended array', length)  #this should be the x comp of particles(?) of one core, yep
                        for n in range(3):
                            self.synthRho_frb[n][nc] = (frbs_cyl_rho[n]['gas','density']).sum()/length
                            self.synthRho_frbmid[n][nc] = (frbs_mid_cyl_rho[n]['gas','density']).sum()/length

                            self.synthField_frb[n][nc] = (frbs_cyl_B[n][B[n]]).sum()/length
                            self.synthField_frbmid[n][nc] = (frbs_mid_cyl_B[n][B[n]]).sum()/length

                # IMAGES: PICKED Z PROJ FOR NOW
                    if 0:  
                        theArray_rho = frbs_cyl_rho[2]['gas','density'].v
                        theArray_midrho = frbs_mid_cyl_rho[2]['gas','density'].v 

                        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                        cmap.set_under('w')
                        minmin = theArray_rho[theArray_rho>0].min()
                        norm = mpl.colors.LogNorm(vmin=minmin,vmax=theArray_rho.max())
                        
                        # TO OVERLAY B FIELD
                        X = np.unique(frbs_cyl_rho[2]['x'])
                        Y = np.unique(frbs_cyl_rho[2]['y'])
                        xx,yy = np.meshgrid(X,Y) 

                        if the_radius == 1/128: 
                            ax[0].imshow(theArray_rho, cmap=cmap,norm=norm)#,shading='nearest')
                            ax[0].streamplot(xx,yy,frbs_cyl_B[2][B[0]],frbs_cyl_B[2][B[1]])
                            print('inside 128')
                        if the_radius == 1/256: 
                            ax[1].imshow(theArray_rho, cmap=cmap,norm=norm)#,shading='nearest')
                            ax[1].streamplot(xx,yy,frbs_cyl_B[2][B[0]],frbs_cyl_B[2][B[1]])
                            print('inside 256')

                            plt.savefig("frbrho_128_256_z_core%d_%s.png"%(core_id,sim))
                            plt.close('all')



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

if 1:  # TO OBTAIN INFLECTION RADIUS
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
DB_frb = {}
DB_mid = {}
DB_mid_frb = {}
for nt,tool in enumerate([scope1]):#,scope2,scope3]): 
    # WHICH CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    core_list = all_cores[2:4]  #DEBUG
    #core_list = all_cores

    # RUN
    tool.qtyRun(nt,rinfs,core_list=core_list) 

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
    the_x = Rho_frb
    the_y = Rho_midfrb
    the_xx = Field_frb
    the_yy = Field_midfrb

    RHO = np.log10(the_x) 
    BLOS = np.log10(abs(the_y))  
    RHOO = np.log10(the_xx) 
    BLOSS = np.log10(abs(the_yy))  
    ok = BLOS > 1

    if 0:
        # WHICH CORES ARE LOW
        core_low = []
        three_core_list = np.concatenate((core_list,core_list,core_list))
        for j,k in enumerate(three_core_list):
            if ok[j] == False:
                core_low.append(k)
        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d BLow_cores: %s \n"%(nt,core_low))
        alphaFile.close()
        print('core_low',core_low)

    # WHAT IS THE ALPHA VALUE
    pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
    alpha = pfit[0]
    BLOS_o = pfit[1]
    
    DB[nt] = []
    DB_frb[nt] = []
    DB_mid[nt] = []
    DB_mid_frb[nt] = []
    DB[nt].append(alpha)

    # ALPHA VALUES & PLOTS 
    if 0:
        fig,ax = plt.subplots(1,1)
        ax.scatter(the_x,the_y,color='r',alpha=0.4)
        ax.scatter(the_x[ok],the_y[ok],color='b',alpha=0.4)
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)
        ax.plot(RHO_X,BLOS_Y,color='k') 
        
        # write it to the records        
        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,DB[nt][0]))
        alphaFile.close()
        print('sims alphas ',DB[nt])

    # PLOTTINGS
    if 0: 
        fig, ax = plt.subplots(1,2)
        ax[0].scatter(the_x,the_y,color='k',alpha=0.5)
        ax[0].axline((0, 0), slope=1, c='k')
        ax[1].scatter(the_xx,the_yy,color='b',alpha=0.5)
        ax[1].axline((0, 0), slope=1, c='k')

        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        label_rho3D =r'$\rho_{3D}$' 
        label_rho2D =r'$\rho_{2D}$' 
        label_rho2Dmid =r'$\rho_{2D,mid}$' 
        label_b3D = r'$B_{3D}$'
        label_b2D = r'$B_{2D}$'
        label_b2Dmid = r'$B_{2D,mid}$'
        label_blos = r'$\left\langle\mid B_i \mid\right\rangle (\mu G)$'
        label_Nlos = r'$\left\langle N_i \right\rangle$'
        ax[0].set_xlabel(label_rho2D)
        ax[0].set_ylabel(label_rho2Dmid)
        #ax[0].set_xlim(1e-1,1e1)
        #ax[0].set_ylim(1e-1,1e1)
        ax[0].set_aspect('equal')

        ax[1].set_xlabel(label_b2D)
        ax[1].set_ylabel(label_b2Dmid)
        #ax[1].set_xlim(1e-1,1e1)
        #ax[1].set_ylim(1e-1,1e1)
        ax[1].set_aspect('equal')

    savename = 'BxyzVsNxyz_cyl128_r1_synth_%s.png'%nt
    savenameD = 'RhoB2DmidvsRhoB2D_128frb_%s.png'%nt 
    #fig.tight_layout() 
    #fig.savefig(savenameD)
    #plt.close(fig)


# EDIT
# call synth_obs_plot.py
# pdb.set_trace()
