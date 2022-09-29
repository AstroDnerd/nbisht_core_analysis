
'''
synthetic observations, version 3
'''

from starter2 import *
import davetools
reload(davetools)
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)


class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

    def qtyRun(self,sim,rinf,core_list=None): 
        print('inside qtyRun')
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:]
        self.synthRho = [np.zeros(len(core_list)) for x in range(3)]
        self.synthRho_mid = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_mid = [np.zeros(len(core_list)) for x in range(3)]
        
        self.synthRho_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.synthRho_frbmid = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_frb = [np.zeros(len(core_list)) for x in range(3)]
        self.synthField_frbmid = [np.zeros(len(core_list)) for x in range(3)]

        # CORES
        for nc,core_id in enumerate(core_list):
            ds = self.this_looper.load(the_frame[0]) 
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # PIECES FOR THE OBJECTS
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_radius = rinf[sim][nc] # OR a fixed 1/128; what is most realistic for a telescope
            the_area= np.pi * (the_radius**2) 
            the_normal = [[1,0,0],[0,1,0],[0,0,1]]

            # MAKE THE OBJECTS:
            xyz = [0,1,2]
            the_cyl = {}
            the_mid_cyl = {}
            for i in range(3):
                the_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=(1,'code_length'))
                the_mid_cyl[xyz[i]] = ds.disk(the_center,the_normal[i],the_radius,height=the_radius) 


            # THE DENSITY & FIELD, ZONE METHOD: 
            B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
            for j in range(3):
                self.synthRho[j][nc] = (the_cyl[j]['density'] * the_cyl[j]['cell_volume']).sum()/the_area
                self.synthField[j][nc] = (the_cyl[j]['density'] * the_cyl[j][B[j]] * the_cyl[j]['cell_volume']).sum()/the_cyl[j]['gas','cell_mass'].sum()

                self.synthRho_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j]['cell_volume']).sum()/the_area
                self.synthField_mid[j][nc] = (the_mid_cyl[j]['density'] * the_mid_cyl[j][B[j]] * the_mid_cyl[j]['cell_volume']).sum()/the_mid_cyl[j]['gas','cell_mass'].sum()

            # PROJECTIONS AND FRBS
            if 0:
                # DENSITY 
                projs_cyl_rho = []
                projs_mid_cyl_rho = []
                frbs_cyl_rho = []
                frbs_mid_cyl_rho = []
                for k in range(3):
                    projs_cyl_rho.append(ds.proj(('gas','density'),k,data_source = the_cyl[k]))
                    projs_mid_cyl_rho.append(ds.proj(('gas','density'),k,data_source = the_mid_cyl[k]))

                    frbs_cyl_rho.append(projs_cyl_rho[k].to_frb(the_radius,[128,128],the_center))
                    frbs_mid_cyl_rho.append(projs_mid_cyl_rho[k].to_frb(the_radius,[128,128],the_center))

                # FIELD 
                projs_cyl_B = []
                projs_mid_cyl_B = []
                frbs_cyl_B = []
                frbs_mid_cyl_B = []
                for m in range(3):
                    projs_cyl_B.append(ds.proj(B[m],m,weight_field =('gas','density'),data_source = the_cyl[m]))
                    projs_mid_cyl_B.append(ds.proj(B[m],m,weight_field =('gas','density'),data_source = the_mid_cyl[m]))

                    frbs_cyl_B.append(projs_cyl_B[m].to_frb(the_radius,[128,128],the_center))
                    frbs_mid_cyl_B.append(projs_mid_cyl_B[m].to_frb(the_radius,[128,128],the_center))

                # THE DENSITY & FIELD, FRB METHOD: 
                length = len(frbs_mid_cyl_rho[0]['gas','density'])
                for n in range(3):
                    self.synthRho_frb[n][nc] = (frbs_cyl_rho[n]['gas','density']).sum()/length
                    self.synthRho_frbmid[n][nc] = (frbs_mid_cyl_rho[n]['gas','density']).sum()/length

                    self.synthField_frb[n][nc] = (frbs_cyl_B[n][B[n]]).sum()/length
                    self.synthField_frbmid[n][nc] = (frbs_mid_cyl_B[n][B[n]]).sum()/length
            #pdb.set_trace()

            # SHOW IMAGES: PICKED Z PROJ FOR NOW (EDIT)
                if 0:  
                    theArray_rho = frbs_cyl_rho[2]['gas','density'].v
                    theArray_midrho = frbs_mid_cyl_rho[2]['gas','density'].v 

                    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                    cmap.set_under('w')
                    minmin = theArray_midrho[theArray_midrho>0].min()
                    norm = mpl.colors.LogNorm(vmin=minmin,vmax=theArray_rho.max()) #what makes most sense for max?

                    plt.close('all')
                    fig, ax = plt.subplots(1,2)
                    ax[0].imshow(theArray_rho, cmap=cmap,norm=norm)#,shading='nearest')
                    ax[1].imshow(theArray_midrho, cmap=cmap,norm=norm)#,shading='nearest')
                    plt.savefig("frb_cylvsmidcyl_core%d_%s.png"%(core_id,sim))




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
for nt,tool in enumerate([scope1,scope2,scope3]): 
    # WHICH CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    #core_list = all_cores[2:4]  #DEBUG
    core_list = all_cores

    tool.qtyRun(nt,rinfs,core_list=core_list) 

    # PICKED DB FOR NOW... (EDIT)
    fig,ax = plt.subplots(1,1)
    pRho = tool.synthRho
    Rho = np.concatenate((pRho[0],pRho[1],pRho[2]))
    pRho_frb = tool.synthRho_frb
    Rho_frb = np.concatenate((pRho_frb[0],pRho_frb[1],pRho_frb[2]))
    pRho_mid = tool.synthRho_mid
    Rho_mid = np.concatenate((pRho_mid[0],pRho_mid[1],pRho_mid[2]))
    pRho_midfrb = tool.synthRho_frbmid
    Rho_midfrb = np.concatenate((pRho_midfrb[0],pRho_midfrb[1],pRho_midfrb[2]))
    
    pField = tool.synthField
    Field = np.concatenate((pField[0],pField[1],pField[2]))
    pField_frb = tool.synthField_frb
    Field_frb = np.concatenate((pField_frb[0],pField_frb[1],pField_frb[2]))
    pField_mid = tool.synthField_mid
    Field_mid = np.concatenate((pField_mid[0],pField_mid[1],pField_mid[2]))
    pField_midfrb = tool.synthField_frbmid
    Field_midfrb = np.concatenate((pField_midfrb[0],pField_midfrb[1],pField_midfrb[2]))
    

    # TO PLOT/OBTAIN ALPHAS
    if 0: 
        DB[nt] = []
        DB_frb[nt] = []
        DB_mid[nt] = []
        DB_mid_frb[nt] = []
        
        # SET THE DESIRED X & Y
        the_x = Rho_frb
        the_y = Field_frb

        RHO = np.log10(the_x) 
        BLOS = np.log10(abs(the_y))
        ok = BLOS > 1

        pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        alpha = pfit[0]
        BLOS_o = pfit[1]
        
        # PICKED DB FOR NOW... (EDIT)
        DB[nt].append(alpha)

        ax.scatter(the_x,the_y,alpha=0.4)
        ax.scatter(the_x[ok],the_y[ok],color='g',alpha=0.4)
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)  #edited  
        ax.plot(RHO_X,BLOS_Y) 

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(r'$\alpha = %f$'%alpha)
        fig.savefig('BxyzVsNxyz_cylfrb128projv2_synth_%s.png'%nt)
        plt.close(fig)

        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,DB[nt][0]))
        alphaFile.close()
        print('sims alphas ',DB[nt])


# PLOT OF THE 2D,3D ALPHAS
if 0:  # EDITTTTTT!!!
    fig,ax = plt.subplots(1,1)
    a3D = [0.456282,0.479462,0.617871]
    a2D_cyl =[0.711603,0.713381,0.731364]
    a2D_reg = [0.520261,0.446243,0.496914]
    ax.plot(a3D, a2D, 'bo',alpha=0.8)
    ax.plot(a2D, a2D_blurr, 'bs',alpha=0.8)
    fig.savefig('alpha_2D_vs_3D')
    plt.close(fig)
