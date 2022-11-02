'''
synthetic observations, version 2
zones: can & boxes, frbs: can & boxes
all get their respective self array to be plotted/compared
'''

from starter2 import *
import davetools
reload(davetools)
import annotate_particles_4_cpy
reload(annotate_particles_4_cpy)
from scipy.ndimage import gaussian_filter


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
        self.synthRhox_can = np.zeros([len(core_list)]) 
        self.synthRhoy_can = np.zeros([len(core_list)])
        self.synthRhoz_can = np.zeros([len(core_list)])
        self.synthRhox_can_frb = np.zeros([len(core_list)]) 
        self.synthRhoy_can_frb = np.zeros([len(core_list)])
        self.synthRhoz_can_frb = np.zeros([len(core_list)])
        self.synthRhox_box = np.zeros([len(core_list)]) 
        self.synthRhoy_box = np.zeros([len(core_list)])
        self.synthRhoz_box = np.zeros([len(core_list)])
        self.synthRhox_box_frb = np.zeros([len(core_list)]) 
        self.synthRhoy_box_frb = np.zeros([len(core_list)])
        self.synthRhoz_box_frb = np.zeros([len(core_list)])

        self.synthRhoz_Gauss = np.zeros([len(core_list)])
        self.synthRhoy_Gauss = np.zeros([len(core_list)])
        self.synthRhox_Gauss = np.zeros([len(core_list)])

        self.synthBx_can = np.zeros([len(core_list)])
        self.synthBy_can = np.zeros([len(core_list)])
        self.synthBz_can = np.zeros([len(core_list)])  
        self.synthBx_can_frb = np.zeros([len(core_list)])
        self.synthBy_can_frb = np.zeros([len(core_list)])
        self.synthBz_can_frb = np.zeros([len(core_list)])  
        self.synthBx_box = np.zeros([len(core_list)])
        self.synthBy_box = np.zeros([len(core_list)])
        self.synthBz_box = np.zeros([len(core_list)])  
        self.synthBx_box_frb = np.zeros([len(core_list)])
        self.synthBy_box_frb = np.zeros([len(core_list)])
        self.synthBz_box_frb = np.zeros([len(core_list)])  

        self.synthBz_Gauss = np.zeros([len(core_list)])
        self.synthBy_Gauss = np.zeros([len(core_list)])
        self.synthBx_Gauss = np.zeros([len(core_list)]) 

        # CORES
        position_dict={}
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)
            ds = self.this_looper.load(the_frame[0])
            
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms

            # THE PIECES FOR THE OBJECT
            all_particles = np.stack([ms.particle_x,ms.particle_y,ms.particle_z])  #stack along 0 axis                
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 

              # FOR THE BOX
            the_radius = rinf[sim][nc]
            the_area= (2*the_radius)**2
            the_left = the_center - the_radius 
            the_right = the_center + the_radius

            the_normalX = [1,0,0] 
            the_normalY = [0,1,0]
            the_normalZ = [0,0,1] 
 
            xax = ds.coordinates.x_axis[2]  #..[proj_axis]  
            yax = ds.coordinates.y_axis[2]
            Rx = all_particles[xax]
            Ry = all_particles[yax] 
            R2d = np.sqrt(Rx**2 + Ry**2)
            #radius = R2d.max()     
            #radius = max([radius,3./128])
            radius = 1/128  #one root grid zone
            area = np.pi * radius**2


            # MAKE THE OBJECT:
              #FOR THE BOX
            the_region = ds.region(the_center,the_left,the_right)
            the_mass = the_region['gas','cell_mass'].sum()  

            the_CylX = ds.disk(the_center,the_normalX,radius,height=(1,'code_length'))
            the_CylY = ds.disk(the_center,the_normalY,radius,height=(1,'code_length'))
            the_CylZ = ds.disk(the_center,the_normalZ,radius,height=(1,'code_length'))
            massX = the_CylX['gas','cell_mass'].sum() 
            #massXb = (the_CylX['gas','density'].mean() * the_CylX['gas','cell_volume']).sum()  #delete later
            massY = the_CylY['gas','cell_mass'].sum() 
            #massYb = (the_CylY['gas','density'].mean() * the_CylY['gas','cell_volume']).sum()
            massZ = the_CylZ['gas','cell_mass'].sum() 
            #massZb = (the_CylZ['gas','density'].mean() * the_CylZ['gas','cell_volume']).sum()

 
            # THE DENSITY
            self.synthRhox_can[nc] =(the_CylX['density'] * the_CylX['cell_volume']).sum()/area  
            self.synthRhoy_can[nc] =(the_CylY['density'] * the_CylY['cell_volume']).sum()/area  
            self.synthRhoz_can[nc] =(the_CylZ['density'] * the_CylZ['cell_volume']).sum()/area  

            self.synthRhox_box[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area 
            self.synthRhoy_box[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area  
            self.synthRhoz_box[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area  

            # THE FIELD
            Bx = 'magnetic_field_x'
            By = 'magnetic_field_y'
            Bz = 'magnetic_field_z'
            self.synthBx_can[nc] = (the_CylX['density'] * the_CylX[Bx] * the_CylX['cell_volume']).sum()/massX
            self.synthBy_can[nc] = (the_CylY['density'] * the_CylY[By] * the_CylY['cell_volume']).sum()/massY
            self.synthBz_can[nc] = (the_CylZ['density'] * the_CylZ[Bz] * the_CylZ['cell_volume']).sum()/massZ
        
            self.synthBx_box[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_x'] * the_region['gas','cell_volume']).sum()/the_mass
            self.synthBy_box[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_y'] * the_region['gas','cell_volume']).sum()/the_mass
            self.synthBz_box[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_z'] * the_region['gas','cell_volume']).sum()/the_mass
            

            # PROJECTIONS
            if 1:
                # ANOTHER PROJ METHOD
                #proj = yt.ProjectionPlot(ds, 2, ('gas','density'),data_source = the_CylZ) 
                #frb = proj.frb

                # DENSITY 
                projX_can = ds.proj(('gas','density'),0,data_source = the_CylX)
                projY_can = ds.proj(('gas','density'),1,data_source = the_CylY)
                projZ_can = ds.proj(('gas','density'),2,data_source = the_CylZ)
                frbX_can = projX_can.to_frb(radius,[128,128],the_center)
                frbY_can = projY_can.to_frb(radius,[128,128],the_center) 
                frbZ_can = projZ_can.to_frb(radius,[128,128],the_center)  #width, rez, center
                
                projX_box = ds.proj(('gas','density'),0,data_source = the_region)
                projY_box = ds.proj(('gas','density'),1,data_source = the_region)
                projZ_box = ds.proj(('gas','density'),2,data_source = the_region)
                frbX_box = projX_box.to_frb(the_radius,[128,128],the_center)
                frbY_box = projY_box.to_frb(the_radius,[128,128],the_center) 
                frbZ_box = projZ_box.to_frb(the_radius,[128,128],the_center)  #width, rez, center

                # FIELD
                projBX_can = ds.proj(('gas','magnetic_field_x'),0,weight_field =('gas','density'),data_source = the_CylX)
                projBY_can = ds.proj(('gas','magnetic_field_y'),1,weight_field =('gas','density'),data_source = the_CylY)
                projBZ_can = ds.proj(('gas','magnetic_field_z'),2,weight_field =('gas','density'),data_source = the_CylZ)
                frbBX_can = projBX_can.to_frb(radius,[128,128],the_center)
                frbBY_can = projBY_can.to_frb(radius,[128,128],the_center) 
                frbBZ_can = projBZ_can.to_frb(radius,[128,128],the_center)

                projBX_box = ds.proj(('gas','magnetic_field_x'),0,weight_field =('gas','density'),data_source = the_region)
                projBY_box = ds.proj(('gas','magnetic_field_y'),1,weight_field =('gas','density'),data_source = the_region)
                projBZ_box = ds.proj(('gas','magnetic_field_z'),2,weight_field =('gas','density'),data_source = the_region)
                frbBX_box = projBX_box.to_frb(the_radius,[128,128],the_center)
                frbBY_box = projBY_box.to_frb(the_radius,[128,128],the_center) 
                frbBZ_box = projBZ_box.to_frb(the_radius,[128,128],the_center)
                
                #pdb.set_trace()
                length = len(frbX_box['gas','density'])

                
                if 0:
                    # APPLY GAUSS BEAM
                    #frb.apply_gauss_beam(nbeam=25,sigma=2.0)  #YT version didn't work, the following takes field and sigma 
                    density_gx = gaussian_filter(frbX_can['gas','density'], 2)
                    density_gy = gaussian_filter(frbY_can['gas','density'], 2)
                    density_gz = gaussian_filter(frbZ_can['gas','density'], 2)
                    b_gx = gaussian_filter(frbBX_can[Bx], 2)
                    b_gy = gaussian_filter(frbBY_can[By], 2)
                    b_gz = gaussian_filter(frbBZ_can[Bz], 2)

                if 1:  # SHOW IMAGES    # BTW: more vim tricks: c i (, d i {
                    plt.close('all')
                    fig, ax = plt.subplots(1,2)
                    ax[0].imshow(frbZ_can['gas','density'].v)  
                    ax[1].imshow(frbZ_box['gas','density'].v)  
                    plt.savefig("frb_canVSbox_core%d_%s.png"%(core_id,sim))

            # FRB SYNTH FIELDS
            if 0:
                self.synthRhox_can_frb[nc] = (np.array(frbX_can['gas','density'])).sum()/length 
                self.synthRhoy_can_frb[nc] = (np.array(frbY_can['gas','density'])).sum()/length 
                self.synthRhoz_can_frb[nc] = (np.array(frbZ_can['gas','density'])).sum()/length  

                self.synthRhox_box_frb[nc] = (np.array(frbX_box['gas','density'])).sum()/length 
                self.synthRhoy_box_frb[nc] = (np.array(frbY_box['gas','density'])).sum()/length 
                self.synthRhoz_box_frb[nc] = (np.array(frbZ_box['gas','density'])).sum()/length  

                self.synthBx_can_frb[nc] = (np.array(frbBX_can[Bx])).sum()/length
                self.synthBy_can_frb[nc] = (np.array(frbBY_can[By])).sum()/length
                self.synthBz_can_frb[nc] = (np.array(frbBZ_can[Bz])).sum()/length

                self.synthBx_box_frb[nc] = (np.array(frbBX_box[Bx])).sum()/length
                self.synthBy_box_frb[nc] = (np.array(frbBY_box[By])).sum()/length
                self.synthBz_box_frb[nc] = (np.array(frbBZ_box[Bz])).sum()/length

            # FRB SYNTH FIELDS: GAUSSIAN FILTERED (EDIT)
            if 0:
                self.synthRhox_Gauss[nc] =(np.array(density_gx)).sum()/length
                self.synthRhoy_Gauss[nc] =(np.array(density_gy)).sum()/length
                self.synthRhoz_Gauss[nc] =(np.array(density_gz)).sum()/length
                self.synthBx_Gauss[nc] =(np.array(b_gx)).sum()/length
                self.synthBy_Gauss[nc] =(np.array(b_gy)).sum()/length
                self.synthBz_Gauss[nc] =(np.array(b_gz)).sum()/length


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

if 1:  # FOR BOXES
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

if 0:  # FOR CANS ONLY
    rinfs = [0,2,3]
    
simnames = ['u601','u602', 'u603']
atf = {}
for nt,tool in enumerate([scope1,scope2,scope3]): 
    # WHICH CORES and pass them
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    #core_list = TL6.loops['u601'].core_by_mode['Binary']
    core_list = all_cores[2:4]  #DEBUG
    #core_list = all_cores

    tool.qtyRun(nt,rinfs,core_list=core_list) 

    fig,ax = plt.subplots(1,1)
    Rhox = tool.synthRhox_can
    Rhoy = tool.synthRhoy_can
    Rhoz = tool.synthRhoz_can
    Rhox_bf = tool.synthRhox_box_frb
    Rhoy_bf = tool.synthRhoy_box_frb
    Rhoz_bf = tool.synthRhoz_box_frb
    Rhoz_Gauss = tool.synthRhoz_Gauss
    Rhoy_Gauss = tool.synthRhoy_Gauss
    Rhox_Gauss = tool.synthRhox_Gauss
    Rho = np.concatenate((Rhox,Rhoy,Rhoz))
    Rho_bf = np.concatenate((Rhox_bf,Rhoy_bf,Rhoz_bf))
    Rho_Gauss = np.concatenate((Rhox_Gauss,Rhoy_Gauss,Rhoz_Gauss))
    Bx = tool.synthBx_can
    By = tool.synthBy_can
    Bz = tool.synthBz_can
    Bx_bf = tool.synthBx_box_frb
    By_bf = tool.synthBy_box_frb
    Bz_bf = tool.synthBz_box_frb
    Bx_Gauss = tool.synthBx_Gauss
    By_Gauss = tool.synthBy_Gauss
    Bz_Gauss = tool.synthBz_Gauss
    Bxyz = np.concatenate((Bx,By,Bz))
    Bxyz_bf = np.concatenate((Bx_bf,By_bf,Bz_bf))
    Bxyz_Gauss = np.concatenate((Bx_Gauss,By_Gauss,Bz_Gauss))

    if 0:
        atf[nt] = []
        #low_cores[nt] = []

        RHO = np.log10(Rho_bf) 
        #print('RHO',RHO)
        BLOS = np.log10(abs(Bxyz_bf))
        #print('BLOS',BLOS)
        ok = BLOS > 1

        pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        #pfit = np.polyfit(RHO,BLOS,1) 
        alpha = pfit[0]
        BLOS_o = pfit[1]  #could use this...

        atf[nt].append(alpha)

        ax.scatter(Rho_bf,Bxyz_bf,alpha=0.4)
        ax.scatter(Rho_bf[ok],Bxyz_bf[ok],color='g',alpha=0.4)
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        #RHO_x = np.linspace(RHO.min(),RHO.max(),num=len(RHO))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)  #edited  
        ax.plot(RHO_X,BLOS_Y) 
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(r'$\alpha = %f$'%alpha)
        print('atf[nt][0]:',atf[nt][0])
        fig.savefig('BxyzVsNxyz_boxFRB_OKfit_synth_%s.png'%nt)
        plt.close(fig)

        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,atf[nt][0]))
        alphaFile.close()
        print('sims alphas ',atf[nt])
        #print('sims alphaless cores ',low_cores[nt][:])
   
# PLOT OF THE 2D,3D ALPHAS
'''
fig,ax = plt.subplots(1,1)
a3D = [0.456282,0.479462,0.617871]
a2D_cyl =[0.711603,0.713381,0.731364]
a2D_reg = [0.520261,0.446243,0.496914]
ax.plot(a3D, a2D, 'bo',alpha=0.8)
ax.plot(a2D, a2D_blurr, 'bs',alpha=0.8)
fig.savefig('alpha_2D_vs_3D')
plt.close(fig)
'''
