
'''
the fixed resolution buffer may help for these purposes...
://yt-project.org/doc/visualizing/manual_plotting.html 

synthetic observations, version 1
first attempt at everything.
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
        self.synthRhox = np.zeros([len(core_list)]) 
        self.synthRhoy = np.zeros([len(core_list)])
        self.synthRhoz = np.zeros([len(core_list)])  #also holds the_region
        self.synthRhoz_Gauss = np.zeros([len(core_list)])
        self.synthRhoy_Gauss = np.zeros([len(core_list)])
        self.synthRhox_Gauss = np.zeros([len(core_list)])

        self.synthBx = np.zeros([len(core_list)])
        self.synthBy = np.zeros([len(core_list)])
        self.synthBz = np.zeros([len(core_list)])   #also holds the_region
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
            #if ms.nparticles < 10:
            #    continue


            # THE PIECES FOR THE OBJECT
            all_particles = np.stack([ms.particle_x,ms.particle_y,ms.particle_z])  #stack along 0 axis                
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            if 0:
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
            if 0:
                the_region = ds.region(the_center,the_left,the_right)
                the_mass = the_region['gas','cell_mass'].sum()  

            the_CylX = ds.disk(the_center,the_normalX,radius,height=(1,'code_length'))
            the_CylY = ds.disk(the_center,the_normalY,radius,height=(1,'code_length'))
            the_CylZ = ds.disk(the_center,the_normalZ,radius,height=(1,'code_length'))
            massX = the_CylX['gas','cell_mass'].sum() 
            massXb = (the_CylX['gas','density'].mean() * the_CylX['gas','cell_volume']).sum()  #mean or average...
            massY = the_CylY['gas','cell_mass'].sum() 
            massYb = (the_CylY['gas','density'].mean() * the_CylY['gas','cell_volume']).sum()
            massZ = the_CylZ['gas','cell_mass'].sum() 
            massZb = (the_CylZ['gas','density'].mean() * the_CylZ['gas','cell_volume']).sum()


            # THE FIELD: cyl or region(make data_source later):  make x, y, z more efficient
            Bx = 'magnetic_field_x'
            By = 'magnetic_field_y'
            Bz = 'magnetic_field_z'
            if 0:
                self.synthBx[nc] = (the_CylX['density'] * the_CylX[Bx] * the_CylX['cell_volume']).sum()/massX
                self.synthBy[nc] = (the_CylY['density'] * the_CylY[By] * the_CylY['cell_volume']).sum()/massY
                self.synthBz[nc] = (the_CylZ['density'] * the_CylZ[Bz] * the_CylZ['cell_volume']).sum()/massZ
                #print('Bz_Zone',self.synthBz[nc])
            if 0:
                self.synthBx[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_x'] * the_region['gas','cell_volume']).sum()/the_mass
                self.synthBy[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_y'] * the_region['gas','cell_volume']).sum()/the_mass
                self.synthBz[nc] = (the_region['gas','density'] * the_region['gas','magnetic_field_z'] * the_region['gas','cell_volume']).sum()/the_mass
 
            # THE DENSITY: cyl or region:  make x, y z more efficient
            if 0:
                self.synthRhox[nc] =(the_CylX['density'] * the_CylX['cell_volume']).sum()/area  
                self.synthRhoy[nc] =(the_CylY['density'] * the_CylY['cell_volume']).sum()/area  
                self.synthRhoz[nc] =(the_CylZ['density'] * the_CylZ['cell_volume']).sum()/area  
                print('RhoZ_Zone',self.synthRhoz[nc])
            if 0:  #... no distinction in the definitions of x,y,z
                self.synthRhox[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area 
                self.synthRhoy[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area  
                self.synthRhoz[nc] =(the_region['gas','density'] * the_region['gas','cell_volume']).sum()/the_area  
            
            if 0: 
                the_CylX = the_region
                the_CylY = the_region
                the_CylZ = the_region

            # PROJECTIONS and adding gaussian beam
            if 0:
                # ANOTHER PROJ METHOD
                #proj = yt.ProjectionPlot(ds, 2, ('gas','density'),data_source = the_CylZ) 
                #frb = proj.frb
                projX = ds.proj(('gas','density'),0,weight_field =('gas','cell_volume'),data_source = the_CylX)
                projY = ds.proj(('gas','density'),1,weight_field =('gas','cell_volume'),data_source = the_CylY)
                projZ = ds.proj(('gas','density'),2,weight_field =('gas','cell_volume'),data_source = the_CylZ)
                frbX = projX.to_frb(radius,[128,128],the_center)
                frbY = projY.to_frb(radius,[128,128],the_center) 
                frbZ = projZ.to_frb(radius,[128,128],the_center)  #width, rez, center

                projBX = ds.proj(('gas','magnetic_field_x'),0,weight_field =('gas','density'),data_source = the_CylX)
                projBY = ds.proj(('gas','magnetic_field_y'),1,weight_field =('gas','density'),data_source = the_CylY)
                projBZ = ds.proj(('gas','magnetic_field_z'),2,weight_field =('gas','density'),data_source = the_CylZ)
                frbBX = projBX.to_frb(radius,[128,128],the_center)
                frbBY = projBY.to_frb(radius,[128,128],the_center) 
                frbBZ = projBZ.to_frb(radius,[128,128],the_center)
                #pdb.set_trace()
                
                length = len(frbX['gas','density'])
                
                if 0:
                    # APPLY GAUSS BEAM
                    #frb.apply_gauss_beam(nbeam=25,sigma=2.0)  #YT version didn't work, the following takes field and sigma 
                    density_gx = gaussian_filter(frbX['gas','density'], 2)
                    density_gy = gaussian_filter(frbY['gas','density'], 2)
                    density_gz = gaussian_filter(frbZ['gas','density'], 2)
                    b_gx = gaussian_filter(frbBX[Bx], 2)
                    b_gy = gaussian_filter(frbBY[By], 2)
                    b_gz = gaussian_filter(frbBZ[Bz], 2)

                if 0: 
                    plt.close('all')
                    fig, ax = plt.subplots(1,2)
                    ax[0].imshow(frbX['gas','density'].v)  
                    ax[1].imshow(frbX['gas','density'].v)  
                    #ax[1].imshow(density_gx)  #vim tricks: c i (, d i {                
                    plt.savefig("frbfilter_regFRBsTest_sig2_core%d_%s.png"%(core_id,sim))

            # NOW FOR THE NON-GAUSSIAN FILTERED DATA...
            if 0:
                self.synthRhox[nc] = (np.array(frbX['gas','density'])).sum()/length 
                self.synthRhoy[nc] = (np.array(frbY['gas','density'])).sum()/length 
                self.synthRhoz[nc] = (np.array(frbZ['gas','density'])).sum()/length  
                print('RhoZ_FRB',self.synthRhoz[nc])
                self.synthBx[nc] = (np.array(frbBX[Bx])).sum()/length
                self.synthBy[nc] = (np.array(frbBY[By])).sum()/length
                self.synthBz[nc] = (np.array(frbBZ[Bz])).sum()/length
                #print('Bz_FRB',self.synthBz[nc])

            # NOW FOR THE GAUSSIAN FILTERED DATA...
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
    #rinf1 = r_inflection.R_INFLECTION(TL6.loops['u601'])
    #rinf_1 = rinf1.run()
if 'scope2' not in dir() or clobber:
    scope2=telescope(TL6.loops['u602'])
    #rinf2 = r_inflection.R_INFLECTION(TL6.loops['u602'])
    #rinf_2 = rinf2.run()
if 'scope3' not in dir() or clobber:
    scope3=telescope(TL6.loops['u603'])
    #rinf3 = r_inflection.R_INFLECTION(TL6.loops['u603'])
    #rinf_3 = rinf3.run()

#rinfs = [rinf_1]#, rinf_2, rinf_3]
rinfs = [0,2,3]
simnames = ['u601','u602', 'u603']

atf = {}
low_cores = {}
for nt,tool in enumerate([scope1,scope2,scope3]): 

    # WHICH CORES and pass them
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    #core_list = TL6.loops['u601'].core_by_mode['Binary']
    core_list = all_cores[2:4]
    #core_list = all_cores

    tool.qtyRun(nt,rinfs,core_list=core_list)  #add rinfs
    #print('after qtyRun')

    fig,ax = plt.subplots(1,1)
    Rhox = tool.synthRhox
    Rhoy = tool.synthRhoy
    Rhoz = tool.synthRhoz
    Rhoz_Gauss = tool.synthRhoz_Gauss
    Rhoy_Gauss = tool.synthRhoy_Gauss
    Rhox_Gauss = tool.synthRhox_Gauss
    Rho = np.concatenate((Rhox,Rhoy,Rhoz))
    Rho_Gauss = np.concatenate((Rhox_Gauss,Rhoy_Gauss,Rhoz_Gauss))
    Bx = tool.synthBx
    By = tool.synthBy
    Bz = tool.synthBz
    Bx_Gauss = tool.synthBx_Gauss
    By_Gauss = tool.synthBy_Gauss
    Bz_Gauss = tool.synthBz_Gauss
    Bxyz = np.concatenate((Bx,By,Bz))
    Bxyz_Gauss = np.concatenate((Bx_Gauss,By_Gauss,Bz_Gauss))

    if 0:
        atf[nt] = []
        low_cores[nt] = []

        RHO = np.log10(Rho) 
        #print('RHO',RHO)
        BLOS = np.log10(abs(Bxyz))
        #print('BLOS',BLOS)
        ok = BLOS > 1

        pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        #pfit = np.polyfit(RHO,BLOS,1) 
        alpha = pfit[0]
        BLOS_o = pfit[1]  #could use this...

        atf[nt].append(alpha)
        #for i in ok: 
        #    print(RHO)
        #    print('ok',ok)
        #    if ok[i] == True:
        #        low_cores[nt].append(nf)

        ax.scatter(Rho,Bxyz,alpha=0.4)
        ax.scatter(Rho[ok],Bxyz[ok],color='g',alpha=0.4)
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        #RHO_x = np.linspace(RHO.min(),RHO.max(),num=len(RHO))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)  #edited  
        ax.plot(RHO_X,BLOS_Y) 
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(r'$\alpha = %f$'%alpha)
        print('atf[nt][0]:',atf[nt][0])
        fig.savefig('BxyzVsNxyzcv_cylMassb_OKfit_synth_%s.png'%nt)
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
