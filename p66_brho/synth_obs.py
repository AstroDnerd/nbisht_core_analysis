
'''
the fixed resolution buffer may help for these purposes...
://yt-project.org/doc/visualizing/manual_plotting.html 
'''

from starter2 import *
import data_locations as dl
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
        self.synthRhoz = np.zeros([len(core_list)])
        self.synthRhoz_Gauss = np.zeros([len(core_list)])

        self.synthBx = np.zeros([len(core_list)])
        self.synthBy = np.zeros([len(core_list)])
        self.synthBz = np.zeros([len(core_list)])
        self.synthBz_Gauss = np.zeros([len(core_list)])

        # CORES
        position_dict={}
        for nc,core_id in enumerate(core_list):
            self.cores_used.append(core_id)
            ds = self.this_looper.load(the_frame[0])

            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue


            # THE PIECES FOR THE OBJECT
            all_particles = np.stack([ms.particle_x,ms.particle_y,ms.particle_z])  #stack along 0 axis                
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_radius = rinf[core_id]  #inflection[sim].rinflection[core_id]  #r_inflection.rinflection_list[core_id]  ;HOWWW

            the_left = the_center - the_radius 
            the_right = the_center + the_radius
            the_normalX = [1,0,0] 
            the_normalY = [0,1,0]
            the_normalZ = [0,0,1] 

           
            # TO OBTAIN THE AREAS
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
            the_region = ds.region(the_center,the_left,the_right)

            the_CylX = ds.disk(the_center,the_normalX,radius,height=(1,'code_length'))
            the_CylY = ds.disk(the_center,the_normalY,radius,height=(1,'code_length'))
            the_CylZ = ds.disk(the_center,the_normalZ,radius,height=(1,'code_length'))
            massX = the_CylX['gas','cell_mass'].sum() 
            massY = the_CylY['gas','cell_mass'].sum() 
            massZ = the_CylZ['gas','cell_mass'].sum() 


            # THE FIELD: cyl or region(make data_source later):  make x, y, z more efficient
            Bx = 'magnetic_field_x'
            By = 'magnetic_field_y'
            Bz = 'magnetic_field_z'
            if 0:
                self.synthBx[nc] = (the_CylX['density'] * the_CylX[Bx] * the_CylX['cell_volume']).sum()/massX
                self.synthBy[nc] = (the_CylY['density'] * the_CylY[By] * the_CylY['cell_volume']).sum()/massY
                self.synthBz[nc] = (the_CylZ['density'] * the_CylZ[Bz] * the_CylZ['cell_volume']).sum()/massZ
                print('what is this Bz',self.synthBz_Gauss[nc])
            #if 0:
                #self.synthBx[nc] = the_region['density']  #edit
 
            # THE DENSITY: cyl or region:  make x, y z more efficient
            if 0:
                self.synthRhox[nc] =(the_CylX['density'] * the_CylX['cell_volume']).sum()/area  
                self.synthRhoy[nc] =(the_CylY['density'] * the_CylY['cell_volume']).sum()/area  
                self.synthRhoz[nc] =(the_CylZ['density'] * the_CylZ['cell_volume']).sum()/area  
                print('what is this rho',self.synthRhoz[nc])

            data_source = the_region
            # PROJECTIONS and adding gaussian beam
            if 1:
                proj = ds.proj(('gas','density'),2,data_source = data_source)  #naming default, cylZ
                frb = proj.to_frb(radius,[128,128],the_center)  #width, rez=128, center
                #proj = yt.ProjectionPlot(ds, 2, ('gas','density'),data_source = the_CylZ) 
                #frb = proj.frb
                
                if 0:
                    # APPLY GAUSS BEAM
                    #frb.apply_gauss_beam(nbeam=25,sigma=2.0)  #YT version didn't work 
                    density_gauss = gaussian_filter(frb['gas','density'], 2)
                    cv_gauss = gaussian_filter(frb['gas','cell_volume'], 2)
                    bz_gauss = gaussian_filter(frb[Bz], 2)

                if 0: 
                    plt.close('all')
                    fig, ax = plt.subplots(1,2)
                    ax[0].imshow(frb['gas','density'].v)  
                    ax[1].imshow(density_gauss)  #vim tricks: c i (, d i {                
                    plt.savefig("frbfilter_theRegionDensityTest_sig2_core%d.png"%core_id)

            # NOW FOR THE GAUSSIAN FILTER DATA...
            if 0:
                self.synthRhoz_Gauss[nc] = (np.array(density_gauss) * np.array(cv_gauss)).sum()/area  
                print('what is this rho_g',self.synthRhoz_Gauss[nc])
                gauss_massZ = np.array(frb['gas','cell_mass']).sum()
                self.synthBz_Gauss[nc] = (np.array(density_gauss) * np.array(bz_gauss) * np.array(cv_gauss)).sum()/gauss_massZ
                print('what is this Bz_g',self.synthBz_Gauss[nc])

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
# AS FOUND IN tools_tracks/r_inflection

if 'inflection' not in dir():
    inflection = {}
    for sim in TL6.loops:
        inflection[sim]=r_inflection.R_INFLECTION(TL6.loops[sim])
        inflection[sim].run()
        print('after .run')


atf = {}
low_cores = {}
for nt,tool in enumerate([scope1,scope2,scope3]): 

    # WHICH CORES and pass them
    all_cores = np.unique(tool.this_looper.tr.core_ids)
    #core_list = TL6.loops['u601'].core_by_mode['Binary']
    #core_list = all_cores[2:3]
    core_list = all_cores

    tool.qtyRun(nt,inflection[nt],core_list=core_list)

    fig,ax = plt.subplots(1,1)
    Rhox = tool.synthRhox
    Rhoy = tool.synthRhoy
    Rhoz = tool.synthRhoz
    Rhoz_Gauss = tool.synthRhoz_Gauss
    Rho = np.concatenate((Rhox,Rhoy,Rhoz))
    Bx = tool.synthBx
    By = tool.synthBy
    Bz = tool.synthBz
    Bz_Gauss = tool.synthBz_Gauss
    Bxyz = np.concatenate((Bx,By,Bz))

    if 1:
        atf[nt] = []
        low_cores[nt] = []

        RHO = np.log10(Rhoz_Gauss)  #previously added ABS...but this shouldn't be the case
        print('RHO',RHO)
        BLOS = np.log10(abs(Bz_Gauss))
        print('BLOS',BLOS)
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

        ax.scatter(Rhoz_Gauss,Bz_Gauss,alpha=0.4)
        ax.scatter(Rhoz_Gauss[ok],Bz_Gauss[ok],color='g',alpha=0.4)
        RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        #RHO_x = np.linspace(RHO.min(),RHO.max(),num=len(RHO))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)  #edited  
        ax.plot(RHO_X,BLOS_Y) 
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(r'$\alpha = %f$'%alpha)
        print('atf[nt][0]:',atf[nt][0])
        fig.savefig('BzVsNcv_Gauss_OKfit_synth_%s.png'%nt)
        plt.close(fig)

        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,atf[nt][0]))
        alphaFile.close()
        print('sims alphas ',atf[nt])
        #print('sims alphaless cores ',low_cores[nt][:])
    

