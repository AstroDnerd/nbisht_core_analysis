
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


class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop
        self.cores_used = []

    def qtyRun(self,sim,core_list=None):
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:]
        self.synthRho = np.zeros([len(core_list)])  #should be in init?
        #self.synthBx = np.zeros([len(core_list)])
        #self.synthBy = np.zeros([len(core_list)])
        self.synthBz = np.zeros([len(core_list)])

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
 
            all_particles = np.stack([ms.particle_x,ms.particle_y,ms.particle_z])  #stack along 0 axis                
            the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
            the_normalX = [1,0,0] 
            the_normalY = [0,1,0]
            the_normalZ = [0,0,1] 

            Bx = 'magnetic_field_x'
            By = 'magnetic_field_y'
            Bz = 'magnetic_field_z'
            xax = ds.coordinates.x_axis[2]  #..[proj_axis]  
            yax = ds.coordinates.y_axis[2]
            #zax = ds.coordinates.z_axis[0]  
            Rx = all_particles[xax]
            Ry = all_particles[yax]
            #Rz = all_particles[zax]
            R2d = np.sqrt(Rx**2 + Ry**2)  #find the other two
            #radius = R2d.max()
            radius = 1/128
            area = np.pi * radius**2
            areaXY = np.pi * radius**2
            #the_radius = max([radius,3./128])
         
            # THE AVERAGED BLOS AND DENSITY TO BE PLOTTED
            the_Cyl = ds.disk(the_center,the_normalZ,radius,height=(1,'code_length'))
            cv = the_Cyl['gas','cell_volume'].sum()
            mass = the_Cyl['gas','cell_mass'].sum() 
            
            self.synthBz[nc] = (the_Cyl['density'] * the_Cyl[Bz] * the_Cyl['cell_volume']).sum()/mass
            #self.synthBx[nc] = (the_Cyl[Bx] * the_Cyl['cell_volume']).sum()/cv
            #self.synthBy[nc] = (the_Cyl[By] * the_Cyl['cell_volume']).sum()/cv
            #self.synthBz[nc] = (the_Cyl[Bz] * the_Cyl['cell_volume']).sum()/cv

            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_mass']).sum()/mass 
            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/cv  

            self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/area  
            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/areaX  
            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/areaY  
            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/areaZ  

            # PROJECTIONS
            if 0:
                proj = ds.proj(('gas','density'),0,data_source = the_Cyl)  #EDIT
                pw = proj.to_pw(center=the_center,origin='domain')  #EDIT
                #pw.annotate_these_particles4(1.0) 
                pw.save('core%d_sim%d'%(core_id,sim))

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

atf = {}
low_cores = {}
for nt,tool in enumerate([scope1,scope2,scope3]):
    # Test/TEST A FEW CORES
    all_cores = np.unique(tool.this_looper.tr.core_ids)

    #core_list = all_cores[2:3]
    core_list = all_cores
    tool.qtyRun(nt,core_list=core_list)

    fig,ax = plt.subplots(1,1)
    Rho = tool.synthRho
    Bz = tool.synthBz

    if 1:
        atf[nt] = []
        low_cores[nt] = []

        RHO = np.log10(Rho)  #previously added ABS...but this shouldn't be the case
        BLOS = np.log10(abs(Bz))
        ok = BLOS > 1
        #pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        pfit = np.polyfit(RHO,BLOS,1) 
        alpha = pfit[0]
        BLOS_o = pfit[1]  #could use this...

        atf[nt].append(alpha)
        #for i in ok: 
        #    print(RHO)
        #    print('ok',ok)
        #    if ok[i] == True:
        #        low_cores[nt].append(nf)

        ax.scatter(Rho,Bz,alpha=0.4)
        ax.scatter(Rho[ok],Bz[ok],color='g',alpha=0.4)
        #RHO_x = np.linspace(RHO[ok].min(),RHO[ok].max(),num=len(RHO[ok]))
        RHO_x = np.linspace(RHO.min(),RHO.max(),num=len(RHO))
        RHO_X = 10 ** RHO_x
        BLOS_Y = 10 ** (alpha*RHO_x + BLOS_o)  #edited  
        ax.plot(RHO_X,BLOS_Y) 
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.savefig('BzVSNCV_noOK_synth_%s.png'%nt)
        plt.close(fig)

        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,atf[nt][0]))
        alphaFile.close()
        print('sims alphas ',atf[nt])
        #print('sims alphaless cores ',low_cores[nt][:])
    

