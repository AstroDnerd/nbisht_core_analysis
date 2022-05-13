
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

        #self.mean_synthBlos = np.empty([0],dtype=float)
        #self.mean_synthRho =np.empty([0],dtype=float)

    def qtyRun(self,sim,core_list=None):
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        # THE FINAL FRAME 
        the_frame = thtr.frames[-1:]
        self.synthRho = np.zeros([len(core_list)])
        self.synthBlos = np.zeros([len(core_list)])

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
            the_normal = [1,0,0]

            Blos = 'magnetic_field_z'
            xax = ds.coordinates.x_axis[0]  #..[proj_axis]  
            yax = ds.coordinates.y_axis[0]
            Rx = all_particles[xax]
            Ry = all_particles[yax]
            R2d = np.sqrt(Rx**2 + Ry**2)
            radius = R2d.max()
            #the_radius = max([radius,3./128])
         
            # THE AVERAGED BLOS AND DENSITY TO BE PLOTTED
            the_Cyl = ds.disk(the_center,the_normal,radius,height=(1,'code_length'))
            cv = the_Cyl['gas','density'].sum()
            mass = the_Cyl['gas','cell_mass'].sum() 
            
            self.synthBlos[nc] = (the_Cyl['density'] * the_Cyl[Blos] * the_Cyl['cell_volume']).sum()/mass

            #self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_mass']).sum()/mass 
            self.synthRho[nc] =(the_Cyl['density'] * the_Cyl['cell_volume']).sum()/cv  

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

    core_list = all_cores[2:3]
    tool.qtyRun(nt,core_list=core_list)

    fig,ax = plt.subplots(1,1)
    Rho = tool.synthRho
    Blos = tool.synthBlos

    if 0:
        atf[nt] = []
        low_cores[nt] = []

        RHO = np.log10(Rho)  #previously added ABS...but this shouldn't be the case
        BLOS = np.log10(abs(Blos))
        ok = BLOS > 1
        pfit = np.polyfit(RHO[ok],BLOS[ok],1) 
        alpha = pfit[0]
        BLOS_o = pfit[1]  #could use this...

        atf[nt].append(alpha)
        #for i in ok: 
        #    print(RHO)
        #    print('ok',ok)
        #    if ok[i] == True:
        #        low_cores[nt].append(nf)

        ax.scatter(RHO,BLOS)
        BLOS_Y = alpha * RHO + BLOS_o 
        ax.plot(RHO,BLOS_Y) 
        fig.savefig('blos_rho_synth_%s.png'%nt)
        plt.close(fig)

        alphaFile = open("p66_brho/alphaRecords.txt",'a')
        alphaFile.write("Sim %d alpha %f \n"%(nt,atf[nt][0]))
        alphaFile.close()
        print('sims alphas ',atf[nt])
        #print('sims alphaless cores ',low_cores[nt][:])
    

