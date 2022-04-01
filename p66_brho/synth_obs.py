
'''
the fixed resolution buffer may help for these purposes...
://yt-project.org/doc/visualizing/manual_plotting.html 
'''

from starter2 import *
import data_locations as dl
import davetools
reload(davetools)


class telescope(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop

        self.mean_synthRho = defaultdict(list)   
        self.mean_synthBlos = defaultdict(list)


    def qtyRun(self,core_list=None):
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()
        tsorted = thtr.times

        self.synthRho = np.zeros([len(core_list),len(tsorted)])
        self.synthBlos = np.zeros([len(core_list),len(tsorted)])

        # CORES
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
 
            the_frames = thtr.frames
            #the_frames = [60]   #test single frames with np.where 
            #for nf,frame in enumerate(the_frames):   #the analogy to doing the following
            for frame in the_frames:
                nf = np.where(thtr.frames == frame)[0][0]
                #print('frame ',frame)
                #print('nf ',nf)
                ds = self.this_looper.load(frame)

                all_particles = np.stack([ms.particle_x,ms.particle_y,ms.particle_z])  #stack along 0 axis
                the_center = ms.mean_center[:,nf]  #OR ms.mean_x...y,z 
                the_normal = [1,0,0]

                Blos = 'magnetic_field_z'
                xax = ds.coordinates.x_axis[0]  #..[proj_axis]  
                yax = ds.coordinates.y_axis[0]
                Rx = all_particles[xax]
                Ry = all_particles[yax]
                R2d = Rx**2 + Ry**2
                the_radius = R2d.max()
             
                # THE AVERAGED BLOS AND DENSITY TO BE PLOTTED
                the_Cyl = ds.disk(the_center,the_normal,the_radius,height=(1,'code_length'))
                cv = the_Cyl['gas','cell_volume'].sum()
                #rho = the_Cyl['gas','density'].sum() 
                
                self.mean_synthBlos[core_id] = (the_Cyl['density'] * the_Cyl[Blos] * the_Cyl['cell_volume']).sum()/cv 
                self.synthBlos[nc,nf] = self.mean_synthBlos[core_id]
                #print('synthBlos[nc,:] ',self.synthBlos[nc,:])
                self.mean_synthRho[core_id]= (the_Cyl['density'] * the_Cyl['cell_volume']).sum()/cv
                self.synthRho[nc,nf] = self.mean_synthRho[core_id]
                #print('synthRho[nc,:] ',self.synthRho[nc,:])



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
for nt,tool in enumerate([scope1]):#,scope2,scope3]):
    #core_list = [27]
    #core_list = core_list 
    tool.qtyRun()

    fig,ax = plt.subplots(1,1)
    Rho = tool.synthRho
    Blos = tool.synthBlos

    atf[nt] = []
    for nf in range(Rho.shape[1]):
        RHO = np.log10(Rho[:,nf])
        #print("RHO ",RHO)
        BLOS = np.log10(abs(Blos[:,nf]))
        #print("BLOS ",BLOS)
        pfit = np.polyfit(RHO,BLOS,1)
        alpha = pfit[0]
        #print('time ',nf)
        #print('ALPHA ',alpha)
        BLOS_o = pfit[1]  #could use this...

        atf[nt].append(alpha)
        ax.scatter(RHO,BLOS)
        #BLOS_Y = alpha * RHO + BLOS_o 
        #ax.plot(RHO,BLOS_Y) 
        fig.savefig('blos_rho_synth_%s.png'%nf)
    print('sims alphas ',atf[nt][:-1])
    

