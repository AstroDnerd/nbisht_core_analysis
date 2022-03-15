
from starter2 import *
import data_locations as dl
import davetools
reload(davetools)


class telescope(): 
    def __init__(self,the_loop):
        self.this_loopers = loop

        self.rhoMean = np.empty([0],dtype=float)  
        self.blosMean = np.empty([0],dtype=float)

    def plotting(the_x, the_y):
        fig ax = plt.subplots(1,1)
        ax.plot(the_x, the_y) 

    def makeCyl():
        # NOTE: needs to catch 3 arguments
        # definte center, normal, radius, height..
        the_Cyl = ds.disk(center, normal, radius, height, fields=None, ds=None)
        
        # the averages, respectively, to be plotted
        BLOS_avg = (the_Cyl[self.blosMean] * the_Cyl[cv] * the_cyl[rho]).sum()/cv
        D_avg= (the_cyl[rho] * the_cyl[cv]).sum()/cv
        # NOTE: pass BLOS_avg and D_avg to plotting() 


    def qtyRun(self,core_list=None):
        thtr = self.this_looper.tr

        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        
        # loop over cores
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.particle_pos(core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            # NOTE: define a circle with ms: reREAD MINI_SCRUBBER(), understand it better!!
           

            b_los = thtr.c([core_id],'magnetic_field_z')
            rho = thtr.c([core_id],'density')
            cv = thtr.c([core_id],'cell_volume')
            cvSum = cv.sum() 

            self.rhoMean[core_id]= (rho*cv).sum(axis=0)/cvSum
            self.blosMean[core_id]= (b_los*cv).sum(axis=0)/cvSum
            # NOTE: reminder of which axis represents what 

        # loop over frames, or do this a different way
        the_frames = thtr.frames
        for nf,frame in enumerate(the_frames):
            ds = self.this_looper.load(frame)
            # NOTE: pass ds, cv, rho to makeCyl()



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
    tool.qtyRun()
# have qtyRun return things and then call the next definitions or embed definitions within qtyRun..


