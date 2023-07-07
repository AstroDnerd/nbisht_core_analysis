
from starter2 import *
import xtra_energy

from scipy.optimize import curve_fit
from scipy import stats
from scipy.ndimage import gaussian_filter
import core_proj_three
reload(core_proj_three)
import other_scrubber
reload(other_scrubber)
#import three_loopers_six as TL
import camera_path

class multipro():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.profiles_gas=defaultdict(list) #{}
        self.profiles_part=defaultdict(list) #{}

    def run(self, core_list=None, frame_list=None, tsing=None, timescale=0, get_particles=False, save_sorts=False):
        self.timescale=timescale
        this_looper=self.this_looper
        thtr=this_looper.tr

        if get_particles:
            suites_to_use=[1]
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin(np.abs( thtr.times/colors.tff-time))
            return index
        Nplots = 4
        if timescale in [2]:
            Ntimes = 4
        else:
            Ntimes = 6
        ext = [extents() for n in range(Nplots+1)]


        # FOR EACH CORE!!
        for core_id in core_list:
            self.profiles_gas[core_id]={}
            self.profiles_part[core_id]={}
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            ms.get_central_at_once(core_id)  #ask?
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')

            if self.timescale==0:
                frame_mask[0]=True
                frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.25 t_{\rm{sing}}$', r'$t=0.5 t_{\rm{sing}}$', r'$t=0.75 t_{\rm{sing}}$',\
                    r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']

            if self.timescale==2:
                frame_mask[0]=True
                #frame_mask[get_time_index(0.25*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.5*tsing.tsing_core[core_id])]=True
                #frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                #half_collapse = 0.5*(tsing.tsing_core[core_id]+tsing.tend_core[core_id])
                #theta=0.5
                #half_collapse = theta*tsing.tsing_core[core_id]+(1-theta)*tsing.tend_core[core_id]
                #frame_mask[get_time_index(half_collapse)]=True
                self.titles=[ r'$t=0$', r'$t=0.5 t_{\rm{sing}}$', r'$t=t_{\rm{sing}}$', r'$t=t_{\rm{sung}}$']

            if sum(frame_mask) == 0:
                pdb.set_trace()
            
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))


            # FOR EACH FRAME IN FRAME LIST!!
            for nframe,frame in enumerate(frame_list):
                self.profiles_gas[core_id][frame]={}
                self.profiles_part[core_id][frame]={}
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                nf = np.where( this_looper.tr.frames == frame)[0][0]


                for suite in suites_to_use:
                    #collector={}
                    collector=defaultdict(list)
                    the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
                    the_normal = [[1,0,0],[0,1,0],[0,0,1]]
                    the_radius_one = 1/128  
                    the_area_one = np.pi * (the_radius_one**2) 
                    the_midcyl_one = {}
                    xyz = [0,1,2]

                    if suite == 0:
                        for i in range(3):  #get one direction right first
                            the_midcyl_one[i] = ds.disk(the_center,the_normal[i],the_radius_one,height=the_radius_one) 

                            dv = the_midcyl_one[i][YT_cell_volume]
                            DD = the_midcyl_one[i][YT_density]
                            BX = the_midcyl_one[i][YT_magnetic_field_x]

                            b_los_mcyl_one = (DD*BX*dv).sum()/(DD*dv).sum()  #note I have cell_mass in other similar definitions. 
                            columnrho_mcyl_one = (DD*dv).sum()/the_area_one

                            collector['B_los_r1'].append(b_los_mcyl_one.v)  
                            collector['N_r1'].append(columnrho_mcyl_one.v)
                            collector['B_los_r1log'].append(nar(np.log10(abs(b_los_mcyl_one.v))))  
                            collector['N_r1log'].append(nar(np.log10(abs(columnrho_mcyl_one.v))))

                    if suite == 1:  
                        mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nframe) 
                        cell_volume = thtr.c([core_id],'cell_volume')[mask,nframe]
                        density = thtr.c([core_id],'density')[mask,nframe]
                        bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                        by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                        bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                        b = [bx, by, bz]
                        for j in range(3):
                            bparticles = (density * b[j] * cell_volume).sum()/(density * cell_volume).sum()  #but a heads up that this is just for 1 axis! 
                            rhoave = (density * cell_volume).sum()/cell_volume.sum()

                            collector['B'].append(bparticles)
                            collector['Rho'].append(rhoave) 
                            collector['Blog'].append(nar(np.log10(abs(bparticles))))
                            collector['Rholog'].append(nar(np.log10(abs(rhoave)))) 


                    if save_sorts:
                        print('maybe for later')
                    if suite==0:
                        self.profiles_gas[core_id][frame]=collector
                    elif suite==1:
                        self.profiles_part[core_id][frame]=collector
        
        if 0:
            sim = this_looper.sim_name 
            hfivename = 'p66_brho/b_gas_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['prof_gas']=self.profiles_gas #TypeError: Object dtype dtype('O') has no native HDF5 equivalent
            Fptr.close()


