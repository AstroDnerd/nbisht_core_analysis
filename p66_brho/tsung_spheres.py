
from starter2 import *

class tsungspheres():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.bmag_sph = defaultdict(list)
        self.rho_sph = defaultdict(list)

    def run(self, core_list=None, frame_list=None, tsing=None, timescale=0, get_particles=False, save_sorts=False):
        self.timescale=timescale
        this_looper=self.this_looper
        thtr=this_looper.tr

        def get_time_index(time):
            index=np.argmin(np.abs( thtr.times/colors.tff-time))
            return index

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        # FOR EACH CORE!!
        for core_id in core_list:

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            ms.get_central_at_once(core_id)  #ask?

            frame_mask = np.zeros_like(thtr.times, dtype='bool')

            tsung = tsing.tend_core[core_id]
            frame_mask[get_time_index(0.125*tsung)]=True  
            frame_mask[get_time_index(0.250*tsung)]=True  
            frame_mask[get_time_index(0.375*tsung)]=True  
            frame_mask[get_time_index(0.5*tsung)]=True  
            frame_mask[get_time_index(0.625*tsung)]=True  
            frame_mask[get_time_index(0.75*tsung)]=True  
            frame_mask[get_time_index(0.875*tsung)]=True  
            frame_mask[get_time_index(tsung)]=True  

            if sum(frame_mask) == 0:
                pdb.set_trace()
            
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))


            # FOR EACH FRAME IN FRAME LIST!!
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                nf = np.where(this_looper.tr.frames == frame)[0][0]

                the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
                the_radius = 1/128  
                the_sphere = ds.sphere(the_center, the_radius)

                dV = the_sphere[YT_cell_volume]
                dM = the_sphere[YT_cell_mass]
                DD = the_sphere[YT_density]
                BB = the_sphere[YT_magnetic_field_strength]

                Bmag_sph = (DD*BB*dV).sum()/dM.sum()  
                Rho_sph = (DD*dV).sum()/dV.sum()

                self.bmag_sph[nframe].append(Bmag_sph.v)  
                self.rho_sph[nframe].append(Rho_sph.v)  

        data_bmagsph = [*self.bmag_sph.values()]
        data_rhosph = [*self.rho_sph.values()]
        if 0:  
            sim = this_looper.sim_name 
            hfivename = 'p66_brho/h5files/brho_sphtsung_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bmag_sph']=data_bmagsph
            Fptr['rho_sph']=data_rhosph

            print('check h5files in p66_brho/h5files')
            Fptr.close()


