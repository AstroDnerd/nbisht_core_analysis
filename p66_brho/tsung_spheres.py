
from starter2 import *
import cfpack as cfp
from cfpack import stop,print

class tsungspheres():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.bmag_sph = defaultdict(list)
        self.rho_sph = defaultdict(list)
        self.tsungtff = [] 

    def run(self, core_list=None, frame_list=None, tsing=None, obs=False, whichdataset=None):
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

            alltime = thtr.times
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

            self.tsungtff.append(thtr.times[get_time_index(tsung)])

            if sum(frame_mask) == 0:
                pdb.set_trace()
            
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list)) 


            # FOR EACH FRAME IN FRAME LIST!!
            if whichdataset == 'allelse':
                for nframe,frame in enumerate(frame_list):
                    ds = this_looper.load(frame)
                    nf = np.where(this_looper.tr.frames == frame)[0][0]

                    the_center = ms.mean_center[:,-1]  #the three coords for the last frame 
                    the_radius = 1/128 
                    the_area = np.pi * (the_radius**2)
                    the_sphere = ds.sphere(the_center, the_radius)

                    dV = the_sphere[YT_cell_volume]
                    dM = the_sphere[YT_cell_mass]
                    DD = the_sphere[YT_density]
                    BB = the_sphere[YT_magnetic_field_strength]
                    BBx = the_sphere[YT_magnetic_field_x]
                    BBy = the_sphere[YT_magnetic_field_y]
                    BBz = the_sphere[YT_magnetic_field_z]
                    B = [BBx, BBy, BBz]

                    if obs == False:
                        Bmag_sph = (DD*BB*dV).sum()/dM.sum()  
                        Rho_sph = (DD*dV).sum()/dV.sum()
                        self.bmag_sph[nframe].append(Bmag_sph.v)  
                        self.rho_sph[nframe].append(Rho_sph.v)  
                    if obs == True:
                        for j in range(3): 
                            Rho_sph = (DD*dV).sum()/the_area
                            Bmag_sph = (DD*B[j]*dV).sum()/dM.sum()  
                            self.bmag_sph[nframe].append(Bmag_sph.v)  
                            self.rho_sph[nframe].append(Rho_sph.v)  

        sim = this_looper.sim_name 
        if whichdataset == 'allelse':  
            data_bmagsph = [*self.bmag_sph.values()]
            data_rhosph = [*self.rho_sph.values()]

            sim = this_looper.sim_name 
            hfivename = 'p66_brho/h5files/blosncol_sphtsung_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bmag_sph']=data_bmagsph
            Fptr['rho_sph']=data_rhosph
            print('check h5files in p66_brho/h5files')
            Fptr.close()

        if whichdataset == 'tsungtimes':  
            hfivename = 'p66_brho/h5files/coreids_tsungtff_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['tsungtff']=self.tsungtff
            print('check h5files in p66_brho/h5files')
            Fptr.close()
        
        return alltime 


