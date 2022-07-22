

from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

class atool():
    def __init__(self,this_looper):
        self.this_looper=this_looper

    def run(self,core_list=None):
        this_looper=self.this_looper

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        name = this_looper.sim_name
        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        times=thtr.times[mask]+0 #the zero makes a copy
        times.shape=times.size,1
        times=times/colors.tff
        G = colors.G
        #gx = thtr.track_dict['grav_x']
        #gy = thtr.track_dict['grav_y']
        #gz = thtr.track_dict['grav_z']
        #GE2 = -1/(8*np.pi)*(gx*gx+gy*gy+gz*gz)
        #ge_min=GE2.min()
        #ge_max=GE2.max()
        self.Pearson_GE_r = np.zeros([len(core_list), len(times)])
        self.P_GE_r = np.zeros([len(core_list), len(times)])
        self.Pearson_rho_r = np.zeros([len(core_list), len(times)])
        self.Alpha_rho_r = np.zeros([len(core_list), len(times)])
        self.Rhonot_rho_r = np.zeros([len(core_list), len(times)])
        self.PeakRho  = np.zeros([len(core_list), len(times)])
        for nc, core_id in enumerate(core_list):
            print('pearson %s %d'%(name,core_id))

                
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            #ms.particle_pos(core_id)

            if True:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4

            rho = ms.density[sl]
            rho = rho[:,mask]

            self.PeakRho[nc,:]=rho.max(axis=0)

            gx = thtr.c([core_id],'grav_x')[sl][:,mask]
            gy = thtr.c([core_id],'grav_y')[sl][:,mask]
            gz = thtr.c([core_id],'grav_z')[sl][:,mask]
            GE2 = 1/(8*np.pi*G)*(gx*gx+gy*gy+gz*gz)

            RRR = ms.r[sl][:,mask]
            for n in range(GE2.shape[1]):
                ok = RRR[:,n]>0
                if ok.sum() <3:
                    continue
                the_x=np.log(RRR[:,n][ok])
                the_y=np.log(GE2[:,n][ok])
                #the_y=rho[:,n]
                r,p=scipy.stats.pearsonr(the_x,the_y)
                self.Pearson_GE_r[nc,n]=r
                self.P_GE_r[nc,n]=p
                the_y=np.log(rho[:,n])
                r,p=scipy.stats.pearsonr(the_x,the_y)
                self.Pearson_rho_r[nc,n]=r


                pfit = np.polyfit( the_x, the_y, 2)
                self.Alpha_rho_r[nc,n]=pfit[0]
                self.Rhonot_rho_r[nc,n]=pfit[1]
                

