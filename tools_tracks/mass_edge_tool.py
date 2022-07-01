from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

from collections import defaultdict

from scipy.interpolate import interp1d
class edge():
    def __init__(self, this_looper):
        self.this_looper=this_looper
        self.edge={}

    def run(self,core_list=None):
        this_looper=self.this_looper
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        frame = this_looper.target_frame
        ds = this_looper.load(frame)
        G = ds['GravitationalConstant']/(4*np.pi)
        xtra_energy.add_energies(ds)
        for core_id in core_list:
            print('Mass Edge %s %d'%(this_looper.sim_name,core_id))

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            c = nar([ms.mean_x[-1], ms.mean_y[-1],ms.mean_z[-1]])
            
            R_SPHERE = 4/128
            rsph = ds.arr(R_SPHERE,'code_length')
            sp = ds.sphere(c,rsph)

            GE = np.abs(sp['grav_energy'])
            dv = np.abs(sp['cell_volume'])
            RR = sp['radius']
            DD = sp['density']

            R_KEEP = R_SPHERE #self.rinflection[core_id]

            #get the zones in the sphere that are
            #within R_KEEP
            ok_fit = np.logical_and(RR  < R_KEEP, GE>0)

            rok=RR[ok_fit].v

            ORDER = np.argsort(rok)
            rho_o = DD[ok_fit][ORDER].v
            ge_o = GE[ok_fit][ORDER].v
            dv_o  = dv[ok_fit][ORDER].v
            rr_o_full  = RR[ok_fit][ORDER].v
            mass_r = (rho_o*dv_o).cumsum()
            enrg_r = (ge_o*dv_o).cumsum()

            rr_o = np.linspace( max([1/2048, rr_o_full.min()]), rr_o_full.max(), 128)
            all_r,all_m=rr_o_full[1:], mass_r[1:]

            my_r = np.linspace(max([1/2048, all_r.min()]),all_r.max(),1024)
            mbins = np.linspace( all_m.min(), all_m.max(), 128)
            thing, mask = np.unique( all_r, return_index=True)
              
            mfunc = interp1d( all_r[mask], all_m[mask])
            my_m = gaussian_filter(mfunc(my_r),2)
            dm=(my_m[1:]- my_m[:-1])
            mm =0.5*(my_m[1:]+my_m[:-1])
            dr=(my_r[1:]-my_r[:-1])
            dm_dr = dm/dr
            rbins = 0.5*(my_r[1:]+my_r[:-1])

            SWITCH=2*rbins*dm_dr/mm 
            findit=np.logical_and( rbins>0.01, SWITCH < 1)
            if findit.sum() > 0:
                ok = np.where(findit)[0][0]
            self.edge[core_id]=my_r[ok-1]
