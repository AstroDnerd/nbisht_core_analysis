
from starter2 import *
import three_loopers_six as TLM
import xtra_energy
this_looper = TLM.loops['u601']

from starter2 import *
import data_locations as dl
from collections import defaultdict

import pcolormesh_helper as pch
reload(pch)
import davetools
reload(davetools)

import colors
plt.close('all')


class vrad():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]

    def run(self,core_list=None,frames=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')

        thtr.sort_time()
        self.times=thtr.times/colors.tff

        if frames is None:
            frames = thtr.frames
        if frames == 'reg':
            times = thtr.times
            #assume the first nonzero time is dt_sim
            dt = times[ times>0][0]
            nframe = times/dt
            outliers = np.round(nframe) - nframe
            frame_mask = np.abs(outliers) < 1e-3
            frames = thtr.frames[frame_mask]
            times = thtr.times[frame_mask]
            self.frames=frames
            self.times=times

        tsorted = thtr.times
        self.core_list=core_list
        self.vrad_std=np.zeros([len(core_list),len(thtr.times)])
        for nc,core_id in enumerate(core_list):
            reload(trackage)
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('go ', core_id)

            this_center=ms.mean_center[:,-1]
            self.cores_used.append(core_id)

            ms.particle_pos(core_id)
            if ms.nparticles <= 1000:
                sl = slice(None)
                nparticles=ms.nparticles
            else:
                sl = slice(None,None,10)
                nparticles=ms.raw_vx[:,0][sl].size
            pos = np.stack([ ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(), ms.particle_z[sl].transpose()])
            vel = np.stack([ ms.raw_vx[sl].transpose(),ms.raw_vy[sl].transpose(), ms.raw_vz[sl].transpose()])
            ms.get_central_at_once(core_id)
            cen_vmag = ms.cen_vmag[sl].transpose()
            rrr = ms.r[sl]+0
            rrr = rrr.transpose()
            rrr[ rrr<1./2048] = 1./2048
            den = ms.density[sl].transpose()
            dv = ms.cell_volume[sl].transpose()
            vrad = ms.vr_rel[sl].transpose()

            vrad_cell_avg = (vrad*dv).sum(axis=1)/dv.sum(axis=1)
            #vrad_cell_avg.shape = (vrad_cell_avg.size,1)
            #vrad_cell_var = np.sqrt(((vrad-vrad_cell_avg)**2*dv).sum(axis=1)/dv.sum(axis=1))
            #vrad_cell_avg.shape=vrad_cell_avg.size
            #vrad_mass_avg = (den*vrad*dv).sum(axis=1)/(den*dv).sum(axis=1)
            if np.isnan(vrad_cell_avg).any():
                pdb.set_trace()
            self.vrad_std[nc,:] = vrad_cell_avg



sim_list = ['u501','u502','u503']
import three_loopers_u500 as TL
import close_tool
if 'nt' not in dir():
    nt={}
    for this_simname in sim_list:
        this_looper=TL.loops[this_simname]
        nt[this_simname] = vrad(this_looper)
        nt[this_simname].run()

if 1:
    fig,ax=plt.subplots(1,3, figsize=(12,8))
    import heat_map
    for ns,this_simname in enumerate(sim_list):
        heat_map.heat_map( nt[this_simname].vrad_std, nt[this_simname].times,ax=ax[ns])
        ax[ns].plot(nt[this_simname].times, [-1]*len(nt[this_simname].times), c='r',linewidth=1)
        ax[ns].plot(nt[this_simname].times, [0]*len(nt[this_simname].times), c='r',linewidth=1)
        axbonk(ax[ns],xlabel=r'$t/t_{ff}$',ylabel=r'\langle V_r\rangle')
    fig.savefig('plots_to_sort/vrad_heat.png')

