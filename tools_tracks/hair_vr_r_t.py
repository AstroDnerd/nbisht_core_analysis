
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


class more_hair():
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
        for core_id in core_list:
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

            name = self.this_looper.sim_name
            
            import colors
            fig,ax=plt.subplots(1,1,figsize=(12,12))
            ax0=ax
            #ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[0][2]
            #ax3=ax[1][0]; ax4=ax[1][1];ax5=ax[1][2]
            #ax1b = ax1.twinx()
            t = thtr.times/colors.tff
            xmin = 1/2048
            xmax = 3.5e-1
            ymin = 1e-2
            ymax = 1e6
            for iframe,frame in enumerate(frames):
                nframe = np.where( thtr.frames == frame)[0][0]
                ax0.clear()
                #ax1.clear()
                #ax1b.clear()
                #ax2.clear()
                #ax3.clear()
                #ax4.clear()
                #ax5.clear()
                c  = [0.5,0.5,0.5,0.5]
                print(pos[0].shape)
                ax0.plot( rrr[:nframe+1], vrad[:nframe+1,:], c=c,linewidth=0.2)
                #ax0.scatter( [t[nframe]]*nparticles, vrad[nframe,:], c='k',s=0.2)
                #ax0.violinplot( vrad[nframe,:], positions=[t[nframe]])
                #axbonk(ax0,xlabel='R',ylabel='Vradial')
                from mpl_toolkits import mplot3d
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                #fig,ax=plt.subplots(1,1)
                for npart in range(pos[0].shape[0]):
                    ax.plot(rrr[npart,:],vrad[npart,:],pos[2][npart,:],c=c)


                #ax.plot( ms.this_x.transpose(), ms.this_y.transpose(), c=[0.5]*4, linewidth=0.2)
                outname='plots_to_sort/hair_3d_vr_%s_c%04d_i%04d'%(name,core_id,iframe)
                fig.savefig(outname)
                print(outname)





sim_list = ['u502']
import three_loopers_u500 as TL
import close_tool
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()
        
if 0:
    for this_simname in sim_list:
        this_looper=TL.loops[this_simname]
        nt = potentool(this_looper)
        nt.run(core_list=[323])
if 0:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run( )
if 0:
    print('wut')
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run(core_list=[76], frames=[10,106])
if 1:
    this_simname = 'u501'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run(core_list=[122], frames=[106])#frames=list(range(0,106,10))+[106])
if 0:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run( frames=[106])

if 0:
    for this_simname in sim_list:
        this_looper=TL.loops[this_simname]
        #get the core with the biggest distance to the next core.
        #Thus, the max of the min for each core.
        d = ct[this_simname].distance_matrix
        d[d==0]=3 #the actual min is itself, which is zero.
        most_isolated_core = ct[this_simname].cores_used[ np.argmax( d.min(axis=0))]
        nt = potentool(this_looper)
        nt.run(core_list=[most_isolated_core])


