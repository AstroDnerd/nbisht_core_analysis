
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

class packing():
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
            if ms.nparticles < 2000:
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
            #pos = np.stack([ ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(), ms.particle_z[sl].transpose()])
            #pos = np.stack([ ms.particle_x[sl].flatten(),ms.particle_y[sl].flatten(), ms.particle_z[sl].flatten()])
            #pos = np.stack([ ms.particle_x[sl].flatten(),ms.particle_y[sl].flatten(), ms.particle_z[sl].flatten()])

            Nx = 128
            Nbins = 128
            ploot = np.zeros([Nbins,len(frames)])
            fig,ax=plt.subplots(1,1)
            rm = rainbow_map(len(frames))
            for nf,frame in enumerate(frames):
                iframe = np.where(thtr.frames==frame)[0][0]
                pos = np.stack([ thtr.c([core_id],'particle_position_x')[:,iframe],
                                 thtr.c([core_id],'particle_position_y')[:,iframe],
                                 thtr.c([core_id],'particle_position_z')[:,iframe]])

                I = (pos*Nx).astype('int')
                II = I[0]+Nx*(I[1]+Nx*I[2])
                U,C = np.unique(II, return_counts=True)
                #ploot[U,nf]=C

                #bins=np.mgrid[0:Ntotal:(Ntotal+1)*1j]

                #plt.clf()
                #plt.hist( II, bins=bins)
                #plt.savefig('plots_to_sort/farts%d.png'%nf)
                #print(nf)
                #plt.clf()
                #plt.scatter(U[:30],C[:30])
                #print(U.size)
                #print(np.unique(U).size)
                #print(U.max())
                #plt.clf()
                bins=np.geomspace(1,8e3,16)
                print(bins)
                #ax.hist(C, color=rm(nf), bins = , histtype='step')
                hist,bins = np.histogram(C, bins=bins)
                bin_cen = 0.5*(bins[1:]+bins[:-1])
                ax.plot(bin_cen,hist,c=rm(nf))


           
            axbonk(ax,xlabel='Nparticles/zone', ylabel='Nzones', xscale='log',yscale='log')
            plt.savefig('plots_to_sort/fartsZZZ.png')
            print('yyooooo')
            if 0:
                vel = np.stack([ ms.raw_vx[sl].transpose(),ms.raw_vy[sl].transpose(), ms.raw_vz[sl].transpose()])
                ms.get_central_at_once(core_id)
                cen_vmag = ms.cen_vmag[sl].transpose()
                rrr = ms.r[sl]+0
                rrr = rrr.transpose()
                rrr[ rrr<1./2048] = 1./2048
                den = ms.density[sl].transpose()
                dv = ms.cell_volume[sl].transpose()

if 1:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = packing(this_looper)
    nt.run(core_list=[2], frames='reg')
