from starter2 import *
import scipy.stats


class radial():
    def __init__(self,loop):
        self.loop=loop

    def run(self,core_list=None, times_tsung=None,tsing=None, axis=0,field='density'):
        sim_name = loop.sim_name
        thtr=loop.tr
        if core_list is None:

            core_list=np.unique(loop.tr.core_ids)

        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        for nc,core_id in enumerate(core_list):
            print('Plot core',core_id)
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=False)
            tsung = tsing.tend_core[core_id]
            if times_tsung == -1:
                times = [tsing.tsing_core[core_id], tsing.tend_core[core_id]]
            else:
                times = nar(times_tsung)*tsung
            print(times)
            nrow=1
            ncol=len(times)
            fig,ax=plt.subplots(nrow, ncol, figsize=(8,4))
            for ntime,time in enumerate(times):
                print('T',time)
                nf = get_time_index(time)
                frame = thtr.frames[nf]
                print('ds')
                ds = loop.load(frame)
                print('cg')
                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                MaxRadius=msR[:,nf].max()
                Radius = max([0.5/colors.length_units_pc, MaxRadius])
                Radius = min([Radius, 0.25])
                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)


                dv = sph[YT_cell_volume]
                RR = sph['radius']
                DD = sph[YT_density]
                ORDER = np.argsort( RR)
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                V_cuml = np.cumsum( dv[ORDER])
                rho_cuml = M_cuml/V_cuml

                rbins = np.geomspace(RR.min(),RR.max(),64)
                mean, bins, counts = scipy.stats.binned_statistic(RR,DD,statistic='mean',bins=rbins)
                rcen=0.5*(rbins[:-1]+bins[1:])

                ax[ntime].plot(rcen,mean,c='r')
                ax[ntime].scatter(sph['radius'],sph['density'])
                ax[ntime].plot(RR[ORDER], rho_cuml,c='k')
                ax[ntime].set(xscale='log',yscale='log', title='%.1f tsung'%(time/tsung))
            fig.savefig('%s/radial_%s_c%04d'%(plot_dir,loop.sim_name,core_id))


                    
                

import track_loader as TL
sims=['u502']
TL.load_tracks(sims)
import tsing
tsing_tool = tsing.get_tsing(TL.tracks)

if 1:
    for sim in sims:
        loop = TL.tracks[sim]
        cores = loop.core_by_mode['Alone'][2:]
        cores = [74]
        III = radial(loop)
        axis=0
        time=[0.1,0.3,0.5,0.7,0.9,1]
        time = [0.5, 0.9, 1]
        #time=[0.3]
        III.run( core_list=cores, times_tsung=time, tsing=tsing_tool[sim],axis=axis)
