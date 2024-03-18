from starter2 import *


class imager():
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
                #sph = ds.sphere(center,rsph)
                left = center - rsph.v
                right = center + rsph.v
                dx = 1/2048
                Nzones = np.round((right-left)/dx)
                region  = ds.covering_grid(4,left,Nzones)
                #region = ds.region(center,left,right)
                print(left)
                print(Nzones)

                fig,ax=plt.subplots(1,1)
                print('imshow')
                proj=np.log10(region[field].sum(axis=axis).v)
                H = [1,2,0][axis]
                V = [2,0,1][axis]
                LL = left[H]*colors.length_units_pc
                RR = right[H]*colors.length_units_pc
                BB = left[V]*colors.length_units_pc
                TT = right[V]*colors.length_units_pc


                LLL = colors.length_units_pc
                circle = plt.Circle( (center[H]*LLL, center[V]*LLL), 0.1,fill=False)
                ax.imshow(proj, origin='lower',interpolation='nearest',extent=[LL,RR,BB,TT])
                ax.add_artist(circle)
                #ax.plot([1.4,2.4],[0,1])
                
                print('save')
                fig.savefig('%s/telescope_%s_c%04d_%d_tsing_ax%d'%(plot_dir,loop.sim_name,core_id,ntime,axis))
                

import track_loader as TL
sims=['u502']
TL.load_tracks(sims)
import tsing
tsing_tool = tsing.get_tsing(TL.tracks)

if 1:
    for sim in sims:
        loop = TL.tracks[sim]
        cores = loop.core_by_mode['Alone'][2:]
        #cores = [74]
        III = imager(loop)
        axis=0
        time=[0.1,0.3,0.5,0.7,0.9,1]
        #time=[0.3]
        III.run( core_list=cores, times_tsung=time, tsing=tsing_tool[sim],axis=axis)
