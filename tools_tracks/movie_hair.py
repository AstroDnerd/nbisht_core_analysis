from starter2 import *

class flow():
    def __init__(self, loop):
        self.loop=loop
        self.cores_used=[]
        
    def run(self,core_list=None, frames=None, axis=2, external_ax=None, external_fig=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.loop.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

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

        rmcore=rainbow_map(len(core_list))
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            ms.particle_pos(core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            self.cores_used.append(core_id)

            axis_to_actually_plot=[axis]
            if external_ax is None:
                fig,ax=plt.subplots(1,1, figsize=(12,8))
            else:
                ax=external_ax
                fig=external_fig
                ax.clear()
            for LOS in axis_to_actually_plot:
                x = [1,0,1][LOS] # Using [1,0,1] and [2,2,0] 
                y = [2,2,0][LOS] # unfolds nicely.
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]

                all_x,all_y,all_z=ms.particle_x,ms.particle_y, ms.particle_z
                all_p = [all_x,all_y,all_z]

                for nff,n2 in enumerate(frames):
                    n1 = np.where(thtr.frames==n2)[0][0]
                    ax.clear()
                    XX,YY= all_p[x].transpose(), all_p[y].transpose()


                    ax.scatter(XX[0,:],YY[0,:], c='k', zorder=2)
                    ax.plot(XX[:n1+1,:],YY[:n1+1,:], c=[0.5]*4, zorder=7, linewidth=0.1)
                    ax.scatter(XX[n1,:],YY[n1,:], c='r', s=1,zorder=1)

                    outname='plots_to_sort/movie_hair_%s_c%04d_%04d.png'%(self.loop.sim_name,core_id,n2)
                    fig.savefig(outname)
                    print('save',outname)

                
if 0:
    import three_loopers_six as TLX
    thing = flow(TLX.loops['u603'])
    core_list=[76]
    thing.run(core_list=core_list)
if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u503'])
    core_list=[76]
    thing.run(core_list=core_list, frames='reg')
if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u503'])
    core_list=[76]
    thing.run(core_list=core_list, frames='reg')
if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u501'])
    core_list=[24]
    thing.run(core_list=core_list, frames='reg')
if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u501'])
    core_list=[24]
    thing.run(core_list=core_list, frames='reg', axis=0)
if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u503'])
    core_list=[16]
    thing.run(core_list=core_list, frames='reg')

if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u503'])
    core_list=[26]
    core_list=[24]
    thing.run(core_list=core_list, frames='reg')

if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u503'])
    core_list=[55]
    #thing.run(core_list=core_list, frames='reg',axis=0)
    thing.run(core_list=core_list, frames='reg',axis=2)

if 0:
    import three_loopers_u500 as TL5
    thing = flow(TL5.loops['u501'])
    core_list=[323]
    #thing.run(core_list=core_list, frames='reg',axis=0)
    thing.run(core_list=core_list, frames='reg',axis=0)
