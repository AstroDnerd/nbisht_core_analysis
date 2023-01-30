from starter2 import *

import find_other_cores
reload(find_other_cores)


import three_loopers_u500 as TL
sim_list = ['u501','u502','u503']
sim_list = ['u502']
if 'mini_scrubbers' not in dir():
    print('mini scrubbers.')
    mini_scrubbers={}
    for sim in sim_list:
        this_looper=TL.loops[sim]
        mini_scrubbers[sim]={}
        all_cores=np.unique(this_looper.tr.core_ids)
        thtr=this_looper.tr

        for core_id in all_cores:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[sim][core_id]=ms

for sim in sim_list:
    reload(find_other_cores)
    this_looper = TL.loops[sim]
    #main_core_list = np.unique(TL.loops[sim].tr.core_ids)
    #main_core_list.sort()
    fig,axes=plt.subplots(2,3,figsize=(12,8))
    fig2,axes2=plt.subplots(1,1)
    fig3,axes3=plt.subplots(1,1)
    ext=extents()
    ext2=extents()
    ext3=extents()
    loost=['Alone', 'Binary','Cluster']
    loost=['Alone']
    for nmode,mode in enumerate(loost):
        main_core_list=TL.loops[sim].core_by_mode[mode]
        #main_core_list=[76,74, 32]
        print('get')
        #main_core_list = np.unique(TL.loops[sim].tr.core_ids)
        min_d,cores_used, distance= find_other_cores.get_all_distances( this_looper, main_core_list, mini_scrubbers[sim])
        print('got')
        no="""
        core_list=TL.loops[sim].core_by_mode[mode]
        min_d[min_d<0]=4000
        min_for_each = min_d.min(axis=1)
        args = np.argsort(min_for_each)

        sorted_mins = min_for_each[args]
        sorted_cores = nar(main_core_list)[args]
        sorted_loc = nar(range(len(main_core_list)))[args]

        for line in min_d:
            lines = line[line<4000]
            lines.sort()
            ext(lines)
            axes[0][nmode].plot( lines)
        nth_distance=distance[:,:,0]
        nth_distance.sort(axis=1)
        nth= nth_distance[:,1:].transpose()
        ext2(nth)
        axes[1][nmode].plot(nth)

        distance1=distance[:,:,0]
        distance1.sort(axis=1)
        distance1=distance1[:,1:]
        distance2=distance[:,:,-1]
        distance2.sort(axis=1)
        distance2=distance1[:,1:]

        the_x,the_y=distance2[:,0], distance2[:,1]
        xlab='t=-1 closest'
        ylab='t=-1 2nd'
        axes2.scatter( the_x,the_y, c='rgb'[nmode])
        axes2.plot( the_x,the_x,c='k')
        ext3(the_x)
        ext3(the_y)

        if 1:
            d1 = distance[:,:,-1]
            min_end = d1.min(axis=1)
            ht = hull_by_frame[sim][0]
            mode_overlap = np.zeros_like(d1)
            total_overlap = np.zeros(len(main_core_list))
            for nc,core_id in enumerate(main_core_list):
                over_index = np.where(ht.cores_used==core_id)[0][0]
                total_overlap[nc] = overlaps[sim][over_index,:,0].sum()


            axes3.scatter(min_end, total_overlap, c='rgb'[nmode])
        """

        #axes[1][nmode].hist( np.log10(first_distance))
#   for nmode in [0,1,2]:
#       axes[0][nmode].set(yscale='log',xlabel='core',ylabel='min ever', ylim=ext.minmax, xlim=[0,10])
#       axes[1][nmode].set(yscale='log',xlabel='core',ylabel='Delta t=0', ylim=ext2.minmax, xlim=[0,10])
#   fig.savefig('plots_to_sort/min_d')
#   axes2.set(xscale='log',xlabel=xlab,yscale='log',ylabel=ylab,xlim=ext3.minmax,ylim=ext3.minmax)
#   fig2.savefig('plots_to_sort/distance_distance')
#   fig3.savefig('plots_to_sort/distance_overlap')


    if 0:
        fig,ax=plt.subplots(1,2)
        ax0=ax[0];ax1=ax[1]
        ax1.plot( min_for_each[args])
        ax1.set(yscale='log')
        ax0.plot(min_for_each[args],marker='*')
        #for nc,core_id in enumerate(main_core_list):
        #    Ncores,b=min_d.shape
        #    D = min_d[nc,:]
        #    #ax.scatter([nc]*(Ncores-1), D[D>0])
        #    #ax0.plot( sorted(D[D>0])[:10])
        ax0.set(yscale='log',xlabel='core',ylabel='min D')
        fig.savefig('plots_to_sort/min_dist')




