from starter2 import *

import find_other_cores
reload(find_other_cores)


#import three_loopers_u500 as TL
#sim_list = ['u501','u502','u503']
import three_loopers_six as TL
sim_list = ['u603']
#sim_list = ['u601','u602','u603']

import anne
anne.make_hulls()
if 'mini_scrubbers' not in dir():
    mini_scrubbers={}

if 1:
    for sim in sim_list:
        if sim in mini_scrubbers:
            continue
        print('mini scrubbers.',sim)
        this_looper=TL.loops[sim]
        mini_scrubbers[sim]={}
        all_cores=np.unique(this_looper.tr.core_ids)
        thtr=this_looper.tr

        for core_id in all_cores:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[sim][core_id]=ms

if 'overlaps' not in dir():
    overlaps={}
    ratios={}

if 1:
    for sim in sim_list:
        ht0=anne.ht[sim]
        loop=TL.loops[sim]
        frame_list=loop.tr.frames
        #frame_list=[0]
        ncores=len(ht0.cores_used)
        nframes=len(frame_list)

        nframes=1
        nframe=0
        if sim not in overlaps:
            overlaps[sim] = np.zeros([ncores,ncores,nframes])
            ratios[sim] = np.zeros([ncores,ncores,nframes])
            htool = anne.ht[sim]
            for nc1, core_id_1 in enumerate(htool.cores_used):
                for nc2, core_id_2 in enumerate(htool.cores_used):
                    if nc1==nc2:
                        continue
                    a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                    #please come back to think about this logic.
                    #david, 2022-09-22
                    #2023-01-06 I think it should be zero of either is zero.
                    #The only issue is full containment and zero volume.
                    ratio = max([a,b])
                    rat= sorted( [a,b])
                    if rat[1] == 0:
                        #ratio = max([a,b])
                        ratio = 0
                    else:
                        ratio=rat[0]/rat[1]
                    ratios[sim][nc1,nc2,nframe] = ratio
                    overlaps[sim][nc1,nc2,nframe] = a


def add_cores_to_group(core_id,group=None,buddies=None,all_cores=None):
    group.append(core_id)
    for core_2 in buddies[core_id]:
        if core_2 in all_cores:
            all_cores.remove( core_2)
            add_cores_to_group(core_2, group=group, buddies=buddies,all_cores=all_cores)

if 1:
    #Mode Grader.
    Delta_Alone = 0.025
    new_mode={}
    new_dict={}
    buddy_list={}
    for sim in sim_list:
        new_mode[sim]={'Alone':[],'Binary':[],'Cluster':[]}
        new_dict[sim]={}
        buddy_list[sim]={}
        this_looper=TL.loops[sim]
        main_core_list=np.unique(this_looper.tr.core_ids)
        min_d,cores_used, distance= find_other_cores.get_all_distances( this_looper, main_core_list, mini_scrubbers[sim])
        distance[ distance==0] = 4000
        dmin = distance.min(axis=2)

        for ncore,core_id in enumerate(main_core_list):
            my_args = np.argsort( dmin[ncore])
            my_min = dmin[ncore][my_args]

            near_cores = my_min < Delta_Alone
            if near_cores.sum() <1:
                my_mode = 'Alone'
            else:
                my_buddies = cores_used[my_args][near_cores]
                buddy_list[sim][core_id] = my_buddies
                my_mode = 'Cluster'
            new_mode[sim][my_mode].append(core_id)

        #all_cores = copy.copy(list(main_core_list))
        A_cores = copy.copy(new_mode[sim]['Alone'])
        #B_cores = copy.copy(new_mode[sim]['Binary'])
        C_cores = set(new_mode[sim]['Cluster'])
        fptr=open('browser_data/CoreMode_%s.tsv'%sim, 'w')
        groups={}
        group_id=0
        new_modes = []
        while len(C_cores):
            groups[group_id]=[]
            core_id=C_cores.pop()
            add_cores_to_group(core_id,groups[group_id], buddy_list[sim],C_cores)
            group_id+=1
        #this sucks, do it better
        B_number=0
        C_number=0
        group_name={}
        for group_id in sorted(groups.keys()):
            if len(groups[group_id])==2:
                group_name[group_id] = "B%d"%B_number
                B_number+=1
            else:
                group_name[group_id] = "C%d"%C_number
                C_number+=1
        for core_id in main_core_list:
            name=None
            for group_id in sorted(groups.keys()):
                if core_id in groups[group_id]:
                    name = group_name[group_id]
            if name is None:
                if core_id not in A_cores:
                    print("serious error.")
                    raise
                name = 'Alone'
            new_modes.append(name)

        print(new_modes)
        fptr=h5py.File("browser_data/core_formation_mode_new_%s.h5"%sim,'w')
        try:
            fptr['modes']=new_modes
            fptr['core_ids']=main_core_list
        except:
            raise
        finally:
            fptr.close()

#        try:
#            for core_id in A_cores:
#                fptr.write( '%s_c%04d %s\n'%(sim,core_id,'Alone'))
#                #all_cores.pop( all_cores.index( core_id))
#                A_cores.pop( A_cores.index( core_id))
#            group_id = 0
#            while len(C_cores):
#                core_id = C_cores.pop(0)
#                print("Group %d"%group_id)
#                groups[group_id]=[core_id]
#
#                for core_2 in buddy_list[sim][core_id]:
#                    groups[group_id].append(core_2)
#                    fptr.write( '%s_c%04d %s\n'%(sim,core_2, 'B%d'%BinaryNumber))
#                    B_cores.pop( B_cores.index(core_2))
#                    #all_cores.pop( all_cores.index(core_2))
#                BinaryNumber +=1
#        except:
#            raise
#        finally:
#            fptr.close()
#


if 0:
    #good distance plot
    for sim in sim_list:
        this_looper=TL.loops[sim]
        main_core_list=np.unique(this_looper.tr.core_ids)
        min_d,cores_used, distance= find_other_cores.get_all_distances( this_looper, main_core_list, mini_scrubbers[sim])
        distance[ distance==0] = 4000
        dmin = distance.min(axis=2)
        d_end = distance[:,:,-1]
        min_end = d_end.min(axis=1)
        fig,axes=plt.subplots(1,3)
        min_of_mins = np.argsort( dmin.min(axis=0))
        counter={}
        ext=extents()
        for core_id in cores_used[min_of_mins]:
            core_index = np.where( cores_used == core_id)[0][0]
            for nmode, mode in enumerate(['Alone', 'Binary','Cluster']):
                mode_list=TL.loops[sim].core_by_mode[mode]
                if core_id in mode_list:
                    counter[mode] = counter.get(mode,0)+1
                    this_mode = mode
                    this_nmode = nmode
                    break
            the_y = nar(sorted(dmin[:,core_index])[:5])
            ncores=the_y.size
            the_x = nar([counter[mode]]*ncores)
            axes[this_nmode].scatter( the_x,the_y)
            ext(the_y)
        for aaa in axes:
            aaa.set(ylim=ext.minmax,yscale='log')
        fig.savefig('plots_to_sort/sorted_distance_%s'%sim)






if 0:
    #Distance and Overlap.
    for sim in sim_list:
        reload(find_other_cores)
        this_looper = TL.loops[sim]
        main_core_list = np.unique(TL.loops[sim].tr.core_ids)
        #main_core_list.sort()
        fig,axes=plt.subplots(2,3,figsize=(12,8))
        fig2,axes2=plt.subplots(1,1)
        fig3,axes3=plt.subplots(1,4, figsize=(15,8))
        a30=axes3[0];a31=axes3[1]; a32=axes3[2];a33=axes3[3]
        ext=extents()
        ext2=extents()
        ext3=extents()
        loost=['Alone', 'Binary','Cluster']
        #loost=['Cluster']
        min_d,cores_used, distance= find_other_cores.get_all_distances( this_looper, main_core_list, mini_scrubbers[sim])
        d_end = distance[:,:,-1]
        min_end = d_end.min(axis=1)
        for nmode,mode in enumerate(loost):
            mode_list=TL.loops[sim].core_by_mode[mode]
            if 1:
                #ht = hull_by_frame[sim][0]
                ht = anne.ht[sim]
                mode_overlap = np.zeros_like(d_end)
                total_overlap = np.zeros(len(mode_list))
                mean_ratio = np.zeros(len(mode_list))
                N_over = np.zeros(len(mode_list))
                min_dist = np.zeros(len(mode_list))
                dist_1 = np.zeros(len(mode_list))
                dist_2 = np.zeros(len(mode_list))
                dist_3 = np.zeros(len(mode_list))
                min_all = np.zeros(len(mode_list))
                for nc,core_id in enumerate(mode_list):
                    over_index = np.where(ht.cores_used==core_id)[0][0]
                    total_overlap[nc] = overlaps[sim][over_index,:,0].sum()
                    N_over[nc] = (overlaps[sim][over_index,:,0]>0).sum()
                    mean_ratio[nc] = ratios[sim][over_index,:,0].mean()
                    min_dist[nc]=sorted(distance[over_index,:,-1])[1]
                    dist_1[nc]=sorted(distance[over_index,:,-1])[2]
                    dist_2[nc]=sorted(distance[over_index,:,-1])[3]
                    dist_3[nc]=sorted(distance[over_index,:,-1])[4]
                    d_all=distance[over_index,:,:].flatten()
                    #min_all[nc]=sorted(d_all[d_all>0])[0]
                    #a31.text(min_dist[nc],mean_ratio[nc],"%s"%core_id)


                #a30.scatter(min_dist, total_overlap, c='rgb'[nmode])
                a30.scatter(min_dist, mean_ratio, c='rgb'[nmode])
                #a32.scatter(min_dist, N_over, c='rgb'[nmode])
                #a33.scatter(min_all, min_dist, c='rgb'[nmode])
                a31.scatter(min_dist/dist_1, dist_1/dist_2,c='rgb'[nmode])
                #a32.scatter( dist_2,dist_2/dist_3,c='rgb'[nmode])
                a32.scatter( min_dist/dist_1, dist_2,c='rgb'[nmode])
                #a33.scatter(min_dist,dist_1,c='rgb'[nmode], alpha=0.5)


            #axes[1][nmode].hist( np.log10(first_distance))
        #for nmode in [0,1,2]:
        #    axes[0][nmode].set(yscale='log',xlabel='core',ylabel='min ever', ylim=ext.minmax, xlim=[0,10])
        #    axes[1][nmode].set(yscale='log',xlabel='core',ylabel='Delta t=0', ylim=ext2.minmax, xlim=[0,10])
        #fig.savefig('plots_to_sort/min_d')
        #axes2.set(xscale='log',xlabel=xlab,yscale='log',ylabel=ylab,xlim=ext3.minmax,ylim=ext3.minmax)
        #fig2.savefig('plots_to_sort/distance_distance')
        a30.set(xlabel=r'$\Delta$',ylabel='mean_ratio')
        #a31.set(xlabel=r'$\Delta$',ylabel='mean_ratio')
        a31.set(xlabel='d0/d1',ylabel='d1/d2')
        fig3.savefig('plots_to_sort/distance_overlap_%s'%sim)


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



