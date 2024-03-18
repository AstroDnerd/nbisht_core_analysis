

from starter2 import *

#import three_loopers_six as TL
#import three_loopers_u500 as TL
import track_loader as TL
import convex_hull_tools as CHT
reload(CHT)
sim_list = ['u601','u602','u603']
sim_list=['u502']
#sim_list=['u502']
TL.load_tracks(sim_list)

if 'hull_by_frame' not in dir():
    #
    # Hull volume vs total cell volume
    #
    hull_by_frame = {}
    #looper_list=[MOD.loops['u401']] #,MOD.loops['u402'],MOD.loops['u403']]
    #looper_list=[MOD.loops['u402']] #,MOD.loops['u402'],MOD.loops['u403']]
    #loopers = MOD.loops
if 1:
    for sim in sim_list:
        if sim in hull_by_frame:
            continue
        loop=TL.loops[sim]
        name = loop.sim_name
        hull_by_frame[name]={}
        frame_list=[0,30] #loop.tr.frames
        #frame_list=[125]
        #frame_list=[0]
        for nframe, frame in enumerate(frame_list):
            hull_by_frame[name][frame]=CHT.hull_tool(loop)
            hull_by_frame[name][frame].make_hulls(frames=[frame])
            hull_by_frame[name][frame].make_overlaps()

if 'overlaps' not in dir():
    overlaps={}
    ratios={}

for sim in sim_list:
    ht0=hull_by_frame[sim][0]
    loop=TL.loops[sim]
    frame_list=loop.tr.frames
    #frame_list=[0]
    ncores=len(ht0.cores_used)
    nframes=len(frame_list)

    if sim not in overlaps:
        overlaps[sim] = np.zeros([ncores,ncores,nframes])
        ratios[sim] = np.zeros([ncores,ncores,nframes])
        for nframe, frame in enumerate(frame_list):
            htool = hull_by_frame[sim][frame]
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

if 0:
    import convex_hull_plot2d as CHP
    reload(CHP)
    import tsing
    reload(tsing)
    if 'tsing_tool' not in dir():
        tsing_tool={}
        for ns,sim in enumerate(sim_list):
            obj=tsing.te_tc(TL.loops[sim])
            tsing_tool[sim]=obj
            tsing_tool[sim].run()
if 0:
    #numbers
    for sim in sim_list:
        loop=TL.loops[sim]
        times=loop.tr.times/colors.tff
        OOO = overlaps[sim]
        Tot = (OOO>0).sum(axis=0)
        for mode in ['Cluster']:
            core_list=loop.core_by_mode[mode]
            tool0 = hull_by_frame[sim][0]
            cores_used = nar(tool0.cores_used)
            times=loop.tr.times/colors.tff
            fig,ax=plt.subplots(1,1)
            for nc1,core_id in enumerate(core_list):
                index = np.where( cores_used == core_id)[0][0]
                any_overlap = OOO[index,:,:] > 0
                line = any_overlap.sum(axis=0)
                ax.plot(times, line)
                print(line)
            fig.savefig('plots_to_sort/Counts_%s_mode_%s'%(loop.sim_name,mode))
            plt.close(fig)
if 0:
    #fractions 
    for sim in sim_list:
        loop=TL.loops[sim]
        times=loop.tr.times/colors.tff
        OOO = overlaps[sim]
        Tot = (OOO>0).sum(axis=0)
        for mode in ['Binary']:
            core_list=loop.core_by_mode[mode]
            #core_list=core_list[:2]
            #core_list=[76]
            tool0 = hull_by_frame[sim][0]
            cores_used = nar(tool0.cores_used)
            times=loop.tr.times/colors.tff
            fig,ax=plt.subplots(1,1)
            for nc1,core_id in enumerate(core_list):
                index = np.where( cores_used == core_id)[0][0]
                my_buddies = OOO[index,:,:].sum(axis=1) > 0
                #my_buddies[index]=True
                other_core_id_list = cores_used[my_buddies]
                print(other_core_id_list)

                for other_core_id in other_core_id_list:
                    other_index = np.where( cores_used == other_core_id)[0][0]
                    line=OOO[index,other_index,:]
                    tnorm = tsing_tool[name].tsing_core[core_id]
                    ax.plot( times, line)
                #print(other_core_id_list)
                for nf,frame in enumerate(loop.tr.frames):
                    continue
                    if nf%10:
                        continue
                    this_over = OOO[index,:,nf] > 0
                    #if frame%10==0:
                    #    print(frame)
                    #    print(cores_used[this_over])

                    CHP.plot_2d(hull_by_frame[sim][frame], core_list=other_core_id_list, 
                                prefix="n%04d_"%frame, axis_to_plot=[-1],accumulate=True,all_plots=False,
                                label_cores=[-1], frames=[frame])
                #fig.savefig('plots_to_sort/Overlaps_%s_c%04d_mode_%s'%(loop.sim_name,core_id,mode))
                #plt.close(fig)
            fig.savefig('plots_to_sort/Overlaps_%s_mode_%s'%(loop.sim_name,mode))
            plt.close(fig)


if 0:
    for sim in sim_list:
        loop=TL.loops[sim]
        times=loop.tr.times/colors.tff
        OOO = overlaps[sim]
        Tot = (OOO>0).sum(axis=0)
        for mode in ['Alone']:
            core_list=loop.core_by_mode['Alone']
            core_list=core_list[:2]
            tool0 = hull_by_frame[sim][0]
            cores_used = nar(tool0.cores_used)
            times=loop.tr.times/colors.tff
            for nc1,core_id in enumerate(core_list):
                fig,ax=plt.subplots(1,1)
                index = np.where( cores_used == core_id)[0][0]
                my_buddies = OOO[index,:,:].sum(axis=1) > 0
                other_core_id_list = cores_used[my_buddies]
                for other_core_id in other_core_id_list:
                    other_index = np.where( cores_used == other_core_id)[0][0]
                    ax.plot( times, OOO[nc1,other_index,:])
                fig.savefig('plots_to_sort/Overlaps_%s_c%04d_mode_%s'%(loop.sim_name,core_id,mode))
                plt.close(fig)

            



if 0:
    for sim in sim_list:
        loop=TL.loops[sim]
        times=loop.tr.times/colors.tff
        OOO = overlaps[sim]
        Tot = (OOO>0).sum(axis=0)
        fig,ax=plt.subplots(1,2)
        args = np.argsort( Tot[:,50])
        if 0:
            ax[0].imshow(Tot[args,:], cmap='Reds')

        if 0:
            ttt = times+0
            ttt.shape=ttt.size,1
            ax[0].plot(ttt,Tot.transpose())



        for n in range(1,30):
            #more_than= ((Tot>=n)*(Tot<n+1)).sum(axis=0)
            more_than= ((Tot>=n)).sum(axis=0)
            ax[1].plot( times, more_than)
        more_than= ((Tot==0)).sum(axis=0)
        ax[1].plot( times, more_than)
        more_than= ((Tot==1)).sum(axis=0)
        ax[1].plot( times, more_than)
        more_than= ((Tot==2)).sum(axis=0)
        ax[1].plot( times, more_than)

        startswith=Tot[:,0] == 17
        ax[0].plot( times, Tot[startswith,:].transpose())




        fig.savefig('plots_to_sort/derp_%s.png'%sim)
if 0:
    for sim in sim_list:
        loop=TL.loops[sim]
        times=loop.tr.times/colors.tff
        OOO = overlaps[sim]
        Tot = (OOO>0).sum(axis=0)
        fig,ax=plt.subplots(1,1)

        for startswith in range(Tot.max()):
            ax.clear()
            if not (Tot[:,0]==startswith).any():
                continue
            xbins = times
            ybinedge = np.arange(50)-0.5
            ybins = 0.5*(ybinedge[1:]+ybinedge[:-1])
            nx = len(xbins) ; ny=len(ybins)
            TheX = np.r_[(ny)*[xbins]].transpose()
            TheY = np.r_[(nx)*[ybins]]
            hist = np.zeros( [xbins.size,ybins.size])
            for ntime, time in enumerate(times):
                thishist,bins = np.histogram(Tot[:,ntime],bins=ybinedge)
                hist[ntime,:]=thishist
            cmap = copy.copy(mpl.cm.get_cmap("viridis"))
            cmap.set_under('w')
            minmin = hist[hist>0].min()
            norm = mpl.colors.LogNorm(vmin=1,vmax=33)
            ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')


            use_these=Tot[:,0] == startswith
            ax.plot( times, Tot[use_these,:].transpose())
            ax.set_title("startswith %d (%d)"%(startswith, use_these.sum()))

            fig.savefig('plots_to_sort/derp_%s_startswith_%03d.png'%(sim,startswith))






    
