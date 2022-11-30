

from starter2 import *

#import three_loopers_six as TL
import three_loopers_u500 as TL
import convex_hull_tools as CHT
reload(CHT)
sim_list = ['u601','u602','u603']
sim_list=['u502']
#sim_list=['u502']

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
        frame_list=loop.tr.frames
        #frame_list=[125]
        #frame_list=[0]
        for nframe, frame in enumerate(frame_list):
            hull_by_frame[name][frame]=CHT.hull_tool(loop)
            hull_by_frame[name][frame].make_hulls(frames=[frame])
            hull_by_frame[name][frame].make_overlaps()

if 'overlaps' not in dir():
    overlaps={}

for sim in sim_list:
    ht0=hull_by_frame[sim][0]
    loop=TL.loops[sim]
    frame_list=loop.tr.frames
    #frame_list=[0]
    ncores=len(ht0.cores_used)
    nframes=len(frame_list)

    if sim not in overlaps:
        overlaps[sim] = np.zeros([ncores,ncores,nframes])
        for nframe, frame in enumerate(frame_list):
            htool = hull_by_frame[sim][frame]
            for nc1, core_id_1 in enumerate(htool.cores_used):
                for nc2, core_id_2 in enumerate(htool.cores_used):
                    if nc1==nc2:
                        continue
                    a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                    #please come back to think about this logic.
                    #david, 2022-09-22
                    ratio = max([a,b])
                    rat= sorted( [a,b])
                    if rat[1] == 0:
                        ratio = max([a,b])
                    else:
                        ratio=rat[0]/rat[1]
                    overlaps[sim][nc1,nc2,nframe] = ratio

if 1:
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



if 1:
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






    
