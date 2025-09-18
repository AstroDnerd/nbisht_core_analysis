
from starter2 import *

import convex_hull_tools as CHT

reload(CHT)
import hair_dryer
reload(hair_dryer)
#import stay_close
#import three_loopers_tenfour as TL4
import three_loopers_six as TL
import close_tool
import convex_hull_tools as CHT
reload(CHT)
sim_list=['u601','u602','u603']
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL.loops[this_simname], ht[this_simname])
        st[this_simname].find()
plt.close('all')



if 1:
    #every pair
    #the good one, don't touch it.
    import means_etc
    for this_simname in sim_list:
        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        ooo = np.zeros( [len(htool.cores_used)]*2) -1
        distance_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        color_matrix = np.zeros_like(distance_matrix)
        neighbor_matrix = np.zeros_like(distance_matrix).astype('int')
        particle_matrix = np.zeros_like(distance_matrix)
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                #if nc2 <= nc1:
                #    continue
                ooo[nc1,nc2] = htool.overlaps[core_id_1][nc2]
                overlap_matrix[nc1,nc2] = min( [htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] ])
                color_matrix[nc1,nc2]=overlap_matrix[nc1,nc2]
                neighbor = 1
                if stool.set_by_core[core_id_1] == stool.set_by_core[core_id_2]:
                    neighbor = 2
                neighbor_matrix[nc1,nc2]=neighbor
                a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                ratio = max([a,b])
                rat= sorted( [a,b])
                if rat[1] == 0:
                    ratio = max([a,b])
                else:
                    ratio=rat[0]/rat[1]

                color_matrix[nc1,nc2]=ratio
        #fig, ax, axtop,axright = means_etc.three_way_bean()
        fig, ax, axtop,axright = multiplots.three_way_bean(figsize=(4,4), left=0.15, width=0.62, bottom=0.11, height=0.62, histdepth=0.02)

        Md = np.triu(ctool.distance_matrix).flatten()
        #PHYSICS UNITS
        Md *= 4.6 #NOW ITS IN PARSECS
        Mo = np.triu(overlap_matrix).flatten()

        ok = (Md > 0)*(Mo >= 0)
        cmap=copy.copy(mpl.cm.get_cmap("viridis"))
        norm=mpl.colors.Normalize(vmin=0,vmax=1)
        color=cmap(norm(color_matrix.flatten()))

        neighbor=(neighbor_matrix).flatten()[ok]
        ax.scatter( Md[ok][ neighbor==1], Mo[ok][ neighbor==1]-0.01,color='r')
        ax.scatter( Md[ok][ neighbor==2], Mo[ok][ neighbor==2],color=color[ok][neighbor==2,:])
        #ax.scatter( Md[ok][ neighbor==2], Mo[ok][ neighbor==2],color=color[neighbor==2])
        #ax.scatter( Md[ok][ neighbor==2], Mo[ok][ neighbor==2],color='k')
        bins = np.geomspace( Md[ok].min(), Md[ok].max(),64)
        axtop.hist( Md[ok][ neighbor == 1], histtype='step',color='r',bins=bins)
        axtop.hist( Md[ok][ neighbor == 2], histtype='step',color='k',bins=bins)
        #axright.hist( Mo[ok][neighbor == 1], histtype='step',color='r',orientation='horizontal',bins=64)
        axright.hist( Mo[ok][neighbor == 2], histtype='step',color='k',orientation='horizontal',bins=64)
        axbonk(ax, xlabel='Distance [pc]', ylabel='overlap',xscale='log',yscale='linear')
        axbonk(axtop, xlabel='', ylabel='N',xscale='log',yscale='log')
        axbonk(axright, xlabel='N', ylabel='',xscale='log',yscale='linear')
        axtop.set_xlim( ax.get_xlim())
        axright.set_ylim( ax.get_ylim())
        axtop.set_xticks([])
        axright.set_yticks([])

        #ax.set_xticks( nar( ax.get_xticks())*4.6)
        print("==========")
        print( ax.get_xticks())



        fig.savefig('plots_to_sort/%s_stay_close.pdf'%(htool.name))
        plt.close('all')


if 0:
    #closest pairs
    import means_etc
    for this_simname in sim_list:
        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        distance_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        color_matrix = np.zeros_like(distance_matrix).astype('int')
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                overlap_matrix[nc1,nc2] = max( [htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] ])
                color = 1
                if stool.set_by_core[core_id_1] == stool.set_by_core[core_id_2]:
                    color = 2
                color_matrix[nc1,nc2]=color
        #fig,ax=plt.subplots(1,1)
        fig, ax, axtop,axright = means_etc.three_way_bean()
        Md1 = ctool.distance_matrix +0
        Mo1 = overlap_matrix

        Md1[ Md1==0] = 5 #the closest core to any one is itself.  So make zero 5 so we can take the min.
        closest = np.argmin(Md1,axis=0)

        Md = Md1[closest, tuple(range(Md1.shape[1]))]
        Mo = Mo1[closest, tuple(range(Md1.shape[1]))]

if 0:

        #Md = ctool.distance_matrix.flatten()
        #Mo = overlap_matrix.flatten()
        #ok = (Md >= 0)*(Mo >= 0)
        ok = slice(None)
        color = (color_matrix[closest, tuple(range(Md1.shape[1]))]).flatten()[ok]
        c = ['grk'[v] for v in color]
        ax.scatter( Md[ok][ color==1], Mo[ok][ color==1]-0.01,color='r')
        ax.scatter( Md[ok][ color==2], Mo[ok][ color==2],color='k')
        bins = np.geomspace( Md[ok].min(), Md[ok].max(),64)
        axtop.hist( Md[ok][ color == 1], histtype='step',color='r',bins=bins)
        axtop.hist( Md[ok][ color == 2], histtype='step',color='k',bins=bins)
        #axright.hist( Mo[ok][color == 1], histtype='step',color='r',orientation='horizontal',bins=64)
        axright.hist( Mo[ok][color == 2], histtype='step',color='k',orientation='horizontal',bins=64)
        axbonk(ax, xlabel='Distance', ylabel='overlap',xscale='log',yscale='linear')
        axbonk(axtop, xlabel=None, ylabel='N',xscale='log',yscale='log')
        axbonk(axright, xlabel='N', ylabel=None,xscale='log',yscale='linear')
        axtop.set_xlim( ax.get_xlim())
        axright.set_ylim( ax.get_ylim())

        fig.savefig('plots_to_sort/%s_stay_closest.png'%(htool.name))
        plt.close('all')





if 0:
    #messing around.
    import means_etc
    fig,ax=plt.subplots(1,3,figsize=(12,8))
    #for aaa in ax:
    #    aaa.set_aspect('equal')
    fig2,ax2=plt.subplots(3,1)
    for ns,this_simname in enumerate(sim_list):
        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]
        ooo = np.zeros( [len(htool.cores_used)]*2) -1
        nnn = np.zeros( [len(htool.cores_used)]*2) -1
        ppp = np.zeros( [len(htool.cores_used)]*2) -1
        color_matrix = np.zeros_like(distance_matrix)
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                ooo[nc1,nc2] = htool.overlaps[core_id_1][nc2]
                nnn[nc1,nc2] = htool.overlap_number[core_id_1][nc2]
                ppp[nc1,nc2] = htool.nparticles[nc1]
                a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                ratio = max([a,b])
                rat= sorted( [a,b])
                if rat[1] == 0:
                    ratio = max([a,b])
                else:
                    ratio=rat[0]/rat[1]

                color_matrix[nc1,nc2]=ratio
        #both = np.concatenate([[ooo],[ooo.transpose()]])
        both = np.concatenate([[nnn],[nnn.transpose()]])
        mmmin = both.min(axis=0)
        mmmax = both.max(axis=0)
        szz = np.concatenate([[ppp],[ppp.transpose()]])
        n_smal = szz.min(axis=0)
        n_big  = szz.max(axis=0)
        top = np.triu(mmmin)
        bot = np.triu(mmmax)
        ok1 = (bot > 0)*(top > -0.25)

        dist = np.triu(ctool.distance_matrix)[ok1]
        ratio1 =  top[ok1]/bot[ok1]
        #the_y = ratio
        #the_y=bot[ok1] #ratio
        #the_y = top[ok1]
        #the_y = np.triu(nnn/ppp.transpose())[ok1]

        ax[ns].scatter(top[ ok1],ratio1)
        axbonk(ax[ns],xlabel='overlap',ylabel='ratio')#,xscale='log')
        #ax[ns].set_yscale('symlog',linthresh=1)
        #ax[ns].set_ylim([0,the_y.max()])

        #derp = nnn>0
        #X = (nnn/ppp)[derp]
        #Y = (nnn.transpose()/ppp)[ derp.transpose()]
        #ax2[ns].scatter(X, Y)
        #ax2[ns].set_yscale('symlog',linthresh=1)

    fig.savefig('plots_to_sort/ug.png')
    #fig2.savefig('plots_to_sort/ug2.png')
