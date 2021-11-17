
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
import colors
plt.close('all')

import convex_hull_tools as CHT
reload(CHT)

#import three_loopers_1tff as tl
#import three_loopers_tenfour as TL4
#MOD = TL4


if 'loopers' not in dir():
    directory = dl.sims['u302']
    savefile =  "u302_long_pos_only.h5"
    print("Loading loop", 'u352', savefile)
    this_looper = looper.core_looper(directory= directory,savefile_only_trackage=savefile)
    looper_list=[this_looper]
    loopers={}
    loopers['u302'] = this_looper

if 'clobber' not in dir():
    clobber=False

if 'hull_by_frame' not in dir():
    #
    # Hull volume vs total cell volume
    #
    hull_by_frame = {}
    #looper_list=[MOD.loops['u401']] #,MOD.loops['u402'],MOD.loops['u403']]
    #looper_list=[MOD.loops['u402']] #,MOD.loops['u402'],MOD.loops['u403']]
    #loopers = MOD.loops
    for loop in looper_list:
        name = loop.out_prefix
        hull_by_frame[name]={}
        for nframe, frame in enumerate(loop.tr.frames):
            hull_by_frame[name][frame]=CHT.hull_tool(loop)
            hull_by_frame[name][frame].make_hulls(frames=[frame])
if 1:
    hullvol = defaultdict(list)
    cellvol = defaultdict(list)
    for loop in looper_list:
        name = loop.out_prefix
        hvol = []
        cvol = []
        for nframe, frame in enumerate(loop.tr.frames):
            hvol.append( hull_by_frame[name][frame].hull_volumes)
            #cvol.append( hull_by_frame[name][frame].cell_volumes)
            cvol.append( hull_by_frame[name][frame].unique_volumes)
        hullvol[name]=nar(hvol).transpose()
        cellvol[name]=nar(cvol).transpose()

if 1:
    #
    # Volumes and ratios by time.  
    #
    times_crossing_01 = []
    for name in hull_by_frame:
        times=loopers[name].tr.times
        fig,ax = plt.subplots(2,4, figsize=(12,8))
        fig2,ax2 = plt.subplots(1,3)
        fig3,ax3 = plt.subplots(1,2)

        H0 = hullvol[name][:,0]
        H0.shape = H0.shape[0], 1
        mean_hv = np.mean(hullvol[name]/H0,axis=0)
        C0 = cellvol[name][:,0]
        C0.shape = C0.shape[0], 1
        mean_cv = np.mean(cellvol[name]/C0,axis=0)

        THRESH = 0.1
        THRESHcv= 0.1
        mean_crossing = np.where(mean_hv < THRESH)[0][0]
        mean_crossing_time = times[mean_crossing]

        for nparticle, vols in enumerate(zip(hullvol[name],cellvol[name])):
            hvol,cvol = vols
            hv1 = hvol/hvol[0]
            cv1 = cvol/cvol[0]

            ax[0][0].plot( times, hvol ,c=[0.5]*3, linewidth=0.1)
            ax[0][1].plot( times, hv1, c=[0.5]*3, linewidth=0.1)
            ax[1][0].plot( times, cvol,c=[0.5]*3, linewidth=0.1)
            ax[1][1].plot( times, cv1,c=[0.5]*3, linewidth=0.1)

            first_crossing = np.where( hv1 < THRESH)[0][0]
            times_crossing_01.append(first_crossing)
            this_crossing_time = loopers[name].tr.times[first_crossing]
            hv2 = hv1#np.interp( times,times*mean_crossing_time/this_crossing_time , hv1)
            hv3 = hv2/mean_hv
            #ax[0][2].plot( times*mean_crossing_time/this_crossing_time, hv2, c=[0.5]*3, linewidth=0.1)
            #ax[0][3].plot( times*mean_crossing_time/this_crossing_time, hv3, c=[0.5]*3, linewidth=0.1)
            ax[0][2].plot( times, hv2, c=[0.5]*3, linewidth=0.1)
            ax[0][3].plot( times, hv3, c=[0.5]*3, linewidth=0.1)


            cross = np.where( cv1 < THRESH)[0][0]
            this_cross = loopers[name].tr.times[cross]
            #this_cross=this_crossing_time
            cv2 = np.interp( times,times*mean_crossing_time/this_cross, cv1)
            cv3 = cv1/mean_cv   
            #ax[1][2].plot( times*mean_crossing_time/this_crossing_time, cv2, c=[0.5]*3, linewidth=0.1)
            #ax[1][3].plot( times*mean_crossing_time/this_crossing_time, cv3, c=[0.5]*3, linewidth=0.1)
            ax[1][2].plot( times, cv2, c=[0.5]*3, linewidth=0.1)
            ax[1][3].plot( times, cv3, c=[0.5]*3, linewidth=0.1)

            #rat = hvol/hvol[:1].mean()/cvol/cvol[:1].mean()
            rat = hvol/cvol
            #rat /= rat[0]
            #rat /= rat[:-2].mean()
            ax2[0].plot( times, rat, c=[0.5]*3,linewidth=0.1)
            ratnorm = (hvol/hvol.mean())/(cvol/cvol.mean())
            ratnorm = rat/rat[:-2].mean()
            ax2[1].plot( times, ratnorm, c=[0.5]*3,linewidth=0.1)
            #ax[2].plot(times, mean_hv/mean_cv,c='k')

        #ax[0].plot(loopers[name].tr.times, mean_hv,c='k',linewidth=4)
        ax[0][1].scatter( mean_crossing_time, THRESH, c='k',s=8)
        axbonk(ax[0][0],yscale='log', title='Hull, raw')#, ylim=[1e-9,1e-1])
        axbonk(ax[0][1],yscale='log', title='Hull/H0')#, ylim=[1e-9,1e-1])
        axbonk(ax[0][2],yscale='log', title='Hull shift')#, ylim=[1e-9,1e-1])
        axbonk(ax[0][3],yscale='log', title='Hull renorm')#, ylim=[1e-9,1e-1])
        axbonk(ax[1][0],yscale='log', title='Cell, raw')#, ylim=[1e-9,1e-1])
        axbonk(ax[1][1],yscale='log', title='Cell/H0')#, ylim=[1e-9,1e-1])
        axbonk(ax[1][2],yscale='log', title='Cell shift')#, ylim=[1e-9,1e-1])
        axbonk(ax[1][3],yscale='log', title='Cell renorm')#, ylim=[1e-9,1e-1])

        axbonk(ax2[0],yscale='log', title='hvol/cvol')
        axbonk(ax2[1],yscale='log',ylim=[0.1,10], title='hvol/cvol/<mean>')
        axbonk(ax2[2],yscale='log',ylim=[0.1,10], title='norm, no shift')

        axbonk(ax3[0],yscale='log')
        axbonk(ax3[1],yscale='log')
        fig.savefig('plots_to_sort/%s_hvol_cvol_time.png'%name)
        fig2.savefig('plots_to_sort/%s_hvol_ratio_time.png'%name)
        plt.close('all')

if 1:
    #
    # Volumes and ratios by time.  
    #
    for name in hull_by_frame:
        fig,ax = plt.subplots(1,1)
        if name != 'u402':
            continue
        h0 = []
        c0 = []
        wtf=[]
        times = hull_by_frame[name][0].this_looper.tr.times
        for nparticle, vols in enumerate(zip(hullvol[name],cellvol[name])):
            hvol,cvol = vols
            oops = hvol == 0
            hvol[ oops] = cvol[oops]
            #hvol /= hvol[0]
            #cvol /= cvol[0]
            ax.plot( mean_hv, mean_cv, c='g', marker='*')
            ax.plot( hvol/hvol[:1].mean(), cvol/cvol[:1].mean(), c=[0.5]*3, linewidth=0.1)#, marker='.')
            #ax.plot( hvol, cvol, c=[0.5]*3, linewidth=0.1)#, marker='.')
            #ax.plot(times+0.01,hvol)
            h0.append(hvol[0])
            c0.append(cvol[0])
        #ax.scatter( h0, c0)
        for n in [0,1,2,3,4,5]:
            dx = (0.5)**n*1./128
            dv = dx**3
            print(1./dx)
            ax.plot( [1e-6,1e-4], [dv]*2, c='r')
            ax.plot(  [dv]*2,[1e-6,1e-4], c='r')
        xlim=None #[1e-7, 1.1]
        ylim=xlim
        #ax.plot(xlim,xlim,c='k')
        ax.plot([1e-8,5e-2],[1e-8,5e-2],c='k')
        axbonk(ax,xscale='log',yscale='log',xlabel='Hull',ylabel='cell',xlim=xlim,ylim=ylim)

        fig.savefig('plots_to_sort/ratio_time_paths.png')


