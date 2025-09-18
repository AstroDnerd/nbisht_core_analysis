
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def counter(this_looper,core_list=None, do_sphere=False, external_ax=None, hist_args={}):
    print('wu',do_sphere)

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    frame_ind=0
    x = this_looper.tr.track_dict['x'][:,frame_ind]
    y = this_looper.tr.track_dict['y'][:,frame_ind]
    z = this_looper.tr.track_dict['z'][:,frame_ind]
    #x = this_looper.tr.track_dict['particle_pos_x'][:,frame_ind]
    #y = this_looper.tr.track_dict['particle_pos_y'][:,frame_ind]
    #z = this_looper.tr.track_dict['particle_pos_z'][:,frame_ind]
    dx=dy=dz=1./128
    i = np.floor( x/dx).astype('int')
    j = np.floor( y/dx).astype('int')
    k = np.floor( z/dx).astype('int')
    III=  np.zeros([128]*3)
    III[(i,j,k)]=1



    if do_sphere:

        i,j,k=np.mgrid[0:128:1,0:128:1,0:128:1]-64
        z = np.sqrt(i**2+j**2+k**2)
        SSS=  np.zeros([128]*3)
        VVV = III.sum()
        R = (3/(4*np.pi)*VVV)**(1./3)
        SSS[z < R]=1
        #plt.clf()
        #plt.imshow(SSS.sum(axis=0))
        #plt.savefig('plots_to_sort/derp3.png')

        V = 4/3* np.pi*(32)**3
        A = np.pi*(32)**2
        print('V',V)
        print('A',A)
        print(SSS.sum())



    def neighbor_counter(III):
        JJJ=  np.zeros(nar(III.shape)+2)
        sl = slice(1,-1)
        JJJ[sl,sl,sl]=III
        #periodic shift
        JJJ[0,:,:]=JJJ[-2,:,:]
        JJJ[:,0,:]=JJJ[:,-2,:]
        JJJ[:,:,0]=JJJ[:,:,-2]
        JJJ[-1,:,:]=JJJ[1,:,:]
        JJJ[:,-1,:]=JJJ[:,1,:]
        JJJ[:,:,-1]=JJJ[:,:,1]

        collector = np.zeros([128]*3)
        nx,ny,nz=JJJ.shape
        for i in  [0,1,2]:
            for j in  [0,1,2]:
                for k in  [0,1,2]:
                    si = slice( i, nx-2+i)
                    sj = slice( j, ny-2+j)
                    sk = slice( k, nz-2+k)
                    collector += JJJ[si,sj,sk]
        collector = collector[III==1]
        return collector


    collector = neighbor_counter(III)
    if external_ax is None:
        fig,ax=plt.subplots(1,1)
    else:
        ax=external_ax
    bins = np.arange(-0.5,28.5)
    oot=ax.hist( collector.flatten(), bins=bins,cumulative=False, histtype='step',**hist_args)
    if do_sphere:
        counter2 = neighbor_counter(SSS)
        ax.hist( counter2, bins=bins,cumulative=False, histtype='step', color=[0.5]*3)
    oot=1
    if external_ax is None:
        axbonk(ax,yscale='log',ylabel='N',xlabel='Nneighbors')
        fig.savefig('plots_to_sort/neighbor_zone_counter_%s.png'%this_looper.sim_name)
    #plt.clf()
    #plt.imshow( JJJ.sum(axis=0))
    #plt.savefig('plots_to_sort/derp.png')
    return oot, collector





sims=['u501', 'u502','u503']
fig,ax=plt.subplots(1,1)
for ns,sim in enumerate(sims):
    do_sphere=False
    if ns==0:
        do_sphere=True

    oot=counter(TL.loops[sim], do_sphere=do_sphere,external_ax=ax,hist_args={'color':colors.color[sim],'label':r'$sim%d$'%(ns+1)})
ax.legend(loc=2)
axbonk(ax,xlabel='N neighbors',ylabel='Nzone',yscale='log')
ax.set_xticks(np.linspace(0,27,10))
fig.savefig('plots_to_sort/surface_to_volume.png')

