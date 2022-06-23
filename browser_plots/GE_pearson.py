
from starter2 import *
from collections import defaultdict
import scipy
import colors

import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

def GE_pearson(this_looper,core_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    name = this_looper.sim_name
    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    G = colors.G
    #gx = thtr.track_dict['grav_x']
    #gy = thtr.track_dict['grav_y']
    #gz = thtr.track_dict['grav_z']
    #GE2 = -1/(8*np.pi)*(gx*gx+gy*gy+gz*gz)
    #ge_min=GE2.min()
    #ge_max=GE2.max()
    PearsonR = np.zeros([len(core_list), len(times)])
    PearsonP = np.zeros([len(core_list), len(times)])
    PearsonRho = np.zeros([len(core_list), len(times)])
    PeakRho  = np.zeros([len(core_list), len(times)])
    for nc, core_id in enumerate(core_list):
        print('GE pearson %s %d'%(name,core_id))

            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
        #ms.particle_pos(core_id)

        if ms.nparticles < 1000:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl]
        rho = rho[:,mask]

        PeakRho[nc,:]=rho.max(axis=0)

        gx = thtr.c([core_id],'grav_x')[sl][:,mask]
        gy = thtr.c([core_id],'grav_y')[sl][:,mask]
        gz = thtr.c([core_id],'grav_z')[sl][:,mask]
        GE2 = 1/(8*np.pi*G)*(gx*gx+gy*gy+gz*gz)

        RRR = ms.r[sl][:,mask]
        for n in range(GE2.shape[1]):
            the_x=np.log(RRR[:,n])
            the_y=np.log(GE2[:,n])
            #the_y=rho[:,n]
            r,p=scipy.stats.pearsonr(the_x,the_y)
            PearsonR[nc,n]=r
            PearsonP[nc,n]=p
            the_y=np.log(rho[:,n])
            r,p=scipy.stats.pearsonr(the_x,the_y)
            PearsonRho[nc,n]=r
            
    if 0:
        fig,ax=plt.subplots(1,2)
        ax[0].plot(times,PearsonR)
        #ax[0].boxplot(PearsonR)
        #ax[1].boxplot(PearsonRho)
        fig.savefig('plots_to_sort/phi_box_%s.png'%name)

    return {'PR':PearsonR, 'PP':PearsonP, 'Prho':PearsonRho, 'T':times, 'PeakRho':PeakRho}



    if 0:
        fig,ax=plt.subplots(1,1)
        ax.plot(times , GE2, c=c, linewidth=0.1)
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$(\nabla \phi)^2/8 pi G$',yscale='log', ylim=[ge_min,ge_max])
        ax2=ax.twinx()
        c=[1.0,0.1,0.1,0.1]
        ax2.plot(times , rho, c=c, linewidth=0.1)
        axbonk(ax2,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log')

        outname='plots_to_sort/%s_GE_t_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)



sims=['u501', 'u502','u503']
if 'stuff' not in dir():
    stuff={}
    for sim in sims:
        core_list = np.unique(TL.loops[sim].tr.core_ids)
        #core_list=core_list[:10]
        stuff[sim] = GE_pearson(TL.loops[sim],core_list=core_list)

if 1:
    for sim in stuff:
        fig,ax=plt.subplots(1,1)
        T = stuff[sim]['T']
        rho=stuff[sim]['PeakRho']
        Rphi=stuff[sim]['PR']
        ax.plot(Rphi.transpose() ,rho.transpose(),c=[0.1]*4)
        axbonk(ax,xlabel='time',ylabel='rho max', yscale='log')
        fig.savefig('plots_to_sort/peak_rho_pearson_phi%s.png'%sim)

if 1:
    for sim in stuff:
        fig,ax=plt.subplots(1,1)
        T = stuff[sim]['T']
        rho=stuff[sim]['PeakRho']
        ax.plot(T,rho.transpose(),c=[0.1]*4)
        axbonk(ax,xlabel='time',ylabel='rho max', yscale='log')
        fig.savefig('plots_to_sort/peak_rho_%s.png'%sim)

if 0:
    for sim in stuff:
        fig,ax=plt.subplots(1,1)
        c=[0.1]*4
        #ax.plot( stuff[sim]['T'], stuff[sim]['PR'].transpose(),c=c)
        #ax.scatter( stuff[sim]['Prho'].transpose(), stuff[sim]['PR'].transpose(),c=c)
        XX,YY= stuff[sim]['Prho'].flatten(), stuff[sim]['PR'].flatten()
        ok = (~np.isnan(XX))*(~np.isnan(YY))
        XX=XX[ok]
        YY=YY[ok]
        xbins = np.linspace( XX.min(), XX.max(), 64)
        ybins = np.linspace( YY.min(), YY.max(), 64)
        hist, xb, yb = np.histogram2d(XX,YY, bins=[xbins,ybins])
        import pcolormesh_helper as pch
        pch.helper(hist,xb,yb,ax=ax)
        fig.savefig('plots_to_sort/RGE_Rrho_%s.png'%sim)

if 1:
    for sim in stuff:
        fig,ax=plt.subplots(1,2)
        Rphi =   stuff[sim]['PR']
        ax[0].boxplot( Rphi )
        ax[0].plot( Rphi.mean(axis=0))
        ax[1].boxplot(  stuff[sim]['Prho'])


        axbonk(ax[0],xlabel='frame',ylabel='Rgrad phi')
        axbonk(ax[1],xlabel='frame',ylabel='R rho')
        fig.savefig('plots_to_sort/Boxes_%s.png'%(sim))


if 0:
    from scipy.ndimage import gaussian_filter
    fig,ax=plt.subplots()
    for sim in stuff:
        Rphi =   stuff[sim]['PR']
        Rrho =   stuff[sim]['Prho']
        ax.plot( gaussian_filter(Rphi.mean(axis=0),1), colors.color[sim] +'--')
        ax.plot( Rrho.mean(axis=0), colors.color[sim])


        axbonk(ax,xlabel='frame',ylabel='Rgrad phi')
        fig.savefig('plots_to_sort/MeanR_%s.png'%(sim))

