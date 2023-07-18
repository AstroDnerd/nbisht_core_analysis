
from starter2 import *

from scipy.ndimage import gaussian_filter

#import three_loopers_u900 as TL
import track_loader as TL

from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')
import pcolormesh_helper as pch
reload(pch)

def splat(array, tcenter, ax, title,bins):
    smooth= ndimage.gaussian_filter1d(array, 2, 0)
    smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
    ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax)
    ax.set_yscale('symlog',linthresh=100)
    ax.set_title(title)
    return ds_x,ds_y,ds_h,ds_dv,ds_p
class dq_dt2():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self,core_list=None,frame_list=None):
        this_looper=self.this_looper
        thtr = this_looper.tr

        if frame_list is None:
            mask = movie_frames.quantized_mask(this_looper)
            times = thtr.times[mask]
            if times[0] == times[1]:
                mask[0]=False
            times = thtr.times[mask]
            frame_list=this_looper.tr.frames[mask]

            mask_m1 = mask[:-1]
            mask_p1 = mask[1:]
            tcenter = 0.5*(times[:-1]+times[1:])
            dt = times[1:]-times[:-1]
            dt2 = times[2:]-times[:-2]
            dt2.shape = dt2.size,1
            tcen2 = 0.5*(times[2:]+times[:-2])

            dt_square = dt+0
            dt_square.shape = dt_square.size,1


        if core_list is None:
            core_list = sorted(np.unique( thtr.core_ids))


        for core_id in core_list:
            print('go',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            if ms.nparticles < 1000:
                sl=slice(None)
                c=[0.5]*4
            else:
                sl = slice(None,None,10)
                #c=[0,0,0,0.1]
                c=[0.1]*4
            rho = ms.density[sl].transpose()
            dv = ms.cell_volume[sl].transpose()[mask,:]
            rho = rho[mask,:]
            Bmag=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
            BP = Bmag**2/2
            B2 = Bmag**2
            divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]



            def extract(arr):
                return arr[sl].transpose()[mask,:]
            bx   =extract(thtr.c([core_id],'magnetic_field_x'))
            by   =extract(thtr.c([core_id],'magnetic_field_y'))
            bz   =extract(thtr.c([core_id],'magnetic_field_z'))
            vx   =extract(thtr.c([core_id],'velocity_x'))
            vy   =extract(thtr.c([core_id],'velocity_y'))
            vz   =extract(thtr.c([core_id],'velocity_z'))
            dxvx =extract(thtr.c([core_id],'dxvx'))
            dxvy =extract(thtr.c([core_id],'dxvy'))
            dxvz =extract(thtr.c([core_id],'dxvz'))
            dyvx =extract(thtr.c([core_id],'dyvx'))
            dyvy =extract(thtr.c([core_id],'dyvy'))
            dyvz =extract(thtr.c([core_id],'dyvz'))
            dzvx =extract(thtr.c([core_id],'dzvx'))
            dzvy =extract(thtr.c([core_id],'dzvy'))
            dzvz =extract(thtr.c([core_id],'dzvz'))
            Sx = bx*dxvx+by*dyvx+bz*dzvx
            Sy = bx*dxvy+by*dyvy+bz*dzvy
            Sz = bx*dxvz+by*dyvz+bz*dzvz
            Stretch= bx*Sx+by*Sy+bz*Sz

            rho_0=rho[0,:]#.mean()
            b_0=BP[0,:]#.mean()

            RB=Stretch/(B2*divv)
            #dumbass=RB[:1,:]
            #print('dumbass',dumbass.shape)
            #RB_Mode = scipy.stats.mode(dumbass,axis=1)[0]
            #RB_Mode = scipy.stats.mode(dumbass)
            RB_Mode = scipy.stats.mode(RB,axis=1)[0]
            RB_Mode = gaussian_filter(RB_Mode,2)
            #RB_Mean = (rho*RB*dv).sum(axis=1)/(rho*dv).sum(axis=1)
            #RB_Mean.shape=RB_Mean.size,1
            RB_Mean = (rho*RB*dv).sum()/(rho*dv).sum()
            RB_Mean_t = (rho*RB*dv).sum(axis=1)/(rho*dv).sum(axis=1)
            Monster = 10
            NoMonsters = np.abs(RB) < Monster
            RB_Mean_10 = (rho*RB*dv*NoMonsters).sum(axis=1)/(rho*dv*NoMonsters).sum(axis=1)
            RB_Mean_10 = gaussian_filter(RB_Mean_10,2)
            RB_Mean_10.shape=RB_Mean_10.size,1
            #print(RB_Mode)

            bins_1 = np.geomspace( 1,1e9,19)
            bins_m1 = -bins_1[::-1]
            bins = nar(list(bins_m1)+list(bins_1))
            bincen = 0.5*(bins[1:]+bins[:-1])

            fig, ax=plt.subplots(4,2,figsize=(8,12))
            ax0=ax[0][0]
            ax1=ax[0][1] 

            ax2=ax[1][0]
            ax3=ax[1][1]

            ax4=ax[2][0]
            ax5=ax[2][1]

            ax6=ax[3][0]
            ax7=ax[3][1]

            if 1: #mess around
                #R take 1
                THE_AX=ax0
                ext_r = extents(RB)
                if 1:
                    #this one shows stuff
                    bins_1 = np.geomspace( np.abs(RB).min(), np.abs(RB).max(),32)
                    bins_1 = np.geomspace( Monster, np.abs(RB).max(),32)
                    bins_2 = np.linspace(-Monster,Monster,5*16)
                    bins_m1 = -bins_1[::-1]
                    bins = nar(list(bins_m1)+list(bins_2)+list(bins_1))
                    bins = np.unique(bins)

                PRB = splat(1-RB,tcenter,THE_AX,'1-RB',bins)
                THE_AX.set_yscale('symlog', linthresh=Monster)

            if 1:
                THE_AX=ax1
                if 1:
                    #this one shows the bulk
                    bins = np.linspace(-Monster,Monster,32)

                #bins=np.geomspace(RB[RB>0].min(),RB.max())
                #bins=np.linspace(RB.min(),RB.max())
                PRB = splat(1-RB,tcenter,THE_AX,'1-RB',bins)
                THE_AX.set_yscale('linear')

            if 1:
                THE_AX=ax1
                THE_AX.plot(times, 1-RB_Mean_10,'r')
                THE_AX.set(ylim=[-Monster,Monster])

            if 1:
                #
                # B and rho
                #
                THIS_AX = ax2

                THIS_AX.plot(times, rho/rho_0, c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, BP/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='B/B0,rho/rho_0')
            if 1:
                #B rho phase
                THIS_AX=ax3
                ext=extents()

                x = (rho/rho_0).flatten()
                y = (B2/b_0).flatten()
                ext(x);ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins]) 
                THIS_AX.set(xscale='log',yscale='log',ylabel='B',xlabel='rho')
                THIS_AX.plot(ext.minmax,ext.minmax,c='k')

                pfit = np.polyfit( np.log10(x), np.log10(y),1)
                #print(pfit, RB_Mean_10)
                ux=np.unique(x)
                lnx=np.log10(ux)
                THIS_AX.plot( ux , 10**(pfit[0]*lnx+pfit[1]),c='orange')

            if 1:
                THIS_AX = ax4

                RHO_FIX=(rho/rho_0)**np.abs((1-RB_Mean_10))

                THIS_AX.plot(times, RHO_FIX, c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, B2/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='Mean RB 10')

                THIS_AX=ax5
                y=B2/b_0
                x = RHO_FIX
                ext(x);ext(y)
                x=x.flatten();y=y.flatten()
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins]) 
                THIS_AX.plot(ext.minmax,ext.minmax)
                THIS_AX.set(xscale='log',yscale='log', ylabel='B',xlabel='rho')



            if 1:
                #check you did the math right.
                x = np.log10((rho/rho_0))
                y = np.log10((B2/b_0))
                time_ax=0; part_ax=1
                xbar = np.mean(x, axis=part_ax)
                ybar = np.mean(y, axis=part_ax)
                xbar.shape =  xbar.size,1
                ybar.shape =  ybar.size,1
                top = ((y - ybar)*(x-xbar)).sum(axis=part_ax)
                bot = ( (x-xbar)**2).sum(axis=part_ax)
                kappa = top/bot
                norm = ybar-kappa*xbar
                #ax3.plot( 10**x, 10**(kappa*x+norm),c='g')
                THIS_AX=ax6
                THIS_AX.hist(kappa, histtype='step', label=r'$\kappa_{c,f}$')

            if 1:
                THIS_AX=ax6
                THIS_AX.hist(1-RB_Mean_10, histtype='step', label=r'$\langle 1- R_B\rangle(t)$')
                THIS_AX.axvline(pfit[0], label=r'$\kappa_c$', color='orange')
                THIS_AX.legend(loc=0)

            if 1:
                THIS_AX=ax7
                THIS_AX.plot(times, kappa, label=r'$\kappa_{c,f}$')
                THIS_AX.plot(times, 1-RB_Mean_10,c='r')
                THIS_AX.legend(loc=0)

                THIS_AX.set(xlabel='time',ylabel='kappa_c')

            if 0:
                #check you did the math right.
                x = np.log10((rho/rho_0))
                y = np.log10((B2/b_0))
                time_ax=0; part_ax=None
                xbar = np.mean(x, axis=part_ax)
                ybar = np.mean(y, axis=part_ax)
                top = ((y - ybar)*(x-xbar)).sum(axis=part_ax)
                bot = ( (x-xbar)**2).sum(axis=part_ax)
                kappa = top/bot
                norm = ybar-kappa*xbar
                ax3.plot( 10**x, 10**(kappa*x+norm),c='g')


            fig.tight_layout()
            fig.savefig('plots_to_sort/Frank_%s_c%04d'%(sim,core_id))

sim_list=['u902']
TL.load_tracks(sim_list)
for sim in sim_list:
    ddd = dq_dt2(TL.loops[sim])
    core_list=None
    #core_list=[7]
    #core_list=[74, 353]
    #core_list=[74]
    #core_list=[112]
    core_list=TL.loops[sim].core_by_mode['Alone']
    #core_list=core_list[:3]
    #core_list=[214, 114]
    P=ddd.run(core_list=core_list)


