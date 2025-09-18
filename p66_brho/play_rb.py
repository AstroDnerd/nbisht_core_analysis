
from starter2 import *

import three_loopers_u900 as TL

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

            fig, ax=plt.subplots(3,6,figsize=(20,12))
            ax0=ax[0][0];
            ax1=ax[0][1] 
            ax2=ax[0][2];
            ax3=ax[0][3]
            ax4=ax[0][4];
            ax5=ax[0][5]

            ax6=ax[1][0];
            ax7=ax[1][1]
            ax8=ax[1][2]
            ax9=ax[1][3]
            ax10=ax[1][4]
            ax11=ax[1][5]

            ax12=ax[2][0]
            ax13=ax[2][1]
            ax14=ax[2][2]
            ax15=ax[2][3]
            ax16=ax[2][4]
            ax17=ax[2][5]

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
            NoMonsters = np.abs(RB) < 10
            RB_Mean_10 = (rho*RB*dv*NoMonsters).sum(axis=1)/(rho*dv*NoMonsters).sum(axis=1)
            RB_Mean_10 = gaussian_filter(RB_Mean_10,2)
            RB_Mean_10.shape=RB_Mean_10.size,1
            #print(RB_Mode)

            if 0:
                THE_AX=ax6
                b = RB[0,:]
                print(RB.shape)
                print(b.shape)
                print('mode t=0',scipy.stats.mode(b))
                print(1-RB_Mode.flatten())
                ax6.plot(b)
                




            bins_1 = np.geomspace( 1,1e9,19)
            bins_m1 = -bins_1[::-1]
            bins = nar(list(bins_m1)+list(bins_1))
            bincen = 0.5*(bins[1:]+bins[:-1])

            if 1: #mess around
                #R take 1
                THE_AX=ax0
                ext_r = extents(RB)
                if 1:
                    #this one shows stuff
                    bins_1 = np.geomspace( np.abs(RB).min(), np.abs(RB).max(),32)
                    bins_1 = np.geomspace( 5, np.abs(RB).max(),32)
                    bins_2 = np.linspace(-5,5,5*16)
                    bins_m1 = -bins_1[::-1]
                    bins = nar(list(bins_m1)+list(bins_2)+list(bins_1))
                    bins = np.unique(bins)

                PRB = splat(1-RB,tcenter,THE_AX,'1-RB',bins)
                THE_AX.set_yscale('symlog', linthresh=1)

            if 1:
                #plot the mode.
                THE_AX=ax1
                THE_AX.plot(times, 1-RB_Mode,'b')
                THE_AX.plot(times, 1-RB_Mean_t,'g')
                THE_AX.plot(times, 1-RB_Mean_10,'r')

                #minrb=(1-RB).min(axis=1)
                #THE_AX.plot(times, minrb)
                THE_AX.set(ylim=[-6,6])

            if 0:
                PeakQ=[]
                for Q in (1-RB):
                    Q.sort()
                    Nbins=1024
                    Npoints = int( Q.size/Nbins)
                    dQ = np.zeros(Nbins)
                    aQ = np.zeros(Nbins)
                    for Ibin in range(Nbins):
                        binI = Ibin*Npoints
                        binIp1 = min([(Ibin+1)*Npoints, Nbins-1])
                        dQ[Ibin] = Q[binIp1]-Q[binI]
                        aQ[Ibin] = 0.5*(Q[binIp1]+Q[binI])
                    PeakQ.append( aQ[np.argmin(dQ)])
                ax0.plot(times, PeakQ,c='purple')
            if 0: #Messy
                #
                # d log B / d ln rho
                #
                THIS_AX=ax6
                smooth_log_b=ndimage.gaussian_filter1d(np.log(B2), 2, 0)/2
                dlogb_dt = (smooth_log_b[2:,:]-smooth_log_b[:-2,:])/dt2
                smooth_rho=ndimage.gaussian_filter1d(np.log(rho), 2, 0)
                drho_dt = (smooth_rho[2:,:]-smooth_rho[:-2,:])/dt2
                dLogs=(dlogb_dt/drho_dt)
                tcen2 = 0.5*(times[:-2]+times[2:])
                #splat( dLogs, tcenter2, ax6, 'a',bins)
                if 0:
                    #get everything
                    bins_1 = np.geomspace( 1,np.abs(dLogs).max(),128)
                    bins_m1 = -bins_1[::-1]
                    bins = nar(list(bins_m1)+list(bins_1))
                if 1:
                    bins = np.linspace(-5,5,32)
                print(dLogs)
                dR_x,dR_y,dR_h,dR_dv,dR_p=heat_map.heat_map( dLogs.transpose(), tcen2, bins=bins, ax=THIS_AX)

                THIS_AX=ax7
                dLogs2 = gaussian_filter(dLogs,1)
                print(dLogs2.shape)
                dR_x,dR_y,dR_h,dR_dv,dR_p=heat_map.heat_map( dLogs2.transpose(), tcen2, bins=bins, ax=THIS_AX)
                MeanDlogs = dLogs2.mean(axis=1)
                THIS_AX.plot(tcen2,MeanDlogs,c='r')
                Filter = np.abs(dLogs)<5
                ModeDlogs = scipy.stats.mode(dLogs2*Filter,axis=1)[0]
                THIS_AX.plot(tcen2,ModeDlogs,c='b')

            if 1:
                THE_AX=ax1
                if 1:
                    #this one shows the bulk
                    bins = np.linspace(-5,5,32)

                #bins=np.geomspace(RB[RB>0].min(),RB.max())
                #bins=np.linspace(RB.min(),RB.max())
                PRB = splat(1-RB,tcenter,THE_AX,'1-RB',bins)
                THE_AX.set_yscale('linear')




            if 0:
                #Mode from CDF; all
                # Doesn't work well
                THE_AX=ax3
                #Q = (1-RB)[0,:]
                PeakQ=[]
                ax0.set(title='2')
                for Q in (1-RB):
                    Q.sort()
                    y = np.arange(Q.size)/Q.size

                    Median = Q[ np.argmin( np.abs(y-0.5))]
                    PeakQ.append(Median)
                    continue

                    #dQ = Q[1:]-Q[:-1]
                    #dy = y[1:]-y[:-1]
                    #Qc = 0.5*(Q[1:]+Q[:-1])
                    
                    Qtest = np.linspace(-10,10,128)
                    #iii = scipy.interpolate.interp( Q,y)
                    iii = np.interp(Qtest, Q,y)
                    iii = gaussian_filter(iii,8)
                    #ax3.plot( Qtest, iii(Qtest), c='k')
                    #ax3.plot( Qtest, iii, c='k')
                    #ax3.set(xlim=[-5,5])

                    dQt = Qtest[1:]-Qtest[:-1]
                    dQc = 0.5*(Qtest[1:]+Qtest[:-1])
                    dyt = iii[1:]-iii[:-1]
                    MaxQ = dQc[np.argmax(dyt/dQt)]
                    PeakQ.append(MaxQ)
                    #ax5.axvline(MaxQ)
                    #ax5.plot(dQc,50*dyt/dQt)

                    shift = Q[np.argmin( np.abs( y-0.5))]
                    #shift=0
                    #ax5.plot( Q-shift, y)
                ax0.plot(times, PeakQ,c='r')
                ax5.set(xlim=[-5,5])


            if 0:
                #CDF; play
                THE_AX=ax3
                Q = (1-RB)[0,:]
                Q.sort()
                y = np.arange(Q.size)/Q.size
                THE_AX.plot(Q,y)
                dQ = Q[1:]-Q[:-1]
                dy = y[1:]-y[:-1]
                Qc = 0.5*(Q[1:]+Q[:-1])
                #THE_AX.set(xlim=[-5,5])
                ax4.plot( Qc, dy/dQ)
                ax5.plot(Qc, dy/dQ)
                ax5.set(xlim=[-5,5])
                
                Qtest = np.linspace(-10,10,1024)
                #iii = scipy.interpolate.interp( Q,y)
                iii = np.interp(Qtest, Q,y)
                iii = gaussian_filter(iii,5)
                #ax3.plot( Qtest, iii(Qtest), c='k')
                ax3.plot( Qtest, iii, c='k')
                ax3.set(xlim=[-5,5])

                dQt = Qtest[1:]-Qtest[:-1]
                dQc = 0.5*(Qtest[1:]+Qtest[:-1])
                dyt = iii[1:]-iii[:-1]
                MaxQ = dQc[np.argmax(dyt/dQt)]
                ax5.axvline(MaxQ)
                ax5.plot(dQc,50*dyt/dQt)


            if 1:
                THE_AX=ax2

                #x = np.abs((RB).flatten())
                #y = np.abs((Stretch/(B2*divv)).flatten())
                #y = np.abs((Stretch).flatten())
                x = np.abs(Stretch).flatten()
                y = np.abs(divv).flatten()
                ext3=extents()
                ext3(x);ext3(y)
                bins =np.geomspace(ext3.minmax[0],ext3.minmax[1])
                pch.simple_phase(x,y,ax=THE_AX, bins=[bins,bins]) 
                THE_AX.plot(ext3.minmax,ext3.minmax)
                THE_AX.set(xscale='log',yscale='log', title='DivV vs S')

            if 1:
                #
                # B and rho
                #
                THIS_AX = ax6

                THIS_AX.plot(times, rho/rho_0, c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, BP/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='B/B0,rho/rho_0')
            if 1:
                #B rho phase
                THIS_AX=ax7
                ext=extents()

                y = (rho/rho_0).flatten()
                x = (B2/b_0).flatten()
                ext(x);ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins]) 
                THIS_AX.set(xscale='log',yscale='log',xlabel='B',ylabel='rho')
                THIS_AX.plot(ext.minmax,ext.minmax,c='k')

            if 1:
                THIS_AX = ax8

                RHO_FIX=(rho/rho_0)**np.abs((1-RB_Mean_10))

                THIS_AX.plot(times, RHO_FIX, c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, B2/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='Mean RB 10')
                THIS_AX=ax9
                x=B2/b_0
                y = RHO_FIX
                ext(x);ext(y)
                x=x.flatten();y=y.flatten()
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins]) 
                THIS_AX.plot(ext.minmax,ext.minmax)
                THIS_AX.set(xscale='log',yscale='log')
            if 0: 
                #B rho phase
                THIS_AX=ax9
                ext=extents()

                MyR=nar(PeakQ)
                MyR.shape=MyR.size,1
                y = ((rho/rho_0)**(1-MyR)).flatten()
                x = (B2/b_0).flatten()
                ext(x);ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins])
                THIS_AX.set(xscale='log',yscale='log',ylabel='rho',xlabel='B', title='Mean 1-R %0.2e'%(1-RB_Mean.mean()))
                THIS_AX.plot(ext.minmax,ext.minmax,c='k')


            fig.savefig('plots_to_sort/R_parts_%s_c%04d'%(sim,core_id))

sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt2(TL.loops[sim])
    core_list=None
    #core_list=[7]
    #core_list=[74, 353]
    core_list=[74]
    #core_list=[112]
    core_list=TL.loops[sim].core_by_mode['Alone']
    #core_list=core_list[:3]
    #core_list=[214, 114]
    P=ddd.run(core_list=core_list)


