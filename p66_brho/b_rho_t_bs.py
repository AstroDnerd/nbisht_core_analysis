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

            fig,ax_square=plt.subplots(2,2)
            ax=ax_square.flatten()

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




            bins_1 = np.geomspace( 1,1e9,19)
            bins_m1 = -bins_1[::-1]
            bins = nar(list(bins_m1)+list(bins_1))
            bincen = 0.5*(bins[1:]+bins[:-1])


            if 1:
                #
                #plot dbdt
                #
                THIS_AX=ax0
                smooth_b=ndimage.gaussian_filter1d(B2, 2, 0)/2
                db_dt = (smooth_b[1:,:]-smooth_b[:-1,:])/dt_square
                db_x,db_y,db_h,db_dv,db_p=heat_map.heat_map( db_dt.transpose(), tcenter, bins=bins, ax=THIS_AX)
                THIS_AX.set_yscale('symlog',linthresh=100)
                THIS_AX.set_title('db/dt')


            if 1:
                #
                #plot stretch, squish.
                #
                Pstretch=splat(Stretch,tcenter,ax1,'Stretch',bins)
                Psquish=splat(-B2/2*divv,tcenter,ax2,'Squish',bins)

                RedChan = np.log10(db_h.transpose())
                GreenChan = np.log10(Psquish[2].transpose())
                BlueChan = np.log10(Pstretch[2].transpose())
                oot=np.stack([RedChan,GreenChan,BlueChan],axis=2)
                #oot=np.stack([moo,mork+pork,mork+pork],axis=2)
                ax3.imshow(oot, origin='lower', extent=[0,1,0,1])

                oot=np.stack([RedChan,GreenChan+BlueChan,GreenChan+BlueChan],axis=2)
                #oot=np.stack([moo,mork+pork,mork+pork],axis=2)
                ax4.imshow(oot, origin='lower', extent=[0,1,0,1])


            if 1:
                #
                #the legend
                #
                THE_AX = ax5
                venn1 = np.zeros([800,800])
                venn2 = np.zeros([800,800])
                venn3 = np.zeros([800,800])
                x,y = np.mgrid[0:1:1/800, 0:1:1/800]-0.5
                r=0.25
                x1,y1 = r*np.cos(2*np.pi/3), r*np.sin(2*np.pi/3)
                ok = (x-x1)**2 + (y-y1)**2 < 1.2*r**2
                
                venn1[ok] = 1
                x2,y2 = r*np.cos(4*np.pi/3), r*np.sin(4*np.pi/3)
                ok = (x-x2)**2 + (y-y2)**2 < 1.2*r**2
                venn2[ok] = 1
                x3,y3 = r*np.cos(6*np.pi/3), r*np.sin(6*np.pi/3)
                ok = (x-x3)**2 + (y-y3)**2 < 1.2*r**2
                venn3[ok] = 1
                venn = np.stack([venn1,venn2,venn3],axis=2)
                THE_AX.imshow(venn,origin='lower', extent=[-0.5,0.5,-0.5,0.5])
                THE_AX.text( y1,x1, 'dBdt', color='white')
                THE_AX.text( y2,x2, 'B2', color='white')
                THE_AX.text( y3,x3, 'Stretch', color='white')


            if 1:
                #R take 1
                THE_AX=ax6
                #shoot, why is this BP?  Should be B2.
                RB=Stretch/(BP*divv)
                print(RB.min(),RB.max())
                if 0:
                    #dumb, use bins from above.
                    #Shows nothing, actually
                    bins=bins
                if 0:
                    #this one shows stuff
                    bins_1 = np.geomspace( 1,1e3,19)
                    bins_m1 = -bins_1[::-1]
                    bins = nar(list(bins_m1)+list(bins_1))
                if 1:
                    #this one shows the bulk
                    bins = np.linspace(-5,5,32)

                #bins=np.geomspace(RB[RB>0].min(),RB.max())
                #bins=np.linspace(RB.min(),RB.max())
                PRB = splat(1-RB,tcenter,THE_AX,'1-RB',bins)
                THE_AX.set_yscale('linear')

            if 0:
                #the full range of RB, just in case.
                #cumulative
                THE_AX=ax7
                QQ = (RB+0).flatten()
                QQ.sort()
                THE_AX.plot(QQ)
                THE_AX.set_title('sorted RB')
                THE_AX.set_yscale('symlog',linthresh=5)
                THE_AX.axhline(0, c=[0.5]*4)
                THE_AX.axhline(2, c=[0.5]*4)
                THE_AX.axhline(-2, c=[0.5]*4)



            if 1:
                #
                # d log B / d ln rho
                #
                THIS_AX=ax7
                smooth_log_b=ndimage.gaussian_filter1d(np.log(B2), 2, 0)/2
                dlogb_dt = (smooth_log_b[2:,:]-smooth_log_b[:-2,:])/dt2
                smooth_rho=ndimage.gaussian_filter1d(np.log(rho), 2, 0)
                drho_dt = (smooth_rho[2:,:]-smooth_rho[:-2,:])/dt2
                dLogs=(dlogb_dt/drho_dt)
                print("%0.2e %0.2e %s"%(dLogs.min(), dLogs.max(),'word'))


                if 0:
                    #get everything
                    bins_1 = np.geomspace( 1,np.abs(dLogs).max(),128)
                    bins_m1 = -bins_1[::-1]
                    bins = nar(list(bins_m1)+list(bins_1))
                if 1:
                    bins = np.linspace(-5,5,32)


                dR_x,dR_y,dR_h,dR_dv,dR_p=heat_map.heat_map( dLogs.transpose(), tcen2, bins=bins, ax=THIS_AX)
                THIS_AX.set_yscale('symlog',linthresh=100)
                THIS_AX.set_title('dlogB/dlogRho')
                THIS_AX.set_yscale('linear')


            if 0:
                #cumulative
                QQQ = dLogs.flatten()
                QQQ.sort()
                ax9.plot(QQQ)
                ax9.axhline(0,c=[0.5]*4)
                ax9.axhline(2,c=[0.5]*4)
                ax9.axhline(1,c=[0.5]*4)
                ax9.axhline(0.5,c=[0.5]*4)
                ax9.set_yscale('symlog',linthresh=1)
                ax9.set_title('dlnb/dlnrho')

            if 0:
                #cumulative
                Q = dlogb_dt.flatten()+0
                Q.sort()
                ax10.plot(Q)
                Q = drho_dt.flatten()+0
                Q.sort()
                ax10.plot(Q)
                ax10.set_title('dlnb and dlnrho')

            if 1:
                #phase diagram to see how right I am.
                THIS_AX = ax8
                x = dLogs.flatten()
                y = 1-RB[1:-1,:].flatten()
                bins = np.linspace(-5,5,128)
                print(x.shape,y.shape)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins])
                arr=[-2,2]
                THIS_AX.plot(arr,arr,c=[0.5]*4)
                THIS_AX.set(ylabel='1-R',xlabel=r'$\frac{d_t \ln |B|}{d_t \ln \rho}$')

            if 0:
                #
                # histograms
                #

                Hsquish = Psquish[2]
                b0=np.arange(0,5,5/16)
                b1=np.geomspace(5,10000,32)
                b2=np.concatenate([b0,b1])
                b3 = np.concatenate([-b2[::-1],b2])

                #ax12.plot(b3)

                rmap = rainbow_map( RB.shape[0])
                print('WAAAAAAAAA',RB.mean())
                for n in [0,50,-1]:
                    #if n%10:
                    #    continue
                    H=1-RB[n,:]
                    print(H.mean())
                    #ax12.hist(H, bins=np.linspace(-10,10,64), histtype='step',color=rmap(n))
                    ax12.hist(H, bins=np.linspace(-5,5,64), histtype='step',color=rmap(n))
                    #H.sort()
                    #v = np.arange(H.size)/H.size
                    #ax12.plot(H,v)
                    #break
                ax12.set(yscale='linear',xscale='linear',title='1-RB')


            if 1:
                #
                # B and rho
                #
                THIS_AX = ax12
                rho_0=rho[0,:].mean()
                b_0=BP[0,:].mean()

                THIS_AX.plot(times, rho/rho_0, c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, BP/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='B/B0,rho/rho_0')
            if 1:
                THIS_AX = ax13
                rho_0=rho[0,:].mean()
                b_0=BP[0,:].mean()
                #RB_Mean = np.abs((rho*RB*dv).sum(axis=0)/(rho*dv).sum(axis=0))
                RB_Mean = (rho*RB*dv).sum()/(rho*dv).sum()
                print("ONE MINUS RB MEAN",1-RB_Mean)

                THIS_AX.plot(times, (rho/rho_0)**np.abs((1-RB_Mean)), c=[0.5]*4, linewidth=0.1)
                THIS_AX.plot(times, BP/b_0, c=[1.0,0.0,0.0,0.5], linewidth=0.1)
                THIS_AX.set(yscale='log',title='B/B0,rho/rho_0')

            if 1:
                #B rho phase
                THIS_AX=ax14
                ext=extents()
                rho_0=rho[0,:].mean()
                b_0=B2[0,:].mean()

                x = (rho/rho_0).flatten()
                y = (B2/b_0).flatten()
                ext(x);ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins])
                THIS_AX.set(xscale='log',yscale='log',xlabel='rho',ylabel='B')
                THIS_AX.plot(ext.minmax,ext.minmax,c='k')
            if 1:
                #B rho phase
                THIS_AX=ax15
                ext=extents()
                rho_0=rho[0,:].mean()
                b_0=B2[0,:].mean()

                x = ((rho/rho_0)**np.abs(1-RB_Mean)).flatten()
                y = (B2/b_0).flatten()
                ext(x);ext(y)
                bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
                pch.simple_phase(x,y,ax=THIS_AX, bins=[bins,bins])
                THIS_AX.set(xscale='log',yscale='log',xlabel='rho',ylabel='B')
                THIS_AX.plot(ext.minmax,ext.minmax,c='k')

            if 0:
                THIS_AX=ax16
                THIS_AX.plot( times, rho)



            #
            # OLD IDEAS
            #
            if 0:
                bins_1 = np.geomspace( 1,1e3,19)
                bins_m1 = -bins_1[::-1]
                bins = nar(list(bins_m1)+list(bins_1))
                goodbins=bins
                bincen = 0.5*(bins[1:]+bins[:-1])


            if 0:
                RB=np.abs(Stretch/(B*divv))
                bins=np.geomspace(RB[RB>0].min(),RB.max(),16)
                output=ax0.hist( RB.flatten(), bins=bins)
                print(output[0])
                #phist(RB.flatten())
            if 0:
                RB=np.abs(B*divv/Stretch)
                bins=np.geomspace(RB[RB>0].min(),RB.max(),16)
                b = RB.flatten()
                b.sort()
                ax0.plot(b)
                ax0.set(yscale='log')
            if 0:
                x =np.abs((B*divv).flatten())
                y =np.abs((Stretch).flatten())

                pch.simple_phase( x,y,ax=ax3,log=True,nBins=128)
                #ax3.plot([1e3,1e9],[1e3,1e9])
                ax3.set(xscale='log',yscale='log', xlabel='B^2 divv',ylabel='BdotS',title='All time')

            if 0:
                #smooth_b=ndimage.gaussian_filter1d(B, 2, 0)
                #db_dt = (smooth_b[1:,:]-smooth_b[:-1,:])/dt_square

                R  = -(Stretch/B*divv)#[0,:]
                ok = R<0
                x=np.abs(R[ok]).flatten()
                y=np.abs(Stretch[ok]).flatten()
                pch.simple_phase(x,y,ax=ax4,log=True,nBins=128)
                ax4.set(xscale='log',yscale='log',xlabel='R',ylabel='stretch<0')
                ok = R>0
                x=np.abs(R[ok]).flatten()
                y=np.abs(Stretch[ok]).flatten()
                pch.simple_phase(x,y,ax=ax2,log=True,nBins=128)
                ax2.set(xscale='log',yscale='log',xlabel='R',ylabel='stretch>0')

            if 0:
                #bins=np.linspace(-100,100,128)
                bins_1 = np.geomspace( 1,1e6,19)
                bins_m1 = -bins_1[::-1]
                bins = nar(list(bins_m1)+list(bins_1))
                goodbins=bins
                bincen = 0.5*(bins[1:]+bins[:-1])
                splat(R,tcenter,ax5,'RB',goodbins)
                #ax5.axhline(-1,c=[0.5]*4)
                #ax5.axhline(-0.5,c=[0.5]*4)
                #ax5.axhline(-2,c=[0.5]*4)
                #ax5.set_yscale('linear')

            if 0:
                bins_1 = np.geomspace( 1,1e6,19)
                bins_m1 = -bins_1[::-1]
                bins = nar(list(bins_m1)+list(bins_1))
                bins=np.linspace(-1,1,128)

                bincen = 0.5*(bins[1:]+bins[:-1])
                j=1/Stretch
                #output=splat(j,tcenter,ax6,'RB',bins)
                #pdb.set_trace()
                ax6.hist( np.log10(np.abs(j.flatten())))
                #ax5.axhline(-1,c=[0.5]*4)
                #ax5.axhline(-0.5,c=[0.5]*4)
                #ax5.axhline(-2,c=[0.5]*4)
                #ax5.set_yscale('linear')

            #pdb.set_trace()


            fig.savefig('plots_to_sort/db_dt_parts_%s_c%04d'%(sim,core_id))




sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt2(TL.loops[sim])
    core_list=None
    #core_list=[7]
    core_list=[74, 112]
    #core_list=[74]
    #core_list=[112]
    #core_list=TL.loops[sim].core_by_mode['Alone']
    P=ddd.run(core_list=core_list)


