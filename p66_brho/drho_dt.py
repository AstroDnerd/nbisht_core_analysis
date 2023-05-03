from starter2 import *

import three_loopers_u900 as TL

from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')

class dq_dt():
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

            dt_square = dt+0
            dt_square.shape = dt_square.size,1


        if core_list is None:
            core_list = sorted(np.unique( thtr.core_ids))

        for core_id in core_list:
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
        rho = rho[mask,:]
        B=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
        B = B**2/2
        divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]

        fig, ax=plt.subplots(3,6,figsize=(20,12))
        ax0=ax[0][0];ax1=ax[0][1] 
        ax2=ax[1][0];ax3=ax[1][1]
        ax4=ax[2][0];ax5=ax[2][1]
        ax6=ax[0][2];
        ax7=ax[1][2]
        ax8=ax[2][2]
        ax9=ax[0][3]
        ax10=ax[1][3]
        ax11=ax[2][3]
        ax12=ax[0][4]
        ax13=ax[1][4]
        ax14=ax[2][4]
        ax15=ax[0][5]
        ax16=ax[1][5]
        ax17=ax[2][5]


        if 0:
            #unsmoothed.  Noisy
            ax0.plot(times , rho, c=c, linewidth=0.1)
            axbonk(ax0,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho$',yscale='log')#, ylim=[rho_min,rho_max])

            c2 = copy.copy(c)
            c2[0]=0
            ax0.plot(times, B, c=c2,linewidth=0.1)
            #axbonk(ax, xlabel=r'$t/t_{ff}$',ylabel='B', yscale='log')


            ##ax2.plot( dt, drho_dt, c=c,linewidth=0.1)
            #ax2.plot(times[1:]-times[:-1])
            #ax2.set( ylim = [dt_square.mean()*(0.99), dt_square.mean()*1.01])
            drho_dt = (rho[1:, :] - rho[:-1,:])/dt_square
            ax2.plot( tcenter, drho_dt, c=c, linewidth=0.1)
            ax2.set_yscale('symlog',linthresh=0.5)
        if 0:
            for frame in nf:
                
                the_x=dsmooth_dt[frame,:]
                ax1.hist( the_x, bins=bins, histtype='step', color=rm(frame), density=True)
            ax1.set_xscale('symlog',linthresh=10)
            ax1.set_yscale('log')

        nf = np.arange(tcenter.size)
        rm = rainbow_map(tcenter.size)

        smooth=ndimage.gaussian_filter1d(rho, 2, 0)
        ax0.plot( times, smooth, c=c,linewidth=0.1)
        ax0.set(yscale='log', ylabel='rho smooth')

        dsmooth_dt = (smooth[1:,:]-smooth[:-1,:])/dt_square
        drho_dt=dsmooth_dt
        ax2.plot( tcenter, dsmooth_dt, c=c,linewidth=0.1)
        ax2.set_yscale('symlog',linthresh=100)
        ax2.set(ylabel='drho/dt')


        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        drdt_x, drdt_y, drdt_h, drdt_dv, drdt_ploot=\
                heat_map.heat_map( dsmooth_dt.transpose(), tcenter, bins=bins, ax=ax4)
        ax4.set_yscale('symlog',linthresh=100)

        twin = ax4.twinx()
        fraction_positive = ( dsmooth_dt > 0).sum(axis=1)/dsmooth_dt.shape[0]
        twin.plot(tcenter, fraction_positive)
        twin.set(ylim=[0,1])

        #
        #
        #

        smooth=ndimage.gaussian_filter1d(B, 2, 0)
        ax1.plot( times, smooth, c=c,linewidth=0.1)
        ax1.set(yscale='log', ylabel='B smooth')

        dsmooth_dt = (smooth[1:,:]-smooth[:-1,:])/dt_square
        dB_dt = dsmooth_dt
        ax3.plot( tcenter, dsmooth_dt, c=c,linewidth=0.1)
        ax3.set_yscale('symlog',linthresh=100)
        ax3.set(ylabel='d/dt')


        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        db_x,db_y,db_h,db_dv,db_p=heat_map.heat_map( dsmooth_dt.transpose(), tcenter, bins=bins, ax=ax5)
        ax5.set_yscale('symlog',linthresh=100)

        twin = ax5.twinx()
        fraction_positive = ( dsmooth_dt > 0).sum(axis=1)/dsmooth_dt.shape[0]
        twin.plot(tcenter, fraction_positive)
        twin.set(ylim=[0,1])

        #
        #
        #
        smooth= ndimage.gaussian_filter1d(divv, 2, 0)
        #ax1.plot( times, smooth, c=c,linewidth=0.1)
        #ax1.set(yscale='log', ylabel='B smooth')

        ax7.plot( times, smooth, c=c,linewidth=0.1)
        ax7.set_yscale('symlog',linthresh=100)
        ax7.set(ylabel='d/dt')


        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        heat_map.heat_map( smooth.transpose(), times, bins=bins, ax=ax8)
        ax8.set_yscale('symlog',linthresh=100)

        twin = ax8.twinx()
        fraction_positive = ( smooth > 0).sum(axis=1)/dsmooth_dt.shape[0]
        twin.plot(times, fraction_positive)
        twin.set(ylim=[0,1])
        #
        #
        #
        smooth= ndimage.gaussian_filter1d(-1*rho*divv, 2, 0)
        smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
        #ax1.plot( times, smooth, c=c,linewidth=0.1)
        #ax1.set(yscale='log', ylabel='B smooth')

        ax10.plot( tcenter, smooth, c=c,linewidth=0.1)
        ax10.set_yscale('symlog',linthresh=100)
        ax10.set(ylabel='d/dt')


        c2=copy.copy(c)
        c2[0]=0
        ax10.plot( tcenter, drho_dt, c=c2, linewidth=0.1)


        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        dr_x,dr_y,dr_h,dr_dv,dr_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax11)
        ax11.set_yscale('symlog',linthresh=100)
        ax11.contour(drdt_x, drdt_y, drdt_h)

        twin = ax11.twinx()
        fraction_positive = ( smooth > 0).sum(axis=1)/dsmooth_dt.shape[0]
        twin.plot(tcenter, fraction_positive)
        twin.set(ylim=[0,1])
        #
        #
        #
        fakeit=np.log10(dr_h.transpose())
        borkus=np.log10(drdt_h.transpose())
        dumb = np.zeros_like(fakeit)
        dumb=fakeit
        oot = np.stack([fakeit,borkus,dumb],axis=2)
        ax9.imshow(oot)

        #
        #
        #
        no="""
        smooth= ndimage.gaussian_filter1d(-1*B*divv, 2, 0)
        smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
        #ax1.plot( times, smooth, c=c,linewidth=0.1)
        #ax1.set(yscale='log', ylabel='B smooth')

        ax13.plot( tcenter, smooth, c=c,linewidth=0.1)
        ax13.set_yscale('symlog',linthresh=100)
        ax13.set(ylabel='d/dt')


        c2=copy.copy(c)
        c2[0]=0
        ax13.plot( tcenter, dB_dt, c=c2, linewidth=0.1)


        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        dbb_x,dbb_y,dbb_h,dbb_dv,dbb_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax14)
        ax14.set_yscale('symlog',linthresh=100)
        #ax14.contour(drdt_x, drdt_y, drdt_h)


        moo = np.log10(db_h.transpose())
        dumb = np.zeros_like(moo)
        mork = np.log10(dbb_h.transpose())
        oot=np.stack([moo,mork,dumb],axis=2)
        ax12.imshow(oot)

        #
        #
        #
        bx = thtr.c([core_id],'magnetic_field_x')
        by = thtr.c([core_id],'magnetic_field_y')
        bz = thtr.c([core_id],'magnetic_field_z')
        vx = thtr.c([core_id],'velocity_x')
        vy = thtr.c([core_id],'velocity_y')
        vz = thtr.c([core_id],'velocity_z')
        dxvx = thtr.c([core_id],'dxvx')
        dxvy = thtr.c([core_id],'dxvy')
        dxvz = thtr.c([core_id],'dxvz')
        dyvx = thtr.c([core_id],'dyvx')
        dyvy = thtr.c([core_id],'dyvy')
        dyvz = thtr.c([core_id],'dyvz')
        dzvx = thtr.c([core_id],'dzvx')
        dzvy = thtr.c([core_id],'dzvy')
        dzvz = thtr.c([core_id],'dzvz')
        Sx = bx*dxvx+by*dyvx+bz*dzvx
        Sy = 0 #bx*dxvy+by*dyvy+bz*dzvy
        Sz = 0 #bx*dxvz+by*dyvz+bz*dzvz
        Stretch= bx*Sx+by*Sy+bz*Sz
        Stretch = Stretch[sl].transpose()[mask,:]

        smooth= ndimage.gaussian_filter1d(-1*Stretch, 2, 0)
        smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
        #ax1.plot( times, smooth, c=c,linewidth=0.1)
        #ax1.set(yscale='log', ylabel='B smooth')

        ax15.plot( tcenter, smooth, c=c,linewidth=0.1)
        ax15.set_yscale('symlog',linthresh=100)
        ax15.set(ylabel='d/dt')

        #c2=copy.copy(c)
        #c2[0]=0
        #ax13.plot( tcenter, dB_dt, c=c2, linewidth=0.1)

        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax16)
        ax16.set_yscale('symlog',linthresh=100)
        #ax14.contour(drdt_x, drdt_y, drdt_h)


        #moo = np.log10(db_h.transpose())
        #dumb = np.zeros_like(moo)
        #mork = np.log10(dbb_h.transpose())
        #oot=np.stack([moo,mork,dumb],axis=2)
        #ax12.imshow(oot)
        fig.savefig('plots_to_sort/b_and_rho_2_c%04d.pdf'%(core_id))
        """

        fig.savefig('plots_to_sort/d_rho_dt_c%04d.pdf'%(core_id))

sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt(TL.loops[sim])
    core_list=[7]
    #core_list=[74]
    ddd.run(core_list=core_list)


