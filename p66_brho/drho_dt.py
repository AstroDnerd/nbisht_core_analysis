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

        #if ms.nparticles < 1000:
        #    sl=slice(None)
        #else:
        #    sl = slice(None,None,10)
        #    #c=[0,0,0,0.1]
        #    c=[0.1]*4
        c=[0.5]*4
        sl=slice(None)
        rho = ms.density[sl].transpose()
        rho = rho[mask,:]
        divv=thtr.c([core_id],'velocity_divergence')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
        #B=thtr.c([core_id],'magnetic_field_strength')[sl].transpose()[mask,:]#/colors.mean_field[this_looper.sim_name]
        #B = B**2/2

        fig, ax=plt.subplots(3,3,figsize=(12,12))
        ax0=ax[0][0];
        ax1=ax[0][1] 
        ax2=ax[0][2];


        ax3=ax[1][0]
        ax4=ax[1][1];
        ax5=ax[1][2]
        ax6=ax[2][0];
        ax7=ax[2][1]
        ax8=ax[2][2]
        #ax9=ax[0][3]
        #ax10=ax[1][3]
        #ax11=ax[2][3]
        #ax12=ax[0][4]
        #ax13=ax[1][4]
        #ax14=ax[2][4]
        #ax15=ax[0][5]
        #ax16=ax[1][5]
        #ax17=ax[2][5]



        rho_smooth=ndimage.gaussian_filter1d(rho, 2, 0)
        ax0.plot( times, rho_smooth, c=c,linewidth=0.1)
        ax0.set(yscale='log', ylabel='rho smooth')

        dsmooth_dt = (rho_smooth[1:,:]-rho_smooth[:-1,:])/dt_square
        ax1.plot( tcenter, dsmooth_dt, c=c,linewidth=0.1)
        ax1.set_yscale('symlog',linthresh=100)
        ax1.set(ylabel='drho/dt, finite')

        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        drdt_x, drdt_y, drdt_h, drdt_dv, drdt_ploot=\
                heat_map.heat_map( dsmooth_dt.transpose(), tcenter, bins=bins, ax=ax2)
        ax2.set_yscale('symlog',linthresh=100)
        ax2.set(ylabel='drho/dt finite')

        #twin = ax4.twinx()
        #fraction_positive = ( dsmooth_dt > 0).sum(axis=1)/dsmooth_dt.shape[0]
        #twin.plot(tcenter, fraction_positive)
        #twin.set(ylim=[0,1])


        #
        #
        #
        divv_smooth= ndimage.gaussian_filter1d(divv, 2, 0)
        ext=extents()
        ext(np.abs(divv_smooth))

        ax3.plot( times, divv_smooth, c=c,linewidth=0.1)
        ax3.set_yscale('symlog',linthresh=100)
        ax3.set(ylabel='div v')


        bins_1 = np.geomspace( 1,ext.minmax[1],19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        heat_map.heat_map( divv_smooth.transpose(), times, bins=bins, ax=ax4)
        ax4.set_yscale('symlog',linthresh=100)
        ax4.set(ylabel='div v')


        #twin = ax8.twinx()
        #fraction_positive = ( smooth > 0).sum(axis=1)/dsmooth_dt.shape[0]
        #twin.plot(times, fraction_positive)
        #twin.set(ylim=[0,1])

        #
        #
        #
        drdt_smooth= ndimage.gaussian_filter1d(-1*rho*divv, 2, 0)
        drdt_smooth_cen = 0.5*(drdt_smooth[1:,:]+drdt_smooth[:-1,:])

        bins_1 = np.geomspace( 1,1e9,19)
        bins_m1 = -bins_1[::-1]
        bins = nar(list(bins_m1)+list(bins_1))

        dr_x,dr_y,dr_h,dr_dv,dr_p=heat_map.heat_map( drdt_smooth.transpose(), tcenter, bins=bins, ax=ax5)
        ax5.set_yscale('symlog',linthresh=100)
        ax5.contour(drdt_x, drdt_y, drdt_h)
        ax5.set(ylabel='-rho divv')

        #ax5.plot( times, drdt_smooth, c=c,linewidth=0.1)
        #ax5.set_yscale('symlog',linthresh=100)
        #ax1.set(yscale='log', ylabel='B smooth')


        #ax7.plot( tcenter, drdt_smooth_cen, c=c,linewidth=0.1)
        #ax7.set_yscale('symlog',linthresh=100)
        #ax7.set(ylabel='d/dt')
        x=np.abs(drdt_smooth_cen.flatten())
        y=np.abs(dsmooth_dt.flatten())
        pch.simple_phase(x,y,log=True, ax=ax6)
        ax6.set(xscale='log',yscale='log', xlabel='-rho divv', ylabel='drho/dt')
        ext=extents()
        ext(x);ext(y)
        ax6.plot(ext.minmax,ext.minmax,c='r')


        #c2=copy.copy(c)
        #c2[0]=0
        #ax10.plot( tcenter, drho_dt, c=c2, linewidth=0.1)



        #twin = ax11.twinx()
        #fraction_positive = ( smooth > 0).sum(axis=1)/dsmooth_dt.shape[0]
        #twin.plot(tcenter, fraction_positive)
        #twin.set(ylim=[0,1])

        log_drhodt = np.zeros_like(drdt_h.transpose())
        ddd = drdt_h.transpose()
        log_drhodt[ ddd>0] = ddd[ddd>0]
        log_rhodivv = np.zeros_like(dr_h.transpose())
        rhodivv = dr_h.transpose()
        log_rhodivv[ rhodivv>0] = rhodivv[rhodivv>0]

        ext_both=extents()
        ext_both(log_rhodivv); ext_both(log_drhodt)
        scale = ext_both.minmax[1]/10


        oot = np.stack([log_drhodt/scale, log_drhodt/scale, log_rhodivv/scale],axis=2)
        print(oot.min(), oot.max())
        ax7.imshow(oot, origin='lower')
        ax7.set_aspect( ax1.get_aspect())


        fig.tight_layout()
        fig.savefig('plots_to_sort/d_rho_dt_c%04d.pdf'%(core_id))

sim_list=['u902']
for sim in sim_list:
    ddd = dq_dt(TL.loops[sim])
    #core_list=[7]
    core_list=[74]
    ddd.run(core_list=core_list)


