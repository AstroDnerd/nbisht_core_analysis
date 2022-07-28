from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
G = colors.G
def why(this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    if frame_list is None:
        frame_list=[this_looper.target_frame]

    rm = rainbow_map(len(frame_list))
    for core_id in core_list:
        fig_radials, ax_radials = plt.subplots(1,1)
        ax_radials.axhline(1, c=[0.5]*4)
        for nframe,frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            nf = np.where( this_looper.tr.frames == frame)[0][0]
            xtra_energy.add_energies(ds)

            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            ms.compute_ge(core_id)
            ms.compute_ke(core_id)
            ms.compute_ke_rel(core_id)
            c = nar([ms.mean_x[nf], ms.mean_y[nf],ms.mean_z[nf]])
            
            rsph = ds.arr(8.0/128,'code_length')
            sp = ds.sphere(c,rsph)


            #Get data arrays

            dv = np.abs(sp[YT_cell_volume])
            RR =sp[YT_radius]
            DD = sp[YT_density]
            EG = np.abs(sp[YT_grav_energy_2])
            EK = np.abs(sp[YT_kinetic_energy])


            y_ext = extents()
            if 0:
                #gets everything
                y_ext(EK)
                y_ext(EG)
            else:
                y_ext(nar([1e-2,1e6]))

            #Histograms

            ebins = np.geomspace( y_ext.minmax[0], y_ext.minmax[1],65)
            rbins = np.geomspace( RR [RR >0].min(), RR .max(),67)

            fig,ax = plt.subplots(2,2)
            ax0=ax[0][0]; ax1=ax[0][1]
            ax2=ax[1][0]; ax3=ax[1][1]
            #ax4=ax[1][2]; ax5=ax[0][2]

            #relative kinetic energy.
            vx = sp[YT_velocity_x].v - ms.mean_vx[nf]
            vy = sp[YT_velocity_y].v - ms.mean_vy[nf]
            vz = sp[YT_velocity_z].v - ms.mean_vz[nf]
            NRK = 0.5*DD*(vx*vx+vy*vy+vz*vz)

            if 0:
                EK_EK, xbek, ybek = np.histogram2d( EK , NRK, bins=[ebins,ebins],weights=dv)
                stuff = pch.helper(EK_EK, xbek, ybek, transpose=False, ax=ax5)
                ax5.plot( ebins, ebins)
                axbonk(ax5,xscale='log',yscale='log')


            ####
            EK=NRK
            #####

            r_cen = 0.5*(rbins[1:]+rbins[:-1]) #we'll need this later.
            EG_hist, xb, yb = np.histogram2d( RR , EG, bins=[rbins,ebins],weights=dv)
            EK_hist, xb, yb = np.histogram2d( RR , EK, bins=[rbins,ebins],weights=dv)

            EG_EK, xbgk, ybgk = np.histogram2d( EG , EK, bins=[ebins,ebins],weights=dv)



            stuff = pch.helper(EG_hist, xb, yb, transpose=False, ax=ax0)
            stuff = pch.helper(EK_hist, xb, yb, transpose=False, ax=ax1)



            axbonk(ax0,xlabel='R',ylabel='EG',xscale='log',yscale='log')
            axbonk(ax1,xlabel='R',ylabel='EK',xscale='log',yscale='log')

            if 1:
                stuff = pch.helper(EG_EK, xbgk, ybgk, transpose=False, ax=ax2)
                ax2.plot( ebins, ebins)
                axbonk(ax2,xscale='log',yscale='log', xlabel='EG',ylabel='EK')

            if 0:
                # Contours.
                # Works ok.
                xc = 0.5*(xb[1:]+xb[:-1])
                yc = 0.5*(yb[1:]+yb[:-1])
                ax2.contour(xc, yc,  np.log10(EG_hist).transpose() , colors='k')
                ax2.contour(xc, yc,  np.log10(EK_hist).transpose() , colors='r')
                #ax2.contour( EK_hist )
                axbonk(ax2, yscale='log',xscale='log')

            #CUMULATIVE

            ORDER = np.argsort( RR)
            RR_cuml = RR[ORDER]
            M_local = DD[ORDER]*dv[ORDER]
            M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
            V_cuml = np.cumsum( dv[ORDER])
            V_local = dv[ORDER]
            EG_cuml = np.cumsum( EG[ORDER]*V_local)/V_cuml
            EK_cuml = np.cumsum( EK[ORDER]*V_local)/V_cuml

            RAT = EK_cuml/EG_cuml
            bound=np.argmin(RAT)
            ax3.scatter( RR_cuml[bound], EK_cuml[bound])
            ax3.text(RR_cuml[bound], EK[bound], "%0.1e"%RAT[bound])
            ax3.hist( EG_cuml, bins = np.geomspace( EG_cuml.min(), EG_cuml.max(), 64) )
            #ax3.plot(  RR_cuml, EG_cuml, c='k')
            #ax3.plot( RR_cuml, EK_cuml, c='r')
            print(EG_cuml)
            axbonk(ax3,xscale='log',yscale='log', ylim=y_ext.minmax)

            ax_radials.plot( RR_cuml, EK_cuml/EG_cuml, c=rm(nframe))


            #particles.
            stuff={'edgecolor':'r','s':30, 'facecolor':'None'}
            ax0.scatter(ms.r[:,nf], np.abs(ms.ge[:,nf]), **stuff)
            ax1.scatter(ms.r[:,nf], np.abs(ms.ke_rel[:,nf]), **stuff)
            ax2.scatter(np.abs(ms.ge[:,nf]),ms.ke_rel[:,nf], **stuff)

            #histograms
            if 0:
                ax4.hist( EG.v, histtype='step', weights=dv.v, bins=ebins, label='Eg', color='k')
                ax4.hist( EK.v, histtype='step', weights=dv.v, bins=ebins, label='Ek', color='r')
                axbonk(ax4,xscale='log')

            fig.suptitle('%s c%04d %s'%(this_looper.sim_name,core_id,str(this_looper.mode_dict[core_id])))
            outname='plots_to_sort/why_%s_c%04d_n%04d.png'%(this_looper.sim_name, core_id,frame)
            fig.savefig(outname)

            plt.close(fig)
            print(outname)
        axbonk(ax_radials,yscale='log',ylim=[1./50,50], ylabel='EK(<)/EG(<)', xscale='log',xlabel='r')
        fig_radials.savefig('plots_to_sort/radials_%s_c%04d.png'%(this_looper.sim_name,core_id))
        plt.close(fig_radials)


import three_loopers_six as TL


import anne
reload(anne)
#anne.make_inflection()
sim_list=['u601','u602','u603']
sim_list=['u602']
if 'ge' not in dir() or True:
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        #core_list=[8]
        core_list=[381]
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        frame_list = [100]

        why( TL.loops[sim],core_list=core_list, do_plots=True, frame_list=frame_list)#, r_inflection=anne.inflection[sim])
