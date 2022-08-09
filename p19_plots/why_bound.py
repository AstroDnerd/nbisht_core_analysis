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
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)
        ms.compute_ke_rel(core_id)
        for nframe,frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            nf = np.where( this_looper.tr.frames == frame)[0][0]
            xtra_energy.add_energies(ds)

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
            y_ext(nar([1e-2,1e6]))


            fig,ax = plt.subplots(1,1)
            ax3=ax

            #relative kinetic energy.
            vx = sp[YT_velocity_x].v - ms.mean_vx[nf]
            vy = sp[YT_velocity_y].v - ms.mean_vy[nf]
            vz = sp[YT_velocity_z].v - ms.mean_vz[nf]
            EK = 0.5*DD*(vx*vx+vy*vy+vz*vz)

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
            #ax3.hist( EG_cuml, bins = np.geomspace( EG_cuml.min(), EG_cuml.max(), 64) )
            ax3.plot(  RR_cuml, EG_cuml, c='k')
            ax3.plot( RR_cuml, EK_cuml, c='r')
            print(EG_cuml)
            axbonk(ax3,xscale='log',yscale='log', ylim=y_ext.minmax)




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
