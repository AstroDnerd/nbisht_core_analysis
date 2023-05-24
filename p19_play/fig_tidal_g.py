from starter2 import *
import xtra_energy

from scipy import stats
import core_proj_three
reload(core_proj_three)
import other_scrubber
reload(other_scrubber)
#import three_loopers_six as TL
import camera_path
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']

if 0:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        core_proj_three.core_proj_multiple(loop,axis_list=[0],core_list=[74],frame_list=[0,10],camera=camera, main_core=74)


class multi_profile():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
    def run(self, core_list=None, frame_list=None, tsing=None):
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        Nplots = 4
        Ntimes = 4
        fig,axes = plt.subplots(Nplots,Ntimes, figsize=(8,12))
        fig.subplots_adjust(hspace=0,wspace=0)
        fig2,axes2 = plt.subplots(2,Ntimes, figsize=(8,12))
        #fig2.subplots_adjust(hspace=0,wspace=0)
        ext = [extents() for n in range(Nplots+1)]
        if Nplots > 4:
            rho_ext = extents()
            for core_id in core_list:
                ms = trackage.mini_scrubber(this_looper.tr,core_id)
                rho_ext(ms.density)
        for core_id in core_list:
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            #frame_mask[0]=True
            #frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tend_core[core_id])]=True
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))
            axes[0][0].set_title(r'$t=0$')
            axes[0][1].set_title(r'$t=0.9 t_{\rm{sing}}$')
            axes[0][2].set_title(r'$t=t_{\rm{sing}}$')
            axes[0][3].set_title(r'$t=t_{\rm{end}}$')

            img_collector=[]
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gravity(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)

                dv = sph[YT_cell_volume]
                RR = sph['radius']
                DD = sph[YT_density]
                ORDER = np.argsort( RR)
                rho_sort = DD[ORDER]
                RR_sort = RR[ORDER]
                dv_sort = dv[ORDER]
                M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
                d2_cuml = np.cumsum( DD[ORDER]**2*dv[ORDER])
                V_cuml = np.cumsum( dv[ORDER])
                EG = sph[YT_grav_energy_2]
                #EK = scrub.ke_rel
                EG_cuml = np.abs(np.cumsum( EG[ORDER]*dv[ORDER]))/V_cuml
                #EK_cuml = np.cumsum( EK[ORDER]*dv[ORDER])/V_cuml

                scrub=other_scrubber.scrubber(sph,do_velocity=False)

                gx = sph[YT_acceleration_x]
                gy = sph[YT_acceleration_y]
                gz = sph[YT_acceleration_z]

                gr = np.abs(scrub.rx_hat*gx+scrub.ry_hat*gy+scrub.rz_hat*gz)
                gt = np.sqrt(gx**2+gy**2+gz**2)
                rr = scrub.r
                r_bins=np.geomspace(1/2048,rr.max(),64)


                total_g_radial, bin_edge, bin_num = stats.binned_statistic( scrub.r, gr*dv, statistic='sum',bins=r_bins)
                total_rho_radial, bin_edge, bin_num = stats.binned_statistic( scrub.r, scrub.density*dv, statistic='sum',bins=r_bins)
                bin_cent = 0.5*(bin_edge[1:]+bin_edge[:-1])

                axes2[0][nframe].plot( bin_cent, total_g_radial)
                axes2[0][nframe].plot( bin_cent, 4*np.pi*colors.G*np.cumsum(total_rho_radial))
                axes2[0][nframe].plot( RR_sort, M_cuml)
                axes2[0][nframe].set(yscale='log',xscale='log')


                args = {'linewidth':0.2, 'c':[0.5]*4}
                ext[-1](RR_sort)
                #axes[0][nframe].plot( RR_sort, rho_sort, c=[0.5]*4)
                axes[0][nframe].plot( RR_sort, M_cuml, **args)
                ext[0](M_cuml)

                #axes[1][nframe].plot(gy)
                #axes[1][nframe].plot(gz)
                #axes[1][nframe].plot( RR_sort, np.sqrt(EG_cuml*8*np.pi*colors.G))
                #axes[1][nframe].plot( RR_sort, GR[ORDER]/V_cuml)

                continue







                if Nplots > 4:
                    bins = np.geomspace(rho_ext.minmax[0], rho_ext.minmax[1],128)
                    rho_to_hist = DD+0
                    rho_to_hist.sort()
                    cuml = np.arange(rho_to_hist.size)/rho_to_hist.size
                    #axes[4][0].plot( rho_to_hist, cuml)
                    axes[4][nframe].hist( DD.v, weights=dv.v, histtype='step', bins=bins, density=True)
                if 'hru' not in dir():
                    hru=0
            #ooo='plots_to_sort/incremental_c%04d'%core_id
            #fig.savefig(ooo)
            #print(ooo)

        #for ax in axes[0]:
        #    ax.set(xscale='log',yscale='log',ylim=ext[0].minmax, ylabel=r'$<\rho>(<r)$', xlim=ext[-1].minmax)
        #for ax in axes[1]:
        #    ax.set(xscale='log',yscale='log')
        for row in axes:
            for ax in row[1:]:
                ax.set(ylabel='', yticks=[])
        #for ax in axes[-1]:
        #    ax.set(xlabel='r')
        fig.savefig('plots_to_sort/radial_profile_%s'%(this_looper.sim_name))
        fig2.savefig('plots_to_sort/r2adial_profile_%s'%(this_looper.sim_name))



import tsing
reload(tsing)
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
if 'multi_proj' not in dir() or True :
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[:3]
        #core_list=core_list[:10]
        core_list=[68]
        #core_list
        #core_list=core_list[10:]
        #core_list=[114]


        #core_list = [114]
        #core_list=[361]
        #core_list=[8]
        #core_list=[381]
        #core_list=[323]
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        #frame_list = [100]

        mp=multi_profile(TL.loops[sim])
        mp.run(core_list=core_list,tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
