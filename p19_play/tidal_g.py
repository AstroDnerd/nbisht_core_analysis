from starter2 import *
import xtra_energy

import core_proj_three
reload(core_proj_three)
import other_scrubber
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


class proj_and_profile():
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
        for core_id in core_list:
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            frame_mask[get_time_index(tsing.tend_core[core_id])]=True
            frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
            frame_mask[0]=True
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))

            fig,axes = plt.subplots(3,len(frame_list))
            img_collector=[]
            ext_rho=extents()
            ext_rho2=extents()
            ext_r=extents()
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_gravity(ds)
                #xtra_energy.add_energies(ds)
                #xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                #center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                center = nar([ms.mean_xm[nf], ms.mean_ym[nf],ms.mean_zm[nf]])
                #center = nar([ms.mean_x[nf], ms.mean_y[nf],ms.mean_z[nf]])
                #msR = ms.rc
                msR = ms.rm[:,nf]
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR.max()
                Radius = max([1.0/128, MaxRadius])
                #print(Radius*128)
                #Radius=1/128
                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)

                dv = sph[YT_cell_volume]
                RR = sph['radius']
                DD = sph[YT_density]


                scrub = other_scrubber.scrubber(sph,do_velocity=False)
                scrub.compute_g_radial()
                Rradial = scrub.r
                GR = scrub.gr

                Rmin=1/2048
                Rradial[Rradial<Rmin]=Rmin

                order = np.argsort(Rradial)

                r_bins = np.geomspace( Rradial.min(),Rradial.max(),128)

                from scipy import stats
                total_g_radial, bin_edges, binnumber=stats.binned_statistic(Rradial, GR, statistic='sum', bins=r_bins)
                bin_cen = 0.5*(bin_edges[1:]+bin_edges[:-1])

                #mass = np.cumsum( DD[order]*dv[order])*4*np.pi*colors.G
                order=np.argsort(msR)
                mass = ( DD[order]*dv[order])*4*np.pi*colors.G

                #axes[1][nframe].plot( bin_cen,np.abs(total_g_radial))
                axes[1][nframe].plot( msR[order],mass[order])
                #ext_rho(np.abs(total_g_radial))
                ext_rho(mass)






            for ax in axes[1]:
                ax.set(xscale='log',yscale='log',ylim=ext_rho.minmax)
            #for ax in axes[2]:
            #    ax.set(xscale='log',yscale='log',ylim=ext_rho2.minmax)
            fig.savefig('plots_to_sort/gravity_metric_%s_c%04d'%(this_looper.sim_name,core_id))



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
        #core_list=[68]
        #core_list = TL.loops[sim].core_by_mode['Alone']
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

        for core_id in core_list:
            frame_list = [tsing_tool[sim].tend_frame[core_id]]
            multi_proj=proj_and_profile(TL.loops[sim])
            multi_proj.run(core_list=[core_id],tsing=tsing_tool[sim], frame_list=frame_list)#, r_inflection=anne.inflection[sim])
