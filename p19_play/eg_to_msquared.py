from starter2 import *
import xtra_energy

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
        Nplots = 1
        Ntimes = 4
        fig,axes = plt.subplots(Nplots,Ntimes, figsize=(12,6))
        fig.subplots_adjust(hspace=0,wspace=0)
        ext=extents()
        ext_r=extents()

        for core_id in core_list:
            fig2,ax2=plt.subplots(1,1)
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            frame_mask[0]=True
            frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(tsing.tend_core[core_id])]=True
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))
            axes[0].set_title(r'$t=0$')
            axes[1].set_title(r'$t=0.9 t_{\rm{sing}}$')
            axes[2].set_title(r'$t=t_{\rm{sing}}$')
            axes[3].set_title(r'$t=t_{\rm{end}}$')

            img_collector=[]
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
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
                V_cuml = np.cumsum( dv[ORDER])
                rho_sort = DD[ORDER]
                dv_sort = dv[ORDER]
                RR_sort = RR[ORDER]
                mass_cuml=np.cumsum(rho_sort*dv_sort)
                grav = sph[YT_grav_energy_2]
                grav_cuml = np.cumsum(grav[ORDER]*dv_sort)
                if 0:
                    vel = []
                    for axis in 'xyz':
                        #vel.append( sp['velocity_%s'%axis][ORDER][:10].mean())
                        #vel.append( sp['velocity_%s'%axis][ORDER].mean())
                        #vel.append( (rho_sort*sp['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/M_cuml)
                        vel.append( (rho_sort*sph['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/(rho_sort*dv_sort)[:30].sum())
                        #vel.append(0)
                    scrub = other_scrubber.scrubber(sph, reference_velocity = vel)


                    vt = scrub.vt_rel
                    vt_cumsum = np.cumsum( rho_sort*vt[ORDER]*dv_sort)/mass_cuml
                    #vt_cumsum = np.cumsum( vt[ORDER]*dv_sort)/V_cuml
                    #pfit = np.polyfit( np.log10(RR_sort),np.log10(vt_cumsum),1)
                    #print(pfit)
                    #axes[nframe].plot( RR_sort, 10**(pfit[0]*np.log10(RR_sort)+pfit[1]))
                    #vt_cumsum = np.cumsum( rho_sort*dv_sort)#/V_cuml
                    args={'color':[0.5]*4}
                    axes[0].plot(RR_sort, vt_cumsum, **args)
                    #axes[nframe].plot(RR_sort, RR_sort**3,c='r')
                    ext(vt_cumsum)
                    ext_r(RR_sort)
                    #rho_to_hist = DD+0
                    #rho_to_hist.sort()
                    #cuml = np.arange(rho_to_hist.size)/rho_to_hist.size
                    #axes[4][0].plot( rho_to_hist, cuml)
                    #axes[nframe].hist( DD.v, weights=dv.v, histtype='step', bins=bins, density=True, color=[0.5]*4)

                    
                    axes[1].plot(RR_sort, mass_cuml/V_cuml)
                    axes[2].plot(RR_sort, MG)
                MG=mass_cuml**2/RR_sort*colors.G
                GG=np.abs(grav_cuml)
                #axes[3].plot(RR_sort, GG)
                ax2.scatter( np.abs(MG),np.abs(GG))
                ax2.plot( np.abs(MG), np.abs(MG))
            ax2.set(xscale='log',yscale='log')
            fig2.savefig('plots_to_sort/EG_M2_%s_c%04d'%(this_looper.sim_name,core_id))
            plt.close(fig2)




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
        core_list=core_list[3:]
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
