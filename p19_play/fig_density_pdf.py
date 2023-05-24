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
        ext = [extents() for n in range(Nplots+1)]

        rho_ext = extents()
        for core_id in core_list:
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            rho_ext(ms.density)


        ds = this_looper.load(0)
        ad = ds.all_data()
        min_rho = ad[YT_density].v.min()
        bins = np.geomspace(min_rho, rho_ext.minmax[1],64)

        for nf,frame in enumerate([0, this_looper.target_frame]):
            ds = this_looper.load(frame)
            ad = ds.all_data()
            pdf, bins = np.histogram( ad[YT_density].v, weights=ad[YT_cell_volume].v,bins=bins, density=True)
            bc = 0.5*(bins[1:]+bins[:-1])
            if nf==0:
                sl = slice(0,1)
            else:
                sl = slice(2,None)
            for ax in axes[sl]:
                ax.plot( bc, pdf, c='k')


        print('word')
        for core_id in core_list:
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
                rho_to_hist = DD+0
                rho_to_hist.sort()
                #cuml = np.arange(rho_to_hist.size)/rho_to_hist.size
                #axes[4][0].plot( rho_to_hist, cuml)
                axes[nframe].hist( DD.v, weights=dv.v, histtype='step', bins=bins, density=True, color=[0.5]*4)

        for ax in axes:
            ax.set(xscale='log',yscale='log',xlabel=r'$\rho$',ylabel='')
        axes[0].set(ylabel=r'$PDF(\rho)$')
        for ax in axes[1:]:
            ax.set(yticks=[])
        fig.savefig('plots_to_sort/fig_density_pdf_%s'%(this_looper.sim_name))



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
        core_list=core_list[:3]
        #core_list=core_list[:10]
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
