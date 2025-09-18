from starter2 import *
import xtra_energy

import core_proj_three
reload(core_proj_three)

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


class simple_multipane_proj():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
    def run(self, core_list=None, frame_list=None):
        this_looper=self.this_looper
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        for core_id in core_list:
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

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


                plot_collector=[]
                ax=0
                field=YT_density
                weight_field=None
                fig,axes=plt.subplots(1,3)
                for ax in [0,1,2]:
                    pw = yt.ProjectionPlot(ds, ax, YT_density, data_source=sph, center=center, origin='window', weight_field=weight_field)
                    frb =  pw.data_source.to_frb( 2*sph.radius,512)
                    axes[ax].imshow( frb[YT_density], norm=mpl.colors.LogNorm())
                    axes[ax].set(xticks=[], yticks=[])
                fig.savefig('plots_to_sort/three_proj_%s_c%04d'%(this_looper.sim_name,core_id))
                plt.close('fig')
                #doesn't work well
                if 0:
                    pw = yt.ProjectionPlot(ds, ax, field, data_source=sph, center=center, origin='window', weight_field=weight_field)
                    plot_collector.append(pw)
                    from mpl_toolkits.axes_grid1 import AxesGrid
                    fig=plt.figure()
                    grid = AxesGrid(
                    fig,(0.075,0.075,0.85,0.85),
                    nrows_ncols=(1, 3),
                    axes_pad=0.1,
                    label_mode="1",
                    share_all=True)
                    axis_list=[0]
                    for axdir in axis_list:
                        #plot = p.plots[field]
                        plot = plot_collector[axdir].plots[YT_density]
                        plot.figure = fig
                        plot.axes = grid[axdir].axes
                        plot.cax = grid.cbar_axes[axdir]
                    for axdir in axis_list:
                        plot_collector[axdir]._setup_plots()
                    fig.set_size_inches(12,4)
                    fig.savefig('plots_to_sort/proj3_%s_c%04d_n%04d'%(this_looper.sim_name,core_id,frame))



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
            multi_proj=simple_multipane_proj(TL.loops[sim])
            multi_proj.run(core_list=[core_id], frame_list=frame_list)#, r_inflection=anne.inflection[sim])
