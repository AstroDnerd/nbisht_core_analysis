from starter2 import *
import xtra_energy
reload(xtra_energy)
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
from collections import defaultdict
from scipy.ndimage import gaussian_filter
import other_scrubber
reload(other_scrubber)
G = colors.G
class radials():
    def __init__(self):
        self.alpha_time=defaultdict(list)
        self.cores_used=[]
    def run(self,this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        if frame_list is None and False:
            mask = movie_frames.quantized_mask(this_looper).flatten()
            ok = np.zeros_like(mask)
            ok[::10] = mask[::10]
            mask=ok ;print('kludge mask')
            times=thtr.times[mask]+0 #the zero makes a copy
            self.times=times+0
            times.shape=times.size,1
            times=times/colors.tff
            frame_list=thtr.frames[mask]
            rm = rainbow_map(len(frame_list))

        y_ext=extents()
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index

        nframes = 4
        nrows = 4
        fig6,ax6=plt.subplots(nrows,nframes)
        ext=[extents() for n in range(nrows)]
        extr=extents()
        for core_id in core_list:
            if 1:
                frame_mask = np.zeros_like(thtr.times, dtype='bool')
                frame_mask[get_time_index(tsing.tend_core[core_id])]=True
                frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
                frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
                frame_mask[0]=True
                times = thtr.times[frame_mask]+0
                times.shape=times.size,1
                frame_list=thtr.frames[frame_mask]
                rm = rainbow_map(len(frame_list))
                #if len(frame_list) != nframes:
                #    raise
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)
            
            r_bins = np.geomspace( 2e-4, 32/128, 32)
            r_cen = 0.5*(r_bins[1:]+r_bins[:-1])
            RV = np.zeros( [len(r_bins)-1, len(frame_list)])
            frame_rmap=rainbow_map(len(frame_list))
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                time = times[nframe,0]

                c = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([8.0/128, MaxRadius])
                rsph = ds.arr(Radius,'code_length')
                sp = ds.sphere(c,rsph)



                #Get data arrays

                dv = sp[YT_cell_volume]
                RR = sp['radius']
                DD = sp[YT_density]
                #EG = sp[YT_grav_energy_2]
                #GDGR=sp[YT_gdotgradrho]

                ORDER = np.argsort( RR)
                RR_sort = RR[ORDER]
                rho_sort = DD[ORDER]
                dv_sort = dv[ORDER]


                M_cuml = np.cumsum( rho_sort*dv_sort)
                d2_cuml = np.cumsum( rho_sort**2*dv_sort)
                V_cuml = np.cumsum( dv[ORDER])
                rho_cuml = M_cuml/V_cuml
                #gdgr_cuml = np.cumsum( GDGR[ORDER]*dv[ORDER])

                if 1:
                    vel = []
                    for axis in 'xyz':
                        #vel.append( sp['velocity_%s'%axis][ORDER][:10].mean())
                        #vel.append( sp['velocity_%s'%axis][ORDER].mean())
                        #vel.append( (rho_sort*sp['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/M_cuml)
                        vel.append( (rho_sort*sp['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/(rho_sort*dv_sort)[:30].sum())
                        #vel.append(0)
                    scrub = other_scrubber.scrubber(sp, reference_velocity = vel)
                    #scrub0 = other_scrubber.scrubber(sp, reference_velocity = [0.0]*3)
                    #scrub.compute_ke_rel()
                linewidth=0.3
                c=frame_rmap(nframe)

                extr(RR_sort)
                if 1:
                    from mpl_toolkits.axes_grid1 import AxesGrid
                    fig=plt.figure()
                    grid = AxesGrid(
                        fig,(0.075,0.075,0.85,0.85),
                        nrows_ncols=(1, 3),
                        axes_pad=0.1,
                        label_mode="1",
                        share_all=True)
                    pw = yt.ProjectionPlot(ds,0,YT_density,data_source=sp, origin='window')
                    plot = pw[YT_density]
                    plot.figure=fig
                    plot.axes=grid[0].axes
                    plot.cax =grid.cbar_axes[0]
                    pw._setup_plots()
                    pw.save('plots_to_sort/temp')
                    fig.savefig('plots_to_sort/ttt')
                    pdb.set_trace()

                args = {'linewidth':linewidth, 'c':[0.5]*4}
                #ax6[0][nframe].plot( RR_sort, rho_cuml, **args)
                ext[0](rho_cuml)

                vr = scrub.vr_rel
                #vr_cumsum = np.cumsum( vr[ORDER]*dv_sort)/V_cuml
                #ax6[1][nframe].plot(RR_sort, vr_cumsum, **args)
                #ext[1](vr_cumsum)

                vt = scrub.vt_rel
                #vt_cumsum = np.cumsum( vt[ORDER]*dv_sort)/V_cuml
                #ax6[2][nframe].plot(RR_sort, vt_cumsum, **args)
                #ext[2](vt_cumsum)

                #velocity reference games
                #vr = scrub0.vr_rel
                #vr_cumsum = np.cumsum( vr[ORDER]*dv_sort)/V_cuml
                #ax6[0][nframe].plot(RR_sort, vr_cumsum, linewidth=linewidth, c=[0.5]*4)
                #ext[0](vr_cumsum)

                #EG=M_cuml**2/RR_sort
                EG=M_cuml**2#/RR_sort
                ax6[3][nframe].plot(RR_sort, EG, **args)
                ext[3](EG)
                EG=rho_sort**2
                ax6[2][nframe].plot(RR_sort, EG, **args)
                ext[2](EG)

                
        if 0:
            for na,aaa in enumerate(ax6[0]):
                aaa.set(xlabel='',xscale='log', yscale='log', ylim=ext[0].minmax, xlim=extr.minmax)
                if na>0:
                    aaa.set(yticks=[])
            for na,aaa in enumerate(ax6[1]):
                aaa.set(xlabel='',xscale='log', ylim=ext[1].minmax, xlim=extr.minmax)
                if na>0:
                    aaa.set(yticks=[])
        if 0:
            for na,aaa in enumerate(ax6[2]):
                aaa.set(xlabel='',xscale='log', ylim=ext[2].minmax, xlim=extr.minmax)
                if na>0:
                    aaa.set(yticks=[])
        if 1:
            for na,aaa in enumerate(ax6[2]):
                aaa.set(xlabel='',xscale='log', ylim=ext[2].minmax, xlim=extr.minmax, yscale='log')
                if na>0:
                    aaa.set(yticks=[])
        for na,aaa in enumerate(ax6[3]):
            aaa.set(xlabel='',xscale='log', ylim=ext[3].minmax, xlim=extr.minmax, yscale='log')
            if na>0:
                aaa.set(yticks=[])

        fig6.savefig('plots_to_sort/radials_%s.png'%(this_looper.sim_name))

            #ax5.plot(times.flatten(), alpha_time[core_id])
        #fig5.savefig('plots_to_sort/several_alpha_%s'%(this_looper.sim_name))
        #ax5[0][0].set(ylabel='Eg', ylim=y_ext.minmax)
        #ax5[1][0].set(ylabel='M^2/R', ylim=y_ext.minmax)
        #ax5[2][0].set(ylabel='EG/M^2/R', ylim=[1e-3,1e3])
        #for aaa in ax5.flatten():
        #    aaa.set(xlabel='R',xscale='log',yscale='log')
        #    aaa.axhline(0.5,c=[0.5]*4)
        #    aaa.axhline(2.0,c=[0.5]*4)
        #ax5[0][1].set(ylim=ax5[0][0].get_ylim())
        #ax5[1][1].set(ylim=ax5[1][0].get_ylim())
        #ax5[2][1].set(ylim=ax5[2][0].get_ylim())
        #fig5.savefig('plots_to_sort/eg_m2_tsing_%s'%(this_looper.sim_name))
#


if 0:
    import three_loopers_six as TL
    sim_list=['u601','u602','u603']
    sim_list=['u602']
import three_loopers_u500 as TL
sim_list=['u501','u502','u503']
sim_list=['u502']

import tsing
reload(tsing)
if 'tsing_tool' not in dir() or True:
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
if 'RV' not in dir() or True:
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        core_list = TL.loops[sim].core_by_mode['Alone']
        core_list=core_list[10:15]
        core_list=[114]


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

        RV=radials()
        RV.run(TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
