from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
def mass_inside(this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    if frame_list is None:
        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        ok = np.zeros_like(mask)
        ok[::10] = mask[::10]
        #mask=ok ;print('kludge mask')
        times=thtr.times[mask]+0 #the zero makes a copy
        times.shape=times.size,1
        times=times/colors.tff
        frame_list=thtr.frames[mask]

    rm = rainbow_map(len(frame_list))
    for core_id in core_list:
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        
        r_bins = np.geomspace( 2e-4, 32/128, 32)
        RV = np.zeros( [len(r_bins)-1, len(frame_list)])
        frame_rmap=rainbow_map(len(frame_list))
        fig6,ax6=plt.subplots(1,1)
        got_tsing=False
        got_tend=False
        for nframe,frame in enumerate(frame_list):
            print("MassIn on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
            ds = this_looper.load(frame)
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

            dv = np.abs(sp[YT_cell_volume])
            RR =sp[YT_radius]
            DD = sp[YT_density]

            ORDER = np.argsort( RR)
            RR_cuml = RR[ORDER]
            M_local = DD[ORDER]*dv[ORDER]
            M_cuml = np.cumsum( DD[ORDER]*dv[ORDER])
            V_cuml = np.cumsum( dv[ORDER])
            V_local = dv[ORDER]

            digitized = np.digitize( RR_cuml, r_bins)
            Mass  =nar([ M_cuml[ digitized == i].max() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])

            ok = Mass>0
            RV[ok,nframe]=Mass[ok]
            r_bins_c = 0.5*(r_bins[1:]+r_bins[:-1])
            #ax6.plot(r_bins_c, Mass, c = frame_rmap(nframe))
            linewidth=0.3
            c=frame_rmap(nframe)
            if tsing:
                if not got_tsing:
                    if time > tsing.tsing_core[core_id]:
                        c = 'g'
                        linewidth=1
                        got_tsing=True
                if not got_tend:
                    if time > tsing.tend_core[core_id]:
                        c = 'r'
                        linewidth=1
                        got_tend=True
                        pfit = np.polyfit( np.log10(RR_cuml[1:]), np.log10(M_cuml[1:]/RR_cuml[1:]**3),1)
                        print(pfit)

            from scipy.ndimage import gaussian_filter

            #ax6.plot(RR_cuml, gaussian_filter(M_cuml,1), c=c,linewidth=linewidth)
            ax6.plot(RR_cuml[1:], M_cuml[1:]/RR_cuml[1:]**3, c=c,linewidth=linewidth)
        ax6.set(xlabel='R',xscale='log',ylabel='M(<r)',yscale='log')
        ax6.set_title('Slope at tend: %0.2f'%pfit[0])
        #ax6.plot( r_bins[1:], 4*np.pi/3*r_bins[1:]**3,'k')
        fig6.savefig('plots_to_sort/mass_interior_vs_radius_%s_c%04d.png'%(this_looper.sim_name, core_id))


        fig,axes=plt.subplots(1,1)
        #fig.subplots_adjust(wspace=0)
        rcen = 0.5*(r_bins[1:]+r_bins[:-1])
        XXX,YYY = np.meshgrid( times.flatten(),rcen)
        ax=axes
        norm = mpl.colors.LogNorm( RV[RV>0].min(), RV.mean())
        cmap=copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')

        plot=ax.pcolormesh( XXX,YYY, RV,  norm=norm, shading='nearest', cmap=cmap)
        fig.colorbar(plot, label='M(<r)')
        ax.set(yscale='log', xlabel='t/tff', ylabel='R [code units]')
        
        if 0:
            rrrr = msR.transpose()
            rrrr = rrrr[mask,:]

            ax.plot(times , rrrr, c=[0.5]*3, linewidth=0.1, alpha=0.5)
        if tsing:
            ax.axvline(  tsing.tsing_core[core_id],c='k')
            #ax1.axvline( tsing.tsing_core[core_id],c='k')
            ax.axvline(  tsing.tend_core[core_id],c='k')
            #ax1.axvline( tsing.tend_core[core_id],c='k')


        fig.savefig('plots_to_sort/mass_color_%s_c%04d.png'%(this_looper.sim_name, core_id))



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

import anne
reload(anne)
#anne.make_inflection()
if 'RV' not in dir() or True:
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[74]
        #core_list = TL.loops[sim].core_by_mode['Alone']
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

        RV=mass_inside( TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
if 0:
    fig,ax=plt.subplots(1,1)
    rcen = 0.5*(RV['rbins'][1:]+RV['rbins'][:-1])
    XXX,YYY = np.meshgrid( RV['times'],rcen)

#norm = mpl.colors.LogNorm( vmin=RV['RV'][RV['RV']>0].min(), vmax=RV['RV'].max())
    norm = mpl.colors.LogNorm( vmin=1e-3, vmax=1e3)
    plot=ax.pcolormesh( XXX,YYY, RV['RV'],  norm=norm, shading='nearest', cmap='seismic')
    fig.colorbar(plot)
    ax.set(yscale='log')
    fig.savefig('plots_to_sort/temp.png')

