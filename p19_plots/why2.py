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
def why(this_looper,core_list=None, do_plots=True, r_inflection=None, frame_list=None, tsing=None):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    if frame_list is None:
        thtr=this_looper.tr
        mask = movie_frames.quantized_mask(this_looper).flatten()
        ok = np.zeros_like(mask)
        ok[::10] = mask[::10]
        #mask=ok
        times=thtr.times[mask]+0 #the zero makes a copy
        times.shape=times.size,1
        times=times/colors.tff
        frame_list=thtr.frames[mask]

    rm = rainbow_map(len(frame_list))
    for core_id in core_list:
        ms = trackage.mini_scrubber(this_looper.tr,core_id)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)
        ms.compute_ke_rel(core_id)
        
        r_bins = np.geomspace( 2e-4, 32/128, 32)
        RV = np.zeros( [len(r_bins)-1, len(frame_list)])
        for nframe,frame in enumerate(frame_list):
            print("%s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
            ds = this_looper.load(frame)
            nf = np.where( this_looper.tr.frames == frame)[0][0]
            xtra_energy.add_energies(ds)

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
            EG = np.abs(sp[YT_grav_energy_2])
            #EK = np.abs(sp[YT_kinetic_energy])


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
            #hist, bins = np.histogram( RR_cuml, bins=r_bins, weights=EG[ORDER])
            #bc = 0.5*(bins[:-1]+bins[1:])
            #ax3.plot( bc, hist, c='c')
            digitized = np.digitize( RR_cuml, r_bins)
            EGV = EG*dv
            EKV = EK*dv
            Rmeans  =nar([ RR_cuml[ digitized == i].mean() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])
            DVmeans =nar([ V_cuml[ digitized == i].mean()  if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])
            EGmeans =nar([ EG_cuml[ digitized == i].mean() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])
            EKmeans =nar([ EK_cuml[ digitized == i].mean() if (digitized==i).any() else -4000 for i in range(1,len(r_bins))])

            
            ok = EKmeans>0
            RV[ok,nframe]=(EGmeans/EKmeans)[ok]


        fig,axes=plt.subplots(1,2)
        fig.subplots_adjust(wspace=0)
        rcen = 0.5*(r_bins[1:]+r_bins[:-1])
        XXX,YYY = np.meshgrid( times.flatten(),rcen)
        ax=axes[0]
        ax1=axes[1]

        #norm = mpl.colors.LogNorm( vmin=RV['RV'][RV['RV']>0].min(), vmax=RV['RV'].max())
        norm = mpl.colors.LogNorm( vmin=1e-3, vmax=1e3)
        plot=ax.pcolormesh( XXX,YYY, RV,  norm=norm, shading='nearest', cmap='seismic')
        plot=ax1.pcolormesh( XXX,YYY, RV,  norm=norm, shading='nearest', cmap='seismic')
        fig.colorbar(plot)
        ax.set(yscale='log')
        ax1.set(yscale='log')
        ax1.set(yticks=[])

        
        rrrr = msR.transpose()
        rrrr = rrrr[mask,:]

        ax.plot(times , rrrr, c=[0.5]*3, linewidth=0.1, alpha=0.5)
        if tsing:
            ax.axvline(  tsing.tsing_core[core_id],c='k')
            #ax1.axvline( tsing.tsing_core[core_id],c='k')
            ax.axvline(  tsing.tend_core[core_id],c='k')
            #ax1.axvline( tsing.tend_core[core_id],c='k')


        ax1.set_ylim( ax.get_ylim())
        fig.savefig('plots_to_sort/why2_%s_c%04d.png'%(this_looper.sim_name, core_id))



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
        core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=[361]
        #core_list=[8]
        #core_list=[381]
        #core_list=[323]
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        #frame_list = [100]

        RV=why( TL.loops[sim],core_list=core_list, do_plots=True, frame_list=None, tsing=tsing_tool[sim])#, r_inflection=anne.inflection[sim])
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

