
from starter2 import *
from collections import defaultdict
import scipy
import colors
import xtra_energy
import camera_path
reload(camera_path)
import hair_dryer
reload(hair_dryer)

import three_loopers_u500 as TL
import movie_frames 

from scipy.ndimage import gaussian_filter

def anatomy(this_looper,core_list=None, do_plots=True, mass=None, dof=None, volume=None,
           annotate_phases=False):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)

    thtr=this_looper.tr
    mask = movie_frames.quantized_mask(this_looper).flatten()
    all_times=thtr.times
    all_frames=thtr.frames
    times=thtr.times[mask]+0 #the zero makes a copy
    times.shape=times.size,1
    times=times/colors.tff
    frames=all_frames[mask]
    rho_all = thtr.track_dict['density']
    rho_min=rho_all.min()
    rho_max=rho_all.max()

    mini_scrubbers={}
    for nc,core_id in enumerate(core_list):
        print('V %s %d'%(this_looper.sim_name,core_id))
            
        ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
        ms.particle_pos(core_id)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)
        ms.compute_ke_rel(core_id)
        mini_scrubbers[core_id]=ms

        if ms.nparticles < 1000 or True:
            sl=slice(None)
            c=[0.5]*4
        else:
            sl = slice(None,None,10)
            #c=[0,0,0,0.1]
            c=[0.1]*4

        rho = ms.density[sl].transpose()#*colors.density_units
        rho = rho[mask,:]

        dv = ms.cell_volume[sl].transpose()[mask,:]
        vv = dv.sum(axis=1)
        vx = ms.rel_vx[sl].transpose()[mask,:]
        vy = ms.rel_vy[sl].transpose()[mask,:]
        vz = ms.rel_vz[sl].transpose()[mask,:]
        v22 = ms.rel_vmag[sl].transpose()[mask,:]

        vr = ms.vr_rel[sl].transpose()[mask,:]
        vt = (ms.vt2_rel[sl].transpose()[mask,:])**0.5

        vrm=vr.mean(axis=1)
        v2 = v22.mean(axis=1)
        vtm=vt.mean(axis=1)

        #
        # Pick a bunch of times to analyze.
        # Especially the start and end of singularity
        #
        UB = gaussian_filter(rho.max(axis=1),1)
        tf = times.flatten()
        dt = tf[1:]-tf[:-1]
        dU = UB[1:]-UB[:-1]
        tc = 0.5*(times[1:]+times[:-1])
        dU = dU[1:-1]
        dt = dt[1:-1]
        tc = tc[1:-1]
        dudt=dU/dt
        thresh = 1e5
        singularity = np.where( dudt >thresh)[0][0]

        if (dudt[singularity:]<=0).any():
            collapse_done = np.where( dudt[singularity:] < 0)[0][0]
            collapse_done += singularity
        else:
            collapse_done=-1

        #ax.axvline( tf[singularity], c='k')
        #ax.axvline( tf[collapse_done], c='k')


        frame_index=[]
        for frac in [0.0, 0.2, 0.5, 0.7]:
            target = frac*frames[singularity]
            argmin = np.argmin( np.abs( frames-target))
            frame_index.append( frames[argmin])

        frame_index = sorted(frame_index)

        frame_index += [frames[singularity], frames[collapse_done]]
        frame_index = nar(frame_index).astype('int')
        

        rmap = rainbow_map(all_times.size)
        rmap = rainbow_map(len(frame_index))
        color_list = [ rmap(frame) for frame in range(len(frame_index))]
        color_list[-2]='g'
        color_list[-1]='r'
        line_list = {frames[singularity]:2, frames[collapse_done]:2}

        #
        # Set up plots
        #


        if 0:
            fig,axes=plt.subplots(3,1, figsize=(12,13))
            ax=axes[0];ax1=axes[1]; ax3=axes[2]
            ax2=ax.twinx()
        elif 0:
            from matplotlib.gridspec import GridSpec
            fig=plt.figure(constrained_layout=True, figsize=(12,12))
            nx = len(frame_index)
            gs = GridSpec(3, nx, figure=fig)
            ax = fig.add_subplot( gs[0,:])
            ax1 = fig.add_subplot( gs[1,:])
            #ax3 = fig.add_subplot( gs[2,:])
            ax3 = [ fig.add_subplot( gs[2,iii], wspace=0) for iii in range(len(frame_index))]
            ax2 = ax.twinx()
        else:
            nx = len(frame_index)
            fig = plt.figure(figsize=(8, 10))
            fig.tight_layout()
            outer_grid = fig.add_gridspec(3, 1)
            #ax = outer_grid[0,0].subgridspec(1,1).subplots()
            #ax2 = ax.twinx()
            ax2 = outer_grid[0,0].subgridspec(1,1).subplots()
            ax = ax2.twinx()
            ax1 = outer_grid[1,0].subgridspec(1,1).subplots()
            ax3 = outer_grid[2,0].subgridspec(1,nx,wspace=0).subplots()

        #density plot
        UN = colors.density_units
        ax.plot(times , rho*UN, c=c, linewidth=0.1)
        ylim=[10,5e9]
        axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel=r'$\rho \rm{[cm^{-3}]}$',yscale='log',
               ylim=ylim)
        ax2.set(xlabel=r'$t/t_{ff}$')

        if annotate_phases:
            ax2.text( 0.150, 9, r"$collection$")
            ax2.text( 0.350, 9, r"$hardening$")
            ax2.text( 0.55, 9, r"$singularity$")
            #ax2.text( 0.9, 9, r"$(mosh)$")
            ax2.text( times[singularity]+0.01,9,r'$t_{\rm{sing}}$')
            ax2.text( times[collapse_done]+0.01,9,r'$t_{\rm{sung}}$')


        #velocity plots
        if 1:
            ok = vrm>0
            ax2.plot(times[ok], vrm[ok], 'r--')
            ax2.plot(times[~ok], np.abs(vrm[~ok]), c='r',label=r'$v_r$')
        ax2.plot(times, vtm, c='c', label=r'$v_t$')
        ax2.plot(times, v2, c='k', label=r'$v$')
        ax2.plot( times, times*0+1, c=[0.5]*4)
        ax2.set( ylim=[0,10], ylabel=r'velocity (preimage)')

        if 1:
            #density CDF
            scale = [colors.density_units*ms.density.min(),colors.density_units*ms.density.max()]
            bins=np.geomspace( scale[0],scale[1],64)
            for nnn,frame in enumerate(frame_index):
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                ax.axvline( all_times[nf]/colors.tff, c=color_list[nnn], linewidth=line_list.get(frame,1))
                rho_to_hist = ms.density[:,nf].flatten()*colors.density_units
                cuml = np.arange(rho_to_hist.size)/rho_to_hist.size
                #ax1.hist( rho_to_hist, histtype='step',color=rmap(n),bins=bins, cumulative=True, density=True)
                ax1.plot( sorted(rho_to_hist), cuml,color=color_list[nnn], linewidth=line_list.get(frame,1))
            axbonk(ax1,xlabel=r'$\rm{Density}\ [\rm{cm^{-3}}]\rm{(preimage)}$',
                   yscale='linear',xscale='log', xlim=scale, ylabel='Density CDF (preimage)')

        #
        # Binding Energy
        #
        camera = camera_path.camera_1(this_looper, 'sphere')
        camera.run([core_id], frames, mini_scrubbers)
        y_ext = extents(nar([.5,2e5]))
        r_ext = extents(colors.length_units_au*nar([1./2048, 0.3]))
        if 1:
            nnn=-1
            for frame in frame_index:
                print('frame',frame)
                nnn+=1

                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                #rsph = ds.arr(8.0/128,'code_length')
                if 0:
                    rsph = max([camera.max_radius[nf], 1/128])
                    center = nar([ms.mean_x[nf], ms.mean_y[nf],ms.mean_z[nf]])
                if 1:
                    ms=mini_scrubbers[core_id]
                    msR = ms.rc
                    msR[ msR<1/2048]=1/2048
                    
                    MaxRadius=msR[:,nf].max()
                    Radius = max([8.0/128, MaxRadius])
                    rsph = ds.arr(Radius,'code_length')
                    center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                sp = ds.sphere(center,rsph)


                #Get data arrays

                dv = np.abs(sp[YT_cell_volume])
                RR =sp[YT_radius]
                DD = sp[YT_density]
                EG = np.abs(sp[YT_grav_energy_2])
                #don't use this, make the relative KE
                #EK = np.abs(sp[YT_kinetic_energy])
                #relative kinetic energy.
                vx = sp[YT_velocity_x].v - ms.mean_vx[nf]
                vy = sp[YT_velocity_y].v - ms.mean_vy[nf]
                vz = sp[YT_velocity_z].v - ms.mean_vz[nf]
                EK = 0.5*DD*(vx*vx+vy*vy+vz*vz)
                print('1')

                if 0:
                    #gas density cdf
                    rho_to_hist=DD.flatten()*colors.density_units
                    rho_to_hist.sort()
                    cdf = np.arange( rho_to_hist.size)/rho_to_hist.size
                    ax1.plot( rho_to_hist, cdf, color=color_list[nnn],
                             linewidth=line_list.get(frame,1), linestyle='--')
                    ax1.set(xlabel='rho',yscale='linear',xscale='log',
                            ylabel='cuml')


                if 1:
                    #RADIAL PLOTS
                    ORDER = np.argsort( RR)
                    V_cuml =  np.cumsum( dv[ORDER])
                    V_local =  dv[ORDER]
                    RR_cuml = RR[ORDER]
                    EG_cuml = np.cumsum( EG[ORDER]*V_local)/V_cuml
                    EK_cuml = np.cumsum( EK[ORDER]*V_local)/V_cuml

                    line=line_list.get(frame,1)
                    label_g=None
                    label_k=None
                    if nnn==0:
                        label_g=r'$\overline{E}_G$'
                        label_k=r'$\overline{E}_K$'
                    ax3[nnn].plot(  RR_cuml*colors.length_units_au, EG_cuml, c=color_list[nnn], linestyle='-', linewidth=line,
                                  label=label_g)
                    ax3[nnn].plot( RR_cuml*colors.length_units_au, EK_cuml,  c=color_list[nnn], linestyle='--', linewidth=line,
                                 label=label_k)
                    if nnn==0:
                        ax3[nnn].legend(loc=2)
                    y_ext(EG_cuml)
                    y_ext(EK_cuml)
                    r_ext(RR_cuml*colors.length_units_au)
                    ax3[nnn].set( xscale='log',yscale='log',xlabel=r'$r\ [\rm{AU}]$', ylabel='')
                    print('2')

                if 0:
                    #PHASE PLOTS
                    xxbins=np.geomspace(5e-3,1e7,128)
                    yybins=np.geomspace(5e-3,1e7,128)
                    #xxbins = np.geomspace(ke.min(),ke.max(),128)
                    #yybins = np.geomspace(ge[ge>0].min(),ge.max(),128)
                    hist, xbins,ybins=np.histogram2d(EK.flatten(),EG.flatten(),bins=[xxbins,yybins])

                    pch.helper(hist,xbins,ybins,ax=ax3[nnn])
                    axbonk(ax3[nnn],xscale='log',yscale='log',xlabel='KE',ylabel='GE')
                    ax3[nnn].plot( xxbins,xxbins,c='k')
                    ax3[nnn].scatter(ms.ke_rel[:,nf],np.abs(ms.ge[:,nf]), edgecolor='r',s=30, facecolor='None')
                    y_ext = extents(xxbins)
                    r_ext = extents(yybins)
                #ax3.hist( EG_cuml)
            for na,aaa in enumerate(ax3):
                print('wut')
                print(y_ext.minmax)
                if na>0:
                    aaa.set(yticks=[],ylabel='')
                else:
                    aaa.set(ylabel='Energy (full sphere)')
                aaa.set( xlim=r_ext.minmax, ylim=y_ext.minmax)







        outname='plots_to_sort/%s_anatomy_c%04d.png'%(this_looper.sim_name,core_id)
        fig.savefig(outname)
        print(outname)
        plt.close(fig)



sims=['u501', 'u502','u503']
import three_loopers_u500 as TL
import mass_tools
if 0:
    if 'mt' not in dir():
        mt={}
    for sim in sims:
        if sim not in mt:
            mt[sim]=mass_tools.mass_tool(TL.loops[sim])
            mt[sim].run()

sims=[ 'u502','u503']
sims=['u502']#, 'u501']
for sim in sims:
    #core_list=[381]
    #core_list={'u501':[323], 'u502':[381]}[sim]
    #core_list={'u501':[323], 'u502':[112]}[sim]
    #core_list=[31,32]
    #core_list=[114]

    core_list=None
    #annotate_phases=False
    annotate_phases=True
    core_list = [114]
    #core_list=[74]
    #core_list = TL.loops[sim].core_by_mode['Alone']
    #core_list=core_list[:5]
    frrt=anatomy(TL.loops[sim], do_plots=True, core_list=core_list, annotate_phases=annotate_phases)#, mass=mt[sim].unique_mass, dof=mt[sim].dof, volume=mt[sim].volume)

