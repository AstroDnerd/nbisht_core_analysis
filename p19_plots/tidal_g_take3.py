from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
import pcolormesh_helper as pch
import other_scrubber
import gravity
reload(gravity)
reload(xtra_energy)


def energy_plots(this_looper, core_list=None, r_inflection=None, r_mass=None):

    sim_name = this_looper.sim_name
    thtr=this_looper.tr
    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    frame=this_looper.target_frame
    frame=0
    ds = this_looper.load(frame)
    xtra_energy.add_energies(ds)
    xtra_energy.add_gravity(ds)
    reload(mountain_top)
    #radius=1e-2
    for core_id in core_list:
        proj_axis=0
        if proj_axis == 1:
            print("ERROR: plotting messed up for axis 1")
            raise
        def from_cg(field):
            ccc = cg[field]#.swapaxes(0,2)
            return ccc.v
        def proj(arr):
            pj = arr.sum(axis=proj_axis)
            avg_pj = pj/arr.shape[proj_axis]
            return avg_pj
        def proj_m(arr,weight):
            pj = (arr*weight).sum(axis=proj_axis)
            wj = (weight).sum(axis=proj_axis)
            return pj/wj
        print("Proj energies %s c%04d"%(sim_name, core_id))
        ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True) 

        if 0:
            r_infl = r_inflection[core_id]
            radius=max([1/128, 1.2*r_infl])
        if 0:
            radius = max([2/128, ms.r[:,-1].max()])
        if 0:
            radius = 2/128
        if 1:
            radius = 1/128

        if 0:
            radius = ds.arr(radius,'code_length')

            peak = this_looper.targets[core_id].peak_location

            left = peak - radius.v
            right = peak + radius.v
            dx = 1/2048
            nzones = np.floor((right-left)/dx).astype('int')
            level=4
        else:
            level=0
            dx = 1/2048
            left = [0.0]*3
            right =[1.0]*3
            center=[0.5]*3
            nzones=[128]*3

        cg = ds.covering_grid(level,left,nzones, num_ghost_zones=1)

        h_axis = ds.coordinates.x_axis[proj_axis]
        v_axis = ds.coordinates.y_axis[proj_axis]
        x_h = proj( from_cg( ('gas','xyz'[v_axis]) )[1:-1,1:-1,1:-1])#/nzones[h_axis]
        x_v = proj( from_cg( ('gas','xyz'[h_axis]) )[1:-1,1:-1,1:-1])#/nzones[v_axis]
        TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))
        X_to_plot = TheY.transpose()
        Y_to_plot = TheX.transpose()

        x_h_f = proj( from_cg( ('gas','xyz'[v_axis]) ))#/nzones[h_axis]
        x_v_f = proj( from_cg( ('gas','xyz'[h_axis]) ))#/nzones[v_axis]
        TheXF, TheYF = np.meshgrid( np.unique( x_h_f), np.unique( x_v_f))
        X_to_plot_f = TheYF.transpose()
        Y_to_plot_f = TheXF.transpose()

        rho = from_cg(YT_density)
        gx = from_cg(YT_acceleration_x)
        gy = from_cg(YT_acceleration_y)
        gz = from_cg(YT_acceleration_z)

        #scrub=other_scrubber.scrubber(cg,do_velocity=False)
        #gr = np.abs(scrub.rx_hat*gx+scrub.ry_hat*gy+scrub.rz_hat*gz)
        #gt = np.sqrt(gx**2+gy**2+gz**2)

        dx=1/128
        dgxdx= ( gx[2:,:,:]-gx[:-2,:,:])/(2*dx)
        dgydy= ( gy[:,2:,:]-gy[:,:-2,:])/(2*dx)
        dgzdz= ( gz[:,:,2:]-gz[:,:,:-2])/(2*dx)

        f1 = (rho[1:-1,1:-1,1:-1])
        f2 = -(dgxdx[:,1:-1,1:-1]+dgydy[1:-1,:,1:-1]+dgzdz[1:-1,1:-1,:])/(4*np.pi*colors.G)
        #f2 = -(dgxdx[:,1:-1,1:-1]+dgydy[1:-1,:,1:-1]+dgzdz[1:-1,1:-1,:])/(colors.G)
        minmin=f2[f2>0].min()
        f2[f2<minmin]=minmin


        if 0:
            fig,axes=plt.subplots(1,2)


            ext=extents()
            PPP = proj(f1)
            ext(PPP)

            vmin,vmax = ext.minmax
            cmap = 'viridis'
            norm = mpl.colors.LogNorm( vmin=vmin, vmax=vmax)
            PPP_to_plot = PPP.transpose()

            ploot=axes[0].pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
            fig.colorbar(ploot, ax=axes[0])

            PPP = proj(f2)
            PPP_to_plot = PPP.transpose()
            ploot=axes[1].pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)

            outname='plots_to_sort/Gradial_%s_%s_c%04d_n%04d.png'%('xyz'[proj_axis],sim_name, core_id, frame)
            fig.savefig(outname)
            plt.close(fig)
            print(outname)

        if 0:
            fig,ax=plt.subplots(1,1)
            pch.simple_phase(f1.flatten(),f2.flatten(),ax=ax,log=True)
            ax.set(xscale='log',yscale='log',xlabel='rho',ylabel='divg')
            ext=extents()
            ext(f1); ext(f2)
            ax.plot( ext.minmax,ext.minmax)
            ax.plot( ext.minmax,nar(ext.minmax)/np.sqrt(colors.G))
            fig.savefig('plots_to_sort/density_divg_%s_c%04d_n%04d'%(sim_name,core_id,frame))
            plt.close(fig)

        if 1:
            #check spectral potential
            #This should be perfect.
            rho = cg[YT_density]
            phi = cg[YT_potential_field]
            import gravity
            ggg = gravity.gravity_real(rho, ds.parameters['GravitationalConstant'])
            ggg.solve()
            ext=extents()
            ext(phi); ext(ggg.phi)
            fig,ax=plt.subplots(1,1)
            pch.simple_phase( phi.v.flatten(), ggg.phi.flatten(), ax=ax)
            print("Phi error", np.abs(phi.v-ggg.phi).sum()/np.abs(ggg.phi).sum())
            ax.plot( ext.minmax,ext.minmax)
            fig.savefig('plots_to_sort/potential_test_%s_c%04d_n%04d'%(sim_name,core_id,frame))
            plt.close(fig)
            return -1

        if 1:
            ggg.get_g()
            #check spectral gravity vs density.
            Np = -(ggg.dxgx+ggg.dygy+ggg.dzgz)/(4*np.pi*colors.G)
            ext=extents()
            P1 = proj(rho)
            P2 = proj(Np).real
            ext(P1)
            #ext(P2)
            vmin,vmax = ext.minmax
            cmap = 'viridis'
            norm = mpl.colors.LogNorm( vmin=vmin, vmax=vmax)
            fig,axes=plt.subplots(1,2)
            ploot=axes[0].pcolormesh( X_to_plot_f,Y_to_plot_f, P1, norm=norm, shading='nearest', cmap=cmap)
            ploot=axes[1].pcolormesh( X_to_plot_f,Y_to_plot_f, P2, norm=norm, shading='nearest', cmap=cmap)
            axes[0].set(title='rho')
            axes[1].set(title='div g')
            fig.savefig('plots_to_sort/density_divg')

        if 1:
            #ggg.get_g()
            #check spectral gravity vs density.
            Np = -(ggg.dxgx+ggg.dygy+ggg.dzgz)/(4*np.pi*colors.G)
            #Np = (ggg.dxgx+ggg.dygy+ggg.dzgz)/-ggg.G
            #Np = ggg.also_also_rho/-ggg.G
            #Np = ggg.also_also_also_rho
            #Np = ggg.also_rho
            #Np = ggg.derp
            print(Np.min())
            minmin=Np[Np>0].min()
            Np[Np<minmin]=minmin
            ext=extents()
            dgxdx_code=dgxdx[:,1:-1,1:-1]
            dgxdx_fft =ggg.dxgx[1:-1,1:-1,1:-1]
            dgydy_code=dgydy[1:-1,:,1:-1]
            dgydy_fft =ggg.dygy[1:-1,1:-1,1:-1]
            dgzdz_code=dgzdz[1:-1,1:-1,:]
            dgzdz_fft =ggg.dxgx[1:-1,1:-1,1:-1]

            grad_code = [dgxdx_code,dgydy_code,dgzdz_code]
            grad_fft = [dgxdx_fft,dgydy_fft,dgzdz_fft]

            n=-1
            fig,axes=plt.subplots(2,3,figsize=(12,8))
            fig2,ax2=plt.subplots(1,3,figsize=(12,8))
            for X1, X2 in zip(grad_code,grad_fft):
                n+=1

            #P1 =np.abs(proj(dgxdx_code)) #projections of accelerations is alwasy zero dummy.
            #P2 =np.abs(proj(dgxdx_fft))
                P1 = X1[:,:,32]
                P2 = X2[:,:,32]
                ext=extents()

                ext(P2)
                ext(P2)
                vmin,vmax = ext.minmax
                cmap = 'viridis'
                norm = mpl.colors.Normalize( vmin=vmin, vmax=vmax)
                ploot=axes[0][n].pcolormesh( X_to_plot,Y_to_plot, P1, norm=norm, shading='nearest', cmap=cmap)
                ploot=axes[1][n].pcolormesh( X_to_plot,Y_to_plot, P2, norm=norm, shading='nearest', cmap=cmap)

                logx,logy=np.log10(np.abs(P1.flatten())),np.log10(np.abs(P2.flatten()))
                ext=extents();ext(logx);ext(logy)
                pch.simple_phase(logx,logy,ax=ax2[n], nBins=128)
                ax2[n].plot(ext.minmax,ext.minmax)
                ax2[n].set_aspect('equal')
            fig2.savefig('plots_to_sort/g_code_fft_phase')
            fig.savefig('plots_to_sort/g_code_fft_image')
            

        if 0:
            #Check accelerations.  Not perfect.
            fig,ax=plt.subplots(1,1)

            ggg.get_g()
            def do(arr):
                return np.log10(np.abs(arr).flatten())
            xxx,yyy=do(gx), do(ggg.gx)
            ext=extents()
            ext(xxx);ext(yyy)

            pch.simple_phase(xxx,yyy ,ax=ax)
            ax.plot(ext.minmax,ext.minmax)
            fig.savefig('plots_to_sort/gx_test_%s_c%04d_n%04d'%(sim_name,core_id,frame))

            fig2,ax2=plt.subplots(1,2)
            cmap='viridis'
            ext=extents()
            #ext(gx)
            Sx = 64
            Sy = slice(None)
            Sz = slice(None)
            Q1 = gy
            Q2 = ggg.gy
            print(Q1[:2,:2,:2])
            print(Q2[:2,:2,:2])
            PPP_to_plot = Q1[Sx,Sy,Sz].transpose()
            PPP2 = Q2[Sx,Sy,Sz].transpose()

            ext(PPP_to_plot)
            #ext(PPP2)
            norm = mpl.colors.Normalize( vmin=ext.minmax[0],vmax=ext.minmax[1])

            ax2[0].pcolormesh( X_to_plot_f, Y_to_plot_f, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
            #ext(PPP2)
            #norm = mpl.colors.Normalize( vmin=ext.minmax[0],vmax=ext.minmax[1])
            ax2[1].pcolormesh( X_to_plot_f, Y_to_plot_f, PPP2, norm=norm, shading='nearest', cmap=cmap)
            #pdb.set_trace()


            fig2.savefig('plots_to_sort/gx_image_%s_c%04d_n%04d'%(sim_name,core_id,frame))

            fig3,ax3=plt.subplots(1,3,figsize=(12,8))
            for DDD in [0,1,2]:
                ax3[DDD].hist( [gx,gy,gz][DDD].flatten(), histtype='step')
                ax3[DDD].hist( [ggg.gx,ggg.gy,ggg.gz][DDD].flatten(), histtype='step')
            fig3.savefig('plots_to_sort/g_hist_%s_c%04d_n%04d'%(sim_name,core_id,frame))

            
            #pch.simple_phase( 


    return cg


if 0:
    fig,ax=plt.subplots(1,3,figsize=(12,4))

    for nsim,sim in enumerate(TL.loops):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        this_looper.read_targets_only(mountain_top_fname)
        for core_id in core_list:
            ax[nsim].scatter( this_looper.targets[core_id].min_density,  this_looper.targets[core_id].peak_density)
        axbonk(ax[nsim],xscale='log',yscale='log', xlabel='min density',ylabel='peak density')

    fig.savefig('plots_to_sort/test.png')
sim_list=['u601','u602','u603']
sim_list=['u601']
#sim_list=['u602','u603']


if 1:
    for nsim,sim in enumerate(sim_list):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)
        infl=None
        massedge=None
        #infl=inflection[sim].rinflection
        #massedge=mass_edge[sim].edge
        #core_list=core_list[13:17]

        #core_list=[323]

        core_list=[74]
        
        output=energy_plots(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)

