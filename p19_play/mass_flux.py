from starter2 import *
import xtra_energy
import three_loopers_u500 as TL
import mountain_top
import xtra_energy
import pcolormesh_helper as pch
reload(xtra_energy)

def mass_flux(this_looper, core_list=None, r_inflection=None, r_mass=None):


    sim_name = this_looper.sim_name
    thtr=this_looper.tr
    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    frame_list=thtr.frames[-2:]
    ds_list = [this_looper.load(f) for f in frame_list]
    frame_id = nar([np.where(thtr.frames==f)[0][0] for f in frame_list])
    times = thtr.times[frame_id]

    #print(times)

    for ds in ds_list:
        xtra_energy.add_energies(ds)
        xtra_energy.add_gravity(ds)
    reload(mountain_top)
    #radius=1e-2
    for core_id in core_list:
        print("Proj energies %s c%04d"%(sim_name, core_id))
        ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True) 
        ms.get_central_at_once(core_id)
        ms.particle_pos(core_id)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)

        if 0:
            r_infl = r_inflection[core_id]
            radius=max([1/128, 1.2*r_infl])
        if 0:
            radius = max([2/128, ms.r[:,-1].max()])
        if 0:
            radius = 2/128
        if 1:
            radius = 1/128

        radius = ds.arr(radius,'code_length')

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

        stuff = defaultdict(list)
        for nf, frame in enumerate(frame_list):
            ds = ds_list[nf]
            #peak = this_looper.targets[core_id].peak_location
            #point = np.column_stack([ ms.mean_xc[frame_id[nf]], ms.mean_yc[frame_id[nf]], ms.mean_zc[frame_id[nf]]])
            point = nar([ms.mean_xc[-1], ms.mean_yc[-1], ms.mean_zc[-1]])

            left =  point - radius.v
            right = point + radius.v
            dx = 1/2048
            nzones = np.floor((right-left)/dx).astype('int')

            cg = ds.covering_grid(4,left,nzones, num_ghost_zones=1)


            if 0:
                fig,ax=plt.subplots(1,1)

                h_axis = ds.coordinates.x_axis[proj_axis]
                v_axis = ds.coordinates.y_axis[proj_axis]
                x_h = proj( from_cg( ('gas','xyz'[v_axis]) ) )#/nzones[h_axis]
                x_v = proj( from_cg( ('gas','xyz'[h_axis]) ) )#/nzones[v_axis]
                TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))

                divv = np.abs(from_cg('velocity_divergence'))
                PPP = divv[divv.shape[0]//2,:,:]#density.sum(axis=proj_axis)
                X_to_plot = TheY.transpose()
                Y_to_plot = TheX.transpose()
                PPP_to_plot = PPP.transpose()

                #vmin,vmax = 1e-3,1e3
                cmap = 'seismic'
                vmin = PPP_to_plot.min()
                vmin = density[density>0].min()
                vmax = PPP_to_plot.max()
                norm = mpl.colors.LogNorm( vmin=vmin, vmax=vmax)

                ploot=ax.pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
                fig.colorbar(ploot, ax=ax)

                outname='plots_to_sort/divv_%s_%s_c%04d_n%04d.png'%('xyz'[proj_axis],sim_name, core_id, frame)
                fig.savefig(outname)
                print(outname)

            S1 = tuple([slice(1,-1),slice(1,-1),slice(1,-1)])
            divv=from_cg('velocity_divergence')
            rho = from_cg(YT_density)
            Fxf = from_cg('velocity_x')*rho
            Fyf = from_cg('velocity_y')*rho
            Fzf = from_cg('velocity_z')*rho
            ip1 = slice(2,None)
            im1 = slice(None,-2)
            iii = slice(1,-1)
            dv = from_cg(YT_cell_volume)
            dx = np.unique( dv**(1./3))
            my_divv = np.zeros_like(Fxf)
            #ij=slice(None);ik=slice(None)
            my_divv[iii,:,:]  = ( Fxf[ip1,:,:]-Fxf[im1,:,:])/(2*dx)
            my_divv[:,iii,:] += ( Fyf[:,ip1,:]-Fyf[:,im1,:])/(2*dx)
            my_divv[:,:,iii] += ( Fzf[:,:,ip1]-Fzf[:,:,im1])/(2*dx)
            #shave=np.zeros_like(vxf)
            #shave[S1]=my_divv[S1]
            #shave[S1]=divv*dx

            i1 = (Fxf[-1,iii,iii]+Fxf[-2,iii,iii] -Fxf[1,iii,iii]- Fxf[0,iii,iii]).sum()/2
            i1 +=(Fyf[iii,-1,iii]+Fyf[iii,-2,iii] -Fyf[iii,1,iii]- Fyf[iii,0,iii]).sum()/2
            i1 +=(Fzf[iii,iii,-1]+Fzf[iii,iii,-2] -Fzf[iii,iii,1]- Fzf[iii,iii,0]).sum()/2

            stuff['mass'].append(cg['cell_mass'].sum())
            stuff['flux'].append(i1*dx**2)
            stuff['drhodt'].append( (my_divv*dv).sum())
        for thing in stuff:
            stuff[thing]=nar(stuff[thing]).flatten()

        #print(stuff)
        M =nar( stuff['mass'])
        dM = M[1:]-M[:-1]
        dt = times[1:]-times[:-1]
        dMdt=dM/dt
        F = -nar(stuff['flux'])
        print('dmdt',dMdt)
        print('f0',F[0])
        print('f1',F[1])
        fh=F.sum()
        print('fh',fh)
        print('wtf',dMdt/fh)
        #print('drhodt',stuff['drhodt'])
        print('should',stuff['flux'] /stuff['drhodt'])



    return cg
sim_list=['u601','u602','u603']
#sim_list=['u601']
#sim_list=['u602','u603']
sim_list=['u502']

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
        core_list=[114]
        
        output=mass_flux(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)

