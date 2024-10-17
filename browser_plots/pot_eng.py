from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
import pcolormesh_helper as pch
reload(xtra_energy)

def energy_plots(this_looper, core_list=None, r_inflection=None, r_mass=None):

    sim_name = this_looper.sim_name
    thtr=this_looper.tr
    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    frame=this_looper.target_frame
    ds = this_looper.load(frame)
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

        stream_name = 't2'
        if 0:
            r_infl = r_inflection[core_id]
            radius=max([1/128, 1.2*r_infl])
        if 1:
            radius = max([2/128, ms.r[:,-1].max()])

        radius = ds.arr(radius,'code_length')

        peak = this_looper.targets[core_id].peak_location
        this_target = this_looper.targets[core_id]
        top1=None
        if 0:
            top1 = mountain_top.top(ds,peak, rhomin = this_target.min_density, peak_id=core_id, radius=radius)
        proj_axis=0
        if proj_axis == 1:
            print("ERROR: plotting messed up for axis 1")
            raise
        h_axis = ds.coordinates.x_axis[proj_axis]
        v_axis = ds.coordinates.y_axis[proj_axis]
        left = peak - radius.v
        right = peak + radius.v
        dx = 1/2048
        nzones = np.floor((right-left)/dx).astype('int')

        cg = ds.covering_grid(4,left,nzones, num_ghost_zones=1)
        #cg = ds.covering_grid(0,[0.0]*3,[128]*3, num_ghost_zones=1)
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

        fig,axes=plt.subplots(2,2, figsize=(22,22))
        axlist=axes.flatten()
        #for ax in axlist:
        #    ax.set_aspect('equal')

        #TheX,TheY,TheU,TheV= x_h, x_v, Q_H, Q_V
        x_h = proj( from_cg( ('gas','xyz'[v_axis]) ) )#/nzones[h_axis]
        x_v = proj( from_cg( ('gas','xyz'[h_axis]) ) )#/nzones[v_axis]

        TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))
        KEP = 'ge_ke'
        field_list=[YT_grav_energy, YT_kinetic_energy, YT_ge_ke]
        #YT_ge_ke = ('gas','ge_ke')
        #def kege(field, data):
        #    out = data.ds.arr(np.zeros_like( data[YT_grav_energy].v), 'dimensionless')
        #    ok = np.abs(data[YT_grav_energy])>0
        #    out[ok]=data[YT_grav_energy][ok].v/data[YT_kinetic_energy][ok].v
        #    return out.v
        #ds.add_field(YT_ge_ke, function=kege, sampling_type='cell')

        proj_ext=extents()
        density = from_cg(YT_density)
        for nf,field in enumerate(field_list):
            if nf > 2:
                continue
            rho= from_cg(field)
            PPP= np.abs(proj(rho))
            proj_ext(PPP[PPP>0])
        for nf,field in enumerate(field_list):
            #field=YT_grav_energy
            if type(field) is tuple:
                field_name = field[1]
            else:
                field_name=field
            if field_name not in ['keprime_ge']:
                rho= np.abs(from_cg(field))
                dv = from_cg(YT_cell_volume)
                if nf != 2:
                    PPP= proj(rho)
                else:
                    PPP = proj_m(rho,density)
            else:
                vx = from_cg(YT_velocity_x) - ms.mean_vx[-1]
                vy = from_cg(YT_velocity_y) - ms.mean_vy[-1]
                vz = from_cg(YT_velocity_z) - ms.mean_vz[-1]
                rho= from_cg(YT_density)
                GE = from_cg(YT_grav_energy)
                p = 0.5*rho*(vx*vx+vy*vy+vz*vz)
                PPP = proj(GE/p)

            if field_name in ['grav_energy']:
                PPP = np.abs(PPP)

            if field_name in ['density', 'grav_energy', 'kinetic_energy_density']:
                norm = mpl.colors.LogNorm(vmin=PPP[PPP>0].min(), vmax=PPP.max())
            else:

                norm = mpl.colors.Normalize(vmin=PPP.min(), vmax=PPP.max())
            if nf < 2:
                vmin, vmax = proj_ext.minmax
                cmap = None
            else:
                vmin,vmax = 1e-3,1e3
                cmap = 'seismic'
            norm = mpl.colors.LogNorm( vmin=vmin, vmax=vmax)

            axlist[nf].set_title(field_name)

            #ax.imshow( PPP,norm=norm)
            X_to_plot = TheY.transpose()
            Y_to_plot = TheX.transpose()
            PPP_to_plot = PPP.transpose()
            ploot=axlist[nf].pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
            fig.colorbar(ploot, ax=axlist[nf])

            if 0:
                circle_x = peak[v_axis]
                circle_y = peak[h_axis]
                circle = plt.Circle( (circle_x, circle_y), r_infl, edgecolor='r', facecolor='None')
                axlist[nf].add_patch(circle)
        if 1:
            #fig7,ax7=plt.subplots(1,1)
            xyz = [ms.particle_x[:,-1], ms.particle_y[:,-1], ms.particle_z[:,-1]]
            axlist[0].scatter( xyz[h_axis], xyz[v_axis], c='k')
            axlist[1].scatter( xyz[h_axis], xyz[v_axis], c='k')
            axlist[2].scatter( xyz[h_axis], xyz[v_axis], c='k')
            #xyz = [ms.this_x[:,-1], ms.this_y[:,-1], ms.this_z[:,-1]]
            #xyz = [ms.particle_x[:,-1], ms.particle_y[:,-1], ms.particle_z[:,-1]]
            #axlist[1].scatter( xyz[h_axis], xyz[v_axis])
            #print(ms.nparticles)
            #fig7.savefig('plots_to_sort/wtf.png')

        if 1:
            ge = np.abs(from_cg(YT_grav_energy))
            ke = np.abs(from_cg(YT_kinetic_energy))
            xxbins = np.geomspace(ke.min(),ke.max(),128)
            yybins = np.geomspace(ge[ge>0].min(),ge.max(),128)
            hist, xbins,ybins=np.histogram2d(ke[ge>0].flatten(),ge[ge>0].flatten(),bins=[xxbins,yybins])
            pch.helper(hist,xbins,ybins,ax=axlist[3])
            axbonk(axlist[3],xscale='log',yscale='log',xlabel='KE',ylabel='GE')
            axlist[3].scatter(ms.ke[:,-1],np.abs(ms.ge[:,-1]), edgecolor='r',s=30, facecolor='None')
            axlist[3].plot(xxbins,xxbins,c='k')
            #ms.compute_ke(core_id)
            #axlist[2].scatter( ms.r[:,-1], np.abs(ms.ke[:,-1]), facecolor='None',edgecolor='r')
            #print(ms.ge[:,-1])
            #axbonk(axlist[2],xlabel='r',ylabel='KE',yscale='log',xscale='log')




            #axbonk(ax, xlabel='xyz'[h_axis], ylabel='xyz'[v_axis])
        outname='plots_to_sort/eng_%s_%s_c%04d_n%04d_%s_%s.png'%('xyz'[proj_axis],sim_name, core_id, frame,field_name, stream_name)
        fig.savefig(outname)
        print(outname)
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
#sim_list=['u601']
#sim_list=['u602','u603']

if 1:
    import three_loopers_six as TL
    import r_inflection
    reload(r_inflection)
    if 'inflection' not in dir():
        inflection = {}
    for sim in sim_list:
        if sim not in inflection:
            inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
            inflection[sim].run()
if 0:
    import three_loopers_six as TL
    import mass_edge_tool
    reload(mass_edge_tool)
    if 'mass_edge' not in dir():
        mass_edge = {}
        for sim in sim_list:
            mass_edge[sim]=mass_edge_tool.edge( TL.loops[sim])
            mass_edge[sim].run()

if 1:
    for nsim,sim in enumerate(sim_list):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)
        infl=None
        massedge=None
        infl=inflection[sim].rinflection
        #massedge=mass_edge[sim].edge
        #core_list=core_list[13:17]
        #core_list=[323]
        
        output=energy_plots(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)

