from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
import pcolormesh_helper as pch
reload(xtra_energy)

def phase_proj(this_looper, core_list=None, r_inflection=None, r_mass=None):

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
            radius = max([6/128, ms.r[:,-1].max()])

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

        field_list=[YT_grav_energy, YT_grav_energy_2]
        #YT_ge_ke = ('gas','ge_ke')
        #def kege(field, data):
        #    out = data.ds.arr(np.zeros_like( data[YT_grav_energy].v), 'dimensionless')
        #    ok = np.abs(data[YT_grav_energy])>0
        #    out[ok]=data[YT_grav_energy][ok].v/data[YT_kinetic_energy][ok].v
        #    return out.v
        #ds.add_field(YT_ge_ke, function=kege, sampling_type='cell')

        proj_ext=extents()
        vx = from_cg(YT_velocity_x) 

        name1 = YT_grav_energy
        name2 = YT_grav_energy_2

        g1 = np.abs(from_cg(name1))
        g2 = np.abs(from_cg(name2))
        p1 = proj(g1)
        p2 = proj(g2)

        fig, axlist = plt.subplots(1,3, figsize=(12,8))
        for aaa in axlist:
            aaa.set_aspect('equal')
        def do_plot(PPP, ax, norm, cmap='jet'):
            X_to_plot = TheY.transpose()
            Y_to_plot = TheX.transpose()
            PPP_to_plot = PPP.transpose()
            ploot=ax.pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
            fig.colorbar(ploot, ax=ax)


        ext = extents()
        ext( g1[g1>0])
        ext( g2[g2>0])
        #pdb.set_trace()


        norm = mpl.colors.LogNorm( vmin=ext.minmax[0], vmax=ext.minmax[1])
        bins = np.geomspace( ext.minmax[0], ext.minmax[1], 64)



        do_plot( p1, axlist[0], norm)
        do_plot( p2, axlist[1], norm)

        hist, xbins, ybins = np.histogram2d( g1.flatten(), g2.flatten(), bins=[bins,bins])

        pch.helper( hist, xbins, ybins, ax=axlist[2])
        axbonk(axlist[2],xscale='log',yscale='log')


        fig.savefig('plots_to_sort/PlotPlotPhase_%s_c%04d.png'%(this_looper.sim_name, core_id))

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

if 0:
    import three_loopers_six as TL
    import r_inflection
    reload(r_inflection)
    if 'inflection' not in dir():
        inflection = {}
    for sim in sim_list:
        if sim not in inflection:
            inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
            inflection[sim].run()

if 1:
    for nsim,sim in enumerate(sim_list):
        if nsim!=1:
            continue
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)
        infl=None
        massedge=None
        infl=None#inflection[sim].rinflection
        #massedge=mass_edge[sim].edge
        core_list=core_list[13:17]
        #core_list=[323]

        
        output=phase_proj(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)

