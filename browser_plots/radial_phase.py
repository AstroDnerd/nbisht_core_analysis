from starter2 import *
import colors
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
reload(xtra_energy)
import pcolormesh_helper as pch

reload(trackage)
def radial_phase(this_looper, core_list=None, r_inflection=None, r_mass=None):

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

        r_infl = r_inflection[core_id]

        radius=max([1/128, r_infl])
        #radius = 1/128
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
        #cg = ds.sphere( peak, radius)
        def from_cg(field):
            ccc = cg[field]#.swapaxes(0,2)
            return ccc.v

        fig,axes=plt.subplots(2,2, figsize=(22,22))
        axlist=axes.flatten()
        #for ax in axlist:
        #    ax.set_aspect('equal')

        ##TheX,TheY,TheU,TheV= x_h, x_v, Q_H, Q_V
        #x_h = proj( from_cg( ('gas','xyz'[h_axis]) ) )#/nzones[h_axis]
        #x_v = proj( from_cg( ('gas','xyz'[v_axis]) ) )#/nzones[v_axis]


        x = from_cg(YT_x)-peak[0]
        y = from_cg(YT_y)-peak[1]
        z = from_cg(YT_z)-peak[2]

        radius = np.sqrt(x*x+y*y+z*z)
        rho=from_cg(YT_density)
        dv=from_cg(YT_cell_volume)
        ge=np.abs(from_cg(YT_grav_energy))
        ke_yt=np.abs(from_cg(YT_kinetic_energy))
        vx=from_cg(YT_velocity_x)
        vy=from_cg(YT_velocity_y)
        vz=from_cg(YT_velocity_z)
        #ke = 0.5*rho*(vx*vx+vy*vy+vz*vz)
        ke=ke_yt

        if 1:
            rbins = np.geomspace(1/2048,radius.max(),128)
            dbins = np.geomspace(rho.min(),rho.max(),127)
            hist, xbins,ybins=np.histogram2d(radius.flatten(),rho.flatten(),bins=[rbins,dbins], weights=dv.flatten())
            pch.helper(hist,xbins,ybins,ax=axlist[0])
            axbonk(axlist[0],xlabel='r',ylabel='rho',yscale='log',xscale='log')
            axlist[0].scatter( ms.r[:,-1], ms.density[:,-1])


            rbins = np.geomspace(1/2048,radius.max(),128)
            GEbins = np.geomspace(ge[ge>0].min(),ge.max(),127)
            hist, xbins,ybins=np.histogram2d(radius[ge>0].flatten(),ge[ge>0].flatten(),bins=[rbins,GEbins], weights=dv[ge>0].flatten())
            pch.helper(hist,xbins,ybins,ax=axlist[1])
            ms.compute_ge(core_id)
            axlist[1].scatter( ms.r[:,-1], np.abs(ms.ge[:,-1]),c='r')
            #print(ms.ge[:,-1])
            axbonk(axlist[1],xlabel='r',ylabel='AbsGE',yscale='log',xscale='log')
        
        rbins = np.geomspace(1/2048,radius.max(),128)
        KEbins = np.geomspace(ke[ke>0].min(),ke.max(),127)
        hist, xbins,ybins=np.histogram2d(radius[ke>0].flatten(),ke[ke>0].flatten(),bins=[rbins,KEbins])
        pch.helper(hist,xbins,ybins,ax=axlist[2])
        ms.compute_ke(core_id)
        axlist[2].scatter( ms.r[:,-1], np.abs(ms.ke[:,-1]), facecolor='None',edgecolor='r')
        #print(ms.ge[:,-1])
        axbonk(axlist[2],xlabel='r',ylabel='KE',yscale='log',xscale='log')

        xxbins = np.geomspace(ke.min(),ke.max(),128)
        yybins = np.geomspace(ke_yt.min(),ke_yt.max(),127)
        hist, xbins,ybins=np.histogram2d(ke.flatten(),ke_yt.flatten(),bins=[xxbins,yybins])
        pch.helper(hist,xbins,ybins,ax=axlist[3])
        ms.compute_ke(core_id)
        #axlist[2].scatter( ms.r[:,-1], np.abs(ms.ke[:,-1]), facecolor='None',edgecolor='r')
        #print(ms.ge[:,-1])
        axbonk(axlist[3],xlabel='r',ylabel='KE',yscale='log',xscale='log')
        #TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))
        #YT_ge_ke = ('gas','ge_ke')
        #KEP = 'ge_ke'
        #field_list=[YT_grav_energy, YT_kinetic_energy, YT_ge_ke, KEP]
        #def kege(field, data):
        #    out = data.ds.arr(np.zeros_like( data[YT_grav_energy].v), 'dimensionless')
        #    ok = np.abs(data[YT_grav_energy])>0
        #    out[ok]=data[YT_grav_energy][ok].v/data[YT_kinetic_energy][ok].v
        #    return out.v
        #ds.add_field(YT_ge_ke, function=kege, sampling_type='cell')




            #axbonk(ax, xlabel='xyz'[h_axis], ylabel='xyz'[v_axis])
        stream_name = 't1'
        field_name='rho_r'
        outname='plots_to_sort/phase_%s_%s_c%04d_n%04d_%s_%s.png'%('xyz'[proj_axis],sim_name, core_id, frame,field_name, stream_name)
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
        core_list=[323]
        
        output=radial_phase(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)
        break
