from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
reload(xtra_energy)

def plot_mountain_top(this_looper, core_list=None, r_inflection=None, r_mass=None):

    sim_name = this_looper.sim_name
    thtr=this_looper.tr
    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    frame=this_looper.target_frame
    ds = this_looper.load(frame)
    if 0:
        #projections.
        #I'm not quite right.
        frame = 0
        ds = yt.load('/data/cb1/Projects/P19_CoreSimulations/new_sims/vel_test/DD%04d/data%04d'%(frame,frame))
        for axis in [0,1,2]:
            proj = ds.proj('density',axis)
            pw=proj.to_pw()
            pw.annotate_velocity()
            pw.save('plots_to_sort/arse')
    xtra_energy.add_energies(ds)
    xtra_energy.add_gravity(ds)
    reload(mountain_top)
    #radius=1e-2
    for core_id in core_list:
        ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True) 
        ms.get_central_at_once(core_id)

        if r_inflection is not None and False:
            radius = r_inflection[core_id]
            radius = ds.arr(radius, 'code_length')
        else:
            radius=1/128
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

        SphereVbar = None
        if 1:
            if r_inflection is not None:
                RRR = r_inflection[core_id]
                sph = ds.sphere( peak, RRR)
                M = (sph[YT_density]*sph[YT_cell_volume]).sum()
                SphereVbar = [( (sph[('gas','velocity_%s'%s)]*sph[YT_cell_mass]).sum()/M ).v for s in 'xyz']
                #SphereVbar = [( (sph[('gas','momentum_%s'%s)]*sph[YT_cell_mass]).sum()/M ).v for s in 'xyz']



        cg = ds.covering_grid(4,left,nzones, num_ghost_zones=1)
        #cg = ds.covering_grid(0,[0.0]*3,[128]*3, num_ghost_zones=1)
        def from_cg(field):
            ccc = cg[field]#.swapaxes(0,2)
            return ccc.v
        def proj(arr):
            return arr.sum(axis=proj_axis)
        fig,ax=plt.subplots(1,1)
        ax.set_aspect('equal')


        #field = ('gas','momentum_x')
        field='density'
        #field=YT_grav_energy
        if type(field) is tuple:
            field_name = field[1]
        else:
            field_name=field
        rho= from_cg(field)
        dv = from_cg(YT_cell_volume)
        PPP= proj(rho)

        if field_name in ['density']:
            norm = mpl.colors.LogNorm(vmin=PPP.min(), vmax=PPP.max())
        else:
            norm = mpl.colors.Normalize(vmin=PPP.min(), vmax=PPP.max())

        x_h = proj( from_cg( ('gas','xyz'[h_axis]) ) )
        x_v = proj( from_cg( ('gas','xyz'[v_axis]) ) )


        "H for horizontal, V for vertical"
        if 0:
            stream_name='momentum'
            stream_h = 'momentum_%s'%'xyz'[h_axis]
            stream_v = 'momentum_%s'%'xyz'[v_axis]
        if 1:
            stream_name='velocity'
            stream_h = 'velocity_%s'%'xyz'[h_axis]
            stream_v = 'velocity_%s'%'xyz'[v_axis]
        if 0:

            stream_name='accel'
            stream_h = 'grav_%s'%'xyz'[h_axis]
            stream_v = 'grav_%s'%'xyz'[v_axis]




        Q_H = from_cg(stream_h)
        Q_V = from_cg(stream_v)
        M = (rho*dv).sum()
        if 0:
            Q_H_mean=0
            Q_V_mean=0
        if 0:
            Q_H_mean = Q_H.mean()
            Q_V_mean = Q_V.mean()
        if 0:
            stream_name += "_whole_mean"
            #works pretty good
            Q_H_mean = (Q_H*dv*rho).sum()/M
            Q_V_mean = (Q_V*dv*rho).sum()/M
        if 0:
            mean_vx = ms.mean_vx[-1]
            mean_vy = ms.mean_vy[-1]
            mean_vz = ms.mean_vz[-1]
            mean = [mean_vx, mean_vy, mean_vz]
            Q_H_mean = mean[ h_axis]
            Q_V_mean = mean[ v_axis]
        if 1:
            stream_name += "_vcentral"
            Q_H_mean = ms.vcentral[h_axis,-1]
            Q_V_mean = ms.vcentral[v_axis,-1]
        if 0:
            LLL = top1.leaf
            Q_H_mean =((LLL['density']*LLL['cell_volume']*LLL[stream_h]).sum()/(LLL['density']*LLL['cell_volume']).sum()).v
            Q_V_mean =((LLL['density']*LLL['cell_volume']*LLL[stream_v]).sum()/(LLL['density']*LLL['cell_volume']).sum()).v
        if 0:
            LLL = top1.leaf
            Q_H_mean =((LLL['cell_volume']*LLL[stream_h]).sum()/(LLL['cell_volume']).sum()).v
            Q_V_mean =((LLL['cell_volume']*LLL[stream_v]).sum()/(LLL['cell_volume']).sum()).v
        if 0:
            Q_H_mean = SphereVbar[ h_axis]
            Q_V_mean = SphereVbar[ v_axis]
        Q_H -= Q_H_mean
        Q_V -= Q_V_mean

        Q_H_p = proj(Q_H)
        Q_V_p = proj(Q_V)



        #TheX,TheY,TheU,TheV= x_h, x_v, Q_H, Q_V
        TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))
        TheU,TheV= Q_H_p, Q_V_p
        if 0:
            print(x_h)
            print('poooo')
            print(TheX.size)
            print(TheY.size)
            print(TheU.size)

        gg = TheU**2+TheV**2
        ok = gg>0


        if 0:
            X_to_plot   = TheX#.transpose()
            Y_to_plot   = TheY#.transpose()
            PPP_to_plot =  PPP#.transpose()
            U_to_plot =   TheU#.transpose()
            V_to_plot =   TheV#.transpose()
        else:
            X_to_plot   = TheY.transpose()
            Y_to_plot   = TheX.transpose()
            PPP_to_plot =  PPP.transpose()
            U_to_plot =   TheU.transpose()
            V_to_plot =   TheV.transpose()



        #ax.imshow( PPP,norm=norm)
        ploot=ax.pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest')
        fig.colorbar(ploot)

        #U_to_plot[ PPP_to_plot == 0] = np.nan
        #V_to_plot[ PPP_to_plot == 0] = np.nan
        ax.streamplot(X_to_plot, Y_to_plot, U_to_plot, V_to_plot, density=2.5)


        axbonk(ax, xlabel='xyz'[h_axis], ylabel='xyz'[v_axis])
        outname='plots_to_sort/proj_%s_%s_c%04d_n%04d_%s_%s.png'%('xyz'[proj_axis],sim_name, core_id, frame,field_name, stream_name)
        fig.savefig(outname)
        print(outname)
#       if 0:
#           acceleration = [YT_grav_x, YT_grav_y, YT_grav_z]
#           pw.annotate_streamlines( acceleration[h_axis], acceleration[v_axis], plot_args={'color':'g'})
#           streamline = '_acceleration'
#       if 0:
#           mean_vx = ms.mean_vx[-1]
#           mean_vy = ms.mean_vy[-1]
#           mean_vz = ms.mean_vz[-1]
#           mean = [mean_vx, mean_vy, mean_vz]
#       if 1:
#           reg = top1.region
#           M = reg['cell_volume'].sum().v
#           mean = [ (reg['cell_volume']*reg['velocity_%s'%s]).sum().v/M for s in 'xyz']
#           velocity = [YT_velocity_x, YT_velocity_y, YT_velocity_z]
#           pw.annotate_streamlines( velocity[h_axis], velocity[v_axis], plot_args={'color':'r'}, norm=[mean[h_axis],mean[v_axis]])
#           #pw.annotate_streamlines( velocity[h_axis], velocity[v_axis], plot_args={'color':'r'}, norm=[0.24,-0.13])
#           streamline = 'velocity'
#           print("MEAN",mean)

#       if 1:
#           fig3,ax3=plt.subplots(1,1)
#           ax3.hist( reg[YT_velocity_x].v, histtype='step',color='r')
#           ax3.hist( reg[YT_velocity_y].v, histtype='step',color='g')
#           ax3.hist( reg[YT_velocity_z].v, histtype='step',color='b')
#           ax3.hist( reg[YT_velocity_magnitude].v, histtype='step',color='k')
#           fig3.savefig('plots_to_sort/vel_hist.png')


#       #pw.annotate_clumps([master_clump]+master_clump.leaves)
#       if top1.leaf['particle_index'].size > 10:
#           p_size = 1
#       else:
#           p_size = 7
#       pw.annotate_these_particles4(1.0,col='r',positions= top1.leaf['particle_position'], p_size=p_size)
#       pw.zoom(0.5/radius.v)
#       pw.set_axes_unit('code_length')

#       pw.annotate_clumps([top1.leaf], plot_args={'color':'y'})


#       if r_inflection is not None:
#           RRR = r_inflection[core_id]
#           pw.annotate_sphere( top1.location, RRR, circle_args={'color':'r'})

#       if r_mass is not None:
#           if core_id in r_mass:
#               RRR = r_mass[core_id]
#               pw.annotate_sphere( top1.location, RRR, circle_args={'color':'b'})
#       print(pw.save('plots_to_sort/mountain_top_%s_c%04d%s'%(this_looper.sim_name, core_id, streamline)))


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
        
        output=plot_mountain_top(this_looper,core_list=None,
                          r_inflection=infl)
                          #r_mass=massedge)
