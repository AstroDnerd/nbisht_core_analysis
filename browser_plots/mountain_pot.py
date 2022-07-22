from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top

def plot_mountain_top(this_looper, core_list=None, r_inflection=None, r_mass=None):

    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    ds = this_looper.load(this_looper.target_frame)
    xtra_energy.add_energies(ds)
    reload(mountain_top)
    #radius=1e-2
    radius=4/128
    radius = ds.arr(radius,'code_length')
    for core_id in core_list:
        peak = this_looper.targets[core_id].peak_location
        this_target = this_looper.targets[core_id]
        top1 = mountain_top.top(ds,peak, rhomin = this_target.min_density, peak_id=core_id, radius=radius)
        proj = ds.proj(YT_grav_energy_2,0,center=top1.location,data_source=top1.region)
        pw = proj.to_pw()
        pw.set_cmap('density','Greys')
        #pw.annotate_clumps([master_clump]+master_clump.leaves)
        if top1.leaf['particle_index'].size > 10:
            p_size = 1
        else:
            p_size = 7
        pw.annotate_these_particles4(1.0,col='r',positions= top1.leaf['particle_position'], p_size=p_size)
        pw.zoom(0.5/radius.v)
        pw.set_axes_unit('code_length')

        pw.annotate_clumps([top1.leaf], plot_args={'color':'y'})


        if r_inflection is not None:
            RRR = r_inflection[core_id]
            pw.annotate_sphere( top1.location, RRR, circle_args={'color':'r'})

        if r_mass is not None:
            if core_id in r_mass:
                RRR = r_mass[core_id]
                pw.annotate_sphere( top1.location, RRR, circle_args={'color':'b'})
        print(pw.save('plots_to_sort/mountain_top_%s_c%04d'%(this_looper.sim_name, core_id)))


if 1:
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
        if nsim != 1:
            continue
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        core_list=core_list[5:6]
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)
        infl=None
        massedge=None
        infl=None #inflection[sim].rinflection
        #massedge=mass_edge[sim].edge
        
        plot_mountain_top(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)
