from starter2 import *
import mountain_top
import track_info
import xtra_energy

def plot_mountain_top(this_looper, core_list=None, r_inflection=None, r_mass=None):

    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    ds = this_looper.load(this_looper.target_frame)
    radius=1e-2
    #radius=4/128
    radius = ds.arr(radius,'code_length')
    for core_id in core_list:
        peak = this_looper.targets[core_id].peak_location
        this_target = this_looper.targets[core_id]
        top1 = mountain_top.top(ds,peak, rhomin = this_target.min_density, peak_id=core_id, radius=radius)
        proj = ds.proj('density',0,center=top1.location,data_source=top1.region)
        pw = proj.to_pw(origin='domain')

        pw.set_cmap('density','Greys')
        #pw.annotate_clumps([master_clump]+master_clump.leaves)
        if top1.leaf[YT_particle_index].size > 10:
            p_size = 1
        else:
            p_size = 7
        pw.annotate_these_particles4(1.0,col='r',positions= top1.leaf[YT_particle_position], p_size=p_size)
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


import track_loader as TL

def image_mountains(trackname):

    TL.load_tracks(trackname)

    this_looper = TL.tracks[trackname]
    core_list=np.unique(this_looper.tr.core_ids)
    mountain_top_fname = track_info.tracks[trackname].mountain_top
    if this_looper.targets is None:
        this_looper.read_targets_only(mountain_top_fname)
    
    plot_mountain_top(this_looper,core_list=core_list)
