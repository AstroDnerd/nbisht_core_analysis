from starter2 import *

import convex_hull_tools as CHT
import convex_hull_plot2d as CHP
reload(CHP)
import three_loopers_u500 as TL
import supersets
sim_list = ['u501','u502','u503']
sim_list = ['u502']
if 'ht' not in dir():
    ht={}
    for this_simname in sim_list:
        if this_simname in ht:
            continue
        print("Analysis: Hulls")
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'st' not in dir():
    st={}
    for sim in sim_list:
        if sim in st:
            continue
        print("Analysis: supersets")
        st[sim] = supersets.superset( TL.loops[sim], ht[sim])
        st[sim].find()


for sim in sim_list:
    for nset,superset in enumerate(st[sim].supersets):

        
        prefix = "isopot_s%04d"%nset

        x_mean = 0
        xs=[]
        for core_id in superset:
            ms = trackage.mini_scrubber(TL.loops[sim].tr, core_id, do_velocity=False)
            c = np.stack([ms.mean_x[0],ms.mean_y[0],ms.mean_z[0]])
            x_mean = c + x_mean
            xs.append(ms.mean_x[0])
        x_mean /= len(superset)
        if x_mean[0] > 1:
            x_mean[0] -= 1
        if x_mean[0] < 0:
            x_mean[0] += 1

        ds = TL.loops[sim].load(0)
        #proj=yt.ProjectionPlot(ds, 0, YT_density)
        proj=yt.SlicePlot(ds, 0,center=[x_mean[0],0.5,0.5],fields=[YT_grav_energy_2])
        #image = proj.data_source.to_frb(1,128)[YT_potential_field]
        image = proj.data_source.to_frb(1,128)[YT_grav_energy_2]
        
        fig,ax=plt.subplots(1,1)
        ax.imshow(image, extent=[0,1,0,1])
        CHP.plot_2d(ht[sim], core_list=superset, accumulate=True,all_plots=False,axis_to_plot=[0],frames=[0], prefix=prefix,
                    external_axis=[ax], plot_particles=False)
        fig.savefig('plots_to_sort/isopot_%s'%prefix)



