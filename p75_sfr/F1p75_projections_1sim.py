
'''
author/original source for p56, figure 1: david collins
co-author/modified for p75: luz jimenez vela
'''
from starter2 import *
import figure_sublabel
reload(figure_sublabel)
from mpl_toolkits.axes_grid1 import AxesGrid
import convex_hull_tools as CHT
import convex_hull_plot2d as CHP
reload(CHP)
reload(CHT)
import colors
import track_loader as TL



# ONE SIM AT A TIME - WORKING
plot_dir = './plots_to_sort'
sim_list=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
          'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260', 'm0270', 'm0280', 'm02100', 'm02110']  

if 'ht' not in dir():
    ht={}
    pt = {}
    stuff = {}
    TL.load_tracks(sim_list)

# LETS DO A TRIO OF ONE SIM AT A TIME TO OBTAIN A MOVIE 
for ns, sim in enumerate(sim_list):
    fig, ax= plt.subplots(1,3, dpi=300, figsize=(12,4))  # rows, columns
    delta = 0.07
    fig.subplots_adjust(right=1-delta, left=delta, top=1-delta, bottom=delta)
    grid = ax.transpose().flatten()
    for aaa in grid:
        aaa.set_aspect('equal')

    # CONVEX HULLING; ax[0]
    if 1:  #maybe delete these switches (3 figs in one) later...
        ht[sim] = CHT.hull_tool(TL.loops[sim])
        full_core_list =np.unique(ht[sim].this_looper.tr.core_ids)
        core_list=None
        color_dict = colors.make_core_cmap(full_core_list)
        CHP.plot_2d(ht[sim], frames=[0], accumulate=True, external_axis=[grid[0]], core_list=core_list,
                    color_dict=color_dict, axis_to_plot=[0], add_jitter=False,  center_image=True)

    # IMAGE CENTROIDS; ax[1]
    if 1:
        import image_centroid
        reload(image_centroid)
        if 'centool' not in dir() or True:
            centool={}
            centool[sim] = image_centroid.image_track(TL.loops[sim])
            centool[sim].run(external_ax=grid[1])

    # PLOT PROJECTIONS; ax[2]
    if 1:
        frame = TL.tracks[sim].target_frame
        ds = TL.tracks[sim].load(frame)

        proj = ds.proj(YT_density,0)
        frb = proj.to_frb(1.0,[2048]*2)
        pt[sim]=frb
        stuff[sim]=[ds,proj] #these need to stay alive, so keep them.

        frb=pt[sim]['density'].transpose()
        dx=1/frb.shape[0]
        xcoord, ycoord = np.mgrid[0:1:dx,0:1:dx]

        grid[2].pcolormesh(xcoord,ycoord,np.log10(frb), cmap='Greys',shading='nearest')
        grid[2].set_xlim(0,1)
        grid[2].set_ylim(0,1)
        figure_sublabel.add_sublabel(grid[2], TL.loops[sim])

    # SAVE THE SIMULATION TRIO FIGURE
    for ggg in grid:
        ggg.set_ylabel('')
        ggg.set_xlabel('')
    ax[0].set_ylabel(r'$y [code\ units]$')
    for ng in [0,1,2]:
        ax[ng].set_xlabel(r'$z [code\ units]$')
        
    print('saving convex hull - image centroid - projection sim trio %s',sim_list[ns])
    fig.savefig('plots_to_sort/grid_trio_%s.png'%sim_list[ns])
    #fig.savefig('plots_to_sort/%d.png'%ns)
    plt.close(fig)



