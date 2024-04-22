
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




# THREE SIMS AT A TIME: WORKING
plot_dir = './plots_to_sort'
proj_array=[]
fig, ax= plt.subplots(3,3, dpi=300, figsize=(12,12))  # rows, columns
# maybe needed, tends to cause problems
# fig.tight_layout()
delta = 0.07
fig.subplots_adjust(right=1-delta, left=delta, top=1-delta, bottom=delta)
grid = ax.transpose().flatten()
for aaa in grid:
    aaa.set_aspect('equal')


# LETS DO A TRIO OF ONE SIM AT A TIME TO OBTAIN A MOVIE 
sim_list=['m0230', 'm0231', 'm0232']
#sim_list=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
#          'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260', 'm0270', 'm0280', 'm02100', 'm02110']  

# CONVEX HULLING
if 1:
    if 'ht' not in dir():
        ht={}
        TL.load_tracks(sim_list)
        for sim in sim_list:
            ht[sim] = CHT.hull_tool(TL.loops[sim])
    for ns,sim in enumerate(sim_list):
        full_core_list =np.unique( ht[sim].this_looper.tr.core_ids)
        core_list=None
        color_dict = colors.make_core_cmap(full_core_list)
        CHP.plot_2d(ht[sim],frames=[0], accumulate=True, 
                    external_axis=[grid[ns]], core_list=core_list,
                    color_dict=color_dict,
                    axis_to_plot=[0],
                    add_jitter=False,  center_image=True)
# IMAGE CENTROIDS
if 1:
    import image_centroid
    reload(image_centroid)
    if 'centool' not in dir() or True:
        centool={}
        for ns,sim in enumerate(sim_list):
            centool[sim] = image_centroid.image_track(TL.loops[sim])
            centool[sim].run(external_ax=grid[ns+3])
# PLOT PROJECTIONS
if 1:
    if 'pt' not in dir():
        pt = {}
        stuff = {}
        for ns,sim_name in enumerate(sim_list):
            frame = TL.tracks[sim_name].target_frame
            ds = TL.tracks[sim_name].load(frame)

            proj = ds.proj(YT_density,0)
            frb = proj.to_frb(1.0,[2048]*2)
            pt[sim_name]=frb
            stuff[sim_name]=[ds,proj] #these need to stay alive, so keep them.

    for ns,sim in enumerate(sim_list):
        frb=pt[sim]['density'].transpose()
        dx=1/frb.shape[0]
        xcoord, ycoord = np.mgrid[0:1:dx,0:1:dx]

        grid[ns+6].pcolormesh(xcoord,ycoord,np.log10(frb), cmap='Greys',shading='nearest')
        grid[ns+6].set_xlim(0,1)
        grid[ns+6].set_ylim(0,1)
        figure_sublabel.add_sublabel(grid[ns+6], TL.loops[sim])

for ggg in grid:
    ggg.set_ylabel('')
    ggg.set_xlabel('')
for ng in [0,1,2]:
    ax[2][ng].set_xlabel(r'$y [code\ units]$')
    ax[ng][0].set_ylabel(r'$z [code\ units]$')
    
print('saving convex hull - image centroid - projection trio')
fig.savefig('plots_to_sort/grid_test_trio')



