
from starter2 import *
import figure_sublabel
reload(figure_sublabel)
from mpl_toolkits.axes_grid1 import AxesGrid
import convex_hull_tools as CHT
import convex_hull_plot2d as CHP
reload(CHP)
reload(CHT)
import colors



import three_loopers_six as TL6

plot_dir = './plots_to_sort'
proj_array=[]
#fig = plt.figure(figsize=(16,16))
fig, ax= plt.subplots(3,3, dpi=300, figsize=(12,12))#,figsize=(12,12),dpi=300)
fig.tight_layout()
delta = 0.07
fig.subplots_adjust(right=1-delta, left=delta, top=1-delta, bottom=delta)
grid = ax.transpose().flatten()
for aaa in grid:
    aaa.set_aspect('equal')

sim_list=['u601','u602','u603']
sim_list=['u603']
if 1:
#import three_loopers_mountain_top as TLM

    if 'ht' not in dir():
        ht={}
        for sim in sim_list:
            ht[sim] = CHT.hull_tool(TL6.loops[sim])
    cores_to_label={'u601':[323], 'u602':[],'u603':[]}
    TL6.loops['u601'].sublabel = figure_sublabel.labs(1.15, -0.2, r'$1a$')
    TL6.loops['u602'].sublabel = figure_sublabel.labs(1.1, -0.1, r'$1d$')
    TL6.loops['u603'].sublabel = figure_sublabel.labs(1.15, -0.2, r'$1g$')
    for ns,sim in enumerate(sim_list):
        #print('KLUDGE core list')
        full_core_list =np.unique( ht[sim].this_looper.tr.core_ids)
        #core_list = np.concatenate([core_list,[323]])
        core_list=None
        color_dict = colors.make_core_cmap(full_core_list)#, cmap = 'tab20', seed = -1)
        CHP.plot_2d(ht[sim],frames=[0], accumulate=True, label_cores=cores_to_label[sim],
                    external_axis=[grid[ns]], core_list=core_list,
                    color_dict=color_dict,
                    axis_to_plot=[0],
                    add_jitter=False,  center_image=True)


if 0:
    import image_centroid
    reload(image_centroid)
    #This is annoyingly manual; subfigure labels
    TL6.loops['u601'].sublabel = figure_sublabel.labs(1.1, -0.18, r'$1b$')
    TL6.loops['u602'].sublabel = figure_sublabel.labs(1.05, -0.1, r'$1e$')
    TL6.loops['u603'].sublabel = figure_sublabel.labs(1.1, -0.2, r'$1h$')
    if 'centool' not in dir() or True:
        centool={}
        for ns,sim in enumerate(sim_list):
            centool[sim] = image_centroid.image_track(TL6.loops[sim])
            centool[sim].run(external_ax=grid[ns+3])

if 0:
    #plot projections
    #TO DO make the xticks not suck
    if 'pt' not in dir():
        pt = {}
        stuff = {}
        for ns,sim_name in enumerate(sim_list):
            directory = dl.sims[sim_name]
            frame = dl.target_frames[sim_name]
            set_name = "%s/DD%04d/data%04d"%(directory,frame,frame)
            ds = yt.load(set_name)
            proj = ds.proj(YT_density,0)
            frb = proj.to_frb(1.0,[2048]*2)
            pt[sim_name]=frb
            stuff[sim_name]=[ds,proj] #these need to stay alive, so keep them.

    TL6.loops['u601'].sublabel = figure_sublabel.labs(0.85,0.1, r'$1c$')
    TL6.loops['u602'].sublabel = figure_sublabel.labs(0.85,0.1, r'$1f$')
    TL6.loops['u603'].sublabel = figure_sublabel.labs(0.85,0.1, r'$1i$')
    for ns,sim in enumerate(sim_list):
        frb=pt[sim]['density'].transpose()
        dx=1/frb.shape[0]
        xcoord, ycoord = np.mgrid[0:1:dx,0:1:dx]

        #grid[ns].imshow(np.log10(frb), cmap='Greys', interpolation='nearest',origin='lower')
        grid[ns+6].pcolormesh(xcoord,ycoord,np.log10(frb), cmap='Greys',shading='nearest')
        grid[ns+6].set_xlim(0,1)
        grid[ns+6].set_ylim(0,1)
        figure_sublabel.add_sublabel( grid[ns+6], TL6.loops[sim])

for ggg in grid:
    ggg.set_ylabel('')
    ggg.set_xlabel('')
for ng in [0,1,2]:
    ax[2][ng].set_xlabel(r'$y [code\ units]$')
    ax[ng][0].set_ylabel(r'$z [code\ units]$')


#for ng, aaa in enumerate(grid):
#    aaa.set_title('1%s'%'abcdefghijklmnopqrstuvwxyz'[ng])
print('Im saving is what Im doing')
fig.savefig('plots_to_sort/grid_test.png')

