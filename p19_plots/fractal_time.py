

from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
import movie_frames
import tools_tracks.fractal as fractal
reload(fractal)
import tsing
import track_loader as TL
sims = ['u502']
TL.load_tracks(sims)
tsing_tool = tsing.get_tsing(TL.loops)

if 'fractals' not in dir():
    for sim in sims:
        this_looper=TL.tracks[sim]
        mask = movie_frames.quantized_mask(this_looper)
        frames = this_looper.tr.frames[mask]
        frames=frames[::10]
        #core_list=[9]
        #frames=[109]
        core_list = this_looper.core_by_mode['Alone']
        core_list=None
        fractals=fractal.make_fractals(this_looper,frames,core_list=core_list)

if 1:
    sim='u502'
    frames = sorted(list(fractals.keys()))
    nc = len(fractals[frames[0]].cores_used)
    V = np.zeros([len(frames),nc])
    fig,ax=plt.subplots(1,1)
    for nf,frame in enumerate(frames):
        V[nf,:]=fractals[frame].fractal_dim
        cores_used=fractals[frame].cores_used

    t = TL.loops[sim].tr.times[mask][::10]/colors.tff
    core_list = this_looper.core_by_mode['Alone']
    #core_list = cores_used
    for nc, core_id in enumerate(core_list):
        nf = np.where(cores_used==core_id)[0][0]

        tsung = tsing_tool[sim].tend_core[core_id]
        #tsung = tsing_tool[sim].tsing_core[core_id]
        Y = V[:,nf]
        #Y /= Y[0]
        ax.plot(t/tsung, Y, c=[0.5]*4)
        ax.set(xlabel=r'$t/t_{\rm{sung}}$', ylabel=r'$D$', xlim=[0,1.2])
        #ax.plot(t, V[:,nc])
        #ax.plot(V-V[0,:])
    fig.savefig('plots_to_sort/fractal.pdf')



if 0:
    rm = rainbow_map(10)
    for name in toolshed:
        dims = []
        #mylooper=looper_list[0]
        mylooper=toolshed[name]['looper']
        for nframe,frame in enumerate(toolshed[name]):
            tool = toolshed[name][frame]
            lab=ft.this_looper.out_prefix + " %s"%frame
            plt.hist(tool.fractal_dim,histtype='step',bins=10,color=rm(nframe),label=lab)
            dims.append(tool.fractal_dim)
        plt.savefig("plots_to_sort/wtf_%d.png"%nf)
        plt.savefig('plots_to_sort/fractal_dist.png')
        fig,ax=plt.subplots(1,1)
        ax.violinplot(dims)
        fig.savefig('plots_to_sort/violin_test.png')
        plt.close(fig)
        fig,axes=plt.subplots(1,2)
        dims=nar(dims)

        cores_used = tool.cores_used
        OverAchievers = []
        for nb,b in enumerate(dims.transpose()):
            core_id = cores_used[nb]
            if (b > 2.5).any():
                OverAchievers.append(core_id)
            nparticles=mylooper.snaps[0][core_id].target_indices.size
            if nparticles < 10:
                continue

            axes[0].plot(mylooper.tr.times[framelist], b/b[0:3].mean())
            axes[1].plot(mylooper.tr.times[framelist], b)
        fig.savefig('plots_to_sort/streams.png')
