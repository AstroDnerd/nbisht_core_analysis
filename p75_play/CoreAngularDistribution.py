import warnings
warnings.filterwarnings("ignore")

import matplotlib.ticker as tck
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=False) # Use LaTeX font

from starter2 import *
import track_loader as TL
import ucore
#reload(ucore)

sim_list=['u502']
result_dir = 'u502_core_analysis/'
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)

def make_dir(dir_path):
        if not os.path.exists(dir_path):
                print("making directory:",dir_path)
                os.makedirs(dir_path)

def CS_Transformation(px,py,pz):
    XsqPlusYsq = px**2 + py**2
    r = np.sqrt(XsqPlusYsq + pz**2)
    theta = np.arctan2(pz,np.sqrt(XsqPlusYsq))
    phi = np.arctan2(py,px)
    return r, theta, phi


def CoreDist_allframes(mon,core_id,prefix='DD', reg_dist = False):
        ms = mon.get_ms(core_id)
        ms.particle_pos(core_id)
        NumOfParticles = len(ms.particle_id)
        t = mon.all_times/colors.tff
        t.shape=t.size,1
        nskip=1   #Don't skip any particle
        rho = ms.density.transpose()[:,::nskip]
        px = ms.particle_x.transpose()[:,::nskip]
        py = ms.particle_y.transpose()[:,::nskip]
        pz = ms.particle_z.transpose()[:,::nskip]

        center = [np.mean(px[-1,:]),np.mean(py[-1,:]),np.mean(pz[-1,:])]

        px_new = px - center[0]
        py_new = py - center[1]
        pz_new = pz - center[2]

        radius, theta, phi = CS_Transformation(px_new,py_new,pz_new)
        timelim = [0,t.max()]
        rholim = np.log10([rho.min(),rho.max()])

        cmap = cm.get_cmap('viridis')
        normalizer = Normalize(rholim[0], rholim[1])
        im = cm.ScalarMappable(norm=normalizer)

        nbins = 32
        xlabels = [r'$-\pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', '$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
        ylabels = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
        for frame in mon.frames[::1]:       #all frames
            nf = mon.get_frame_index(frame)
            fig,ax=plt.subplots(1,1)
            if reg_dist == False:
                ax.hist2d(theta[nf,:],phi[nf,:],weights=rho[nf,:], bins = [nbins//2,nbins], cmap=cmap, norm=normalizer)
                ax.set_xlim([-np.pi, np.pi])
                ax.set_ylim([0, np.pi])
                ax.xaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
                ax.xaxis.set_major_locator(MultipleLocator(base=np.pi/2))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.0g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
                ax.yaxis.set_major_locator(MultipleLocator(base=np.pi/4))
            else:
                ax.hist2d(px_new[nf,:],py_new[nf,:],weights=rho[nf,:], bins = [nbins//2,nbins], cmap=cmap, norm=normalizer)
                ax.set_xlim([px_new.min(), px_new.max()])
                ax.set_ylim([py_new.min(), py_new.max()])
            plt.title("Number of Tracers: "+str(NumOfParticles))
            fig.savefig('%s/%sf%04d'%(plot_dir,prefix,nf), papertype = 'a4', dpi = 200)
            plt.close(fig)

def CoreDist_CoreWise(mon,core_id,prefix='DD'):
        ms = mon.get_ms(core_id)
        ms.particle_pos(core_id)
        NumOfParticles = len(ms.particle_id)
        t = mon.all_times/colors.tff
        t.shape=t.size,1
        nskip=1   #Don't skip any particle
        rho = ms.density.transpose()[:,::nskip]
        px = ms.particle_x.transpose()[:,::nskip]
        py = ms.particle_y.transpose()[:,::nskip]
        pz = ms.particle_z.transpose()[:,::nskip]

        center = [np.mean(px[-1,:]),np.mean(py[-1,:]),np.mean(pz[-1,:])]

        px_new = px - center[0]
        py_new = py - center[1]
        pz_new = pz - center[2]

        radius, theta, phi = CS_Transformation(px_new,py_new,pz_new)
        timelim = [0,t.max()]
        rholim = np.log10([rho.min(),rho.max()])

        cmap = cm.get_cmap('viridis')
        normalizer = Normalize(rholim[0], rholim[1])
        im = cm.ScalarMappable(norm=normalizer)
        important_frames = [0,10,20,40,50,60,90,100,117] #9 frames, 3 each for each phase?, collection, hardening, singularity
        fig,ax=plt.subplots(3,3)
        plt.suptitle("Number of Tracers: "+str(NumOfParticles))

        nbins = 32
        xlabels = [r'$-\pi$', r'$-3\pi/4$', r'$-\pi/2$', r'$-\pi/4$', '$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
        xlabels = [r'$-\pi$', r'$-\pi/2$', '$0$', r'$\pi/2$', r'$\pi$']
        ylabels = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']



        for frame_ind in range(len(important_frames)):
            i = frame_ind//3
            j = frame_ind%3
            frame = important_frames[frame_ind]
            nf = mon.get_frame_index(frame)
            ax[i,j].hist2d(theta[nf,:],phi[nf,:],weights=rho[nf,:], bins = [nbins//2,nbins], cmap=cmap, norm=normalizer)
            ax[i,j].set_xlim([-np.pi, np.pi])
            ax[i,j].set_ylim([0, np.pi])
            ax[i,j].set_xticks(np.arange(-(np.pi+0.01),(np.pi+0.01) , np.pi/2))
            ax[i,j].set_xticklabels(xlabels)
            ax[i,j].set_yticks(np.arange(0,(np.pi+0.01) , np.pi/4))
            ax[i,j].set_yticklabels(ylabels)
            if j == 1 or j == 2:
                ax[i,j].get_yaxis().set_visible(False)
            if i == 0 or i == 1:
                ax[i,j].get_xaxis().set_visible(False)
        plt.subplots_adjust(hspace=0.1, wspace=0.1)
        fig.text(0.5, 0.04, 'Theta (Plane wrt z where particles exist)', ha='center', va='center')
        fig.text(0.06, 0.5, 'Phi (Angle inside plane where particles exist)', ha='center', va='center', rotation='vertical')
        fig.colorbar(im, ax=ax.ravel().tolist(), cmap = cmap)
        fig.savefig('%s/%s'%(plot_dir,prefix), papertype = 'a4', dpi = 500)
        plt.close(fig)


if 'things' not in dir():
        things={}

for sim in sim_list:
        mon = monster.closet[sim]
        all_cores=np.unique(TL.loops[sim].tr.core_ids)
        core_list = [114]#list(all_cores) 
        for core_id in core_list:
            ms = mon.get_ms(core_id)
            #XY Plot Framewise
            make_dir(plot_dir+'/'+result_dir+'CartesianDistribution/AllFrames/'+'Core_'+str(core_id)+'/')
            CoreDist_allframes(mon,core_id,prefix=result_dir+'CartesianDistribution/AllFrames/'+'Core_'+str(core_id)+'/DD_', reg_dist = True)
            #Theta-Phi Plot Framewise
            make_dir(plot_dir+'/'+result_dir+'AngularDistribution/AllFrames/'+'Core_'+str(core_id)+'/')
            CoreDist_allframes(mon,core_id,prefix=result_dir+'AngularDistribution/AllFrames/'+'Core_'+str(core_id)+'/DD_')
            #Theta-Phi Plot CoreWise
            make_dir(plot_dir+'/'+result_dir+'AngularDistribution/CoreWise/')
            CoreDist_CoreWise(mon,core_id,prefix=result_dir+'AngularDistribution/CoreWise/'+'Core_'+str(core_id).zfill(3))
