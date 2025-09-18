
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)

file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
#file_list=glob.glob('/home/dcollins/scratch/Paper19/all_frames/track_all_frames_*.h5')
#file_list = ['./sphere2.h5']
#file_list = ['./u16_sphere2_t2.h5']
#out_prefix = 'u16_sphere2_t2'
#file_list = ['./u17_test.h5']
#out_prefix = 'u17_linear'
if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):

        trw.load_loop(this_looper,fname)
        print(nfile/len(file_list))
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

asort =  np.argsort(thtr.times)
n0=asort[0]
nmax=asort[-1]
frame=nmax
tsorted = thtr.times[asort]
thtr = this_looper.tr
stuff={}
for nf, ni in enumerate(asort):
    print(nf)
    rm = rainbow_map(len(all_cores))
    
    core_list = np.unique(thtr.core_ids)
    frame = thtr.frames[ni]

    if 1:
        #time plots

        for nc,core_id in enumerate(core_list[::5]):
            plt.close('all')
            fig,ax=plt.subplots(1,1)
            ax.step(prof0['the_x'],prof0['the_y']*prof0['the_x'])
            ax.step(prof1['the_x'],prof1['the_y']*prof1['the_x'])
            ax.set_title("core %d nparticles %d t=%0.3f n=%d"%(core_id, density[:,ni].size, tsorted[nf], frame))
            #ms = trackage.mini_scrubber(thtr,core_id)
            #tmap=rainbow_map(ms.ntimes)
            density = thtr.c([core_id],'density')
            cell_volume = thtr.c([core_id],'cell_volume')
            vals,bins = np.histogram(density[:,ni],weights=cell_volume[:,ni]*density[:,ni])
            #vals,bins = np.histogram(density[:,ni],weights=cell_volume[:,ni])
            stuff[frame]={'vals':vals,'bins':bins}
            bc = 0.5*(bins[1:]+bins[:-1])
            ax.step(bc,vals,c=[0.5]*4)
            #norm = mpl.colors.Normalize()
            #norm.autoscale( np.log10(density[:,n0]))
            #cmap = mpl.cm.jet
            #color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
        #vals,bins = np.histogram(density.flatten(),weights=cell_volume.flatten())
        #ax.step(bc,vals,c='k')

            #axbonk(ax,xlabel=r'$\rho$',ylabel=r'$P(\rho)$',xscale='log',yscale='log')
            axbonk(ax,xlabel=r'$\rho$',ylabel=r'$M(\rho)$',xscale='log',yscale='log')
            ax.set_xlim(5e-3,1e8)
            ax.set_ylim(1e-9,0.2)
            fig.savefig('/home/dcollins4096/PigPen/PDFs_mass_c%04d_n%04d.png'%(core_id,nf))
            #fig.savefig('/home/dcollins4096/PigPen/PDFs_allcores_n%04d.pdf'%nf)
