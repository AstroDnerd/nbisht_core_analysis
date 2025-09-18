
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)

file_list=glob.glob('%s/rack_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)
#file_list=glob.glob('/home/dcollins/scratch/Paper19/all_frames/track_all_frames_*.h5')
#file_list = ['./sphere2.h5']
#file_list = ['./u16_sphere2_t2.h5']
#out_prefix = 'u16_sphere2_t2'
#file_list = ['./u17_test.h5']
#out_prefix = 'u17_linear'
file_list = ['mag_work_tmp_2_c0034.h5']
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
all_cores = np.unique(thtr.core_ids)

fig,ax=plt.subplots(1,1)

rm = davetools.rainbow_map(len(all_cores))
asort =  np.argsort(thtr.times)
n0=asort[0]
tsorted = thtr.times[asort]
fig,ax=plt.subplots(1,1)
for nc,core_id in enumerate(all_cores):
    density = thtr.c([core_id],'density')

    for npart in list(range(ms.nparticles))[::10]:
        c = color_map.to_rgba(density[npart,n0])
        ax.plot( tsorted, thtr['mag_work'][npart,asort],c=c,linewidth=.1)#linestyle=':')
