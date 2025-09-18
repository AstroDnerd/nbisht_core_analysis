
from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='test'
file_list = file_list[:2]

import tools_tracks.alpha_properties as ap 
this_looper=ap.this_looper
thtr=this_looper.tr
thtr.sort_time()

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))

core_list = [14] #all_cores
fig_many, ax_many = plt.subplots(1,1)
for core_id in core_list:
    ms = trackage.mini_scrubber(thtr,core_id)
    tmap=rainbow_map(ms.ntimes)
    asort =  np.argsort(thtr.times)
    delta=0.1

    rrr = ms.r[:,asort]
    dv = rrr[:,1:] - rrr[:,:-1]
    mask = dv[:,1]*dv[:,2] < 0
    for it,nt in enumerate(asort):
        if it < 2:
            continue
        nt0 = asort[it-1]
        ntm2 = asort[it-2]
        d2v = (rrr[:,nt] - rrr[:,nt0])*(rrr[:,nt0]-rrr[:,ntm2])
        mask = d2v<0

        print("NNNN",(d2v < 0).sum())


        ax_many.clear()
        ax_many.plot([0,1,1,0,0],[0,0,1,1,0])
        this_frame = sorted(thtr.frames)[nt]
        cen=this_looper.snaps[this_frame][core_id].R_centroid
        cenx = cen[0]; ceny=cen[1]
        ax_many.scatter(ms.this_x[mask,nt],ms.this_y[mask,nt],s=0.1)
        ax_many.scatter( cenx,ceny,c='k',s=3)
        ax_many.set_xlim(-delta,1+delta); ax_many.set_ylim(-delta,1+delta)
        title='t=%0.5f'%thtr.times[nt]
        ax_many.set_title(title)
        ax_many.set_aspect('equal')
        outname = '%s/xy_t_c%04d_n%04d.png'%(dl.output_directory,core_id,it)
        fig_many.savefig(outname)
        print("Wrote "+outname)

