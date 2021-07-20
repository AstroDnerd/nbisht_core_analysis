from starter2 import *


radius = 1e-2
if 'den' not in dir():
    den = {}
    for this_simname in ['u201','u202','u203']:
        densities = []
        fptr = h5py.File(dl.peak_list[this_simname])
        frame = dl.target_frames[this_simname]
        ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
        try:
            peaks = fptr['peaks'][()]
        except:
            raise
        finally:
            fptr.close()
        print(peaks)
        for peak_id,center in enumerate(peaks):
            this_peak = ds.arr(center, 'code_length')
            print("%s %d/%d"%(this_simname, peak_id, len(peaks)))
            test_main = ds.sphere(this_peak, radius)
            peak_density = get_density(this_peak,test_main)
            densities.append(peak_density)
            if len(peak_density) == 0:
                raise
        den[this_simname]=densities

fig, ax=plt.subplots(1,1)
for nplot,this_simname in enumerate(['u201','u202','u203']):
    density = den[this_simname]
    cdf = np.arange( len(density))/len(density)
    ax.plot(sorted(density),cdf,c='rgb'[nplot])
axbonk(ax,xscale='log',yscale='log',xlabel='rho',ylabel='cdf')
fig.savefig('plots_to_sort/cdf.png')
