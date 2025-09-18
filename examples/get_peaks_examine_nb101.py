from starter2 import *
from yt.data_objects.level_sets.clump_handling import \
            Clump, \
            find_clumps
from collections import defaultdict

def get_density(peak, test_main):
    X = test_main['x'].to('code_length')
    Y = test_main['y'].to('code_length')
    Z = test_main['z'].to('code_length')
    tolerance = 0.5/2048
    mask = ( np.abs(X - peak[0])<tolerance)*( np.abs(Y-peak[1])<tolerance)*(np.abs(Z-peak[2])<tolerance)
    if mask.sum() != 1:
        pdb.set_trace()
    return test_main['density'][mask]
sim_list = ['nb101']

if 0:
    #
    # Image all peaks and nearby particles.
    #
    density_list=defaultdict(list)
    nzones_list=defaultdict(list)
    peak_id_list=defaultdict(list)
    radius_list=defaultdict(list)
    nparticles_list = defaultdict(list)
    for this_simname in  sim_list:
        frame = dl.target_frames[this_simname]
        ds_name="%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame)
        ds = yt.load(ds_name)
        radius = ds.arr(1e-2,'code_length')
        peak_name =dl.peak_list[this_simname]
        fptr = h5py.File(peak_name,'r')
        print("OPENING %s and %s"%(ds_name, peak_name))
        try:
            peaks = fptr['peaks'][()]
        except:
            raise
        finally:
            fptr.close()
        #peaks = nar([ np.array([0.89331055, 0.1159668 , 0.4440918 ])])
        for peak_id,center_o in enumerate(peaks):
            center = ds.arr(center_o, 'code_length')
            test_main = ds.sphere(center, radius)
            peak_density = get_density(center,test_main)
            proj = ds.proj('density',1,center=center,data_source=test_main)
            pw = proj.to_pw()
            pw.set_cmap('density','Greys')

            #pw.annotate_clumps([master_clump]+master_clump.leaves)
            #pw.annotate_these_particles2(1.0,col='r',positions= test_main['particle_position'])
            pw.zoom(0.5/radius.v)
            pw.set_axes_unit('code_length')
            mountain_top = Clump(test_main,('gas','density'))

            min_density = 1
            step = 1.2
            nzones=(mountain_top['density']>min_density).sum()
            #find_clumps(mountain_top, min_density,peak_density,step)
            min_density = peak_density**(3./4)
            print("DENSITIES peak %0.2e min %0.2e"%(peak_density, min_density))
            mountain_top.find_children( min_density, peak_density)
            for leaf in mountain_top.leaves:
                min_radius=leaf['radius'].min()
                if min_radius < 1/2048:
                    break
            max_radius = leaf['radius'].max()

            pw.annotate_clumps([leaf])
            pw.annotate_title("Peak id %d N particles %d in leaf %d rho %0.2e"%(peak_id,test_main['particle_index'].size, leaf['particle_index'].size, peak_density ))
            pw.save('~/cdbreak_desktop/nikhilb_home/results/plots_to_sort/%s_peak_p%04d_%04d'%(this_simname,peak_id,frame))
            #            proj = ds.proj(field,ax,center=center, data_source = sph) 

            density_list[this_simname].append(peak_density)
            nzones_list[this_simname].append( leaf['radius'].size)
            peak_id_list[this_simname].append(peak_id)
            radius_list[this_simname].append(max_radius)
            nparticles_list[this_simname].append( leaf['particle_index'].size)



if 1:
    density_list=defaultdict(list)
    nzones_list=defaultdict(list)
    peak_id_list=defaultdict(list)
    radius_list=defaultdict(list)
    nparticles_list = defaultdict(list)
    for this_simname in sim_list:
        frame = dl.target_frames[this_simname]
        fig,ax=plt.subplots(2,2)
        ax0 = ax[0][0]; ax1=ax[0][1]
        ax2 = ax[1][0]
        ax0.scatter( density_list[this_simname], nzones_list[this_simname])
        axbonk(ax0,xscale='log',yscale='log',xlabel='Peak Density', ylabel = 'Nzones')
        ax1.scatter( density_list[this_simname], radius_list[this_simname])
        axbonk(ax1,xscale='log',yscale='log',xlabel='peak density', ylabel = 'radius')
        ax2.scatter( radius_list[this_simname], nzones_list[this_simname])
        axbonk(ax2,xscale='log',yscale='log', xlabel='Radius', ylabel='Nzones')

        fig.savefig('~/cdbreak_desktop/nikhilb_home/results/plots_to_sort/%s_density_nzones_mountain_top_%04d.png'%(this_simname,frame))


