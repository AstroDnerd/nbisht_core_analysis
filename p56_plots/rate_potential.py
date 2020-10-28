from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

plt.close('all')


file_list=glob.glob('/scratch1/dcollins/Paper19/Datasets/primitive_with_potential_n0000/*h5')
output_prefix='potential'
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    
if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

file_list_2=glob.glob('%s/*h5'%dl.sixteen_frame)
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    
if 'this_looper_2' not in dir():
    this_looper_2=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list_2):
        this_looper_2.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr2 = this_looper_2.tr
    thtr2.sort_time()

if 1:
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

    if 'rho_extents' not in dir():
        rho_extents=davetools.extents()
        r_extents=davetools.extents()
        for nc,core_id in enumerate(all_cores):
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles == 1:
                continue
            density = thtr.c([core_id],'density')
            rho_extents(density)
            r_extents(ms.r)

    dx = 1./2048
    nx = 1./dx
    frame_list = thtr.frames
    plt.close('all')

    frame_list=[0]
    if 1:
        fig1,ax1 = plt.subplots(1,1)
        collapse_times=[]
        potential_mean=[]
        potential_std=[]
        for nc,core_id in enumerate(core_list):
            density_all = thtr2.c([core_id],'density')
            if density_all[:,0].size <= 3:
                continue

            if (density_all>1e3).any():
                first_collapse_index = np.where( (density_all>1e3).sum(axis=0) >0 )[0][0] 
                t_collapse = thtr2.times[first_collapse_index]
                rho_col = density_all[:,first_collapse_index].max()
            else:
                t_collapse=-1
                rho_col = -1
            collapse_times.append(t_collapse)

            for nf,frame in enumerate(frame_list):
                x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]#or whatever the number of zones is
                y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
                z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
                density = thtr.c([core_id],'density')[:,nf]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
                index = x + nx*(y * nx*z)
                ar = np.argsort(index)
                rs = np.argsort(ar)
                isorted=index[ar]
                mask = np.ones_like(density,dtype='bool')
                mask[1:] = isorted[1:]-isorted[:-1] != 0
                mask2 = mask[ rs]
                mass = (density[mask2]*cell_volume[mask2]).sum()
                potential = thtr.c([core_id],'PotentialField')[:,nf]

                potential_mean.append(   (potential*cell_volume).sum()/cell_volume.sum())
                potential_std.append(   np.sqrt((potential**2*cell_volume).sum()/cell_volume.sum()))

                vals, bin_edge = np.histogram(potential)
                bin_center = 0.5*(bin_edge[1:]+bin_edge[:-1])
                ax1.plot(bin_center,vals,c=[0.5]*4)

        outname='plots_to_sort/potential_hists_n0000.pdf'
        axbonk(ax1,xlabel=r'$\phi$',ylabel=r'$V(\phi)$',xscale='linear',yscale='log')
        fig1.savefig(outname)
        print(outname)

    fig2,ax2l = plt.subplots(1,2)
    ax2_0 = ax2l[0]
    ax2_1 = ax2l[1]
    fig3,ax3=plt.subplots(1,1)

    collapse_times=nar(collapse_times)
    potential_mean = nar(potential_mean)
    potential_std  = nar(potential_std)
    ok = collapse_times>0
    ax2_0.scatter( collapse_times[ok],potential_mean[ok] )
    ax2_1.scatter( collapse_times[ok],potential_std[ok] )

    ax3.errorbar(  collapse_times[ok],potential_mean[ok], yerr=potential_std[ok],fmt='o')
    axbonk(ax3,xlabel=r'$t_c$', ylabel=r'$\phi$')
    fig3.savefig('plots_to_sort/potential_err.png')


    outname='plots_to_sort/potential_tc_n0000.pdf'
    axbonk(ax2_0,xlabel=r'$t_c$', ylabel=r'$\langle \phi \rangle$')
    axbonk(ax2_1,xlabel=r'$t_c$', ylabel=r'$\sqrt{\langle \phi^2 \rangle}$')
    fig2.savefig(outname)
    print(outname)
