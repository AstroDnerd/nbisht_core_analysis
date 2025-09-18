
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

core_list = all_cores
fig_many, ax_many = plt.subplots(1,1)

asort =  np.argsort(thtr.times)
n0=asort[0]
tsorted = thtr.times[asort]
#fig,ax=plt.subplots(1,1)
fig,axes=plt.subplots(1,2)
ax=axes[0]; ax1=axes[1]
rext = extents()
vext = extents()
for nc,core_id in enumerate([14]):
    ax.clear()
    ax1.clear()
    ms = trackage.mini_scrubber(thtr,core_id)
    density = thtr.c([core_id],'density')
    #if field in thtr.track_dict.keys():
    #    the_field = thtr.c([core_id],field)
    #elif field == 'radius':
    the_field = ms.r
    tmap=rainbow_map(ms.ntimes)
    norm = mpl.colors.Normalize()
    norm.autoscale( np.log10(density[:,n0]))
    cmap = mpl.cm.jet
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

    for npart in list(range(ms.nparticles))[::10]:
        #c = color_map.to_rgba(np.log10(density[npart,n0]))
        c = color_map.to_rgba(density[npart,n0])
        #ax.plot( tsorted, density[npart,asort],c='k',linestyle=':',marker='*')
        pos = the_field[npart,asort]
        tcen = 0.5*(tsorted[:-1]+tsorted[1:])
        dt = tsorted[1:]-tsorted[:-1]
        dpos = pos[1:]-pos[:-1]
        drdt = dpos/dt
        if drdt.max() * drdt.min() > 0:
            continue

        ok = dt!=0
        ax.plot( tsorted, pos,c=c,linewidth=.1)#linestyle=':')
        ax1.plot( tcen[ok], drdt[ok],c=c,linewidth=.1)#linestyle=':')
        rext(pos)
        vext(drdt[ok])
    axbonk(ax,xscale='linear',yscale='linear',xlabel='t',ylabel=r'$r$')
    ax1.set_yscale('symlog',linthreshy=1)
    ax1.set_ylim(vext)
    ax.set_ylim(rext)
    oname = "%s/%s_%s_time_c%04d"%(dl.output_directory,out_prefix,'radius',core_id)
    fig.savefig(oname)
    print(oname)
#
# Other undocumented models for what r(t) should be.
#
    #theta = np.arange(0,2*np.pi,0.1)
    #G = 1
    #rho0=10
    #test_r = ms.r[:,n0][30]
    #M = 4./3*np.pi*rho0*test_r**2
    #Eabs= G*M/test_r
    #B=G*M/(2*Eabs)**1.5
    #t = 10*B*(theta-np.sin(theta))
    #r0 = test_r# ms.r[:,n0]
    #guess_r = r0*(1-np.sin(theta))
    #ax.plot(t,guess_r)

    #err= np.exp(np.log(density).std(axis=0)[asort])
    #ax.plot(tsorted, density.mean(axis=0)[asort],c='k')
    #ax.errorbar(tsorted, density.mean(axis=0)[asort],c='k',yerr=err)
    #ax.plot(tsorted, density.mean(axis=0),c='k')

    #t0 = thtr.times[asort][0]
    #t1 = thtr.times[asort][-1]
    #rho0 =1.1 #10 # np.mean(density[:,asort[0]])
    #rho1 = density.max() # np.mean(density[:,asort[-1]])
    #alpha = 1.8
    #tc =t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
    #G=1 #np.pi*4#1620./(4*np.pi)
    #tff_global = np.sqrt(3*np.pi/(32*G*1))
    #tff_local = np.sqrt(3*np.pi/(32*G*rho0))
    ##tc=tff_local #kludge
    #rhot = rho0*(1-(tsorted/tc)**2)**-alpha
    #rho_c = 3*np.pi/(32*G*tc**2)

    #ok = np.isnan(rhot)==False
    #ax.plot( tsorted[ok], rhot[ok], c='r',label=r'$tc/tff = %0.2e$'%(tc/tff_local))
    #ax.text( tsorted[0], rho1, r'$tc = %0.2e \rho_c = %0.2e$'%(tc,rho_c))
    #ax.text( tsorted[0], 0.5*rho1, 
    #ax.legend(loc=0)
#           for i,n in enumerate(asort):
#               timeline=plt.plot( [tsorted[i]]*2,[1,1e8],c=[0.5]*3,linewidth=0.1)
#               timetext=plt.text( tsorted[i], 1e8, 'n=%d'%thtr.frames[n])
#               outname = 'image_tracks/rho_t_fit2_c%04d_s%04d.png'%(core_id,i)
#               plt.yscale('log')
#               plt.savefig(outname)
#               timeline[0].remove()
#               timetext.remove()

#               print('saved '+outname)

#               frame = thtr.frames[n]
#               frame = i #you suck.
#               basedir = "/home/dcollins/RESEARCH2/Paper19_47_overlap/0000_density_tracks/x/ALL"
#               core_dir = "%s/c%04d"%(basedir,core_id)
#               image = "%s/test_full_particles_c%04d_n%04d_Projection_x_density.png"%(core_dir,core_id,frame)
#               img1 = mpimg.imread(image) 
#               img2 = mpimg.imread(outname) 
#               both = np.zeros( [ max([img1.shape[0],img2.shape[0]]), img1.shape[1]+img2.shape[1],4])
#               both[ 0:img1.shape[0] , 0:img1.shape[1]]=img1
#               both[ 0:img2.shape[0] , img1.shape[1]:]=img2
#               oname2 = "image_tracks/density_2_c%04d_n%04d"%(core_id,frame)
#               print(oname2)
#               plt.imsave(oname2,both)

