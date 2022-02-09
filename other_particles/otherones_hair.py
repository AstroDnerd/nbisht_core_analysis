
from starter2 import *
reload(trackage)
import convex_hull_tools as CHT
import matplotlib.colors as colors

reload(CHT)
import hair_dryer
reload(hair_dryer)
import stay_close
import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']
#sim_list=['u402']
import supersets
reload(supersets)
import hair_dryer
reload(hair_dryer)
import colors



class images():
    def __init__(self,other_looper,first_looper):
        self.other_looper=other_looper
        self.first_looper=first_looper
        self.cores_used=[]
        self.name = self.other_looper.sim_name

    def run(self, output_prefix='NAME',core_list=None, frames=None,
            verbose=False, los_list=[-1], external_axis=None):
        print('w3',core_list)
        dx=1./2048
        nx = 1./dx
        thtr = self.other_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')
        if frames is None:
            frames = thtr.frames

        thtr.sort_time()

        tsorted = thtr.times/colors.tff
        self.core_list=core_list
        mini_scrubbers = {}
        other_looper=self.other_looper
        if verbose:
            print("MS 1")
        print('w4')
        print(core_list)
        
        for core_id in core_list:
            ms=  trackage.mini_scrubber(other_looper.tr,core_id, do_velocity=False)
            #ms.make_floats(core_id)

            for nt,frame in [[0,0]]:#[[0,0],[-1,118]]: #kludge needs to be updates to work on frame only
                fig,axes=plt.subplots(2,2, figsize=(12,12))
                ax=axes[0][0]; ax1=axes[0][1]
                ax2=axes[1][0];ax3=axes[1][1]
                ax.set_aspect('equal')

                delta = 0.1
                ax.plot([0,1,1,0,0], [0,0,1,1,0], c=[0.5]*3)

                mask=slice(None)

                if 1:
                    #kludge for bad shift
                    too_low_mask = ms.this_z[:,nt] < -0.3
                    core_mask = ms.trk.core_ids == core_id
                    core_particle_ids = ms.trk.particle_ids[core_mask]
                    wrong_particles = core_particle_ids[ too_low_mask]

                    FirstWrong = np.where(too_low_mask)[0][0]
                    WP = wrong_particles[0]
                    print('WP',WP)
                    d1=other_looper.tr.p([WP],'density')
                    d2=self.first_looper.big_loop.tr.p([WP],'density')
                    print(d1)
                    print(d2)
                    s1=other_looper.tr.p([WP],'shift_z')
                    s2=self.first_looper.big_loop.tr.p([WP],'shift_z')
                    z1=other_looper.tr.p([WP],'z')
                    z2=self.first_looper.big_loop.tr.p([WP],'z')
                    r1=other_looper.tr.p([WP],'test_rho')
                    r2=self.first_looper.big_loop.tr.p([WP],'test_rho')





                    #print("Did I get the particles?", np.abs(core_particle_ids - ms.particle_ids).sum())

                    #print( ms.shift_z[too_low_mask,:])
                    pdb.set_trace()
                    #for npart, part in enumerate(wrong_particles):
                    #    print('shiftj',part,other_looper.tr.p([part],'shift_x'))
                    #    print('shiftk',part,self.first_looper.big_loop.tr.p([part],'shift_x'))


                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
                print("THIS", this_x.shape)
                this_p = [this_x,this_y,this_z]



                x=2;y=0

                x_min, x_max, y_min, y_max = [1,0,1,0]

                x_min = min([x_min,this_p[x].min(), -delta])
                x_max = max([x_max,this_p[x].max(), 1+delta])
                y_min = min([y_min,this_p[y].min(), -delta])
                y_max = max([y_max,this_p[y].max(), 1+delta])

                xx,yy=np.mgrid[0:1:1./128,0:1:1./128]+0.5*1./128
                ii = np.floor((nar(this_p)*128)).astype('int')
                imin = ii.min(axis=1)
                ii[0]-=imin[0]
                ii[1]-=imin[1]
                ii[2]-=imin[2]
                image = np.zeros([ii.max()+1]*3)
                image[ii[0], ii[1], ii[2]]=1
                ax.scatter( this_p[x], this_p[y])
                #ax.scatter( ii[x], ii[y])
                #ax.scatter(ii[x].flatten(), ii[y].flatten())
                if 0:
                    TheZ = image.sum(axis=1)
                    norm = mpl.colors.LogNorm(vmin=1,vmax=TheZ.max())
                    cmap=copy.copy(mpl.cm.get_cmap("viridis"))
                    cmap.set_under('w')
                    #ax.pcolormesh(xx,yy,TheZ, norm=norm,cmap=cmap)
                    ax.imshow(TheZ, cmap=cmap, norm=norm, origin='lower',interpolation='nearest')


                ms2 = trackage.mini_scrubber(self.first_looper.tr, core_id, do_velocity=False)
                first_x,first_y,first_z=ms2.this_x[:,nt],ms2.this_y[:,nt], ms2.this_z[:,nt]
                first_p = np.stack([first_x,first_y,first_z]) #*128

                print("WWW", first_p[x].shape)
                #ax.scatter( first_p[x]-imin[x], first_p[y]-imin[y], s=100)
                ax.scatter( first_p[x], first_p[y], s=100)


                other_density = other_looper.tr.c([core_id],'density')[:,nt]
                first_density = self.first_looper.tr.c([core_id],'density')[:,nt]

                other_x = ms.this_x[:,nt]-ms2.mean_x[nt]
                other_y = ms.this_y[:,nt]-ms2.mean_y[nt]
                other_z = ms.this_z[:,nt]-ms2.mean_z[nt]
                other_r = np.sqrt(other_x**2+other_y**2+other_z**2)
                r_ext = extents()
                d_ext = extents()
                r_ext(other_r); d_ext(other_density)
                r_ext(ms2.r[:,nt]); d_ext(first_density)
                ax1.scatter( other_r,other_density, c='r')
                ax3.scatter( ms2.r[:,nt], first_density, c='g')
                axbonk(ax1,xscale='log',yscale='log', xlim=r_ext.minmax,ylim=d_ext.minmax)
                axbonk(ax3,xscale='log',yscale='log', xlim=r_ext.minmax,ylim=d_ext.minmax)



                outname='plots_to_sort/image_%s_c%04d_%04d.png'%(output_prefix,core_id,frame)
                fig.savefig(outname)
                print(outname)


if 0:
    print('w1')
    if 'loopb' not in dir():
        savefile="otherones_b002_temp.h5"
        loopb = looper.core_looper(directory= TL4.loops['u402'].directory,savefile_only_trackage=savefile)
    print('w1')
#hd = hair_dryer.hair_time(loopb)
#color_dict={1:'r'}
#fig,ax=plt.subplots(1,1)
#hd.run(color_dict=color_dict,external_axis=ax)
#fig.savefig('plots_to_sort/hair2.png')

    IM = images(loop, TL4.loops['u402'])
    IM.run(frames=[0,118])
