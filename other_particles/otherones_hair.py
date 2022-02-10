
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
        
        for core_id in core_list:
            print('Otherones core %d'%core_id)
            ms=  trackage.mini_scrubber(other_looper.tr,core_id, do_velocity=False)
            #ms.make_floats(core_id)

            Ncol = len(frames)
            frame_index = [ np.where( other_looper.tr.frames == frame)[0][0] for frame in frames]
            frame_str = "_%04d"*Ncol%tuple(frames)
            nrows=len(frame_index)
            counter=-1
            fig,axes=plt.subplots(2,nrows, figsize=(12,12))
            for nt,frame in zip(frame_index, frames):
                counter+=1


                ax = axes[counter][0]
                ax1 = axes[counter][1]
                ax.set_aspect('equal')


                delta = 0.1

                mask = slice(None)
                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
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
                #ax.scatter( this_p[x], this_p[y])
                #ax.scatter( ii[x], ii[y])
                #ax.scatter(ii[x].flatten(), ii[y].flatten())
                if 1:
                    TheZ = image.sum(axis=1)
                    norm = mpl.colors.LogNorm(vmin=1,vmax=TheZ.max())
                    cmap=copy.copy(mpl.cm.get_cmap("viridis"))
                    cmap.set_under('w')
                    #ax.pcolormesh(xx,yy,TheZ, norm=norm,cmap=cmap)
                    ax.imshow(TheZ, cmap=cmap, norm=norm, origin='lower',interpolation='nearest')

                box_x, box_y = [0,1,1,0,0], [0,0,1,1,0]
                box_x, box_y = nar(box_x)*128-imin[x], nar(box_y)*128-imin[y]
                ax.plot(box_x,box_y, c=[0.5]*3)

                ms2 = trackage.mini_scrubber(self.first_looper.tr, core_id, do_velocity=False)

                first_p = np.stack([ms2.this_x[:,nt],ms2.this_y[:,nt], ms2.this_z[:,nt]])

                ax.scatter( first_p[x]*128-imin[x], first_p[y]*128-imin[y], s=100, facecolors='none',edgecolors='k')
                #ax.scatter( first_p[x], first_p[y], s=100)

                xticks = np.mgrid[-imin[x]:128-imin[x]:11j]
                labs = ["%0.1f"%val for val in (xticks+imin[x])/128]
                ax.set_xticks(xticks)
                ax.set_xticklabels(labs)

                yticks = np.mgrid[-imin[y]:128-imin[y]:11j]
                labs = ["%0.1f"%val for val in (yticks+imin[y])/128]
                ax.set_yticks(yticks)
                ax.set_yticklabels(labs)
                ax.set_xlabel(r'$%s$'%('xyz'[x]))
                ax.set_ylabel(r'$%s$'%('xyz'[y]))

                other_density = other_looper.tr.c([core_id],'density')[:,nt]
                first_density = self.first_looper.tr.c([core_id],'density')[:,nt]

                other_x = ms.this_x[:,nt]-ms2.mean_x[nt]
                other_y = ms.this_y[:,nt]-ms2.mean_y[nt]
                other_z = ms.this_z[:,nt]-ms2.mean_z[nt]
                other_r = np.sqrt(other_x**2+other_y**2+other_z**2)
                rmin = 1/2048
                rmax=other_r.max()
                other_r[ other_r < rmin] = rmin
                rrr=ms2.r[:,nt]+0
                rrr[rrr<rmin]=rmin
                r_ext = extents()
                d_ext = extents()
                r_ext(other_r); d_ext(other_density)
                r_ext(rrr); d_ext(first_density)

                import pcolormesh_helper as pch
                reload(pch)
                rho_bins = np.geomspace( other_density.min(), other_density.max(), 64)
                r_bins = np.geomspace( *r_ext.minmax)
                hist, xb, yb = np.histogram2d(other_r,other_density, bins=[r_bins,rho_bins])
                pch.helper(hist,xb,yb,ax=ax1, cmap_name='Reds')

                ax1.scatter( rrr, first_density, s=100, facecolors='none',edgecolors='k')
                axbonk(ax1,xscale='log',yscale='log', xlim=r_ext.minmax,ylim=d_ext.minmax, xlabel=r'$r$', ylabel=r'$\rho$')



            outname='plots_to_sort/otherones_%s_c%04d_%s.png'%(output_prefix,core_id,frame_str)
            fig.savefig(outname)
            plt.close(fig)
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
