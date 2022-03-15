
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})
import time
from scipy.spatial import Delaunay
def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
	from https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
    also see this
    https://stackoverflow.com/questions/64310174/in-scipy-spatial-delaunay-what-does-find-simplex-method-return
    """
    #from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)


    #t0 = time.time()
    simplex = hull.find_simplex(p)
    #t1 = time.time()
    #print('woot %f'%(t1-t0))

    good = simplex >=0
    return good

class hull_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.name = this_looper.sim_name
        self.cell_volumes=[]
        self.unique_volumes=[]
        self.hull_volumes=[]
        self.hulls={}
        self.points_3d={}
        self.cores_used=[]
        self.nparticles=[]
    def make_overlaps(self):
        self.overlaps=defaultdict(list)
        self.overlap_number=defaultdict(list)
        #self.make_hulls(frames=[frame])
        for core_1 in self.cores_used:
            print("overlap li,", self.name, core_1)
            for core_2 in self.cores_used:
                fraction,number = self.check_hull_overlap(core_1,core_2)
                mult=1
                if core_1 == core_2:
                    mult=-1
                self.overlaps[core_1].append(mult*fraction)
                self.overlap_number[core_1].append(mult*number)
    def check_hull_overlap(self,core_1, core_2, do_plots=False):
        hull_1 =  self.hulls[core_1]
        hull_2 =  self.hulls[core_2]
        if hull_1 is None or hull_2 is None:
            return 0
        vert_1 = self.points_3d[core_1][hull_1.vertices,:]
        vert_2 = self.points_3d[core_2][hull_2.vertices,:]
        points_1 = self.points_3d[core_1]
        points_2 = self.points_3d[core_2]

        in_1_2 = in_hull(points_1, vert_2)
        number =  in_1_2.sum()
        fraction =  number/points_1.shape[0]

        return fraction, number

    def make_hulls(self,do_3d_plots=False,core_list=None,frames=[0]):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))
        if frames is None:
            frames = thtr.frames

        if core_list is None:
            core_list = all_cores
        if core_list == "short":
            core_list = all_cores[:10]
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            if ms.r.shape[0] <= 10:
                continue
            print("hull on %s n%04d c%04d"%(self.this_looper.sim_name, frames[0], core_id))
            self.cores_used.append(core_id)
            self.nparticles.append( ms.nparticles)
            asort =  np.argsort(thtr.times)
            delta=0.1

            mask = slice(None)
            for it,frame in enumerate(frames):#asort):
                nt= np.where( nar(thtr.frames) == frame)[0][0]
                mask2 = ms.compute_unique_mask(core_id, dx=1./2048,frame=nt)

                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
                this_p = [this_x,this_y,this_z]

                self.points_3d[core_id] = np.array(list(zip(this_x,this_y,this_z)))

                N_unique = min([np.unique(arr).size for arr in (this_x,this_y,this_z)])
                if N_unique > 2:
                    hull_3d = ConvexHull(self.points_3d[core_id])
                    self.hulls[core_id]=hull_3d
                    self.hull_volumes.append(hull_3d.volume)
                else:
                    self.hulls[core_id] = None
                    self.hull_volumes.append(0)
                self.cell_volumes.append( thtr.c([core_id],'cell_volume')[mask,nt].sum())

                self.unique_volumes.append( thtr.c([core_id],'cell_volume')[mask2,nt].sum())

                if do_3d_plots:
                    plt.clf()
                    plt.scatter( this_x, this_y)
                    vert_x = self.points_3d[core_id][hull_3d.vertices,0]
                    vert_y = self.points_3d[core_id][hull_3d.vertices,1]
                    plt.plot(vert_x,vert_y,c='r')
                    outname="plots_to_sort/%s_hull3d_c%04d_n%04d"%(this_looper.out_prefix,core_id,frame)
                    plt.savefig(outname)
                    print(outname)

class rainbow_trout():
    def __init__(self,n=None, vmin=None,vmax=None, cmap='jet'):
        norm = mpl.colors.Normalize()
        if vmin is not None and vmax is not None:
            norm.autoscale([vmin,vmax])
        else:
            norm.autoscale(np.arange(n))
        #cmap = mpl.cm.jet
        self.color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    def __call__(self,val,n_fields=0):
        this_value = self.color_map.to_rgba(val)
        if n_fields > 0:
            this_value = [this_value]*n_fields
        return this_value

def plot_watershed(htool,core_list=None,accumulate=False,frames=[0],all_plots=False, label_cores=[],prefix="",
            color_dict=None,axis_to_plot=[0]):

    thtr = htool.this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    rm = rainbow_map(len(all_cores))
    if frames is None:
        frames = thtr.frames
    if core_list is None:
        core_list = all_cores

    if -1 in axis_to_plot:
        fig_many, ax = plt.subplots(2,2,figsize=(8,8))
        ax=ax.flatten()
    else:
        fig_many, ax = plt.subplots(1,1,figsize=(8,8))
    x_min, x_max, y_min, y_max = [1,0,1,0]
    import colors
    if color_dict is  None:
        color_dict = colors.make_core_cmap( core_list)

    ms_dict={}
    for ncore,core_id in enumerate(core_list):
        ms = trackage.mini_scrubber(thtr,core_id)
        ms.particle_pos(core_id)
        ms_dict[core_id]=ms

    vx1 = thtr.c(core_list[:1], 'velocity_x')[:,0]
    vy1 = thtr.c(core_list[:1], 'velocity_y')[:,0]
    vz1 = thtr.c(core_list[:1], 'velocity_z')[:,0]
    vx2 = thtr.c(core_list[1:], 'velocity_x')[:,0]
    vy2 = thtr.c(core_list[1:], 'velocity_y')[:,0]
    vz2 = thtr.c(core_list[1:], 'velocity_z')[:,0]
    gx1 = thtr.c(core_list[:1], 'grav_x')[:,0]
    gy1 = thtr.c(core_list[:1], 'grav_y')[:,0]
    gz1 = thtr.c(core_list[:1], 'grav_z')[:,0]
    gx2 = thtr.c(core_list[1:], 'grav_x')[:,0]
    gy2 = thtr.c(core_list[1:], 'grav_y')[:,0]
    gz2 = thtr.c(core_list[1:], 'grav_z')[:,0]
    dvx1 = thtr.c(core_list[:1], 'velocity_x')[:,2] -thtr.c(core_list[:1], 'velocity_x')[:,0]
    dvy1 = thtr.c(core_list[:1], 'velocity_y')[:,2] -thtr.c(core_list[:1], 'velocity_y')[:,0]
    dvz1 = thtr.c(core_list[:1], 'velocity_z')[:,2] -thtr.c(core_list[:1], 'velocity_z')[:,0]
    dvx2 = thtr.c(core_list[1:], 'velocity_x')[:,2] -thtr.c(core_list[1:], 'velocity_x')[:,0]
    dvy2 = thtr.c(core_list[1:], 'velocity_y')[:,2] -thtr.c(core_list[1:], 'velocity_y')[:,0]
    dvz2 = thtr.c(core_list[1:], 'velocity_z')[:,2] -thtr.c(core_list[1:], 'velocity_z')[:,0]
    if 0:
        VM1 = np.sqrt(vx1**2+vy1**2+vz1**2)
        VM2 = np.sqrt(vx2**2+vy2**2+vz2**2)
        vx1bar=np.mean(vx1)
        vy1bar=np.mean(vy1)
        vz1bar=np.mean(vz1)
        vmbar = np.sqrt(vx1bar**2+vy1bar**2+vz1bar**2)
        cosine1 =(vx1*vx1bar+vy1*vy1bar+vz1*vz1bar)/(VM1*vmbar)
        cosine2 =(vx2*vx1bar+vy2*vy1bar+vz2*vz1bar)/(VM2*vmbar)
    if 0:
        DVM1 = np.sqrt(dvx1**2+dvy1**2+dvz1**2)
        DVM2 = np.sqrt(dvx2**2+dvy2**2+dvz2**2)
        dvx1bar=np.mean(dvx1)
        dvy1bar=np.mean(dvy1)
        dvz1bar=np.mean(dvz1)
        dvmbar = np.sqrt(dvx1bar**2+dvy1bar**2+dvz1bar**2)
        cosine1 =(dvx1*dvx1bar+dvy1*dvy1bar+dvz1*dvz1bar)/(DVM1*dvmbar)
        cosine2 =(dvx2*dvx1bar+dvy2*dvy1bar+dvz2*dvz1bar)/(DVM2*dvmbar)
    if 1:
        GM1 = np.sqrt(gx1**2+gy1**2+gz1**2)
        GM2 = np.sqrt(gx2**2+gy2**2+gz2**2)
        gx1bar=np.mean(gx1)
        gy1bar=np.mean(gy1)
        gz1bar=np.mean(gz1)
        g1bar = np.sqrt(gx1bar**2+gy1bar**2+gz1bar**2)
        gx2bar=np.mean(gx2)
        gy2bar=np.mean(gy2)
        gz2bar=np.mean(gz2)
        g2bar = np.sqrt(gx2bar**2+gy2bar**2+gz2bar**2)
        cosine1 =(gx1*gx1bar+gy1*gy1bar+gz1*gz1bar)/(GM1*g1bar)
        cosine2 =(gx2*gx2bar+gy2*gy2bar+gz2*gz2bar)/(GM2*g2bar)
    colors=[cosine1,cosine2]
    rtmap = rainbow_trout(vmin=-1,vmax=1)
    for it,frame in enumerate(frames):#asort):
        nt= np.where( nar(thtr.frames) == frame)[0][0]
        if -1 in axis_to_plot:
            for aaa in ax:
                aaa.clear()
        else:
                ax.clear()
        ax[-1].hist(cosine1, label='c%04d'%core_list[0],histtype='step')
        ax[-1].hist(cosine2, label='c%04d'%core_list[1],histtype='step')
        ax[-1].legend(loc=0)
        ax[-1].set_yscale('log')
        for ncore,core_id in enumerate(core_list):
            ms = ms_dict[core_id]
            if ms.r.shape[0] <= 4:
                continue
            print('plot core %s %d'%(htool.this_looper.out_prefix,core_id))
            delta=0.1

            mask = slice(None)

            if not accumulate:
                if -1 in axis_to_plot:
                    for aaa in ax:
                        aaa.clear()
                        aaa.set_aspect('equal')
                        aaa.plot([0,1,1,0,0],[0,0,1,1,0])
                else:
                    ax.clear()
                    ax.set_aspect('equal')
                    ax.plot([0,1,1,0,0],[0,0,1,1,0])
            #this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
            this_x,this_y,this_z=ms.particle_x[mask,nt],ms.particle_y[mask,nt], ms.particle_z[mask,nt]

            all_x,all_y,all_z=ms.particle_x,ms.particle_y, ms.particle_z
            all_p = [all_x,all_y,all_z]

            do_hull = True

            if np.unique(this_x).size < 4  or\
               np.unique(this_y).size < 4  or\
               np.unique(this_z).size < 4 :
                print("Not enough degrees of freedom")
                do_hull = False

            this_p = [this_x,this_y,this_z]


            if -1 in axis_to_plot:
                axis_to_actually_plot = [0,1,2]
            else:
                axis_to_actually_plot=axis_to_plot
            for LOS in axis_to_actually_plot:
                x = [1,0,0][LOS]
                y = [2,2,1][LOS]
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]


                if -1 in axis_to_plot:
                    this_ax = ax[LOS]
                else:
                    this_ax = ax
                n_particles = len(this_p[0])

                #
                # the plot
                #
                this_ax.scatter(this_p[x], this_p[y],s=2, c=colors[ncore])

                #streaks.  Use with caution.
                #this_ax.plot(all_p[x].transpose(), all_p[y].transpose(), c=color_dict[core_id], linewidth=.1)


                if do_hull:
                    points_2d = np.array(list(zip(this_p[x],this_p[y])))
                    hull_2d = ConvexHull(points_2d)
                    vert_x = points_2d[hull_2d.vertices,0]
                    vert_y = points_2d[hull_2d.vertices,1]
                    vert_x = np.concatenate([vert_x,vert_x[0:1]])
                    vert_y = np.concatenate([vert_y,vert_y[0:1]])
                    this_ax.plot(vert_x, vert_y, 'k', linewidth=0.3)

                if core_id in label_cores or -1 in label_cores:
                    this_ax.text( this_p[x].max(), this_p[y].max(), r'$%s$'%core_id)
                x_min = min([x_min,this_p[x].min(), -delta])
                x_max = max([x_max,this_p[x].max(), 1+delta])
                y_min = min([y_min,this_p[y].min(), -delta])
                y_max = max([y_max,this_p[y].max(), 1+delta])

                this_ax.plot([0,1,1,0,0], [0,0,1,1,0], c=[0.5]*3)


                axbonk(this_ax,xlabel=xlab,ylabel=ylab,xlim=[x_min,x_max],ylim=[y_min,y_max])
                #ax_many.set_title(title)
            cumltext=""
            if accumulate:
                cumltext="%04d"%ncore
            if all_plots:
                outname = '%s/%s_hull_3d_t_%sc%04d_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,cumltext,core_id,frame)
                fig_many.savefig(outname)
                print("Wrote "+outname)
        if accumulate:
            outname = '%s/%s_hull_3d_t_%scXXXX_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,prefix,frame)
            fig_many.savefig(outname)
            print("Wrote "+outname)




def image_overlap(self,core_1, core_2, do_plots=False):
    hull_1 =  self.hulls[core_1]
    hull_2 =  self.hulls[core_2]
    vert_1 = self.points_3d[core_1][hull_1.vertices,:]
    vert_2 = self.points_3d[core_2][hull_2.vertices,:]
    points_1 = self.points_3d[core_1]
    points_2 = self.points_3d[core_2]

    in_1_2 = in_hull(points_1, vert_2)
    in_2_1 = in_hull(points_2, vert_1)
    fraction =  in_1_2.sum()/points_1.shape[0]

    print(in_1_2.shape)
    fig_many, ax_many = plt.subplots(2,2,figsize=(8,8))
    ax_all = [ax_many[0][0], ax_many[0][1], ax_many[1][0]]
    for ax in ax_all:
        ax.clear()
        ax.set_aspect('equal')
        #ax.plot([0,1,1,0,0],[0,0,1,1,0])
    c1 = nar(['b']*points_1.shape[0])
    c2 = nar(['g']*points_2.shape[0])
    c1[in_1_2<=0] = 'r'
    c2[in_2_1<=0] = 'm'
    for LOS in [0,1,2]:
        x = [0,2,0][LOS]
        y = [1,1,2][LOS]
        xlab='xyz'[x]
        ylab='yzx'[y]

        dx = 1./128/2
        scatter_x = np.random.random(points_1.shape[0])*dx
        scatter_y = np.random.random(points_1.shape[0])*dx
        ax_all[LOS].scatter(points_1[:,x]+scatter_x, points_1[:,y]+scatter_y,s=0.1,c=c1)
        scatter_x = np.random.random(points_2.shape[0])*dx
        scatter_y = np.random.random(points_2.shape[0])*dx
        ax_all[LOS].scatter(points_2[:,x]+scatter_x, points_2[:,y]+scatter_y,s=0.1,c=c2)

    outname="plots_to_sort/overlap_test_%s_c%04d_c%04d.pdf"%(self.this_looper.out_prefix,core_1,core_2)
    fig_many.savefig(outname)
    print(outname)
    plt.close(fig)
    return fraction

def get_overlapping_cores(self,core_id):
    with_overlap = nar(self.overlaps[core_id]) > 0
    used_with = nar(self.cores_used)[with_overlap]
    overlap_with = nar(self.overlaps[core_id])[with_overlap]
    argsort = np.argsort(overlap_with)
    return overlap_with[argsort], used_with[argsort]
