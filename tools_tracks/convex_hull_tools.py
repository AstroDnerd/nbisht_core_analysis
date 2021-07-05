
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}
color.update({'u201':'r','u202':'g','u203':'b'})

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

    simplex = hull.find_simplex(p)

    good = simplex >=0
    return good

class hull_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cell_volumes=[]
        self.hull_volumes=[]
        self.hulls={}
        self.points_3d={}
        self.cores_used=[]
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
        fraction =  in_1_2.sum()/points_1.shape[0]

        if do_plots:
            print(in_1_2.shape)
            fig_many, ax_many = plt.subplots(2,2,figsize=(8,8))
            ax_all = [ax_many[0][0], ax_many[0][1], ax_many[1][0]]
            for ax in ax_all:
                ax.clear()
                ax.set_aspect('equal')
                #ax.plot([0,1,1,0,0],[0,0,1,1,0])
            c1 = nar(['k']*points_1.shape[0])
            c1[in_1_2<=0] = 'r'
            for LOS in [0,1,2]:
                x = [0,2,0][LOS]
                y = [1,1,2][LOS]
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]

                ax_all[LOS].scatter(points_1[:,x], points_1[:,y],s=0.1,c=c1)

            outname="plots_to_sort/overlap_test_%s_c%04d_c%04d.png"%(self.this_looper.out_prefix,core_1,core_2)
            fig_many.savefig(outname)
            print(outname)
            plt.close(fig)
        return fraction

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
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.r.shape[0] <= 10:
                continue
            print("hull on ", core_id)
            self.cores_used.append(core_id)
            asort =  np.argsort(thtr.times)
            delta=0.1

            mask = slice(None)
            for it,frame in enumerate(frames):#asort):
                nt= np.where( nar(thtr.frames) == frame)[0][0]

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

                if do_3d_plots:
                    plt.clf()
                    plt.scatter( this_x, this_y)
                    vert_x = self.points_3d[core_id][hull_3d.vertices,0]
                    vert_y = self.points_3d[core_id][hull_3d.vertices,1]
                    plt.plot(vert_x,vert_y,c='r')
                    outname="plots_to_sort/%s_hull3d_c%04d_n%04d"%(this_looper.out_prefix,core_id,frame)
                    plt.savefig(outname)
                    print(outname)

def plot_2d_full(htool,core_list=None,accumulate=False,frames=[0], all_plots=False):

    thtr = htool.this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    rm = rainbow_map(len(all_cores))
    if frames is None:
        frames = thtr.frames
    if core_list is None:
        core_list = all_cores
    fig_many, ax_many = plt.subplots(2,2,figsize=(8,8))
    ax_all = [ax_many[0][0], ax_many[0][1], ax_many[1][0]]
    ax4 = ax_many[1][1]
    for ncore,core_id in enumerate(core_list):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.r.shape[0] <= 4:
            continue
        delta=0.1

        mask = slice(None)
        for it,frame in enumerate(frames):#asort):
            nt= np.where( nar(thtr.frames) == frame)[0][0]

            if not accumulate:
                for ax in ax_all:
                    ax.clear()
                    ax.set_aspect('equal')
                    ax.plot([0,1,1,0,0],[0,0,1,1,0])
                ax4.clear()
            this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]

            do_hull = True
            if np.unique(this_x).size < 4  or\
               np.unique(this_y).size < 4  or\
               np.unique(this_z).size < 4 :
                print("Not enough degrees of freedom")
                ax4.text(0,0,"Not enought DOF")
                do_hull = False

            this_p = [this_x,this_y,this_z]


            for LOS in [0,1,2]:
                print("WARNING these LOS are wrong")
                x = [0,2,0][LOS]
                y = [1,1,2][LOS]
                xlab='xyz'[x]
                ylab='yzx'[y]



                ax_all[LOS].scatter(this_p[x], this_p[y],s=0.1)

                if do_hull:
                    points_2d = np.array(list(zip(this_p[x],this_p[y])))
                    hull_2d = ConvexHull(points_2d)
                    vert_x = points_2d[hull_2d.vertices,0]
                    vert_y = points_2d[hull_2d.vertices,1]
                    vert_x = np.concatenate([vert_x,vert_x[0:1]])
                    vert_y = np.concatenate([vert_y,vert_y[0:1]])
                    ax_all[LOS].plot(vert_x, vert_y, 'k')

                x_min = min([this_p[x].min(), -delta])
                x_max = max([this_p[x].max(), 1+delta])
                y_min = min([this_p[y].min(), -delta])
                y_max = max([this_p[y].max(), 1+delta])


                axbonk(ax_all[LOS],xlabel=xlab,ylabel=ylab,xlim=[x_min,x_max],ylim=[y_min,y_max])
                #ax_many.set_title(title)
            cumltext=""
            if accumulate:
                cumltext="%04d"%ncore
            if all_plots:
                raise
                outname = '%s/%s_hull_3d_t_%sc%04d_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,cumltext,core_id,frame)
                fig_many.savefig(outname)
                print("Wrote "+outname)
    if accumulate:
        outname = '%s/%s_hull_3d_t_cXXXX_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,frame)
        fig_many.savefig(outname)
        print("Wrote "+outname)


def plot_2d(htool,core_list=None,accumulate=False,frames=[0],all_plots=False, label_cores=[]):

    thtr = htool.this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    rm = rainbow_map(len(all_cores))
    if frames is None:
        frames = thtr.frames
    if core_list is None:
        core_list = all_cores
    fig_many, ax = plt.subplots(1,1,figsize=(8,8))
    x_min, x_max, y_min, y_max = [1,0,1,0]
    for ncore,core_id in enumerate(core_list):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.r.shape[0] <= 4:
            continue
        delta=0.1

        mask = slice(None)
        for it,frame in enumerate(frames):#asort):
            nt= np.where( nar(thtr.frames) == frame)[0][0]

            if not accumulate:
                ax.clear()
                ax.set_aspect('equal')
                ax.plot([0,1,1,0,0],[0,0,1,1,0])
            this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]

            do_hull = True
            if np.unique(this_x).size < 4  or\
               np.unique(this_y).size < 4  or\
               np.unique(this_z).size < 4 :
                print("Not enough degrees of freedom")
                do_hull = False

            this_p = [this_x,this_y,this_z]


            for LOS in [0]:#,1,2]:
                x = [1,0,0][LOS]
                y = [2,2,1][LOS]
                xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
                ylab=r'$%s \rm(code\ length)$'%'xyz'[y]
                this_ax = ax
                this_ax.scatter(this_p[x], this_p[y],s=2)

                if do_hull:
                    points_2d = np.array(list(zip(this_p[x],this_p[y])))
                    hull_2d = ConvexHull(points_2d)
                    vert_x = points_2d[hull_2d.vertices,0]
                    vert_y = points_2d[hull_2d.vertices,1]
                    vert_x = np.concatenate([vert_x,vert_x[0:1]])
                    vert_y = np.concatenate([vert_y,vert_y[0:1]])
                    this_ax.plot(vert_x, vert_y, 'k')

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
        outname = '%s/%s_hull_3d_t_cXXXX_n%04d.png'%(dl.output_directory,htool.this_looper.out_prefix,frames[0])
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
