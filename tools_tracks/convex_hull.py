
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')
color={'u05':'r','u10':'g','u11':'b'}

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
                xlab='xyz'[x]
                ylab='yzx'[y]

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
            if ms.r.shape[0] <= 4:
                continue
            print("hull on ", core_id)
            self.cores_used.append(core_id)
            asort =  np.argsort(thtr.times)
            delta=0.1

            mask = slice(None)
            for it,nt in enumerate(frames):#asort):

                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
                this_p = [this_x,this_y,this_z]

                self.points_3d[core_id] = np.array(list(zip(this_x,this_y,this_z)))
                hull_3d = ConvexHull(self.points_3d[core_id])
                self.hulls[core_id]=hull_3d
                self.hull_volumes.append(hull_3d.volume)
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


    def plot_2d(self,core_list=None,accumulate=False,frames=[0]):

        thtr = self.this_looper.tr
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
                outname = '%s/%s_hull_3d_t_%sc%04d_n%04d.png'%(dl.output_directory,self.this_looper.out_prefix,cumltext,core_id,frame)
                fig_many.savefig(outname)
                print("Wrote "+outname)

import three_loopers as tl
if 'ht1' not in dir() or clobber: 
    ht1 = hull_tool(tl.looper1) 
    tl.looper1.frame_list=[0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 125]

#    ht1.make_hulls()
#    ht1.plot_2d(frames=None)#,core_list=[10,11])
if 'ht2' not in dir() or clobber:
    ht2 = hull_tool(tl.looper2)
#    ht2.plot_2d(frames=None)
if 'ht3' not in dir() or clobber:
    ht3 = hull_tool(tl.looper3)
#    ht3.plot_2d(frames=None)

#make a lot of plots
#for htool in [ht1, ht2, ht3]:
#    htool.plot_2d(frames=None)

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

#image_overlap(ht1,10,11)

if len(ht1.cores_used) == 0:
    for htool in [ht1, ht2, ht3]:
        htool.overlaps=defaultdict(list)
        htool.make_hulls()
        for core_1 in htool.cores_used:
            print("overlap li,", core_1)
            for core_2 in htool.cores_used:
                result = htool.check_hull_overlap(core_1,core_2)
                if core_1 == core_2:
                    result = -result
                htool.overlaps[core_1].append(result)
def get_overlapping_cores(self,core_id):
    with_overlap = nar(self.overlaps[core_id]) > 0
    used_with = nar(self.cores_used)[with_overlap]
    overlap_with = nar(self.overlaps[core_id])[with_overlap]
    argsort = np.argsort(overlap_with)
    return overlap_with[argsort], used_with[argsort]

fractions,cores=get_overlapping_cores(ht3,185)
catman = np.concatenate
cores = catman([cores,[185]])[::-1]
ht3b = hull_tool(tl.looper3)
ht3b.plot_2d(core_list = cores, accumulate=True)


if 0:
    fig3,ax3=plt.subplots(1,1)
    for htool in [ht1, ht2, ht3]:
        c=color[ htool.this_looper.out_prefix]
        next_fraction=[]

        no_overlap=0
        all_overlap = 0
        for core_1 in htool.cores_used:
            this_over =  nar(htool.overlaps[core_1])
            next_fraction.append(this_over.max())
            if this_over[this_over>=0].sum() < 1e-32:
                no_overlap  += 1
            all_overlap += (this_over > 0.99).any()
        print( "all ", all_overlap)
        print("no overlap",no_overlap)
        ax3.hist( next_fraction, histtype='step',color=c,label="%s"%htool.this_looper.out_prefix)
        ax3.scatter([0],[no_overlap],c=c,marker="*")#,s=1.)
        ax3.scatter([1],[all_overlap],c=c,marker="*")#,s=1.)
        ax3.legend(loc=2)
        axbonk(ax3,ylabel=r'$N_{\rm{overlap}}$',xlabel=r'$f_{\rm{overlap}}$',ylim=[0,70])
        #fig.savefig('plots_to_sort/%s_overlaps.png'%htool.this_looper.out_prefix)
    fig3.savefig('plots_to_sort/next_overlap_dist.png')
            


if 0:
    import p56_plots.density_AC as AC
    reload(AC)
    if 'a1' not in dir():
        a1 = AC.ac_thing('u05'); a1.plot()
        a2 = AC.ac_thing('u10'); a2.plot()
        a3 = AC.ac_thing('u11'); a3.plot()
        acs={'u05':a1,'u10':a2,'u11':a3}

if 0:
    fig,ax=plt.subplots(1,1)
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        c=color[ ht.this_looper.out_prefix]
        hull_lengths = nar(ht.hull_volumes)**(1./3)
        vals, bins = np.histogram(hull_lengths)
        bc = 0.5*(bins[1:]+bins[:-1])
        ax.plot(bc,vals/vals.sum(),color=c,label=ht.this_looper.out_prefix,linestyle="--")
        axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$',ylim=[0,0.6])

        ac = acs[ht.this_looper.out_prefix]
        ax.plot(ac.binned[1],ac.binned[2], c=c)
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/hull_lengths.pdf')
    plt.close('fig')

if 0:
    fig,ax = plt.subplots(1,1)
    def n_neighbors_above_f(self,fraction):
        N_neighbors=[]
        for core_id in self.cores_used:
            N_neighbors.append( (self.overlaps[core_id]>=fraction).sum())
        return N_neighbors

    Nn = 100
    Nf = 10
    N_nei_f = np.zeros([Nf, Nn])
    for htool in [ht1, ht2, ht3]:
        for iFr,fr in enumerate(np.arange(.1,1,.1)):
            print("do ", htool.this_looper.out_prefix, " fr",fr)
            N_nei = nar(n_neighbors_above_f(htool,fr))
            print( "Max neighbors", max(N_nei))
            for iN in range(max(N_nei)+1):
                N_nei_f[iFr,iN] = (N_nei >= iN).sum()





    

if 0:
    fig,ax=plt.subplots(1,1)
    if 'ext_hull' not in dir():
        ext_hull=extents()
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        odd = nar(ht.hull_volumes) < nar(ht.cell_volumes)
        not_odd = nar(ht.hull_volumes) >= nar(ht.cell_volumes)
        ax.scatter(ht.hull_volumes,ht.cell_volumes,c='k')
        ax.scatter(nar(ht.hull_volumes)[odd],nar(ht.cell_volumes)[odd],c='r')
        ax.set_aspect('equal')
        #axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=[0,0.07],ylim=[0,0.07])
        ext_hull(nar(ht.hull_volumes))
        ext_hull(nar(ht.cell_volumes))
        ax.plot(ext_hull.minmax, ext_hull.minmax,c='g')
        axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',
               xlim=ext_hull.minmax,ylim=ext_hull.minmax, xscale='log',yscale='log')
        fig.savefig('plots_to_sort/%s_volume_fractions.pdf'%ht.this_looper.out_prefix)

if 0:
    fig,ax=plt.subplots(1,1)
    #ax.hist(hull_lengths,histtype='step',color='k',normed=True)
    vals, bins = np.histogram(hull_lengths)
    bc = 0.5*(bins[1:]+bins[:-1])
    db = (bins[1:]-bins[:-1])
    ax.plot(bc,vals/vals.sum())
    ax.plot(AC.binned[1],AC.binned[2]/AC.binned[2][0],c='r')
    rect=patches.Rectangle((0,0),AC.L,AC.ACb[0],facecolor=[0.8]*3)
    ax.add_patch(rect)
    fig.savefig('plots_to_sort/%s_sizes.pdf'%this_simname)


