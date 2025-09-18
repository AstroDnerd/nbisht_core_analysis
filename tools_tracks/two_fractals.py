


from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(trackage)
reload(dl)
plt.close('all')

dx=1./128
class two_images():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.fractal_dim=[]

    def run(self,core_list=None,nf=0):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))
        if core_list is None:
            print("Must pass in a core_list of two cores.")
        if len(core_list) != 2:
            print("Must pass in a core_list of two cores.")

        msd = {}
        ix = {}
        iy = {}
        iz = {}
        image={}
        print(core_list)
        mins=defaultdict(list)
        for core_id in core_list:
            msd[core_id] = trackage.mini_scrubber(thtr,core_id)
            ms = msd[core_id]

            ix[core_id] = (np.floor(ms.this_x/dx)[:,nf] ).astype('int')
            iy[core_id] = (np.floor(ms.this_y/dx)[:,nf] ).astype('int')
            iz[core_id] = (np.floor(ms.this_z/dx)[:,nf] ).astype('int')
            mins['x'].append(ix[core_id].min())
            mins['y'].append(iy[core_id].min())
            mins['z'].append(iz[core_id].min())
        for core_id in core_list:
            ms = msd[core_id]

            ix[core_id] = (np.floor(ms.this_x/dx)[:,nf] ).astype('int')
            iy[core_id] = (np.floor(ms.this_y/dx)[:,nf] ).astype('int')
            iz[core_id] = (np.floor(ms.this_z/dx)[:,nf] ).astype('int')
            ix[core_id] -= min( mins['x'])
            iy[core_id] -= min( mins['y'])
            iz[core_id] -= min( mins['z'])

            #ix[core_id] -= ix[core_id].min()
            #iy[core_id] -= iy[core_id].min()
            #iz[core_id] -= iz[core_id].min()

            #density = thtr.c([core_id],'density')[:,nf]
            #cell_volume = thtr.c([core_id],'cell_volume')[:,nf]

            image[core_id] = np.zeros([128,128,128])
            #mask = (x < 128)*(y<128)*(z<128)
            #x = x[ mask]
            #y = y[ mask]
            #z = z[ mask]
            #image = np.zeros( image_size)
            x,y,z = np.mgrid[0:1:dx, 0:1:dx, 0:1:dx]
            image[core_id][(ix[core_id],iy[core_id],iz[core_id])]=1
        fig,axes = plt.subplots(2,2,figsize=(12,12))
        axlist=axes.flatten()
        for nax,axis in enumerate([0,1,2]):
            x = [1,0,1][axis]
            y = [2,2,0][axis]
            ax=axlist[nax]
            xlab=r'$%s \rm(code\ length)$'%'xyz'[x]
            ylab=r'$%s \rm(code\ length)$'%'xyz'[y]
            axbonk(ax,xlabel=xlab,ylabel=ylab,xlim=[0,128],ylim=[0,128])
       
            image_2d_1 = image[core_list[0]].sum(axis=axis)
            image_2d_2 = image[core_list[1]].sum(axis=axis)
            if axis != 2:
                image_2d_1=image_2d_1.transpose()
                image_2d_2=image_2d_2.transpose()

            if 0:
                image_2d_1[image_2d_1>0]=1
                image_2d_2[image_2d_2>0]=1
                total_max = 1
                total_min = 0
            elif 0:
                #image_2d_1[image_2d_1>0]=1
                #image_2d_2[image_2d_2>0]=1
                total_max = max( [image_2d_1.max(),image_2d_2.max()])
                total_min = -2*total_max
            else:
                total_min = 1
                total_max = 2 #min( [image_2d_1.max(),image_2d_2.max()])
                #total_max = min( [total_max, image_2d_1.max()*2, image_2d_2.max()*2])
            norm1 = mpl.colors.Normalize(vmin=total_min,vmax=total_max)
            norm2 = mpl.colors.Normalize(vmin=total_min,vmax=total_max)
            n1=norm1( image_2d_1)
            n2=norm2( image_2d_2)
            n3=norm2( image_2d_1)*0
            cmap1 = copy.copy(mpl.cm.get_cmap("Reds"))
            cmap1.set_under('w')
            cmap2 = copy.copy(mpl.cm.get_cmap("Blues"))
            cmap2.set_under('w')
            c1 = cmap1(n1)
            c2 = cmap2(n2)
            final = (c1+c2)/2

            ax.imshow( final, interpolation='nearest',origin='lower')
            #plt.pcolormesh( x,y, z)

        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        norm = mpl.colors.Normalize(vmin=1,vmax = image_2d_1.max())
        THIS = n1
        axes[1][1].imshow( THIS, interpolation='nearest',cmap=cmap,norm=norm,origin='lower')
        fig.savefig('plots_to_sort/%s_two_fractal_c%04d_c%04d.png'%(self.this_looper.out_prefix,core_list[0],core_list[1]))


import three_loopers_tenfour as TL4
sim_list=['u401','u402','u403']
import analysis_block as ab


if 'ft' not in dir() or clobber:
    ft = {}

for sim_name in sim_list[:1]:
    ftool=two_images(TL4.loops['u401'])
    htool = ab.ht[sim_name]

    overlaps = np.zeros( [len(htool.cores_used)]*2) -1
    overlap_min = np.zeros( [len(htool.cores_used)]*2) -1
    for nc1,core_id_1 in enumerate(htool.cores_used):
        for nc2,core_id_2 in enumerate(htool.cores_used):
            overlaps[nc1,nc2] = htool.overlaps[core_id_1][nc2]
            overlap_min[nc1,nc2] = min([htool.overlaps[core_id_1][nc2],
                                         htool.overlaps[core_id_2][nc1]])



    if 1:
        #max overlap for 89
        core_1 = 89
        cores_used = nar(ab.ht[sim_name].cores_used)
        id_1 = np.where( cores_used == core_1)[0][0]
        id_2 = np.argmax(overlap_min[id_1,:])
        core_2 = ab.ht[sim_name].cores_used[id_2]

        ftool.run(core_list=[core_1,core_2])

    if 0:
        #apple banana for everyone
        for core_1 in ab.ht[sim_name].cores_used:
            cores_used = nar(ab.ht[sim_name].cores_used)
            id_1 = np.where( cores_used == core_1)[0][0]

            cores_used=nar(ab.ht[sim_name].cores_used)
            apples =  cores_used[(overlaps[id_1,:] == 0)*(overlaps[:,id_1]>0)]
            bananas = cores_used[(overlaps[id_1,:] > 0)*(overlaps[:,id_1]==0)]
            both = np.concatenate([apples,bananas])
            for core_2 in both:
                ftool.run(core_list=[core_1,core_2])
    if 0:
        ftool.run(core_list=[276,263])



