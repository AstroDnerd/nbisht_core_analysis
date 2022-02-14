


from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')

#from https://gist.github.com/viveksck/1110dfca01e4ec2c608515f0d5a5b1d1
def fractal_dimension(Z, threshold=0.9):
    do_plots = False

    # Only for 2d image
    assert(len(Z.shape) == 3)

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k,ax):
        S=np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0)
        S=np.add.reduceat(S, np.arange(0, Z.shape[1], k), axis=1)
        S=np.add.reduceat(S, np.arange(0, Z.shape[2], k), axis=2)
        #S = np.add.reduceat(
        #    np.add.reduceat(
        #    np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
        #                       np.arange(0, Z.shape[1], k), axis=1),
        #                       np.arange(0, Z.shape[1], k), axis=2)

        ## the original did this.  We want full boxes as well.
        ## We count non-empty (0) and non-full boxes (k*k)
        #return len(np.where((S > 0) & (S < k*k))[0])
        if do_plots:
            ax.imshow(np.sum(S,axis=1))
        return len(np.where((S > 0)  )[0])


    # Transform Z into a binary array
    Z = (Z > threshold)


    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)

    #don't repeat sizes that are too big.
    image_length= np.max( nar(np.where(Z)))
    image_length = 2**np.round(np.log2(image_length))
    sizes = sizes[ sizes <= image_length]

    # Actual box counting with decreasing size
    counts = []
    if do_plots:
        fig, ax = plt.subplots(1,len(sizes), figsize=(4*len(sizes),4))
    for ns,size in enumerate(sizes):
        counts.append(boxcount(Z, size,aaa))
        if do_plots:
            aaa=ax[ns]
            aaa.set_title("%d %d"%(size,counts[-1]))
    if do_plots:
        fig.savefig('plots_to_sort/tmp.png')

    # Fit the successive log(sizes) with log (counts)

    sizes=nar(sizes)
    counts=nar(counts)
    ok = counts > 1
    if ok.sum():
        coeffs = np.polyfit(np.log(sizes[ok]), np.log(counts[ok]), 1)
    else:
        coeffs = [0,0]
    #plt.clf()
    #plt.plot(np.log(sizes), np.log(counts),c='k',marker='x')
    #plt.plot(np.log(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
    #nfractal=len(glob.glob("plots_to_sort/nfractal*")) #plt.savefig("plots_to_sort/nfractal_%s.png"%nfractal)
    return -coeffs[0], {'counts':counts,'sizes':sizes, 'coeffs':coeffs}

dx=1./128
class fractal_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.fractal_dim=[]

    def run(self,core_list=None,nf=0):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))
        if core_list is None:
            core_list = all_cores
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.r.shape[0] <= 4:
                continue
            print("fractal dim on ", core_id)
            self.cores_used.append(core_id)

            #x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]#or whatever the number of zones is
            #y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
            #z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
            x = (np.floor(ms.this_x/dx)[:,nf] ).astype('int')
            y = (np.floor(ms.this_y/dx)[:,nf] ).astype('int')
            z = (np.floor(ms.this_z/dx)[:,nf] ).astype('int')
            #x,y,z = np.mgrid[0:128:1, 0:128:1, 37:38:1].astype('int')
            density = thtr.c([core_id],'density')[:,nf]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nf]

            x = x-x.min()
            y = y-y.min()
            z = z-z.min()

            image_size=[x.max()+1,y.max()+1,z.max()+1]
            image = np.zeros([128,128,128])
            mask = (x < 128)*(y<128)*(z<128)
            x = x[ mask]
            y = y[ mask]
            z = z[ mask]
            image[(x,y,z)]=1
            #image = np.zeros( image_size)
            if 0:
                print("image size",core_id,image_size)
                ix,iy,iz=np.mgrid[0:128:1,0:128:1,0:128:1]
                #ok = ( (ix**2+iy**2+iz**2) < 56**2 )*( (ix**2+iy**2+iz**2) > 54**2)
                ok = ( (iy**2+iz**2) < 56**2 )*( (iy**2+iz**2) > 54**2)
                ok += (iy==0)*(iz<56)
                ok += (iz==0)*(iy<56)

                image=image*0
                image[ok]=1
            def fill_in(A):
                oot=np.ones_like(A)
                for dim in [0,1,2]:#[0,1,2]:
                    P = np.cumsum(A,axis=dim)
                    E = np.sum(A,axis=dim)
                    oot *= P*(E-P)
                    print("OOT",oot.sum())

                return oot
            cmap=copy.copy(mpl.cm.get_cmap("viridis"))
            cmap.set_under('w')
            fig7,ax7=plt.subplots(1,2)
            norm=mpl.colors.Normalize(1,20)
            ax7[0].imshow(np.sum(image,axis=0), origin='lower',interpolation='nearest',cmap=cmap,norm=norm)
            B =  fill_in(image)  + image
            print("B sum",B.sum())
            #pdb.set_trace()
            ax7[1].imshow(np.sum(B,axis=0), origin='lower',interpolation='nearest',cmap=cmap,norm=norm)
            fig7.savefig('plots_to_sort/fill_in_c%04d.png'%core_id)
            plt.close(fig7)
                        
                
            print("Total points", x.size, image.sum())
            stuff = fractal_dimension(image)
            dimension = stuff[0]
            self.fractal_dim.append(dimension)
            print("fractal dim c%04d %0.2f"%(core_id,self.fractal_dim[-1]))
            if 1:
                plt.clf()
                cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                cmap.set_under('w')
                image_2d = image.sum(axis=0)
                norm = mpl.colors.Normalize(vmin=1,vmax=image_2d.max())
                plt.imshow( image_2d,cmap=cmap,norm=norm)

                plt.savefig('plots_to_sort/%s_fractal_c%04d.png'%(self.this_looper.out_prefix,core_id))

                counts = stuff[1]['counts']
                sizes  = stuff[1]['sizes']
                coeffs  = stuff[1]['coeffs']
                plt.clf()
                plt.plot(np.log10(sizes), np.log(counts),c='k',marker='x')
                plt.plot(np.log10(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
                plt.savefig('plots_to_sort/%s_counts_c%04d.png'%(self.this_looper.out_prefix,core_id))


#import three_loopers_1tff as tl
#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4




if 'toolshed' not in dir():
    toolshed = {}

framelist=[0]#,1,2,3,4,5,6,7,8]
if 1:
    looper_list = TL4.loops
    
    for simname in looper_list:
        looper=looper_list[simname]
        name =  looper.sim_name

        if name not in toolshed:
            toolshed[ name ] = {}
            toolshed[ name ]['looper'] = looper
        for nframe, frame in enumerate(framelist):
            if nframe not in framelist:
                continue
            if frame in toolshed[ name]:
                continue
            toolshed[name][frame] = fractal_tool( looper )
            toolshed[name][frame].run(nf=nframe)

if 1:
    tool = fractal_tool( TL4.loops['u403'])
    tool.run(core_list=[0])
if 0:
    fig,ax=plt.subplots(1,1)
    for nsim,simname in enumerate(toolshed):
        print(nsim,simname)
        ax.hist(toolshed[simname][0].fractal_dim, histtype='step',color=colors.color[simname])
    outname="%s/fractal_dims.png"%"plots_to_sort"
    fig.savefig(outname)

if 0:
    for nf,ft in enumerate([ft1,ft2,ft3]):
        name=ft.this_looper.out_prefix
        plt.hist(ft.fractal_dim,histtype='step',bins=10,color=colors.color[name],label=name)
        #plt.savefig("plots_to_sort/wtf_%d.png"%nf)

    plt.legend(loc=0)
    plt.savefig('plots_to_sort/fractal_dist.png')

if 0:
    rm = rainbow_map(10)
    for name in toolshed:
        dims = []
        #mylooper=looper_list[0]
        mylooper=toolshed[name]['looper']
        for nframe,frame in enumerate(toolshed[name]):
            tool = toolshed[name][frame]
            lab=ft.this_looper.out_prefix + " %s"%frame
            plt.hist(tool.fractal_dim,histtype='step',bins=10,color=rm(nframe),label=lab)
            dims.append(tool.fractal_dim)
        plt.savefig("plots_to_sort/wtf_%d.png"%nf)
        plt.savefig('plots_to_sort/fractal_dist.png')
        fig,ax=plt.subplots(1,1)
        ax.violinplot(dims)
        fig.savefig('plots_to_sort/violin_test.png')
        plt.close(fig)
        fig,axes=plt.subplots(1,2)
        dims=nar(dims)

        cores_used = tool.cores_used
        OverAchievers = []
        for nb,b in enumerate(dims.transpose()):
            core_id = cores_used[nb]
            if (b > 2.5).any():
                OverAchievers.append(core_id)
            nparticles=mylooper.snaps[0][core_id].target_indices.size
            if nparticles < 10:
                continue

            axes[0].plot(mylooper.tr.times[framelist], b/b[0:3].mean())
            axes[1].plot(mylooper.tr.times[framelist], b)
        fig.savefig('plots_to_sort/streams.png')

#FFF=fractal_dimension(image)

#Z=image
#k=4
#b=np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0)
#b=np.add.reduceat(b, np.arange(0, Z.shape[1], k), axis=1)
#b=np.add.reduceat(b, np.arange(0, Z.shape[2], k), axis=2)
