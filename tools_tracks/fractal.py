


from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')

#from https://gist.github.com/viveksck/1110dfca01e4ec2c608515f0d5a5b1d1
def fractal_dimension(Z, threshold=0.9):

    # Only for 2d image
    assert(len(Z.shape) == 3)

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
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
        return len(np.where((S > 0)  )[0])


    print("S2 ",Z.sum())
    # Transform Z into a binary array
    Z = (Z > threshold)
    print("S3 ",Z.sum())

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)

    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    # Fit the successive log(sizes) with log (counts)

    sizes=nar(sizes)
    counts=nar(counts)
    ok = counts > 1
    coeffs = np.polyfit(np.log(sizes[ok]), np.log(counts[ok]), 1)
    #plt.clf()
    #plt.plot(np.log(sizes), np.log(counts),c='k',marker='x')
    #plt.plot(np.log(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
    #nfractal=len(glob.glob("plots_to_sort/nfractal*"))
    #plt.savefig("plots_to_sort/nfractal_%s.png"%nfractal)
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
            density = thtr.c([core_id],'density')[:,nf]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nf]

            x = x-x.min()
            y = y-y.min()
            z = z-z.min()

            image_size=[x.max()+1,y.max()+1,z.max()+1]
            image = np.zeros([128,128,128])
            #image = np.zeros( image_size)
            print("image size",core_id,image_size)
            image[(x,y,z)]=1
            print("Total points", x.size, image.sum())
            stuff = fractal_dimension(image)
            dimension = stuff[0]
            self.fractal_dim.append(dimension)
            print("fractal dim c%04d %0.2f"%(core_id,self.fractal_dim[-1]))
            if 1:
                plt.clf()
                plt.imshow( image.sum(axis=0))

                plt.savefig('plots_to_sort/%s_fractal_c%04d.png'%(self.this_looper.out_prefix,core_id))

                counts = stuff[1]['counts']
                sizes  = stuff[1]['sizes']
                coeffs  = stuff[1]['coeffs']
                plt.clf()
                plt.plot(np.log(sizes), np.log(counts),c='k',marker='x')
                plt.plot(np.log(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
                plt.savefig('plots_to_sort/%s_counts_c%04d.png'%(self.this_looper.out_prefix,core_id))




if 'ft1' not in dir() or clobber:
    ft1=fractal_tool(tl.looper1)
    ft1.run()
if 'ft2' not in dir() or clobber:
    ft2=fractal_tool(tl.looper2)
    ft2.run()
if 'ft3' not in dir() or clobber:
    ft3=fractal_tool(tl.looper3)
    ft3.run()

plt.clf()
color={'u05':'r','u10':'g','u11':'b', 'u201':'r', 'u202':'g','u203':'b'}

for nf,ft in enumerate([ft1,ft2,ft3]):
    name=ft.this_looper.out_prefix
    plt.hist(ft.fractal_dim,histtype='step',bins=10,color=color[name],label=name)
    plt.savefig("plots_to_sort/wtf_%d.png"%nf)

plt.legend(loc=0)
plt.savefig('plots_to_sort/fractal_dist.png')
#FFF=fractal_dimension(image)

#Z=image
#k=4
#b=np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0)
#b=np.add.reduceat(b, np.arange(0, Z.shape[1], k), axis=1)
#b=np.add.reduceat(b, np.arange(0, Z.shape[2], k), axis=2)
