
from starter2 import *
import matplotlib.image as mpimg
from collections import defaultdict

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')

G = 1620/(4*np.pi)
class neighbors():
    def __init__(self,this_looper):
        self.this_looper=this_looper

        self.core_neighbors={}
        self.cores_used=[]
        self.n_neighbors={}
        self.mass={}  #mass for each core
        self.velocity = {}  #mean velocity for each core
        self.r = {}    #positions

        self.distance=defaultdict(list)

        ### temporary arrays for debugging
        self.KE=[]
        self.GE=[]
        self.R12=[]
        self.PAIR = []
        self.MS = {}
        self.DIST=[]

    def make_nd_arrays(self):
        self.mass_arr = np.array([self.mass[core_id] for core_id in self.cores_used])
        self.velocity_arr = np.array([self.velocity[core_id] for core_id in self.cores_used])
        self.neighbor_arr = np.array([self.n_neighbors[core_id] for core_id in self.cores_used])
        self.POS = np.array([ np.sqrt((self.r[core_id]*self.r[core_id]).sum()) for core_id in self.cores_used])
        self.RRR = np.array([ self.r[core_id] for core_id in self.cores_used])

        self.NEXT_DISTANCE = np.array( [sorted(self.distance[core_id]) for core_id in self.cores_used])
    def run(self,core_list=None):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        thtr.sort_time()
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            self.MS[core_id]=ms
            if ms.nparticles < 3:
                continue
            self.cores_used.append(core_id)

            nf = -1 #last frame
            dx = 1./2048
            nx = 1./dx

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
            self.mass[core_id]= (density[mask2]*cell_volume[mask2]).sum() 
            v = [ np.mean( ms.raw_vx[:,nf]),  np.mean( ms.raw_vy[:,nf]),  np.mean( ms.raw_vz[:,nf])]
            self.velocity[core_id]=np.array(v)
            self.r[core_id]= np.array([ms.mean_x[nf], ms.mean_y[nf], ms.mean_z[nf]])

        n_cores = len(self.cores_used)

        for i1 in range(0,n_cores-1):
            c1 = self.cores_used[i1] 
            for i2 in range(i1+1,n_cores):
                c2 = self.cores_used[i2] 
                reduced_mass = self.mass[c1]*self.mass[c2]/(self.mass[c1]+self.mass[c2])
                delta_v = self.velocity[c1] - self.velocity[c2]
                KE = 0.5*reduced_mass*(delta_v*delta_v).sum()
                r21 = np.sqrt(((self.r[c1] - self.r[c2])**2).sum())
                GE = G*self.mass[c1]*self.mass[c2]/r21

                self.distance[c1].append(r21)
                self.distance[c2].append(r21)


                if c1 not in self.core_neighbors:
                    self.core_neighbors[c1] = []
                if c2 not in self.core_neighbors:
                    self.core_neighbors[c2] = []
                if KE < GE:
                    self.core_neighbors[c1].append(c2)
                    self.core_neighbors[c2].append(c1)

                self.KE.append(KE)
                self.R12.append(r21)
                self.GE.append(GE)
                self.PAIR.append( [c1,c2])
                self.DIST.append(r21)
        self.KE  =nar(self.KE)
        self.R12 =nar(self.R12)
        self.GE  =nar(self.GE)
        self.PAIR=nar(self.PAIR)

        for core_id in self.cores_used:
            self.n_neighbors[core_id] = len(self.core_neighbors[core_id])
        self.make_nd_arrays()

        self.mean_distance=(self.NEXT_DISTANCE[:,0]+
                      self.NEXT_DISTANCE[:,1]+
                      self.NEXT_DISTANCE[:,2]+
                      self.NEXT_DISTANCE[:,3])/4
        #c=mean_distance
        #norm=mpl.colors.LogNorm(vmin=dmin, vmax=dmax)

        self.mean_density=4/(4*np.pi/3*self.mean_distance**3)


if 'do_all_plots' not in dir():
    do_all_plots = False


if do_all_plots:
    import three_loopers as tl
    reload(tl)
    if 'clobber' not in dir():
        clobber=False
    if 'neighbor_tool1' not in dir() or clobber:
        neighbor_tool1=neighbors(tl.looper1)
        neighbor_tool1.run()
    if 'neighbor_tool2' not in dir() or clobber:
        neighbor_tool2=neighbors(tl.looper2)
        neighbor_tool2.run()
    if 'neighbor_tool3' not in dir() or clobber:
        neighbor_tool3=neighbors(tl.looper3)
        neighbor_tool3.run()


    def plots(tool, fig, ax):
        ax1=ax[0][0]
        ax2=ax[0][1]
        ax3=ax[1][0]
        ax4=ax[1][1]

        minmin = min([min(tool.KE),min(tool.GE)])
        maxmax = max([max(tool.KE),max(tool.GE)])
        ax1.plot([minmin,maxmax],[minmin,maxmax],c=[0.5]*3)
        ax1.scatter( tool.KE, tool.GE,c='k',s=0.1)
        axbonk(ax1,xscale='log',yscale='log',xlabel='KE',ylabel='GE',xlim=[minmin,maxmax],ylim=[minmin,maxmax])


        ax2.scatter( tool.cores_used, tool.neighbor_arr)
        axbonk(ax2,xscale='linear',yscale='linear',xlabel='core id',ylabel='GE')

        #ax3.hist( np.log10(tool.NEXT_DISTANCE[0]))
        ax3.hist( np.log10(tool.NEXT_DISTANCE[:,0]))
        axbonk(ax3, xlabel='MinDistanceToNextCore',ylabel='N')
        #ax4.hist( np.log10(tool.NEXT_DISTANCE[:,1]), histtype='step',label=1)
        ax4.hist( np.log10(tool.NEXT_DISTANCE[:,0]), histtype='step',label=0)
        ax4.hist( np.log10(tool.NEXT_DISTANCE[:,1]), histtype='step',label=1)
        ax4.hist( np.log10(tool.NEXT_DISTANCE[:,3]), histtype='step',label=3)
        ax4.legend(loc=0)
        axbonk(ax4, xlabel='MinDistanceToThirdCore',ylabel='N')
        #ax3.hist(tool.DIST)
        #ax4.hist(tool.POS) 
        #ax4.hist(np.log10(tool.mass_arr))
        #ax3.hist(np.log10(tool.KE))
        #ax4.hist(np.log10(tool.GE))

        outname='plots_to_sort/%s_neighbors.pdf'%tool.this_looper.out_prefix
        fig.savefig(outname)
        print("saved",outname)
        plt.close(fig)
#fig,ax=plt.subplots(2,2)
#plots(neighbor_tool1,fig,ax)
#fig,ax=plt.subplots(2,2)
#plots(neighbor_tool2,fig,ax)
#fig,ax=plt.subplots(2,2)
#plots(neighbor_tool3,fig,ax)


    def neighbor_image(tool):
        #c=tool.NEXT_DISTANCE[:,1]
        #norm=mpl.colors.LogNorm(vmin=dmin, vmax=dmax)

        mean_distance=(tool.NEXT_DISTANCE[:,0]+
                      tool.NEXT_DISTANCE[:,1]+
                      tool.NEXT_DISTANCE[:,2]+
                      tool.NEXT_DISTANCE[:,3])/4
        #c=mean_distance
        #norm=mpl.colors.LogNorm(vmin=dmin, vmax=dmax)

        mean_density=4/(4*np.pi/3*mean_distance**3)
        c=mean_density
        
        c[c>1e4]=c.max()

        tool.c=c
        norm=mpl.colors.LogNorm(vmin=c.min(), vmax=c.max())
        fig,ax=plt.subplots(1,1)
        ax.hist(np.log10(mean_density),histtype='step')
        axbonk(ax,xlabel=r'$\log_{10}\rho_5$')
        fig.savefig('plots_to_sort/%s_density_dist.pdf'%(tool.this_looper.out_prefix))
        plt.close('fig')

        for line in [0,1,2]:
            fig,ax=plt.subplots(1,1)
            delta = 0.1
            ax.plot([0,1,1,0,0],[0,0,1,1,0])

            ax.set_xlim(-delta,1+delta); ax.set_ylim(-delta,1+delta)

            rm = rainbow_map(np.log10(np.sqrt(3)))
            dmin = 1./2048
            dmax=np.sqrt(3)
            x_coord=[1,2,0][line]
            y_coord=[2,0,1][line]
            x_lab=['y','z','x'][line]
            y_lab=['z','x','y'][line]



            p=ax.scatter( tool.RRR[:,x_coord], tool.RRR[:,y_coord],norm=norm, c=c)
            axbonk(ax,xlabel=x_lab,ylabel=y_lab)
            cb=plt.colorbar(p)
            fig.savefig('plots_to_sort/%s_neighbor_image_%s.pdf'%(tool.this_looper.out_prefix,line))
            plt.close('fig')


    neighbor_image(neighbor_tool1)
    neighbor_image(neighbor_tool2)
    neighbor_image(neighbor_tool3)
