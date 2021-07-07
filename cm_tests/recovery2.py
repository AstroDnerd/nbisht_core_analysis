

from starter2 import *
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D  

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')

def two_point(ms,density,cell_volume,i1,i2,i3,frame):
    x0 = ms.this_x[:,frame]
    y0 = ms.this_y[:,frame]
    z0 = ms.this_z[:,frame]
    x1 = ms.this_x[:,frame+1]
    y1 = ms.this_y[:,frame+1]
    z1 = ms.this_z[:,frame+1]

    a1 = [x0[i1],y0[i1],z0[i1]]
    a2 = [x0[i2],y0[i2],z0[i2]]
    a3 = [x0[i3],y0[i3],z0[i3]]
    dA = np.column_stack([a1,a2,a3])
    rho_a1 = density[i1,frame]
    rho_a2 = density[i2,frame]
    rho_a3 = density[i3,frame]
    rho_a = (rho_a1+rho_a2+rho_a3)/3
    b1 = [x1[i1],y1[i1],z1[i1]]
    b2 = [x1[i2],y1[i2],z1[i2]]
    b3 = [x1[i3],y1[i3],z1[i3]]
    dB = np.column_stack([b1,b2,b3])
    rho_b1 = density[i1,frame+1]
    rho_b2 = density[i2,frame+1]
    rho_b3 = density[i3,frame+1]
    rho_b = (rho_b1+rho_b2+rho_b3)/3
    #the transpose switches left and right multiply
    rho_a = (cell_volume[:,frame]*density[:,frame]).sum()/cell_volume[:,frame].sum()
    rho_b = (cell_volume[:,frame+1]*density[:,frame+1]).sum()/cell_volume[:,frame+1].sum()
    det_dA = np.linalg.det(dA)
    if np.abs(det_dA ) > 1e-16:
        F = np.linalg.solve(dA.T,dB.T).T
    else:
        F = None
    return F, rho_a, rho_b


def all_simplices(d,nvert):
        return [i for i,s in enumerate(d.simplices) if nvert in s]
class lagrange():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.density_ratio=[]
        self.det = []

    def run(self,do_all_plots=True,core_list=None,frame1=10, frame2=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()
        if frame2 is None:
            frame2 = frame1+1

        tsorted = thtr.times
        self.core_list=core_list
        fig,ax=plt.subplots(1,1)


        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, projection='3d')
        from scipy.spatial import Delaunay
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            rm = rainbow_map(ms.nparticles)
            if ms.nparticles < 3:
                continue
            print('go ', core_id, ms.nparticles)
            self.cores_used.append(core_id)
            indices = list(range(ms.nparticles))
            density=thtr.c([core_id],'density')
            density1 = density[:,frame1]
            density2 = density[:,frame2]
            cell_volume=thtr.c([core_id],'cell_volume')
            cell_volume1 = cell_volume[:,frame1]
            cell_volume2 = cell_volume[:,frame2]
            x1=self.this_looper.snaps[frame1][core_id].pos
                #half edge structure
            self.tri = Delaunay(x1)


            x2=self.this_looper.snaps[frame2][core_id].pos
            print(x1.shape)
            rm = rainbow_map(4)

            ds1 = self.this_looper.load(frame1)
            domain_dimensions=ds1.domain_dimensions
            cg1 = ds1.covering_grid(0,[0.0]*3,domain_dimensions)
            density_proj_1 = cg1['density'].sum(axis=2)/domain_dimensions[2]
            projx = cg1['x'].sum(axis=2)/domain_dimensions[2]
            projy = cg1['y'].sum(axis=2)/domain_dimensions[2]

            ds2 = self.this_looper.load(frame2)
            domain_dimensions=ds2.domain_dimensions
            cg2 = ds2.covering_grid(0,[0.0]*3,domain_dimensions)
            density_proj_2 = cg2['density'].sum(axis=2)/domain_dimensions[2]
            projx = cg2['x'].sum(axis=2)/domain_dimensions[2]
            projy = cg2['y'].sum(axis=2)/domain_dimensions[2]

            figp, axesp = plt.subplots(2,2, figsize=(8,8))
            axp=axesp[0]


            pcolormesh_args = {'cmap':'gray','shading':'nearest'}
            axp[0].pcolormesh(projx,projy,density_proj_1,**pcolormesh_args)
            axp[0].set_aspect('equal')
            axp[1].pcolormesh(projx,projy,density_proj_1,**pcolormesh_args)
            axp[1].set_aspect('equal')

            axesp[1][0].pcolormesh(projx,projy,density_proj_2,**pcolormesh_args)
            axesp[1][0].set_aspect('equal')
            ratio_proj = (cg2['density']/cg1['density']).sum(axis=2)/domain_dimensions[2]
            axesp[1][1].pcolormesh(projx,projy,ratio_proj,**pcolormesh_args)
            axesp[1][1].set_aspect('equal')

            self.ijk_1 = (x1*domain_dimensions).astype('int')
            self.ijk_2 = (x2*domain_dimensions).astype('int')

            self.cell_jump = ((self.ijk_1-self.ijk_2)**2).sum(axis=1)
            self.colors = np.array(['y']*self.cell_jump.size)
            self.colors[ self.cell_jump > 0] = 'r'



            for pi in indices:
                print('p', pi, ms.nparticles)
                #print(all_simplices(self.tri,pi))
                #Loop over simplices.  Each simplex is a triangle.
                #Average FTF.
                x1i = x1[pi,:]
                x2i = x2[pi,:]
                my_simplices = all_simplices(self.tri,pi)
                massB = 0; massA = 0; detFTF = 0
                volumeB = 0; volumeA = 0;
                for n_simplex,my_simplex_id in enumerate(my_simplices):
                    my_simplex= self.tri.simplices[my_simplex_id]
                    other_points=my_simplex[my_simplex != pi]
                    dA = x1[other_points]-x1i
                    dB = x2[other_points]-x2i
                    F = np.linalg.solve(dA.T,dB.T).T
                    #The right Cauchy–Green deformation tensor
                    C = F.T@F
                    detFTF += np.linalg.det(C)
                    massA += (cell_volume1[my_simplex]*density1[my_simplex]).sum()
                    massB += (cell_volume2[my_simplex]*density2[my_simplex]).sum()
                    volumeA += cell_volume1[my_simplex].sum()
                    volumeB += cell_volume2[my_simplex].sum()
                #densityA = density1[pi]
                #densityB = density2[pi]
                densityA = massA/volumeA
                densityB = massB/volumeB
                detFTF /= len(my_simplices)
                detF = np.sqrt(detFTF)
                self.density_ratio.append(densityB/densityA)
                self.det.append(detF)
                full_metric = self.density_ratio[-1]*self.det[-1]
                if True: #full_metric <= 1.5:
                    print('   plot')
                    axp[0].plot([x1i[0],x2i[0]],[x1i[1],x2i[1]], c=rm(1/self.det[-1]),marker='o')
                    axp[1].plot([x1i[0],x2i[0]],[x1i[1],x2i[1]], c=rm(self.density_ratio[-1]),marker='o')
                    axp[0].set_title('1/det F')
                    axp[1].set_title('rho2/rho2')
                if False:
                    for n_simplex,my_simplex_id in enumerate(my_simplices):
                        my_simplex= self.tri.simplices[my_simplex_id]
                        other_points=my_simplex[my_simplex != pi]
                        dA = x1[other_points]-x1i
                        dB = x2[other_points]-x2i
                        for i in [0,1,2]:
                            c = rm( self.density_ratio[-1]*self.det[-1])
                            ax.plot([x1i[0],x1i[0]+dA[i,0]],[x1i[1],x1i[1]+dA[i,1]],c=c)
                            #ax.plot([x2i[0],x2i[0]+dB[i,0]],[x2i[1],x2i[1]+dB[i,1]],c=)
                        if np.abs(self.det[-1]) > 1e-16:
                            ax2.scatter(x1i[0],x1i[1],np.abs(1/self.det[-1]),marker='*')


        figp.savefig('plots_to_sort/%s_proj1.png'%prefix)
        fig.savefig('plots_to_sort/%s_volumes.png'%self.this_looper.out_prefix)
        plt.close(fig)

        ax2.set_xlabel('x')
        fig2.savefig('plots_to_sort/%s_dets.png'%self.this_looper.out_prefix)
        self.density_ratio=np.array(self.density_ratio)
        self.det = np.array(self.det)






import three_loopers as tl
import looper_u14 as u14

if 0:
    tool1=lagrange(tl.looper1)
    tool1.run()#core_list=[10,11])

if 1:
    tool1 = lagrange(u14.looper14)
    tool1.run(frame1=5,frame2=6)

if 1:
    fig,ax=plt.subplots(1,1)
    prefix=tool1.this_looper.out_prefix
    ok = tool1.density_ratio > 0
    ok = ok * (np.abs(tool1.det) > 1e-16)
    the_y = np.abs(1./tool1.det[ok])
    the_x = tool1.density_ratio[ok]
    ax.scatter(the_x,the_y,c=tool1.colors)
    minmin = min([the_y.min(),the_x.min()])
    maxmax = max([the_y.max(),the_x.max()])
    minmax=[minmin,maxmax]
    axbonk(ax,xscale='linear',yscale='linear',xlabel='rho2/rho1',ylabel='det(F)',xlim=minmax,ylim=minmax)
    ax.set_xlim( minmin,maxmax)
    ax.set_ylim( minmin,maxmax)
    ax.plot([minmin,maxmax],[minmin,maxmax],c=[0.5]*4)
    #ax.set_xlim([1,10])
    fig.savefig('plots_to_sort/%s_lagrangian.png'%prefix)
    ax.clear()
    ax.hist(the_y,label='1/detF x',histtype='step',bins=16)
    ax.hist(the_x,label='rho2/rho1 x',histtype='step',bins=16)
    ax.legend(loc=1)
    fig.savefig('plots_to_sort/%s_lagrang_hists.png'%prefix)
    ax.clear()
    #ax.hist(np.abs(1./tool1.det[ok])/(tool1.density_ratio[ok]), histtype='step',bins=16)
    ax.hist(np.abs(tool1.det[ok])*(tool1.density_ratio[ok]), histtype='step',bins=16)
    axbonk(ax,xlabel=r'$\rho_2/\rho_1 det(F)$',ylabel=r'$N$')
    fig.savefig('plots_to_sort/%s_the_hist.png'%prefix)
