from starter2 import *
import three_loopers_u500 as TL
from collections import defaultdict
import tsing
reload(tsing)
plt.close('all')
class TEKE():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.tsing_list=[]
        self.mode_list=[]
        self.quan_list=[defaultdict(list), defaultdict(list)]
    def run(self, core_list=None,tsing=None,frame=0):
        this_looper=self.this_looper
        thtr=this_looper.tr

        ds = this_looper.load(frame)
        ad=ds.all_data()
        D = ad[YT_density].v
        TE = D*np.log(D)+1/np.e
        XX = TE
        YY= ad[YT_kinetic_energy].v

        
        #xbins_1 = np.linspace( -1/np.e,1/np.e,64)
        ##xbins_2 = np.geomspace(1/np.e,TE.max(),128)
        #dx =  xbins_1[1]-xbins_1[0]
        #xbins_2 = np.exp(np.arange(np.log(1/np.e),np.log(TE.max()),dx))
        #xbins = np.unique( np.concatenate([xbins_1,xbins_2]))
        xbins = np.geomspace(TE.min(),TE.max(),128)
        ybins = np.geomspace( YY.min(), YY.max(), 128)

        fig,ax=plt.subplots(1,1)
        ax.hist(TE, bins=xbins)
        print(xbins)
        #ax.set(yscale='log')
        #ax.set_xscale('symlog',linthresh=1/np.e)
        ax.set_xscale('log')
        ax.axvline(1/np.e)
        ax.set_yscale('log')

        fig.savefig('plots_to_sort/teke1_%s'%this_looper.sim_name)

        hist, xb, yb = np.histogram2d(XX,YY, bins=[xbins,ybins], density=True)
        import pcolormesh_helper as pch
        reload(pch)
        from scipy.ndimage import gaussian_filter as gf
        def dosmo(arr,n):
            #return arr
            return gf(arr,1)

        fig,ax=plt.subplots(1,3,figsize=(12,8))
        ax0=ax[0];ax1=ax[1]; ax2=ax[2]
        pch.helper(dosmo(hist,1),xb,yb,ax=ax0)
        #ax.set_xscale('symlog',linthresh=1/np.e)
        ax0.set_xscale('log')
        ax0.set_yscale('log')

        rho = thtr.track_dict['density'][:,0]
        vx = thtr.track_dict['velocity_x'][:,0]
        vy = thtr.track_dict['velocity_y'][:,0]
        vz = thtr.track_dict['velocity_z'][:,0]

        TEp = rho*np.log(rho)+1/np.e
        KEp = 0.5*(rho*(vx**2+vy**2+vz**2))
        #ax1.scatter(TEp,KEp,s=0.1,alpha=0.1)
        h2,x2,y2 = np.histogram2d( TEp, KEp, bins=[xbins,ybins], density=True)
        pch.helper(dosmo(h2,1),x2,y2,ax=ax1)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        ratio = np.zeros_like(h2)
        ok=hist>0
        ratio[ok]=dosmo(h2,1)[ok]/dosmo(hist,1)[ok]
        out=pch.helper(ratio,x2,y2,ax=ax2)#,zlim=[1e-1,1])
        fig.colorbar(out['plot'])
        ax2.set_xscale('log')
        ax2.set_yscale('log')

        #ax0.set(xlabel='TE+1/e',ylabel='KE',title='k

        MachPart = np.sqrt(vx**2+vy**2+vz**2).mean()
        MachCode = ad['velocity_magnitude'].mean()
        print(MachPart, MachCode)




        fig.savefig('plots_to_sort/TEKE_%s.png'%this_looper.sim_name)






sim_list=['u501']#,'u502','u503']

if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

for sim in sim_list:
    tool = TEKE(TL.loops[sim])
    tool.run()
