from starter2 import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
import tools_spectra.spectra_tools as st
reload(st)
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0
twopi = np.pi*2

class do_ac():
    def __init__(self,rho):
        self.rho = rho
        self.rhohat = np.fft.fftn(self.rho)
        self.rho2 = self.rhohat*np.conj(self.rhohat)
        self.AC1 = np.fft.ifftn(self.rho2)
        self.ACc = self.AC1 #np.fft.fftshift(self.AC1)
        self.ACc = np.fft.fftshift(self.AC1)
        self.AC = np.real(self.ACc)/self.rho.size
        self.ACunit = self.AC/self.rho.sum()
        
class make_sf():
    def __init__(self,this_looper,frame):
        dx = 1./128
        self.bins = np.arange(0.5*dx,1,dx)
        db = self.bins[1:]-self.bins[:-1]
        x0 = -0.5
        self.rall=np.mgrid[0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx]
        self.this_looper=this_looper
        self.ds = self.this_looper.load(frame)
        left=[0.0]*3
        resolution = self.ds['TopGridDimensions'] 
        self.cg=self.ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

        self.vx = self.cg['velocity_x'].v
        self.vy = self.cg['velocity_y'].v
        self.vz = self.cg['velocity_z'].v

        self.sigma_vx = np.std(self.vx)
        self.sigma_vy = np.std(self.vy)
        self.sigma_vz = np.std(self.vz)

        self.ac_x = do_ac(self.vx)
        self.ac_y = do_ac(self.vy)
        self.ac_z = do_ac(self.vz)
        self.ACx =self.ac_x.AC
        self.ACy =self.ac_y.AC
        self.ACz =self.ac_z.AC

        print('binning')
        norm = 1 #128**3/twopi
        self.v2a = (self.ACx+self.ACy+self.ACz)/norm
        self.sigma_3d = self.sigma_vx + self.sigma_vy + self.sigma_vz
    def bin_not_kludged(self):

        self.v2 = 2*(self.sigma_3d-2*(self.v2a))
        #v2n = v2/128**3*twopi
        self.binned2=rb.rb( self.rall, self.v2,bins=self.bins)
        #a23.plot( binned2[1], binned2[2])
    def bin_kludged(self):

        k1=1
        k2 = 1/5/k1
        self.v2 = (2*(self.sigma_3d*k1)-2*(self.v2a*k2))/4
        #v2n = v2/128**3*twopi
        self.binned2=rb.rb( self.rall, self.v2,bins=self.bins)
        #a23.plot( binned2[1], binned2[2])
    def plot(self):
        dx = 1./128
        bins = np.arange(0.5*dx,1,dx)
        fig,ax=plt.subplots(2,2)
        ax[0][0].imshow( self.ACx.sum(axis=0))
        ax[1][0].imshow( self.ACy.sum(axis=0))
        ax[0][1].imshow( self.ACz.sum(axis=0))
        v2abin=rb.rb( self.rall, self.v2a,bins=bins)
        ax[1][1].plot(v2abin[1],v2abin[2])
        self.v2abin=v2abin
        fig.savefig('plots_to_sort/ac.png')
        plt.close(fig)

