from starter2 import *
import data_locations as dl
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0

this_simname = 'u10'
#this_simname = 'u05'
class ac_thing():
    def __init__(self,this_looper, field='density'):

        self.field=field
        self.this_simname = this_looper.sim_name
        self.directory = dl.sims[self.this_simname]
        frame=0
        self.ds = this_looper.load(frame)
        self.prefix="%s_n%04d"%(this_simname,frame)

        left=[0.0]*3
        resolution = self.ds['TopGridDimensions'] 
        self.cg=self.ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

        self.rho = self.cg[field].v -1 #[:40,:40,:40]

        self.rhohat = np.fft.fftn(self.rho)
        self.rho2 = self.rhohat*np.conj(self.rhohat)
        self.AC1 = np.fft.ifftn(self.rho2)
        self.ACc = np.fft.fftshift(self.AC1)
        self.AC = np.real(self.ACc)
        self.ACB = self.AC
        self.AC = self.AC/self.AC.max()

        #fig2.savefig("plots_to_sort/%s_proj.png"%self.this_simname)
        dx = 1./self.rho.shape[0]
        x,y,z = np.mgrid[0:1:dx,0:1:dx,0:1:dx]
        xcen, ycen, zcen= 0.5, 0.5, 0.5,
        self.rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)
        self.binned=rb.rb2( self.rmag, self.AC.real)#,bins=bins)

        self.ACb = self.binned[2]
        self.db = self.binned[0][1:] - self.binned[0][:-1]
        self.L = np.sum(self.ACb*self.db)/self.ACb[0]

    def plot(self):

        fig2,axes2=plt.subplots(2,2,figsize=figsize)
        a20,a21=axes2[0]
        a22,a23=axes2[1]


        if len(self.rho.shape) != 2:
            im1 = self.rho.sum(axis=axis)
            im2 = self.AC.sum(axis=axis)

        else:
            im1 = self.rho
            im2 = self.AC

        a20.imshow(im1, cmap='Greys' )
        a21.imshow(im2, cmap='Greys')
        a22.plot( self.binned[1],self.binned[2],c='r',label='binned ACfft')

        rect=patches.Rectangle((0,0),self.L,self.ACb[0],facecolor=[0.8]*3)
        a22.add_patch(rect)

        axbonk(a22,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
               #yscale='log',xscale='log')
               yscale='linear',xscale='linear')
        a23.legend(loc=0)

        fig2.savefig('plots_to_sort/%s_%s_AC.pdf'%(self.prefix, self.field))
        print('saved')

#a1 = ac_thing('u05'); a1.plot()
#a2 = ac_thing('u10'); a2.plot()
#a3 = ac_thing('u11'); a3.plot()
if 0:
    r1 = np.arange(31)/32 #rmag[33,33:]
    anne = 2.0*(np.arccos(r1)- r1*np.sqrt(1.0-r1**2))/np.pi;
    a22.plot(anne,c='g')

if 0:
    """fit to power laws.  Clearly this is not a great thing to do."""
    popt, pcov = curve_fit(powerlaw, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,1,-2])
    fit_rho0, fit_r0, fit_alpha = popt

    popt2, pcov2 = curve_fit(powerlaw2, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,-2])
    fit_r02, fit_alpha2 = popt2

    rrr = binned2[0]
    acac = binned2[1]
    a24.plot( rrr, 10**powerlaw(rrr, fit_rho0, fit_r0, fit_alpha),label='powerlaw')
    a24.plot( rrr, 10**powerlaw2(rrr, fit_r02, fit_alpha2),label='tweak2')
    a24.plot( rrr, 10**powerlaw2(rrr, L, fit_alpha2),label='tweak2')


if 0:
    def make_k_freqs(nk,real=False, d=1):
        ny = nk
        nz = nk
        if real:
            nz = nk//2+1

        k_freq = np.zeros([3,nk,ny,nz])
        k1=np.fft.fftfreq(nk,d=d)
        #kx, ky, kz = np.meshgrid(k1,k1,k1)
        #k_freq[0,...]=kx
        #k_freq[1,...]=ky
        #k_freq[2,...]=kz
        x = np.repeat(k1,nk*nk)
        x.shape = (nk,nk,nk)
        y = np.repeat(k1,nk*nk)
        y.shape = (nk,nk,nk)
        y=y.swapaxes(0,1)
        z = np.repeat(k1,nk*nk)
        z.shape = (nk,nk,nk)
        z=z.swapaxes(0,2)
        if real:
            x = x[:,:,:nz]
            y = y[:,:,:nz]
            z = z[:,:,:nz]
        k_freq[0,...]=x
        k_freq[1,...]=y
        k_freq[2,...]=z
        return k_freq



    if 'AC3d' not in dir() or True:
        print('Correlate')
        AC3d=scipy.signal.correlate(rho,rho,mode='same',method='fft')
        AC3d=np.roll(AC3d, AC3d.shape[0]//2,axis=0)
        AC3d=np.roll(AC3d, AC3d.shape[1]//2,axis=1)
        AC3d=np.roll(AC3d, AC3d.shape[2]//2,axis=2)
        print('rolled')
              

    print( 'fft')
    rhohat = np.fft.ifftn(rho)
    rhohatdag = rhohat.conj()
    rho_prod = rhohat*rhohatdag
    print( 'fft2')
    AC3dft = np.fft.fftn(rho_prod)

    k_array = make_k_freqs( AC3dft.shape[0],real=False)
    kmag = np.sqrt(k_array[0,...]**2 + k_array[1,...]**2 + k_array[2,...]**2)
    ktt=kmag.flatten()


    a21.imshow( (AC3d.real).sum(axis=axis) , cmap='Greys')
    proj = AC3dft.real.sum(axis=axis)
    a23.imshow(proj,cmap='Greys')
