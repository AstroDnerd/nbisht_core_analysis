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
if 1:
    #sims
    if 1:   
        """This is the one"""
        directory = dl.sims[this_simname]
        frame=0
        ds = yt.load("%s/DD%04d/data%04d"%(directory,frame,frame))
        prefix="%s_n%04d"%(this_simname,frame)
    if 0:   
        """dead"""
        directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u05'
    if 0:
        """a test"""
        directory = "/scratch1/dcollins/Paper19/u21_sphere_large"
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u21'
    if 0:
        """a smaller test test"""
        directory = "/scratch1/dcollins/Paper19/u22_sphere_128"
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u22'

    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

    rho = cg['density'].v -1 #[:40,:40,:40]
    size = rho.shape[0]
    dx = 1./size
    x,y,z = np.mgrid[0:size:1,0:size:1,0:size:1]
    xcen, ycen, zcen= size/2.,size/2., size/2.
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)

    dx = 1./size
    x,y,z = np.mgrid[0:1:dx,0:1:dx,0:1:dx]
    xcen, ycen, zcen= 0.5, 0.5, 0.5,
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)

R_sphere = 0.25
if 0:
    prefix='sphere_128'
    size=128
    dx = 1./size
    x,y,z = np.mgrid[0:size:1,0:size:1,0:size:1]

    xcen, ycen, zcen= size/2.,size/2., size/2.
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)
    rho = np.zeros_like(rmag)
    rho[ rmag < x.shape[0]/4] = 1.

if 0:
    prefix='sphere'
    size=64
    dx = 1./size
    x,y,z = np.mgrid[0:size:1,0:size:1,0:size:1]

    xcen, ycen, zcen= size/2.,size/2., size/2.
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)
    rho = np.zeros_like(rmag)
    rho[ rmag < x.shape[0]/4] = 1.


if 0:
    prefix='circle_manual'
    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    xcen, ycen = size/2.,size/2.
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2)
    rho = np.zeros_like(rmag)
    rho[ rmag < x.shape[0]/4] = 1.

if 0:
    prefix='circle'
    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    xcen, ycen = size/2.,size/2.
    dv = np.ones_like(x)*dx**3
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2)
    rho = np.array(Image.open("circ.tif")).astype('float')

if 1:
    #
    # Do the AC. 
    #
    rhohat = np.fft.fftn(rho)
    rho2 = rhohat*np.conj(rhohat)
    AC1 = np.fft.ifftn(rho2)
    ACc = np.fft.fftshift(AC1)
    AC = np.real(ACc)
    ACB = AC
    AC = AC/AC.max()

if 1:
    fig2,axes2=plt.subplots(2,2,figsize=figsize)
    a20,a21=axes2[0]
    a22,a23=axes2[1]


    if len(rho.shape) != 2:
        im1 = rho.sum(axis=axis)
        im2 = AC.sum(axis=axis)

    else:
        im1 = rho
        im2 = AC
    a20.imshow(im1, cmap='Greys' )
    a21.imshow(im2, cmap='Greys')


if 1:    
    binned=rb.rb2( rmag, AC.real)#,bins=bins)
    a22.plot( binned[1],binned[2],c='r',label='binned ACfft')

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

if 1:
    """actual correlation length"""

    ACb = binned[2]
    db = binned[0][1:] - binned[0][:-1]
    L = np.sum(ACb*db)/ACb[0]

if 1:
    """do we want the rectangle?"""
    rect=patches.Rectangle((0,0),L,ACb[0],facecolor=[0.8]*3)
    a22.add_patch(rect)

if 1:
    axbonk(a22,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
           #yscale='log',xscale='log')
           yscale='linear',xscale='linear')
    #a24.set_yscale('log')
    #x0a24.set_xscale('log')
    

if 0:
    a24.legend(loc=0)

fig2.savefig('plots_to_sort/%s_density_AC.pdf'%prefix)
print('saved')

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
