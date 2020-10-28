from starter2 import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
import spectra_tools as st
reload(st)
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0
twopi = np.pi*2

if 'vx' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

    #rho1 = cg['density'].v #-1 #[:40,:40,:40]
    #rho2=rho1
    def make_sphere(size):
        dx = 1./size
        x,y,z = np.mgrid[0:size:1,0:size:1,0:size:1]
        xcen, ycen, zcen= size/2.,size/2., size/2.
        rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2+ (z-zcen)**2)
        rho = np.zeros_like(rmag)
        rho[ rmag < x.shape[0]/4] = 9.
        pdb.set_trace()
        return rho
    #vx = make_sphere(128)
    #vy = np.zeros_like(vx)
    #vz = np.zeros_like(vx)
    vx = cg['velocity_x'].v
    vy = cg['velocity_y'].v
    vz = cg['velocity_z'].v

    sigma_vx = np.std(vx)
    sigma_vy = np.std(vy)
    sigma_vz = np.std(vz)

if 1:
    fig2,axes2=plt.subplots(2,2,figsize=figsize)
    a20,a21=axes2[0]
    a22,a23=axes2[1]

    fig3,ax2=plt.subplots(1,1,figsize=figsize)

class do_ac():
    def __init__(self,rho):
        self.rho = rho
        self.rhohat = np.fft.fftn(self.rho)
        self.rho2 = self.rhohat*np.conj(self.rhohat)
        self.AC1 = np.fft.ifftn(self.rho2)
        self.ACc = self.AC1 #np.fft.fftshift(self.AC1)
        self.ACc = np.fft.fftshift(self.AC1)
        self.AC = np.real(self.ACc)
        self.ACunit = self.AC/self.rho.sum()

if 1:
    ac_x = do_ac(vx)
    ac_y = do_ac(vy)
    ac_z = do_ac(vz)
    ACx = ac_x.AC
    ACy = ac_y.AC
    ACz = ac_z.AC

if 0:
#   if 'ACx' not in dir():
#       print('Correlate')
#       ACx=scipy.signal.correlate(vx,vx,mode='same',method='fft')
#       ACy=scipy.signal.correlate(vy,vy,mode='same',method='fft')
#       ACz=scipy.signal.correlate(vz,vz,mode='same',method='fft')


    a20.imshow( (ACx.real).sum(axis=axis) , cmap='Greys')
    a21.imshow( (ACy.real).sum(axis=axis) , cmap='Greys')
    #a22.imshow( (ACz.real).sum(axis=axis) , cmap='Greys')

#
# check normalizations of rho/rhohat and rhohat^2/ AC with parseval's thm.
#
if 0:
    sigma_vxhat = (ac_x.rhohat*ac_x.rhohat.conj()).real.sum()/ac_x.rhohat.size/twopi
    sigma_vxhatn = twopi*sigma_vxhat/128**3
    sigma_vx_m = (ac_x.rho**2).sum()/ac_x.rho.size
    parseval_error = 1-(sigma_vxhatn)/sigma_vx
    print("Sum_rhohat %0.1e sum_rho %0.1e error %0.2e"%( sigma_vxhatn,sigma_vx_m,parseval_error ))

if 0:
    sigma_v2h = (ac_x.rho2*ac_x.rho2.conj()).real.sum()/ac_x.rhohat.size/twopi
    sigma_v2hn = twopi*sigma_v2h/128**3
    sigma_AC = (ac_x.AC1*ac_x.AC1.conj()).sum()/ac_x.AC1.size
    parseval_error = 1-(sigma_v2hn)/sigma_AC
    print("Sum_rhohat %0.1e sum_rho %0.1e error %0.2e"%( sigma_v2hn,sigma_AC,parseval_error ))

if 1:    
    dx = 1./128
    bins = np.arange(0.5*dx,1,dx)
    db = bins[1:]-bins[:-1]
    x0 = -0.5
    rall=np.mgrid[0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx]

if 0:
    print('binning')
    norm = 128**3/twopi
    v2a = (ACx+ACy+ACz)/norm
    sigma_3d = sigma_vx + sigma_vy + sigma_vz
    v2 = 2*sigma_3d-2*v2a
    #v2n = v2/128**3*twopi
    binned2=rb.rb( rall, v2,bins=bins)
    a23.plot( binned2[1], binned2[2])

if 0:
    vxhat = np.fft.fftn(vx)/vx.size
    vyhat = np.fft.fftn(vy)/vx.size
    vzhat = np.fft.fftn(vz)/vx.size
    power = vxhat*vxhat.conj() + vyhat*vyhat.conj() + vzhat*vzhat.conj()
    kspace, power_1d, ff, ksize = st.shell_average_raw(power)
    power_1d = power_1d.real#/ksize
    a22.plot(kspace, power_1d)
    kolmog = kspace**(-5./3)
    kolmog *= power_1d[2]/kolmog[2]
    a22.plot(kspace, kolmog,c=[0.5]*4)
    axbonk(a22, xlabel='k',ylabel='P_v', xscale='log',yscale='log')

if 1:


    hathat = np.fft.ifft(power_1d)
    a23.plot(np.linspace(0,1,hathat.size), sigma_3d - hathat,label='horse')
    a23.legend(loc=0)




    fig2.savefig('plots_to_sort/vel.png')




if 0:
    #a24.plot( binned[0],binned[1],c='r',label='binned ACfft')
    a24.plot( binned2[0],binned2[1],c='g',label=r'$AC_\rho(r)$')

    ok = binned2[1]>0
    #fits = np.polyfit(np.log10(binned2[0][ok]), np.log10(binned2[1][ok]),2)
    #a24.plot( binned2[0][ok], binned2[0][ok]**fits[-1])

if 0:
    #fit to power laws.  Clearly this is not a great thing to do.
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
    #actual correlation length

    AC = binned2[1]
    L = np.sum(AC*db)/AC[0]

if 0:
    #do we want the rectangle?
    rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.2]*3)
    a24.add_patch(rect)

if 0:
    axbonk(a24,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
           #yscale='log',xscale='log')
           yscale='linear',xscale='linear')
    #a24.set_yscale('log')
    #x0a24.set_xscale('log')
    

if 0:
    a24.legend(loc=0)

fig2.savefig('plots_to_sort/p56_convolve_4.png')
fig3.savefig('plots_to_sort/p56_convolve_5.png')
print('saved')

