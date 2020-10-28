import spectra_tools as st
reload(st)
from scipy import fftpack

def make_spectra(directory,frame):
    """This makes 3d power spectra of velocity, acceleration, 
    magnetic field, and density.  Can be very slow.
    FFTs are stored in DD????.products"""
    oober = st.short_oober(directory, frame=frame)
    #st.MakeVelocitySpectra(oober,frame)
    #st.MakeAccelSpectra(oober,frame)
    #st.MakeMagneticSpectra(oober,frame)
    st.MakeDensitySpectra(oober,frame)

def read_spectra(directory,frame):
    """read 3d spectra"""
    spec={}
    #spec['vspec']=dpy( "%s/DD%04d.products/power_velocity.h5"%(directory,frame) , ['k','power'])
    spec['dspec']=dpy( "%s/DD%04d.products/power_density.h5"%(directory,frame) , ['k','power'])
    #spec['hspec']=dpy( "%s/DD%04d.products/power_magnetic.h5"%(directory,frame) , ['k','power'])
    return spec

frame =0
#directory = "/scratch1/dcollins/Paper19/u20_sphere_g1600/"
directory = "/scratch1/dcollins/Paper19/u21_sphere_large"
#directory = "/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential"
#ax_rho_x.plot(the_x[ok], this_dumb_line[ok],'g:')
#make_spectra(directory,frame)
#spec=read_spectra("/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential",frame)
spec=read_spectra(directory,frame)
plt.clf()
plt.close('all')
fig,ax=plt.subplots(2,2)
ax_v=ax[0][0]
ax_h=ax[0][1]
ax_rho_x=ax[1][0]
ax_rho_k=ax[1][1]

if 0:
    the_k=spec['vspec'][0]
    ax_v.plot( the_k, np.abs(spec['vspec'][1]))
    k0 = the_k[5]
    P0 = spec['vspec'][1][5].real

if 0:
    Kolmogorov_ish = P0*(the_k[5:10]/k0)**(-5./3)
    ax_v.plot( the_k[5:10], Kolmogorov_ish)
    axbonk(ax_v,xscale='log',yscale='log',xlabel=r'$k$',ylabel='$P_v(k)$')

if 0:
    the_k=spec['hspec'][0]
    ax_h.plot( the_k, np.abs(spec['hspec'][1]))
    axbonk(ax_h,xscale='log',yscale='log',xlabel='$k$',ylabel='$P_H(k)$')

if 1:
    the_k=spec['dspec'][0]
    ax_rho_k.plot( the_k, np.abs(spec['dspec'][1]))
    k0 = the_k[5]
    #P0 = spec['dspec'][1][5]
    #Kolmogorov_ish = P0*(the_k[5:10]/k0)**(-5./3)
    #ax_v.plot( the_k[5:10], Kolmogorov_ish)
    axbonk(ax_rho_k,xscale='log',yscale='log',xlabel='$k$',ylabel=r'$P_\rho(k)$')

if 1:
    AC = np.fft.irfft(spec['dspec'][1])* 128/(np.pi*2)
    the_x=np.fft.fftfreq(len(AC)) 
    #the_x = np.linspace(0,1,AC.size)
    #AC = AC[:the_k.size]
    dL = the_x[1:]-the_x[:-1]
    ACcen = np.abs(0.5*(AC[1:]+AC[:-1]))
    Length = (dL*ACcen).sum()/ACcen[0]
    print("LENGTH",Length)

    rect = patches.Rectangle((0,0),Length,ACcen[0],linewidth=1,edgecolor='r',facecolor=[1,0,0,0.2])

    # Add the patch to the Axes
    ax_rho_x.add_patch(rect)
    #the_x = np.linspace(0,1,AC.size)
    this_dumb_line = np.abs(AC)
    ok = the_x >= 0
    ax_rho_x.plot(the_x[ok], this_dumb_line[ok])#,'g',marker='*')
    ax_rho_x.plot(the_x[ok], AC.real[ok],'g--')
    ax_rho_x.plot(the_x[ok], AC.imag[ok],'g:')
    #ax_rho_p.plot(the_x, AC.real,'g:')
    #ax_rho_p.plot( the_k, AC.imag,'r')


    R = 0.25
    #h= the_x
    #AC_fit = np.pi/3*h**2*(3*R-h)
    #v = np.pi/3*4*R**3

    h = (R - the_x/2)
    AC_fit = 100*2*np.pi/3*h**2*(3*R-h)
    ax_rho_x.plot(the_x[ok], AC_fit[ok],'r:')
    axbonk(ax_rho_x,xscale='linear',yscale='linear',xlabel='$r$',
           ylabel=r'$\langle \rho(x)\rho(x+r)\rangle$')

    if 0:
        ax_rho_p.clear()
        from scipy.optimize import curve_fit
        def density_ac(r, A, r0, gamma):
            return A*(1+(r/r0)**gamma)
        popt, pcov = curve_fit(density_ac, the_x, AC, p0=[1,1,-2])
        po2 = AC[0]/2, 1, -1
        the_x = np.linspace(0,0.25,100)
        ax_rho_p.plot( the_x, density_ac(the_x, po2[0],po2[1],po2[2]),c='k')


    fig.savefig('plots_to_sort/spectra_sphere_n%04d.pdf'%frame)

