import spectra_tools as st
def make_spectra(directory,frame):
    """This makes 3d power spectra of velocity, acceleration, 
    magnetic field, and density.  Can be very slow.
    FFTs are stored in DD????.products"""
    oober = st.short_oober(directory, frame=frame)
    st.MakeVelocitySpectra(oober,frame)
    #st.MakeAccelSpectra(oober,frame)
    st.MakeMagneticSpectra(oober,frame)
    st.MakeDensitySpectra(oober,frame)

def read_spectra(directory,frame):
    """read 3d spectra"""
    spec={}
    spec['vspec']=dpy( "%s/DD%04d.products/power_velocity.h5"%(directory,frame) , ['k','power'])
    spec['dspec']=dpy( "%s/DD%04d.products/power_density.h5"%(directory,frame) , ['k','power'])
    spec['hspec']=dpy( "%s/DD%04d.products/power_magnetic.h5"%(directory,frame) , ['k','power'])
    return spec

#make_spectra("/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential",0)
spec=read_spectra("/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential",0)
plt.clf()
plt.close('all')
fig,ax=plt.subplots(2,2)
ax0=ax[0][0]
ax1=ax[1][0]
ax2=ax[0][1]
ax3=ax[1][1]

the_k=spec['vspec'][0]
ax0.plot( the_k, np.abs(spec['vspec'][1]))
k0 = the_k[5]
P0 = spec['vspec'][1][5].real
if 1:
    Kolmogorov_ish = P0*(the_k[5:10]/k0)**(-5./3)
    ax0.plot( the_k[5:10], Kolmogorov_ish)
    axbonk(ax0,xscale='log',yscale='log',xlabel=r'$k$',ylabel='$P_v(k)$')
    the_k=spec['hspec'][0]
    ax3.plot( the_k, np.abs(spec['hspec'][1]))
    axbonk(ax3,xscale='log',yscale='log',xlabel='$k$',ylabel='$P_H(k)$')

    the_k=spec['dspec'][0]
    ax1.plot( the_k, np.abs(spec['dspec'][1]))
    k0 = the_k[5]
    #P0 = spec['dspec'][1][5]
    #Kolmogorov_ish = P0*(the_k[5:10]/k0)**(-5./3)
    #ax0.plot( the_k[5:10], Kolmogorov_ish)
    axbonk(ax1,xscale='log',yscale='log',xlabel='$k$',ylabel=r'$P_\rho(k)$')

    AC = np.fft.irfft(spec['dspec'][1])*128/(np.pi*2)
    the_k=spec['dspec'][0]
    AC = AC[:the_k.size]
    dL = the_k[1:]-the_k[:-1]
    ACcen = np.abs(0.5*(AC[1:]+AC[:-1]))
    Length = (dL*ACcen).sum()/ACcen[0]
    print("LENGTH",Length)

    rect = patches.Rectangle((0,0),Length,ACcen[0],linewidth=1,edgecolor='r',facecolor='r')

    # Add the patch to the Axes
    ax2.add_patch(rect)
    the_x = np.linspace(0,1,AC.size)
    this_dumb_line = np.abs(AC)
    ax2.plot(the_x, this_dumb_line,'g:')
    #ax2.plot(the_x, AC.real,'g:')
    #ax2.plot( the_k, AC.imag,'r')
    axbonk(ax2,xscale='log',yscale='linear',xlabel='$r$',ylabel=r'$\langle \rho(x)\rho(x+r)\rangle$')


    fig.savefig('plots_to_sort/spectra.pdf')

    plt.close('fig')
