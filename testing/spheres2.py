"""
Plots both the code and what comes out of a live run.
"""
from starter2 import *
import xtra_energy

def read_output(fname,prefix):
    fptr=open(fname)
    lines=fptr.readlines()
    fptr.close()
    output = []
    for line in lines:
        if line.startswith(prefix):
            float_line = nar(line.split()[1:], dtype='float')
            output.append(float_line)
    output = nar(output,dtype='float')
    return output


frame=1
dire = '/data/cb1/Projects/P19_CoreSimulations/new_sims/u56c_little_big_horse/'
DUMP = dire+'/dump'
fname = dire+'/DD%04d/data%04d'%(frame,frame)
ds = yt.load(fname)
xtra_energy.add_energies(ds)
cg = ds.covering_grid(0,ds.domain_left_edge, ds.domain_dimensions,fields = ['density', YT_grav_energy])

if 0:
    #projections.
    plt.clf()

    rho = read_output(DUMP,"INITIALR")
    rho.shape=(16,16,18)
    plt.clf()
    ploot=plt.imshow( rho.sum(axis=0))
    plt.colorbar(ploot)
    plt.title('density, live')
    plt.savefig('plots_to_sort/grav_1.png')

    plt.clf()
    ploot=plt.imshow( cg['density'].sum(axis=0))
    plt.colorbar(ploot)
    plt.title('density, disk')
    plt.savefig('plots_to_sort/grav_2.png')

if 0:
    #The fft from the live run, 
    #then ifft to make rho again.
    plt.clf()
    #rhohat = np.fromfile('/home/dcollins4096/PigPen/derp4.txt', sep=" ")
    rhohat = read_output(DUMP,"OUTREGION")
    rhohatc = rhohat[:,0] + 1j*rhohat[:,1]
    rhohatc.shape = (16,16,9)
    ploot=plt.imshow( np.abs(rhohatc).sum(axis=0))
    plt.colorbar(ploot)
    plt.title( 'rhohat, live')
    plt.savefig('plots_to_sort/grav_3.png')

    plt.clf()
    rho_inv = np.fft.irfftn(rhohatc)
    ploot=plt.imshow( rho_inv.sum(axis=0))
    plt.colorbar(ploot)
    plt.title( 'rho, fft from live')
    plt.savefig('plots_to_sort/grav_4.png')

if 0:
    #FFT works
    rhohat_me = np.fft.rfftn(rho)
    plt.clf()
    ploot=plt.imshow( np.abs(rhohat_me).sum(axis=0))
    plt.colorbar(ploot)
    plt.title('fft rho')
    plt.savefig('plots_to_sort/grav_5.png')

if 1:
    #phi hat live
    phihat_2 = read_output(DUMP,"PHIHAT")
    phihat = phihat_2[:,0]+1j*phihat_2[:,1]

    phihat.shape = (16,16,9)
    plt.clf()
    ploot=plt.imshow( np.abs(phihat).sum(axis=0))
    plt.colorbar(ploot)
    plt.title( 'PhiHat live')
    plt.savefig('plots_to_sort/grav_6.png')

if 0:
    #enzo green funciton
    green = read_output(DUMP,"GREEEENf")
    green.shape = (16,16,9)
    plt.clf()
    plt.title('green live')
    ploot=plt.imshow( np.abs(green).sum(axis=0))
    plt.colorbar(ploot)
    plt.savefig('plots_to_sort/grav_7.png')

if 1:
    #green me
    nxy = rhohat_me.shape
    khat = np.mgrid[-nxy[0]/2:nxy[0]/2:1, -nxy[1]/2:nxy[1]/2:1, 0:nxy[2]:1]*np.pi*2
    k2=khat[0,...]**2+khat[1,...]**2+khat[2,...]**2
    k2 = np.fft.fftshift(k2,axes=0)
    k2 = np.fft.fftshift(k2,axes=1)
    k2inv = np.zeros_like(k2)
    k2inv[ k2>0] = 1/k2[k2>0]
    k2inv *= -ds.parameters['GravitationalConstant']
    plt.clf()
    ploot=plt.imshow( np.abs(k2inv).sum(axis=0))
    plt.colorbar(ploot)
    plt.title( 'green manual ')
    plt.savefig('plots_to_sort/grav_8.png')

if 0:
    #test that I understand the greens function
    green1 = read_output(DUMP,"GREEN1")
    gx = green1[:,0]*2*np.pi
    gy = green1[:,1]*2*np.pi
    gz = green1[:,2]*2*np.pi
    gg = green1[:,3]
    g2 = gx**2+gy**2+gz**2
    plt.clf()
    #plt.scatter( g2, -1/gg)
    plt.scatter( gg, -gg/k2inv.flatten())
    plt.savefig('plots_to_sort/grav_temp.png')

if 1:
    #phihat me
    phihat_me = rhohat_me*k2inv
    fig,ax=plt.subplots(1,2)
    ploot=ax[0].imshow( np.abs(phihat_me).sum(axis=0))
    ax[0].set_title('phihat me')
    fig.colorbar(ploot,ax=ax[0])
    ploot=ax[1].imshow( np.abs(phihat).sum(axis=0))
    ax[1].set_title('phihat live')
    fig.colorbar(ploot,ax=ax[1])
    fig.savefig('plots_to_sort/grav_9.png')

if 1:
    phi_mi = np.fft.irfftn(phihat_me)
    fig,ax=plt.subplots(1,3, figsize=(18,4))
    ploot=ax[0].imshow( np.abs(phi_mi).sum(axis=0))
    ax[0].set_title('phi me')
    fig.colorbar(ploot,ax=ax[0])

    phi_disk = cg[YT_potential_field]
    ploot=ax[1].imshow( np.abs(phi_disk).sum(axis=0))
    ax[1].set_title('phi disk')
    fig.colorbar(ploot,ax=ax[1])

    ax[2].scatter( phi_disk, phi_mi)

    fig.savefig('plots_to_sort/grav_10.png')
