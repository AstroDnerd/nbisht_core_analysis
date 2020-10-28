from starter2 import *
import scipy.signal
import matplotlib.patches as patches
from PIL import Image
plt.close('all')

if 1:

    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    x,y = np.mgrid[0:size:1 , 0:size:1]
    xcen, ycen = 0.5,0.5
    xcen, ycen = size/2.,size/2.
    dv = np.ones_like(x)*dx**3
    rmag = np.sqrt((x-xcen)**2 + (y-ycen)**2)


    #rho = np.array(Image.open("circ.tif")).astype('float')
    rho = np.zeros_like(rmag)
    radius=size/4
    rho[ rmag<radius] = 1


    #AC=scipy.signal.correlate(rho,rho,mode='same',method='fft')
    rhohat = np.fft.fftn(rho)
    rho2 = rhohat*np.conj(rhohat)
    AC1 = np.fft.ifftn(rho2)
    ACc = np.fft.fftshift(AC1)
    AC = np.real(ACc)
    ACB = AC
    AC = AC/AC.max()

    fig,axes = plt.subplots(2,2)
    ax1=axes[0,0]
    ax2=axes[0,1]
    ax3=axes[1,0]
    ax4=axes[1,1]

    if len(rho.shape)== 2:
        im1 = rho
        im2 = AC
    else:
        im1=rho.sum(axis=1)
        im2=AC.sum(axis=1)

    ax1.imshow(im1)
    ax2.imshow(im2)

    #ax1.plot(rmag[33,33:],  AC[33,33:],c='k')
    ACstripe = AC[33,33:]
    ax3.plot(  ACstripe,c='k',label='slice')
    r1 = np.arange(31)/32 #rmag[33,33:]
    y = 2.0*(np.arccos(r1)- r1*np.sqrt(1.0-r1**2))/np.pi;
    #y = radius*(radius*np.arccos(r1/radius)- r1*np.sqrt(1.0-r1**2/radius**2))/np.pi/128
    ax3.plot( y,c='g',label='ann')

    import radial_binner as rb
    ACbins,ACcen,ACvals = rb.rb2(rmag, AC,nbins=33)
    ax3.plot(ACcen,ACvals,c='r', label='binned')
    dr = ACbins[1:]-ACbins[:-1]

    L = np.sum(ACvals*dr)/ACvals[0]
    rect=patches.Rectangle((0,0),L,ACvals[0],facecolor=[0.2]*3)
    ax3.add_patch(rect)
    print("L = %0.2f R = %0.2f"%(L,radius))

    bins,cen,vals = rb.rb2(rmag, rho/rho.max())
    ax3.plot(cen,vals,c='m')

    ax3.legend(loc=0)
    fig.savefig('plots_to_sort/AC1.png')

if 0:
    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    rmag = np.sqrt((x-0.5)**2 + (x-0.5)**2)
    rho = np.zeros_like(x)
    R_sph =  0.25
    rho[ rmag < R_sph] = 1




    fig,ax1 = plt.subplots(1,1)
    #ax1.plot(rmag, AC/AC.size,c='k')
    ax1.plot(cen,bins,c='r')
    ax1.plot(np.linspace(0,1,31),AC[33,33:]/AC.max(),c='g')
    d = cen
    d[ d>R_sph]=R_sph
    theta = 2*np.arccos(d/R_sph)

    A = 2*R_sph**2/2*(theta-np.sin(theta))
    A = 2*R_sph**2*(np.arccos(cen) - cen*np.sqrt(1-cen**2))
    ax1.plot(cen,A,c='g')
    fig.savefig('plots_to_sort/ac_test.png')
    print('wrote plots_to_sort/ac_test.png')

