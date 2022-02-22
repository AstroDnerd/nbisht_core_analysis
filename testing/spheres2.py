
from starter2 import *

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




if 0:
    Greens = np.fromfile('/home/dcollins4096/PigPen/derp.txt', sep=" ")
    Greens.shape = int(Greens.size/3), 3
    plt.clf()
    plt.scatter(Greens[:,1], -1/Greens[:,0]/16**3)
    plt.savefig('plots_to_sort/green1.png')
if 0:
    Greens = np.fromfile('/home/dcollins4096/PigPen/derp2.txt', sep=" ")
    Greens.shape = int(Greens.size/3), 3
    plt.clf()
    plt.scatter(Greens[:,0], -1./Greens[:,1])
    plt.savefig('plots_to_sort/green2.png')
if 1:
    plt.clf()
    #rhohat = np.fromfile('/home/dcollins4096/PigPen/derp4.txt', sep=" ")
    rhohat = read_output("/scratch1/dcollins/Paper19/u24_point_periodic/dump","OUTREGION")

    rhohat.shape = int(rhohat.size/3), 3
    rhohatc = rhohat[:,1] + 1j*rhohat[:,2]
    rhohatc.shape = (16,16,9)
    plt.imshow( np.abs(rhohatc).sum(axis=0))
    #rhohatc=np.fft.ifftshift(rhohatc)
    rho = np.fft.irfftn(rhohatc)
    #rho = np.fft.ifftshift(rho)
    plt.savefig('plots_to_sort/green3.png')
    plt.clf()
    plt.imshow( rho.sum(axis=0))
    plt.savefig('plots_to_sort/green4.png')

    plt.clf()
    poot=plt.imshow(np.abs(rhohatc).sum(axis=0))
    plt.colorbar(poot)
    plt.savefig('plots_to_sort/green6.png')


if 1:
    plt.clf()
    #rhohat = np.fromfile('/home/dcollins4096/PigPen/derp4.txt', sep=" ")
    rhohat = read_output("/scratch1/dcollins/Paper19/u24_point_periodic/dump","INITIALREG")
    rho = rhohat[:,1]
    rho.shape=(16,16,18)
    plt.clf()
    plt.imshow( rho.sum(axis=0))
    plt.savefig('plots_to_sort/green5.png')

    rho = rho[:,:,:-2]

    rhohatb = np.fft.rfftn(rho)
    plt.clf()
    poot=plt.imshow( np.abs(rhohatb).sum(axis=0))
    plt.colorbar(poot)
    plt.savefig('plots_to_sort/green7.png')

    G = read_output("/scratch1/dcollins/Paper19/u24_point_periodic/dump","GREEEEN")
    G.shape = (16,16,9)

    phihat = G*rhohatc
    phi = np.fft.irfftn(phihat)


