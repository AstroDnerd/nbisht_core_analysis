
'''
a simple signal to check
'''

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import numpy as np
import pdb


x,y = np.mgrid[0:1:1/16,0:1:1/16]
B = np.sin(2*np.pi*x) + np.sin(8*np.pi*y)
C = np.fft.fftn(B)
A = np.abs(C)**2 
A = A.flatten()

npix = B.shape[0]
kfreq = np.fft.fftfreq(npix) * npix
kfreq2D = np.meshgrid(kfreq, kfreq)
knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
knrm = knrm.flatten()
kbins = np.arange(0.5, npix//2+1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])

Abins, _, _ = stats.binned_statistic(knrm, A, \
                                     statistic = "mean", bins = kbins)
Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)


# PLOTTINGS SETUP
if 0: 
    fig,ax = plt.subplots(1,2)
    fig.subplots_adjust(wspace=0, hspace=0)

    ax[0].imshow(B)
    ax[1].loglog(kvals,Abins)
    plt.xlabel('$k$')
    plt.ylabel('$P(k)$')
    plt.tight_layout()
    plt.savefig('asignal.png')
    pdb.set_trace()
# TEST PASEVALLS 
if 1:
    dx = 1/B.size
    dk = 1/C.size  #same as dx
    fx_sq = ((B**2)).sum()  #ratio is 1 without dx here... 
                            #from numpy dft def "n" from "1/n" cancels the "dx", x being 1/n
    fk_sq = ((np.abs(C)**2) * dk).sum()
    ratio = fx_sq/fk_sq
    pdb.set_trace()


