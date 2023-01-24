

'''
from the cloud tutorial
'''

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pdb

image = mpimg.imread("clouds.png")
#pdb.set_trace()

npix = image.shape[0]
#image.shape[0] == image.shape[1]

fourier_image = np.fft.fftn(image)
fourier_amplitudes = np.abs(fourier_image)**2

kfreq = np.fft.fftfreq(npix) * npix
kfreq2D = np.meshgrid(kfreq, kfreq)
knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)

knrm = knrm.flatten()
fourier_amplitudes = fourier_amplitudes.flatten()

kbins = np.arange(0.5, npix//2+1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])

Abins, _, _ = stats.binned_statistic(knrm, fourier_amplitudes, \
                                     statistic = "mean", bins = kbins)

Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)


# PLOTTINGS SETUP
if 0: 
    fig,ax = plt.subplots(1,2)
    fig.subplots_adjust(wspace=0, hspace=0)

    ax[0].imshow(image)
    ax[1].loglog(kvals,Abins)
    plt.xlabel('$k$')
    plt.ylabel('$P(k)$')
    plt.tight_layout()
    plt.savefig('greyscale_powerspec.png')

if 0:  #this works but is not right according to equations 
    image_sq = ((image **2) * image.size).sum() 
    ftt_sq = (np.abs(fourier_image) **2).sum()
    ratio = image_sq/ftt_sq


