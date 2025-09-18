from starter2 import *


def quantized_mask(this_looper):
    #Takes only those separated by dt
    #where dt is the first nonzero frame size.
    #essentially assume that the frames are mostly separated
    #by dt with some outliers.
    times = this_looper.tr.times
    dt = times[ times>0][0]
    nframe = times/dt
    outliers = np.round(nframe) - nframe
    frame_mask = np.abs(outliers) < 1e-3
    return frame_mask

if 0:
    #import three_loopers_u500 as TL5
    loop = TL5.loops['u503']
    times = loop.tr.times
    plt.clf()
    a=times/times[2]
    plt.plot(a-np.round(a), marker='*')
    plt.savefig('plots_to_sort/times.png')
