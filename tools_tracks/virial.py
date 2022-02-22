


# have not started this one yet.

from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
from collections import defaultdict
reload(dl)
reload(trackage)
plt.close('all')


import alpha_tools
reload(alpha_tools)


if 'do_all_plots' not in dir():
    do_all_plots = False

#import three_loopers as TL
#import three_loopers_1tff as TL
#import three_loopers_mountain_top as TLM
import three_loopers_tenfour as TL4
import sf2
frame=0
sim_list = ['u401']
if 'Lsubtool' not in dir():
