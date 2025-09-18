
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

from collections import defaultdict
import r_inflection
reload(r_inflection)
sim_list=['u601']
import three_loopers_six as TL
for sim in sim_list:
    inflection[sim]=r_inflection.R_INFLECTION( TL.loops[sim])
inflection[sim].run( core_list=[1], do_plots=True)
