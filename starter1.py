#
# The typical python packages.
#

import pdb
import sys
import re
import copy
import h5py
import os
import copy
import warnings
import time
import glob
from importlib import reload
#import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colorbar as cb
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab

import astropy.io.fits as pyfits
import platform
ver=platform.python_version()
python_version=  int(ver[0])
if python_version == 3:
    from importlib import reload
import numpy as np
import math
nar = np.array
from scipy.stats import kurtosis
from scipy.stats import skew
from collections import defaultdict
#ef=execfile

x_dict = [1,0,0]
y_dict = [2,2,1]
import IPython
terminal = IPython.get_ipython()
terminal.pt_app.auto_suggest = None
