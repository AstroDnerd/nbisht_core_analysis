
'''
mach number vs time
'''
from starter2 import *


import three_loopers_six as TL6

frame = 60
ds = TL6.loops['u601'].load(frame)
ad = ds.all_data()



# ok, but get a MACH NUMBER for ALL FRAMES
#for nf,frame in enumerate(thtr.frames):

