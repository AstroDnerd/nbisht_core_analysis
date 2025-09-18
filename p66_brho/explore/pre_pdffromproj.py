
'''
making a pdf out of the projected density

'''
from starter2 import *


#class laser_nif():
#    def __init__(self,this_looper):
#        self.this_looper = this_looper



# MAIN
import three_loopers_six as TL6
#if 'laser_nif01' not in dir():
#    laser_nif01 = laser_nif(TL6.loops['u601'])

frame = 60
#ds = laser_nif01.load(frame)
# NOTE: even shorter:
ds = TL6.loops['u601'].load(frame)

proj = ds.ProjectionPlot(ds,'z',('gas','density'))
#proj.save()



