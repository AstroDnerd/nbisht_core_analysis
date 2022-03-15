
'''
mach number vs time
'''
from starter2 import *


class machVsTime():
    def __init__(self,this_looper):
        self.this_looper=this_looper

    def run(self,nt):
        def machNum(ad):
            # with sound speed = 1, becuase pressure is set to 1...
            #machNum_xyz = (np.sqrt(ad['velocity_x']**2 + ad['velocity_y']**2 + ad['velocity_z']**2) * ad['cell_volume']).sum()/ad['cell_volume'].sum()
            # should be the same as... 
            machNum = (ad['velocity_magnitude']*ad['cell_volume']).sum()/ad['cell_volume'].sum()
            print('sim ',nt)
            print('mach num ',machNum)    
            machFile.write("Sim %s MachNum %f, \n"%(nt,machNum)) 

        machFile = open("machRecords.txt",'a')
        thtr = self.this_looper.tr
        for nf,frame in enumerate(thtr.frames):
            ds = self.this_looper.load(frame)
            ad = ds.all_data() 
            machNum(ad)
        machFile.close()

# MAIN
import three_loopers_six as TL6
if 'machTime1' not in dir():
    machTime1 = machVsTime(TL6.loops['u601'])
if 'machTime2' not in dir():
    machTime2 = machVsTime(TL6.loops['u602'])
if 'machTime3' not in dir():
    machTime3 = machVsTime(TL6.loops['u603'])

for nt,tool in enumerate([machTime1,machTime2,machTime3]):
    tool.run(nt)



