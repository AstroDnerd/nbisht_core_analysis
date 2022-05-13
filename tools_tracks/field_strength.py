from starter2 import *
from collections import defaultdict
import three_loopers_u500 as TL


class rotoboto():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.q0=defaultdict(list)
        self.q1=defaultdict(list)
        self.q2=defaultdict(list)

    def run(self,core_list=None):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            bx = thtr.c([core_id],'magnetic_field_x')
            by = thtr.c([core_id],'magnetic_field_y')
            bz = thtr.c([core_id],'magnetic_field_z')
            vx = thtr.c([core_id],'velocity_x')
            vy = thtr.c([core_id],'velocity_y')
            vz = thtr.c([core_id],'velocity_z')

            bt = np.sqrt(bx*bx+by*by+bz*bz)
            vt = np.sqrt(vx*vx+vy*vy+vz*vz)

            #0,1,2 are aligned with the field and velocity.
            #x,y,z are aligned with the grid
            h0x = bx/bt
            h0y = by/bt
            h0z = bz/bt
            
            vhatx = vx/vt
            vhaty = vy/vt
            vhatz = vz/vt

            h2x = h0y*vhatz - h0z*vhaty
            h2y = h0z*vhatx - h0x*vhatz
            h2z = h0x*vhaty - h0y*vhatx

            h1x = h2y*h0z - h2z*h0y
            h1y = h2z*h0x - h2x*h0z
            h1z = h2x*h0y - h2y*h0x

            B0 = h0x*bx+h0y*by+h0z*bz
            B1 = h1x*bx+h1y*by+h1z*bz
            B2 = h2x*bx+h2y*by+h2z*bz

            print(B1/bt)
            #print(B1)
            #print(B2)




for sim in TL.loops:
    this_looper=TL.loops[sim]
    rb = rotoboto(this_looper)
    rb.run(core_list=[323])
    break

