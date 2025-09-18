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

            #1,2,3 are aligned with the field and velocity.
            #x,y,z are aligned with the grid
            h1x = bx/bt
            h1y = by/bt
            h1z = bz/bt
            
            vhatx = vx/vt
            vhaty = vy/vt
            vhatz = vz/vt

            h3x = h1y*vhatz - h1z*vhaty
            h3y = h1z*vhatx - h1x*vhatz
            h3z = h1x*vhaty - h1y*vhatx

            h2x = h3y*h1z - h3z*h1y
            h2y = h3z*h1x - h3x*h1z
            h2z = h3x*h1y - h3y*h1x

            B1 = h1x*bx+h1y*by+h1z*bz
            B2 = h2x*bx+h2y*by+h2z*bz
            B3 = h3x*bx+h3y*by+h3z*bz
            #V1 = h1x*vx+h1y*vy+h1z*vz
            #V2 = h2x*vx+h2y*vy+h2z*vz
            #V3 = h3x*vx+h3y*vy+h3z*vz





for sim in TL.loops:
    this_looper=TL.loops[sim]
    rb = rotoboto(this_looper)
    rb.run(core_list=[323])
    break

