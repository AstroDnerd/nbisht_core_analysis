
'''
making a pdf out of the projected density

'''
from starter2 import *
import data_locations as dl


class laser_nif():
    def __init__(self,this_looper):
        self.this_looper = this_looper

    def laserRun(self):
        thtr = self.this_looper.tr
        theframes = thtr.frames

        #for later
        #for nf,frame in enumerate(theframes):
        #    ds = self.this_looper.load(frame)

        #for now; 
        ds = self.this_looper.load(theframes[7])
        #for now; ProjectionPlot
        if 0:
            proj = yt.ProjectionPlot(ds,'z',('gas','density'))
            proj.save()
        #for now; ds.proj()
        if 1:
            midpt = ds.find_max('density')[1]  #OR _,c = ds.find_max('gas','density')
            width = 1
            rez = [256,256]

            proj = ds.proj(('gas','density'),2)  #there's proj.profile(), and "most obv" proj.to_frb() that may work
            frb = proj.to_frb(width,rez,center=midpt)

            if 0:
                # Explorations with dir(frb)
                #plt.plot(frb['gas','density'],frb['gas','cell_volume'],c='k')  # um, yeah NO
                frb.save_as_dataset()  # saves a .h5 file that could be used for profile
                # NOTE: re explore dir(proj.profile) 

            if 0:
                plt.imshow(np.log10(np.array(frb['gas','density'])),interpolation='nearest',origin='lower')
                plt.savefig('test_forpdf_B.png')

        #to explore/debug:
        if 1:
            breakpoint()  #to do dir(proj)


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
    
if 'laser_nif01' not in dir():
    laser_nif01 = laser_nif(TL6.loops['u601'])


for nt,tool in enumerate([laser_nif01]):
    tool.laserRun()

