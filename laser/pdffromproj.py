
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
        ad = ds.all_data()
        #for now; ProjectionPlot
        if 0:
            proj = yt.ProjectionPlot(ds,'z',('gas','density'))
            proj.save()
        #for now; ds.proj()
        if 1:
            midpt = ds.find_max('density')[1]  #OR _,c = ds.find_max('gas','density')
            width = 1
            rez = [256,256]

            proj_d = ds.proj(('gas','density'),2)  #there's proj.profile(), and "most obv" proj.to_frb() that may work
            proj_cv = ds.proj(('gas','cell_volume'),2)
            frb_d = proj_d.to_frb(width,rez,center=midpt)
            frb_cv = proj_cv.to_frb(width,rez,center=midpt)  #there's no apparent difference in doing this!

            if 0:
                # Explorations with dir(frb)
                #plt.plot(frb['gas','density'],frb['gas','cell_volume'],c='k')  # um, yeah NO
                frb.save_as_dataset()  # saves a .h5 file that could be used for profile
                # NOTE: re explore dir(proj.profile) 

            if 0:
                plt.imshow(np.log10(np.array(frb['gas','density'])),interpolation='nearest',origin='lower')
                plt.savefig('test_forpdf_B.png')
            if 1: #idea
                the_x = np.log10(np.array(frb_d['gas','density']))
                #the_bins = # will find out 
                the_weight = np.log10(np.array(frb_cv['gas','cell_volume']))  #or do another frb for cell_volume?
            
                the_array, xbins = np.histogram(the_x, weights = the_weight, density=True) 
                bin_centers = 0.5*(xbins[1:]+xbins[:-1])
                plot_x = bin_centers
                plot_y = the_array
                plt.plot(plot_x,plot_y, c='k')     
                plt.savefig('test_forpdf_10_diffFRBS.png')




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

