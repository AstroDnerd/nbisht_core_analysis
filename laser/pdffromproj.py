
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
            midpt = ds.find_max('density')[1]  #OR _,c = ds.find_max('gas','density') ...does not work  #understand the [1]
            midpt_cv = ds.find_max('cell_volume')[1]  #no needd...
            midpt_vz = ds.find_max('velocity_z')[1]  #no needd...
            width = 1
            rez = [256,256]

            proj_d = ds.proj(('gas','density'),2)  #there's proj.profile(), and "most obv" proj.to_frb() that may work
            proj_vz = ds.proj(('gas','velocity_z'),2)
            proj_cv = ds.proj(('gas','cell_volume'),2)       #there's no apparent difference in doing this!
            frb_d = proj_d.to_frb(width,rez,center=midpt)
            frb_vz = proj_vz.to_frb(width,rez,center=midpt)
            frb_cv = proj_cv.to_frb(width,rez,center=midpt_cv)  #and therefore this

            if 0:
                # Explorations with dir(frb)
                #plt.plot(frb['gas','density'],frb['gas','cell_volume'],c='k')  # um, yeah NO
                frb.save_as_dataset()  # saves a .h5 file that could be used for profile
                # NOTE: re explore dir(proj.profile) 

            if 0:
                frb = frb_d
                plt.imshow(np.log10(np.array(frb['gas','density'])),interpolation='nearest',origin='lower')
                plt.savefig('test_pdfRho.png')
            if 1: #idea
                frb = frb_d
                the_x = np.log10(np.array(frb['gas','density']))  #simply inputing 'velocity_z' does not work, TRY YTOrthoRay NEXT!! 
                #the_bins = # will find out 
                the_weight = np.log10(np.array(frb['gas','cell_volume']))  #or do another frb for cell_volume?  NOPE
            
                the_array, xbins = np.histogram(the_x, weights = None, density=True) 
                bin_centers = 0.5*(xbins[1:]+xbins[:-1])
                plot_x = bin_centers
                plot_y = the_array
                plt.plot(plot_x,plot_y, c='k')     
                plt.savefig('pdf_frame7.png')


        #to explore/debug:
        if 0:
            breakpoint()  #do dir(proj) ..etc


# MAIN
import three_loopers_six as TL6
if 'clobber' not in dir():
    clobber=True
    
if 'laser_nif01' not in dir():
    laser_nif01 = laser_nif(TL6.loops['u601'])


for nt,tool in enumerate([laser_nif01]):
    tool.laserRun()

