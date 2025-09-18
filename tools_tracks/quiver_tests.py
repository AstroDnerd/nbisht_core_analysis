
from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')


class trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]



    def run(self,do_all_plots=True,core_list=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        fig,ax=plt.subplots(1,1)
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 3:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0
            ax.clear()
            rmin = 1./2048
            for nt, time in enumerate(tsorted):
                density = thtr.c([core_id],'density')[:,nt]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nt]

                this_r = ms.r[:,nt]
                this_r[ this_r < rmin] = rmin
                asort = np.argsort(this_r)
                unsort = np.argsort(asort)
                rsort = this_r[asort]
                vr = ms.vr_rel[:,nt]
                vrs = vr[asort]
                dv = cell_volume[asort]
                vmean = np.sqrt(np.cumsum(vrs**2*dv)/np.cumsum(dv))
                if vmean[0]>1:
                    plt.clf()
                    #plt.quiver( ms.raw_x[:,nt], ms.raw_y[:,nt], ms.vr_x[:,nt],ms.vr_y[:,nt])
                    #plt.quiver( ms.this_x[:,nt], ms.this_y[:,nt], ms.vr_x[:,nt],ms.vr_y[:,nt])
                    #plt.quiver( ms.this_x[:,nt], ms.this_y[:,nt], ms.rx_hat[:,nt],ms.ry_hat[:,nt])
                    
                    #plt.quiver( ms.this_x[:,nt], ms.this_y[:,nt], ms.rx_rel[:,nt],ms.ry_rel[:,nt])
                    x0, y0=ms.mean_x[nt], ms.mean_y[nt]
                    plt.scatter( x0,y0, marker="*",c='g')
                    counter=0
                    #This does not work.
                    def dumb_angle_plotter(x,y,qx,qy):
                        qmag = np.sqrt(qx**2+qy**2)
                        #theta = np.arccos(qx/qmag)
                        theta = np.arctan2(qy,qx)
                        dx = qmag*np.cos(theta)
                        dy = qmag*np.sin(theta)
                        print(dx)
                        x0 = x-0.5*dx
                        x1 = x+0.5*dx
                        y0 = y-0.5*dy
                        y1 = y+0.5*dy
                        plt.plot([x0,x1],[y0,y1],c='r')
                    for tx,ty,tz in zip(ms.this_x[:,nt],ms.this_y[:,nt], ms.this_z[:,nt]):
                        counter += 1
                        plt.plot( [x0, tx], [y0,ty],c='k',lw=0.1)
                        dx = tx-x0; dy = ty-y0
                        r = np.sqrt(dx**2+dy**2)
                        theta = np.arccos(dx/r)
                        #plt.quiver( tx, ty, tx-x0, ty-y0)
                        dumb_angle_plotter(tx,ty, tx-x0, ty-y0)

                        #Py = r*np.sin(theta)
                        #Px = r*np.cos(theta)
                        #x0=0.01
                        #y0=0.015
                        #x1=x0+Px
                        #y1=y0+Py
                        #plt.plot([x0, x1],[y0,y1],c='k',lw=0.1)
                        #plt.quiver( x1, y1, x1-x0, y1-y0)
                        ##for theta in np.arange(0,2*np.pi,0.1):
                        #    r = np.random.random()*0.1+0.2
                    #for theta in np.arange(0,2*np.pi,0.1):
                    #    r = 0.1#np.random.random()*0.1+0.2
                    #    Py = r*np.sin(theta)
                    #    Px = r*np.cos(theta)
                    #    Py2 =2*r*np.sin(theta)
                    #    Px2 =2*r*np.cos(theta)
                    #    x0=0.01
                    #    y0=0.015
                    #    x1=x0+Px
                    #    y1=y0+Py
                    #    x2=x0+Px2
                    #    y2=y0+Py2
                    #    plt.plot([x0, x2],[y0,y2],c='k',lw=0.1)
                    #    #plt.quiver( x1, y1, x1-x0, y1-y0)
                    #    plt.quiver( x1, y1, Px, Py)
                    plt.savefig("%s/test.png"%dl.output_directory)

                    return -1
                ax.plot(rsort,vmean, c=tmap(nt))
            axbonk(ax,xlabel='r',ylabel=r'$\sigma_v(r)$',xscale='log',xlim=[rmin,0.4], ylim=[0,10])
            #ax.set_yscale('symlog',linthresh=10)

            outname = "%s/%s_vmeans_c%04d.png"%(dl.output_directory, self.this_looper.out_prefix, core_id)
            print(outname)
            fig.savefig(outname)


        plt.close(fig)





if 'do_all_plots' not in dir():
    do_all_plots = False


import three_loopers as TL

run1 = trial(TL.looper1)
run1.run()#core_list=[10,11,31])
