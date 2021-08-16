
from starter2 import *
from collections import defaultdict
import scipy
import colors

class lyapunov_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.name = self.this_looper.out_prefix
        self.slopes = []

    def run(self,core_list=None,do_plots=True, frame=0):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores
        if hasattr(core_list,'v'):
            core_list=core_list.v #needs to not have unit.
            core_list=core_list.astype('int')

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        tmap = rainbow_map( len(thtr.frames))
        ds = self.this_looper.load(frame)

        frame_ind = frame
        if frame > 0:
            print("Need to fix the frame index")
            return

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            ms.particle_pos(core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue

            self.cores_used.append(core_id)
            times = thtr.times
            r = ms.r
            this_x = ms.this_x 
            this_y = ms.this_y
            this_z = ms.this_z
            king = np.argmin(r[:,0])

            dx = this_x - this_x[king]
            dy = this_y - this_y[king]
            dz = this_z - this_z[king]
            r = np.sqrt( dx**2+dy**2+dz**2)
            #fig,ax=plt.subplots(1,1)

            #x0 = dx[:,0]
            #tf = times[-1]
            #for myr in r:
            #    ok = myr>0
            #    ax.plot(times[ok], np.log10(myr[ok]))
            #for xonly in dx:
                #XXX = np.abs(xonly)
                #pfit = np.polyfit( times,XXX,1)
                ##x2 = XXX[0] - XXX[0]/tf*times
                #x3 = pfit[0]*times+pfit[1]
                #self.slopes.append(pfit[1])
                #ax.plot(times,x3)
                #ax.plot( np.abs(XXX- x2))

            #outname = 'plots_to_sort/%s_lyapunov_c%04d.png'%(self.name, core_id)
            #fig.savefig(outname)
            #print(outname)
            #plt.close(fig)
            #plt.clf()
            #plt.hist(self.slopes,histtype='step')
            #plt.savefig('plots_to_sort/%s_slopes_c%04d.png'%(self.name,core_id))

            #plt.clf()
            #for ip in range(dx.shape[0]):
            #    dy3 = ms.particle_y[ip,:]-ms.mean_yc
            #    dz3 = ms.particle_z[ip,:]-ms.mean_zc
            #    plt.plot(   dy3, dz3, c=[0.5]*3, linewidth=0.1)
            #    plt.scatter(dy3[0], dz3[0], c='k', s=0.3)
            #plt.savefig('plots_to_sort/%s_hair_c%04d.pdf'%(self.name,core_id))

            if 1:
                #fig,ax=plt.subplots(2,2,figsize=(8,8))
                fig,ax=plt.subplots(1,1)
                #ax.set_aspect('equal')
                endy,endz=ms.particle_y[:,-1].mean(), ms.particle_z[:,-1].mean()
                dy2 = ms.particle_y -endy
                dz2 = ms.particle_z -endz

                SL = slice(None)
                for ip in range(dx.shape[0])[SL]:
                    ax.plot(   ms.particle_y[ip,:], ms.particle_z[ip,:], c=[0.5]*3, linewidth=0.1)

                    continue
                    dy3 = ms.particle_y[ip,:]-ms.mean_yc
                    dz3 = ms.particle_z[ip,:]-ms.mean_zc
                    ax[0][1].plot(   dy3, dz3, c=[0.5]*3, linewidth=0.1)
                    ax[0][1].scatter(dy3[0], dz3[0], c='k', s=0.3)
                    #ax[0].plot(   dy2[ip,:], dz2[ip,:], c=[0.5]*3, linewidth=0.1)
                    #ax[0].scatter(dy2[ip,0], dz2[ip,0], c='k', s=0.3)

                    theta = np.arctan2(dz3,dy3)
                    r = np.sqrt(dy3**2+dz3**2)
                    #r = np.arange(0,5,0.01)
                    #theta=2*np.pi*r
                    log_r = r**0.5# np.log10(10*r +1 )
                    the_x = log_r*np.cos(theta)
                    the_y = log_r*np.sin(theta)
                    #print(np.abs(the_x-dy3).sum())
                    ax[1][0].plot( the_x, the_y, c=[0.5]*3, linewidth=0.1)
                    ax[1][0].scatter( the_x[0], the_y[0], c=[[0.5]*3], s=0.3)
                    ax[1][0].scatter( the_x[-1], the_y[-1], c='g', s=0.3, marker='*')
                ax.scatter(ms.particle_y[:,0],  ms.particle_z[:,0], c='k', s=0.3)
                ax.scatter(ms.particle_y[:,-1], ms.particle_z[:,-1], c='r', s=0.3)
                ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))


                outname = 'plots_to_sort/%s_blowing_hair_c%04d.png'%(self.name,core_id)
                fig.savefig(outname)
                print(outname)
                plt.close(fig)

            #fig,ax=plt.subplots(1,1)
            #for ip in range(dx.shape[0]):
            #    dy3 = ms.particle_y[ip,:]-ms.mean_yc
            #    dz3 = ms.particle_z[ip,:]-ms.mean_zc
            #    theta = np.arctan2(dz3,dy3)
            #    r = np.sqrt(dy3**2+dz3**2)
            #    #r = np.arange(0,5,0.01)
            #    #theta=2*np.pi*r
            #    log_r = r**0.5# np.log10(10*r +1 )
            #    the_x = log_r*np.cos(theta)
            #    the_y = log_r*np.sin(theta)
            #    #print(np.abs(the_x-dy3).sum())
            #    ax.plot( the_x, the_y, c=[0.5]*3, linewidth=0.1)
            #    ax.scatter( the_x[0], the_y[0], c=[[0.5]*3], s=0.3)
            #    ax.scatter( the_x[-1], the_y[-1], c='g', s=0.3, marker='*')
            #    #ax.plot( dy3,dz3,  c=[0.5]*3, linewidth=0.1)
            #    #ax[0].plot(   dy3, dz3, c=[0.5]*3, linewidth=0.1)
            #    #ax[0].scatter(dy3[0], dz3[0], c='k', s=0.3)
            #fig.savefig('plots_to_sort/%s_logradius_2_c%04d.pdf'%(self.name,core_id))

            #fig,ax=plt.subplots(1,subplot_kw={'projection': 'polar'})
            #for ip in [31,32]:#range(dx.shape[0]):
            #    dy3 = ms.particle_y[ip,:]-ms.mean_yc
            #    dz3 = ms.particle_z[ip,:]-ms.mean_zc
            #    theta = np.arctan2(dy3,dz3)
            #    r = np.sqrt(dy3**2+dz3**2)
            #    ax.plot(theta,r,c=[0.5]*3, linewidth=0.1)
            #    #ax[0].plot(   dy3, dz3, c=[0.5]*3, linewidth=0.1)
            #    #ax[0].scatter(dy3[0], dz3[0], c='k', s=0.3)
            #ax.set_rlim(1e-2,0.1)
            #ax.set_rscale('log')
            #ax.set_rlim(1e-2,0.1)

            #fig.savefig('plots_to_sort/%s_logradius_c%04d.pdf'%(self.name,core_id))
            #print('b2')


import three_loopers_mountain_top as TLM
sim_list=['u301','u302','u303']
#sim_list=['u302','u303']
if 'lylist' not in dir() or True:
    lylist={}
    for this_simname in sim_list:
        lylist[this_simname]= lyapunov_tool( TLM.loops[this_simname])

#    for this_simname in  sim_list:
#        lylist[this_simname].run( )#core_list=[10,11])#core_list=[10])

    lylist['u301'].run( core_list = [0,8,27,32,37,44,84,275])
    lylist['u301'].run( core_list = [323])
    lylist['u302'].run( core_list = [30,32])
    lylist['u303'].run( core_list = [233,235])


