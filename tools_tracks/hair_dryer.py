from starter2 import *
class hair_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.name = self.this_looper.out_prefix
        self.slopes = []

    def run(self,core_list=None,do_plots=True, frame=0, colors=None, newplots=True, name=''):
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

        fig,ax=plt.subplots(1,1, figsize=(12,8))
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue
            self.cores_used.append(core_id)
            if 1:
                if newplots:
                    ax.clear()
                c=[0.5]*3
                if colors is not None:
                    c = colors[core_id]

                skip=1
                if ms.nparticles > 1000:
                    skip = 10
                for ip in range(0,ms.nparticles,skip):
                    ax.plot(   ms.particle_y[ip,:], ms.particle_z[ip,:], c=c, linewidth=0.1)
                ax.scatter(ms.particle_y[:,0],  ms.particle_z[:,0], c=[c]*ms.nparticles, s=0.3)
                ax.scatter(ms.particle_y[:,-1], ms.particle_z[:,-1], c='r', s=0.3)
                ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))
                outname = 'plots_to_sort/%s_blowing_hair_%sc%04d.png'%(self.name,name,core_id)
                fig.savefig(outname)
                print(outname)
            if 0:
                cy = np.tile( ms.particle_y.mean(axis=0), (ms.nparticles,1))
                cz = np.tile( ms.particle_z.mean(axis=0), (ms.nparticles,1))
                dy2 = ms.particle_y - cy
                dz2 = ms.particle_z - cz

                if 0:
                    theta = np.arctan2(dz2,dy2)
                    r = np.sqrt(dy2**2+dz2**2)

                    the_x = r**0.5*np.cos(theta)
                    the_y = r**0.5*np.sin(theta)

                    dy2 = the_x
                    dz2 = the_y



                fig,ax=plt.subplots(1,1, figsize=(12,8))
                for ip in range(ms.nparticles):
                    ax.plot(   dy2[ip,:], dz2[ip,:], c=[0.5]*3, linewidth=0.1)
                    #ax.plot(   ms.particle_y[ip,:], ms.particle_z[ip,:], c=[0.5]*3, linewidth=0.1)
                ax.scatter(dy2[:,0],  dz2[:,0], c='k', s=0.3)
                ax.scatter(dy2[:,-1], dz2[:,-1], c='r', s=0.3)
                ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))

            if 0:
                #fig,ax=plt.subplots(2,2,figsize=(8,8))
                fig,ax=plt.subplots(1,1, figsize=(12,8))
                #ax.set_aspect('equal')
                centroid_y = ms.mean_yc
                centroid_y.shape = (1,centroid_y.shape[0])
                centroid_z = ms.mean_zc
                centroid_z.shape = (1,centroid_z.shape[0])
                cy = np.tile( ms.particle_y.mean(axis=0), (ms.nparticles,1))
                cz = np.tile( ms.particle_z.mean(axis=0), (ms.nparticles,1))

                dy2 = ms.particle_y -cy
                dz2 = ms.particle_z -cz
                dy2 = ms.particle_y
                dz2 = ms.particle_z

                SL = slice(None)
                for ip in range(ms.nparticles):
                    ax.plot(   dy2, dz2, c=[0.5]*3, linewidth=0.1)
                    continue
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
                ax.scatter(dy2[:,0],  dz2[:,0],  c='k', s=0.3)
                ax.scatter(dy2[:,-1], dz2[:,-1], c='r', s=0.3)
                ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))


                outname = 'plots_to_sort/%s_with_hair_c%04d.png'%(self.name,core_id)
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
