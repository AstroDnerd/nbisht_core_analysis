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
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue
            self.cores_used.append(core_id)
            if 0:
                fig,ax=plt.subplots(1,1, figsize=(12,8))
                for ip in range(ms.nparticles):
                    ax.plot(   ms.particle_y[ip,:], ms.particle_z[ip,:], c=[0.5]*3, linewidth=0.1)
                ax.scatter(ms.particle_y[:,0],  ms.particle_z[:,0], c='k', s=0.3)
                ax.scatter(ms.particle_y[:,-1], ms.particle_z[:,-1], c='r', s=0.3)
                ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))
                outname = 'plots_to_sort/%s_blowing_hair_c%04d.png'%(self.name,core_id)
                fig.savefig(outname)
                print(outname)
                plt.close(fig)
            if 1:
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


            if 1:
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


sim_list=['u301','u302','u303']
#sim_list=['u302','u303']
if 'lylist' not in dir() or True:
    import three_loopers_mountain_top as TLM
    lylist={}
    for this_simname in sim_list:
        lylist[this_simname]= lyapunov_tool( TLM.loops[this_simname])

#    for this_simname in  sim_list:
#        lylist[this_simname].run( )#core_list=[10,11])#core_list=[10])

    #lylist['u301'].run( core_list = [0,8,27,32,37,44,84,275])
    #lylist['u301'].run( core_list = [323])
    #lylist['u301'].run( core_list = [24])
    #lylist['u302'].run( core_list = [30,32])
    #lylist['u303'].run( core_list = [233,235])
    #lylist['u303'].run( core_list = [184])
    lylist['u303'].run( core_list = [186])


if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()

    bad = np.where(thtr.track_dict['density'] <= 0)
    bad_pids = thtr.particle_ids[bad[0]]
    for bad_id in bad_pids:
        print(" strip bad particle ", bad_id)
        n_bad_densities= (thtr.p([bad_id],'density')  <= 0).sum()
        if n_bad_densities == 0:
            print("No bad densities.")
            raise
        particle_index = np.where( thtr.particle_ids == bad_id)
        for field in thtr.track_dict.keys():
            arr = thtr.track_dict[field]
            smaller = np.delete( arr, particle_index,axis=0)
            thtr.track_dict[field]=smaller
        thtr.particle_ids =  np.delete(thtr.particle_ids, particle_index)
        thtr.core_ids =  np.delete(thtr.core_ids, particle_index)




stillbad = np.where(thtr.track_dict['density'] <= 0)
print("STILL BAD", stillbad)

if 0:
    long_tool = lyapunov_tool( set_looper)
    long_tool.run( core_list=[24])#,24,184])


