from starter2 import *
from collections import defaultdict
import scipy
import colors

class hair_3d_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.name = self.this_looper.out_prefix
        self.slopes = []

    def run(self,core_list=None,do_plots=True, frame=0, colors=None, newplots=True, rotato=False):
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

        frame_ind = np.where(thtr.frames == frame)[0][0]

        #fig,ax=plt.subplots(1,1, projection='3d')
        fig = plt.figure(figsize=(12,12))
        ax = plt.axes(projection='3d')
        fig2,ax2=plt.subplots(1,1)
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)
            self.ms = ms

            if ms.nparticles < 10:
                continue
            self.cores_used.append(core_id)
            if newplots:
                ax.clear()
                c=[0.5]*3
            if colors is not None:
                c = colors[core_id]
            ax.plot3D(ms.mean_xc[:],   ms.mean_yc[:], ms.mean_zc[:], c=c)
            ax2.plot(   ms.mean_yc[:], ms.mean_zc[:], c=c)
            ax.scatter3D(ms.mean_xc[-1], ms.mean_yc[-1], ms.mean_zc[-1], c=[c],s=15)
            ax.text(ms.mean_xc[-1], ms.mean_yc[-1], ms.mean_zc[-1],"%d"%core_id, color=c)
            print("N part", core_id, ms.nparticles)
            size = None
            if ms.nparticles > 1000:
                size = 0.5
            else:
                size = 15
            print("WWWWW",len(ms.particle_x[0,:]))
            ax.scatter3D(ms.particle_x[:,frame_ind],   ms.particle_y[:,frame_ind], ms.particle_z[:,frame_ind], c=[c]*ms.nparticles,s=size)

            #for ip in range(ms.nparticles):
            #    ax.plot3D(ms.particle_x[ip,:],   ms.particle_y[ip,:], ms.particle_z[ip,:], c=c, linewidth=0.1)
            #    ax.scatter(ms.particle_y[:,0],  ms.particle_z[:,0], c=[c]*ms.nparticles, s=0.3)
            #    ax.scatter(ms.particle_y[:,-1], ms.particle_z[:,-1], c='r', s=0.3)
            #    ax.set_title(r'$\rm{%s}\ \rm{core}\ %d$'%(self.name, core_id))
            if newplots:
                outname = 'plots_to_sort/%s_blowing_hair_c%04d.png'%(self.name,core_id)
                fig.savefig(outname)
                print(outname)
        if not newplots:
            outname = 'plots_to_sort/%s_hair_3D_multi_n%04d.png'%(self.name, frame)
            ax.view_init(20,440)
            fig.savefig(outname)
            print(outname)
            fig2.savefig('plots_to_sort/%s_harish_2d.pdf'%self.name)
        if rotato:
            Nrev = 4
            angles = range(0,Nrev*360,10)
            theta_max = 30
            theta_min = 0
            Nangles = len(angles)/2
            for na,angle in enumerate(angles):
                if na < Nangles:
                    theta = theta_min + na * (theta_max-theta_min)/Nangles
                else:
                    theta = theta_max + (na-Nangles) * (theta_min-theta_max)/Nangles
                print(theta)
                ax.view_init(theta,angle)
                outname = "plots_to_sort/%s_hair_3d_centroid_multi_rotat_a%04d.png"%(self.name,na)
                ax.set_title("theta %0.1f phi %0.1f"%(theta,angle))
                fig.savefig(outname)
                print(outname)


if 1:
    import three_loopers_mountain_top as TLM
    tool = hair_3d_tool( TLM.loops['u302'])
    for ns, superset in enumerate(supersets):
        tool.name = 'u401_S%02d'%ns
        core_list = list(superset)
        color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
        rotato = True
        frame_list=[0]
        for frame in frame_list:
            tool.run( core_list = core_list, newplots = False, colors = color_dict,rotato=rotato, frame=frame)

"""
sim_list=['u301','u302','u303']
#sim_list=['u302','u303']
if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()
if 'lylist' not in dir() or True:
    import three_loopers_tenfour as TL4 
    tool = lyapunov_tool( TL4.['u401'])


    frame_list=tool.this_looper.frame_list
    frame_list=[0]
    rotato = True

    for frame in frame_list:
        tool.run( core_list = core_list, newplots = False, colors = color_dict,rotato=rotato, frame=frame)





if 0:
    import three_loopers_mountain_top as TLM
    lylist={}
    #for this_simname in sim_list:
    #    lylist[this_simname]= lyapunov_tool( TLM.loops[this_simname])
    lylist['u301'] = lyapunov_tool( set_looper)

#    for this_simname in  sim_list:
#        lylist[this_simname].run( )#core_list=[10,11])#core_list=[10])

    #lylist['u301'].run( core_list = [0,8,27,32,37,44,84,275])
    #lylist['u301'].run( core_list = [323])
    #lylist['u301'].run( core_list = [24])
    #lylist['u302'].run( core_list = [30,32])
    #lylist['u303'].run( core_list = [233,235])
    #lylist['u303'].run( core_list = [184])
    #lylist['u303'].run( core_list = [186])

    #core_list=[12,13,14,15,16,17,18,19,21] 
    target_core=31
    core_list = nar(ht1.cores_used)[ nar(ht1.overlaps[target_core   ])>0] 
    core_list = np.append(core_list,[target_core])
    popper  = set([124,108, 66,65,62,61])
    core_set = set(core_list)
    core_list = list(core_set.difference(popper))
    reload(colors)
    color_dict = colors.make_core_cmap(core_list, cmap = 'tab20c', seed = -1)
    #lylist['u301'].run( core_list = core_list, newplots = False, colors = color_dict,rotato=False, frame=30)
    frame_list=lylist['u301'].this_looper.frame_list
    frame_list=[0]
    rotato = True
    for frame in frame_list:
        lylist['u301'].run( core_list = core_list, newplots = False, colors = color_dict,rotato=rotato, frame=frame)


"""
