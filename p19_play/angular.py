
from starter2 import *
import other_scrubber
reload(other_scrubber)


class angular():
    def __init__(self,loop):
        self.loop=loop

    def run(self,core_list=None, times_tsung=None,tsing_tool=None, rfix=None, do_plots=True):
        sim_name = self.loop.sim_name
        thtr=self.loop.tr
        if core_list is None:
            core_list=np.unique(loop.tr.core_ids)

        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        self.pjx=[]
        self.pjy=[]
        self.pjz=[]
        self.gjx=[]
        self.gjy=[]
        self.gjz=[]
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms.momenta()

            tsung = tsing_tool.tend_core[core_id]
            tsing = tsing_tool.tsing_core[core_id]

            if do_plots:
                fig,axes=plt.subplots(2,2,figsize=(12,12))
                ax0=axes[0][0];ax1=axes[0][1]
                ax2=axes[1][0];ax3=axes[1][1]
            for ntime,time in enumerate([tsing,tsung]):
                nf = get_time_index(time)
                frame = thtr.frames[nf]
                ds = self.loop.load(frame)
                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                msR = ms.rc
                MaxRadius=msR[:,nf].max()
                #Radius = max([0.5/colors.length_units_pc, MaxRadius])
                Radius = MaxRadius
                print("RADIUS",Radius, "MAX", MaxRadius)
                if rfix is not None:
                    Radius=rfix

                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)
                if 1:
                    proj = ds.proj('density',0,center=center,data_source=sph)
                    pw = proj.to_pw()
                    pw.save('%s/angular_core_%s_c%04d_n%d'%(plot_dir,self.loop.sim_name,core_id,ntime))


                vel = []
                if 1:
                    for axis in 'xyz':
                        vel.append( (sph['velocity_%s'%axis]*sph['cell_mass']).sum()/sph['cell_mass'].sum())
                        #vel.append( sp['velocity_%s'%axis][ORDER][:10].mean())
                        #vel.append( sp['velocity_%s'%axis][ORDER].mean())
                        #vel.append( (rho_sort*sp['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/M_cuml)
                        #vel.append( (rho_sort*sph['velocity_%s'%axis][ORDER]*dv_sort)[:30].sum()/(rho_sort*dv_sort)[:30].sum())
                        #vel.append(0)
                scrub = other_scrubber.scrubber(sph, reference_velocity = vel)
                scrub.momenta()

                self.pjx.append((ms.mass*ms.angular_momentum_rel_x)[:,nf].sum()/ms.mass[:,nf].sum())
                self.pjy.append((ms.mass*ms.angular_momentum_rel_y)[:,nf].sum()/ms.mass[:,nf].sum())
                self.pjz.append((ms.mass*ms.angular_momentum_rel_z)[:,nf].sum()/ms.mass[:,nf].sum())
                self.gjx.append((scrub.mass*scrub.angular_momentum_rel_x).sum()/scrub.mass.sum())
                self.gjy.append((scrub.mass*scrub.angular_momentum_rel_y).sum()/scrub.mass.sum())
                self.gjz.append((scrub.mass*scrub.angular_momentum_rel_z).sum()/scrub.mass.sum())
                pjtheta=ms.j_hat_theta[:,nf]
                pjphi=ms.j_hat_phi[:,nf]
                gjtheta = scrub.j_hat_theta
                gjphi = scrub.j_hat_phi
                if do_plots:
                    axes[ntime][0].hist( pjtheta )
                    axes[ntime][0].hist( pjphi )
                    axes[ntime][1].hist( gjtheta)
                    axes[ntime][1].hist( gjphi)
            if do_plots:
                fig.savefig('%s/angular_singsung_%s_c%04d.pdf'%(plot_dir,self.loop.sim_name,core_id))
                plt.close(fig)

                








                


import track_loader as TL
sims=['u502']
TL.load_tracks(sims)
import tsing
tsing_tool = tsing.get_tsing(TL.tracks)

if 'thing' not in dir():
    for sim in sims:
        core_list=[74]
        core_list=TL.loops[sim].core_by_mode['Alone']
        track=TL.tracks[sim]
        thing=angular(track)
        thing.run(tsing_tool=tsing_tool[sim],core_list=core_list)

if 1:
    pjx = nar(thing.pjx)
    pjy = nar(thing.pjy)
    pjz = nar(thing.pjz)
    pjm = np.sqrt(pjx**2+pjy**2+pjz**2)
    gjx = nar(thing.gjx)
    gjy = nar(thing.gjy)
    gjz = nar(thing.gjz)
    gjm = np.sqrt(gjx**2+gjy**2+gjz**2)
    p_theta = np.arctan2(pjy,pjx)
    p_phi   = np.arccos(pjz/pjm)
    g_theta = np.arctan2(gjy,gjx)
    g_phi   = np.arccos(gjz/gjm)

    fig,axes=plt.subplots(2,2,figsize=(12,12))
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]

    kwargs={'cumulative':True,'histtype':'step'}
    def cuml(arr,ax):
        x = sorted(arr)
        y = np.arange(0,1,1/len(x))

        ax.plot(x,y)
    cuml( p_theta[::2] , ax0)
    cuml( p_theta[1::2], ax1)
    cuml( p_phi[::2]   , ax0)
    cuml( p_phi[1::2]  , ax1)
    cuml( g_theta[::2] , ax2)
    cuml( g_theta[1::2], ax3)
    cuml( g_phi[::2]   , ax2)
    cuml( g_phi[1::2]  , ax3)
    fig.savefig('%s/angular_by_core'%plot_dir)



