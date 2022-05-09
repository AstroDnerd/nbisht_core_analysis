
from starter2 import *
import three_loopers_six as TLM
import xtra_energy
this_looper = TLM.loops['u601']

from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

import colors
plt.close('all')


class potentool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
    def run(self,core_list=None,frames=None):
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

        if frames is None:
            frames = thtr.frames
        if frames == 'reg':
            times = thtr.times
            #assume the first nonzero time is dt_sim
            dt = times[ times>0][0]
            nframe = times/dt
            outliers = np.round(nframe) - nframe
            frame_mask = np.abs(outliers) < 1e-3
            frames = thtr.frames[frame_mask]
            times = thtr.times[frame_mask]
            self.frames=frames
            self.times=times

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('go ', core_id)

            this_center=ms.mean_center[:,-1]
            self.cores_used.append(core_id)
            frame=self.this_looper.target_frame

            ds = self.this_looper.load(frame)
            xtra_energy.add_energies(ds)

            ms.particle_pos(core_id)
            pos = np.stack([ ms.particle_x.transpose(),ms.particle_y.transpose(), ms.particle_z.transpose()])

            sphere = ds.sphere(center= this_center, radius=0.25)
            name = self.this_looper.sim_name
            if 0:
                proj = ds.proj(YT_grav_energy,0,data_source=sphere)
                pw=proj.to_pw(center=this_center,origin='domain')
                pw.set_cmap(YT_grav_energy,'Greys')
                pw.zoom(4)
                pw.annotate_sphere( this_center, ds.arr(7e-3,'code_length'))
                pw.save('plots_to_sort/isohol_%s_n%04d_'%(name, frame))

            if 0:
                proj = ds.proj('density',0,data_source=sphere)
                pw=proj.to_pw(center=this_center,origin='domain')
                pw.zoom(4)
                pw.set_cmap('density','Greys')
                pw.annotate_sphere( this_center, ds.arr(7e-3,'code_length'))
                pw.save('plots_to_sort/isohol_%s_n%04d_'%(name, frame))

            if 0:
                plot = yt.PhasePlot( sphere, YT_density, YT_cell_volume, [YT_cell_mass], weight_field=None)
                plot.save('plots_to_sort/isorho_%s_n%04d'%(name, frame))
                pdb.set_trace()
            if 0:
                fig,ax=plt.subplots(1,1)
                ax.scatter( sphere[YT_radius], sphere[YT_density])
                ax.scatter( ms.r[:,-1], ms.density[:,-1], facecolor='None', edgecolor='r')
                axbonk(ax,xscale='log',yscale='log')
                fig.savefig('plots_to_sort/density_radius_%s_c%04d.png'%(name,core_id))
            if 0:
                fig,ax=plt.subplots(1,1)
                ax.scatter( sphere[YT_radius], -sphere[YT_grav_energy])
                axbonk(ax,xscale='log',yscale='log',xlabel='r',ylabel='|grad phi|')
                fig.savefig('plots_to_sort/potential_radius_%s_c%04d.png'%(name,core_id))
            if 0:
                fig,ax=plt.subplots(1,1)
                ax.scatter( sphere[YT_radius], sphere[YT_kinetic_energy])
                axbonk(ax,xscale='log',yscale='log',xlabel='r',ylabel='|grad phi|')
                fig.savefig('plots_to_sort/KE_radius.png_%s_c%04d.png'%(name,core_id))
            if 0:
                fig,ax=plt.subplots(1,1)
                ok = sphere[YT_density]>300
                vx = sphere[YT_velocity_x]
                vy = sphere[YT_velocity_y]
                vz = sphere[YT_velocity_z]
                dv = sphere[YT_cell_volume]
                rho= sphere[YT_density]
                m = (rho[ok]*dv[ok]).sum()
                vx_bar = (vx*dv*rho)[ok].sum()/m
                vy_bar = (vy*dv*rho)[ok].sum()/m
                vz_bar = (vz*dv*rho)[ok].sum()/m
                vbar = 0.5*( (vx-vx_bar)**2+(vy-vy_bar)**2+(vz-vz_bar)**2)
                ax.scatter( sphere[YT_radius], vbar)
                axbonk(ax,xscale='log',yscale='log',xlabel='r',ylabel='KE')
                fig.savefig('plots_to_sort/Velocity_radius_%s_c%04d.png'%(name,core_id))

            if 0:
                fig,ax=plt.subplots(1,1)
                rrr = ms.r+0
                rrr[ rrr<1./2048] = 1./2048
                rrr = rrr.transpose()
                den = ms.density.transpose()
                fig,ax=plt.subplots(1,1)
                xmin = 1/2048
                xmax = 3.5e-1
                ymin = 1e-2
                ymax = 1e6
                for frame in frames:
                    nframe = np.where( thtr.frames == frame)
                    ax.clear()
                    ax.plot( rrr[:nframe,:], den[:nframe,:] , c=[0.5]*4, linewidth=0.2)
                    #ax.plot( ms.this_x.transpose(), ms.this_y.transpose(), c=[0.5]*4, linewidth=0.2)
                    axbonk(ax,xlabel='r',ylabel='rho',xscale='log',yscale='log', xlim=[xmin,xmax],ylim=[ymin,ymax])
                    outname='plots_to_sort/rho-r-t_%s_c%04d_i%04d'%(name,core_id,nframe)
                    fig.savefig(outname)
                    print(outname)
            if 0:
                import colors
                fig,ax=plt.subplots(1,1)
                t = thtr.times/colors.tff
                rrr = ms.r+0
                rrr[ rrr<1./2048] = 1./2048
                rrr = rrr.transpose()
                den = ms.density.transpose()
                fig,ax=plt.subplots(1,1)
                xmin = 1/2048
                xmax = 3.5e-1
                ymin = 1e-2
                ymax = 1e6
                ax.plot( t, den, c=[0.5]*4, linewidth=0.2)
                axbonk(ax,xlabel=r'$t/t_{ff}$', ylabel='rho',yscale='log')
                fig.savefig('plots_to_sort/rho-t_%s_c%04d'%(name,core_id))
            if 1:
                import colors
                fig,ax=plt.subplots(2,2,figsize=(12,12))
                ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[1][0]
                rrr = ms.r+0
                rrr[ rrr<1./2048] = 1./2048
                rrr = rrr.transpose()
                t = thtr.times/colors.tff
                den = ms.density.transpose()
                xmin = 1/2048
                xmax = 3.5e-1
                ymin = 1e-2
                ymax = 1e6
                for iframe,frame in enumerate(frames):
                    nframe = np.where( thtr.frames == frame)[0][0]
                    ax0.clear()
                    ax1.clear()
                    ax0.plot( rrr[:nframe+1,:], den[:nframe+1,:] , c=[0.5]*4, linewidth=0.2)
                    ax0.scatter( rrr[nframe,:], den[nframe,:], c='k',s=0.2)

                    ax1.plot( t[:nframe+1], den[:nframe+1,:], c=[0.5]*4, linewidth=0.2)
                    ax1.scatter( [t[nframe]]*ms.nparticles, den[nframe,:], c='k',s=0.2)

                    ax2.plot( pos[0][:nframe+1,:], pos[1][:nframe+1,:], c=[0.5]*4, linewidth=0.2)
                    ax2.scatter( pos[0][nframe+1,:], pos[1][nframe+1,:], c='k', s=0.2)
                    #ax.plot( ms.this_x.transpose(), ms.this_y.transpose(), c=[0.5]*4, linewidth=0.2)
                    axbonk(ax0,xlabel='r',ylabel='rho',xscale='log',yscale='log', xlim=[xmin,xmax],ylim=[ymin,ymax])
                    outname='plots_to_sort/so_much_hair_%s_c%04d_i%04d'%(name,core_id,iframe)
                    fig.savefig(outname)
                    print(outname)





sim_list = ['u502']
import three_loopers_u500 as TL
import close_tool
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()
        
if 0:
    for this_simname in sim_list:
        this_looper=TL.loops[this_simname]
        nt = potentool(this_looper)
        nt.run(core_list=[323])
if 1:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = potentool(this_looper)
    nt.run(core_list=[76], frames=[20])

if 0:
    for this_simname in sim_list:
        this_looper=TL.loops[this_simname]
        #get the core with the biggest distance to the next core.
        #Thus, the max of the min for each core.
        d = ct[this_simname].distance_matrix
        d[d==0]=3 #the actual min is itself, which is zero.
        most_isolated_core = ct[this_simname].cores_used[ np.argmax( d.min(axis=0))]
        nt = potentool(this_looper)
        nt.run(core_list=[most_isolated_core])


