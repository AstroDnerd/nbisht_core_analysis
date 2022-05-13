
from starter2 import *
import three_loopers_six as TLM
import xtra_energy
this_looper = TLM.loops['u601']

from starter2 import *
import data_locations as dl
from collections import defaultdict

import pcolormesh_helper as pch
reload(pch)
import davetools
reload(davetools)

import colors
plt.close('all')


class more_hair():
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
            reload(trackage)
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('go ', core_id)

            this_center=ms.mean_center[:,-1]
            self.cores_used.append(core_id)

            ms.particle_pos(core_id)
            if ms.nparticles <= 1000:
                sl = slice(None)
                nparticles=ms.nparticles
            else:
                sl = slice(None,None,10)
                nparticles=ms.raw_vx[:,0][sl].size
            pos = np.stack([ ms.particle_x[sl].transpose(),ms.particle_y[sl].transpose(), ms.particle_z[sl].transpose()])
            vel = np.stack([ ms.raw_vx[sl].transpose(),ms.raw_vy[sl].transpose(), ms.raw_vz[sl].transpose()])
            ms.get_central_at_once(core_id)
            cen_vmag = ms.cen_vmag[sl].transpose()
            rrr = ms.r[sl]+0
            rrr = rrr.transpose()
            rrr[ rrr<1./2048] = 1./2048
            den = ms.density[sl].transpose()
            dv = ms.cell_volume[sl].transpose()
            vrad = ms.vr_rel[sl].transpose()

            vrad_cell_avg = (vrad*dv).sum(axis=1)/dv.sum(axis=1)
            vrad_cell_avg.shape = (vrad_cell_avg.size,1)
            vrad_cell_var = np.sqrt(((vrad-vrad_cell_avg)**2*dv).sum(axis=1)/dv.sum(axis=1))
            vrad_cell_avg.shape=vrad_cell_avg.size
            vrad_mass_avg = (den*vrad*dv).sum(axis=1)/(den*dv).sum(axis=1)
            print('ffff',vrad_mass_avg.shape)
            print('www',vrad.shape)

            name = self.this_looper.sim_name
            
            import colors
            fig,ax=plt.subplots(2,3,figsize=(12,12))
            ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[0][2]
            ax3=ax[1][0]; ax4=ax[1][1];ax5=ax[1][2]
            ax1b = ax1.twinx()
            t = thtr.times/colors.tff
            xmin = 1/2048
            xmax = 3.5e-1
            ymin = 1e-2
            ymax = 1e6
            for iframe,frame in enumerate(frames):
                nframe = np.where( thtr.frames == frame)[0][0]
                print("nframe",nframe,thtr.frames.size)
                ax0.clear()
                ax1.clear()
                ax1b.clear()
                ax2.clear()
                ax3.clear()
                ax4.clear()
                ax5.clear()
                c  = [0.5,0.5,0.5,0.5]
                print(pos[0].shape)

                ax0.plot( pos[0][:nframe+1,:], pos[1][:nframe+1,:], c=c, linewidth=0.2)
                ax0.scatter( pos[0][nframe,:], pos[1][nframe,:], c='k', s=0.2)

                #xbins = np.linspace( pos[0].min(),pos[0].max(),129)
                #ybins = np.linspace( pos[1].min(),pos[1].max(),129)
                #hist, xb, yb = np.histogram2d( pos[0].flatten(), pos[1].flatten(), bins=[xbins,ybins],weights=dv.flatten())
                ##pch.helper(hist.transpose(), yb.transpose(), xb.transpose(), ax=ax0)
                #pch.helper(hist, xb, yb, ax=ax0, transpose=False)
                axbonk(ax0,xlabel='x',ylabel='y', xlim=[pos[0].min(),pos[0].max()], ylim=[ pos[1].min(), pos[1].max()])

                #dbins = np.geomspace( den[d>0].min(), den.max(),64)
                #bbins = np.geomspace( b[b>0].min(), b.max(),64)

                if 1:
                    ax1.plot( t[:nframe+1], den[:nframe+1,:], c=c, linewidth=0.2)
                    ax1.scatter( [t[nframe]]*nparticles, den[nframe,:], c='k',s=0.2)
                    ax1b.violinplot( np.log10(den[nframe,:]), positions=[t[nframe]])
                    axbonk(ax1,xlabel='t',ylabel='rho',yscale='log',ylim=[ymin,ymax])
                    axbonk(ax1b,xlabel='',ylabel='',yscale='linear',ylim=[np.log10(ymin),np.log10(ymax)])

                    ax2.plot( t[:nframe+1], vrad[:nframe+1,:], c=c,linewidth=0.2)
                    ax2.scatter( [t[nframe]]*nparticles, vrad[nframe,:], c='k',s=0.2)
                    ax2.violinplot( vrad[nframe,:], positions=[t[nframe]])
                    axbonk(ax2,xlabel='t/tff',ylabel='Vradial')
                    ax2.plot(t,vrad_cell_avg,c='r',label='cell')
                    ax2.plot(t,vrad_cell_avg+vrad_cell_var,c='r')
                    ax2.plot(t,vrad_cell_avg-vrad_cell_var,c='r')
                    ax2.plot(t,vrad_mass_avg,c='g',label='mass')
                    ax2.legend(loc=0)


                    #density-radius
                    #ax2.plot( rrr[:nframe+1,:], den[:nframe+1,:] , c=c, linewidth=0.2)
                    #ax2.scatter( rrr[nframe,:], den[nframe,:], c='k',s=0.2)
                    #axbonk(ax2,xlabel='r',ylabel='rho',xscale='log',yscale='log', xlim=[xmin,xmax],ylim=[ymin,ymax])

                    #ax3.plot( t[:nframe+1], cen_vmag[:nframe+1,:], c=[0.5]*4,linewidth=0.2)
                    #ax3.scatter( [t[nframe]]*nparticles, cen_vmag[nframe,:], c='k',s=0.2)
                    #ax3.violinplot( cen_vmag[nframe,:], positions=[t[nframe]])

                    #ax4.plot( vel[0][:nframe+1,:], vel[1][:nframe+1,:], c=[0.5]*4, linewidth=0.2)
                    #ax4.scatter( vel[0][nframe+1,:], vel[1][nframe+1,:], c='k', s=0.2)
                    #axbonk(ax4,xlim=[vel[0].min(),vel[0].max()], ylim=[vel[1].min(),vel[1].max()])

                    for aaa,bbb in enumerate([ax3,ax4,ax5]):
                        bbb.plot( t[:nframe+1], vel[aaa][:nframe+1,:], c=c,linewidth=0.2)
                        bbb.scatter( [t[nframe]]*nparticles, vel[aaa][nframe,:], c='k',s=0.2)
                        bbb.violinplot( vel[aaa][nframe,:], positions=[t[nframe]])
                        axbonk(bbb,xlabel='t/tff',ylabel='V'+'xyz'[aaa])
                


                #ax.plot( ms.this_x.transpose(), ms.this_y.transpose(), c=[0.5]*4, linewidth=0.2)
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
if 0:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run( )
if 0:
    print('wut')
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run(core_list=[76], frames=[10,106])
if 0:
    this_simname = 'u501'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run(core_list=[122], frames=[125])#frames=list(range(0,106,10))+[106])
if 1:
    this_simname = 'u501'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run(core_list=None, frames=[125])#frames=list(range(0,106,10))+[106])
if 0:
    this_simname = 'u503'
    this_looper=TL.loops[this_simname]
    nt = more_hair(this_looper)
    nt.run( frames=[106])

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


