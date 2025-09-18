
from starter2 import *
from collections import defaultdict
import scipy
import colors
import xtra_energy
import hair_dryer
reload(hair_dryer)
import pcolormesh_helper as pch

#import three_loopers_six as TL
import three_loopers_u500 as TL
import movie_frames 

class energy_caul():
    def __init__(self,this_looper):
        self.this_looper=this_looper
    def run(self,core_list=None, do_plots=True, mass=None, dof=None, volume=None, frames='movie'):

        this_looper=self.this_looper
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        thtr=this_looper.tr
        all_times=thtr.times
        all_frames=thtr.frames
        if frames in ['movie', 'some']:
            mask = movie_frames.quantized_mask(this_looper).flatten()
            if frames == 'some':
                mask[:-10]=False


            times=thtr.times[mask]
            times.shape=times.size,1
            times=times/colors.tff
            frames=all_frames[mask]
        else:
            frame_index = nar([np.where( all_frames==frame)[0][0] for frame in frames])
            times = all_times[frame_index]

        rho_all = thtr.track_dict['density']
        rho_min=rho_all.min()
        rho_max=rho_all.max()

        mini_scrubbers={}
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
            mini_scrubbers[core_id] = ms
            ms.make_floats(core_id)
            #ms.particle_pos(core_id)
            ms.compute_ge(core_id)
            #ms.compute_ke(core_id)
            ms.compute_ke_rel(core_id)

        import camera_path
        reload(camera_path)
        #Run the camera here will all cores if you're making a multiple study.
        #camera = camera_path.camera_1(this_looper, 'sphere')
        #camera.run(core_list, frames, mini_scrubbers)



        for nc,core_id in enumerate(core_list):
            camera = camera_path.camera_1(this_looper, 'sphere')
            camera.run([core_id], frames, mini_scrubbers)
            print('V %s %d'%(this_looper.sim_name,core_id))
                
            ms = mini_scrubbers[core_id]

            #iframe is the frame position in frames
            #nf is the frame position in all_frames
            for iframe,frame in enumerate(frames):

                ds = this_looper.load(frame)
                xtra_energy.add_energies(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]
                rsph = camera.max_radius[iframe]
                rsph = max([ rsph, 1/128])
                center = camera.all_center[:,iframe]
                sp = ds.sphere(center,rsph)


                #Get data arrays

                dv = np.abs(sp[YT_cell_volume])
                RR =sp[YT_radius]
                DD = sp[YT_density]
                EG = np.abs(sp[YT_grav_energy_2])
                #don't use this, make the relative KE
                #EK = np.abs(sp[YT_kinetic_energy])
                #relative kinetic energy.
                vx = sp[YT_velocity_x].v - ms.mean_vx[nf]
                vy = sp[YT_velocity_y].v - ms.mean_vy[nf]
                vz = sp[YT_velocity_z].v - ms.mean_vz[nf]
                EK = 0.5*DD*(vx*vx+vy*vy+vz*vz)

                xxbins=np.geomspace(5e-3,1e7,128)
                yybins=np.geomspace(5e-3,1e7,128)
                #xxbins = np.geomspace(ke.min(),ke.max(),128)
                #yybins = np.geomspace(ge[ge>0].min(),ge.max(),128)
                hist, xbins,ybins=np.histogram2d(EK.flatten(),EG.flatten(),bins=[xxbins,yybins])

                fig,ax=plt.subplots(1,1)
                pch.helper(hist,xbins,ybins,ax=ax)
                axbonk(ax,xscale='log',yscale='log',xlabel='KE',ylabel='GE')
                ax.plot( xxbins,xxbins,c='k')
                ax.scatter(ms.ke_rel[:,nf],np.abs(ms.ge[:,nf]), edgecolor='r',s=30, facecolor='None')

                #ORDER = np.argsort( RR)
                #V_cuml =  np.cumsum( dv[ORDER])
                #V_local =  dv[ORDER]
                #RR_cuml = RR[ORDER]
                #EG_cuml = np.cumsum( EG[ORDER]*V_local)/V_cuml
                #EK_cuml = np.cumsum( EK[ORDER]*V_local)/V_cuml

                #line=line_list.get(frame,1)
                ##ax3[nnn].plot(  RR_cuml, EG_cuml, c=color_list[nnn], linestyle='-', linewidth=line)
                ##ax3[nnn].plot( RR_cuml, EK_cuml,  c=color_list[nnn], linestyle='--', linewidth=line)
                #y_ext(EG_cuml)
                #y_ext(EK_cuml)
                #r_ext(RR_cuml)









                outname='plots_to_sort/%s_energy_phase_c%04d_n%04d.png'%(this_looper.sim_name,core_id,frame)
                fig.savefig(outname)
                print(outname)
                plt.close(fig)



sims=['u501', 'u502','u503']
#import three_loopers_u500 as TL

sims=['u501', 'u502','u503']
sims=['u502']#, 'u501']
if 'caul_tool' not in dir():
    caul_tool={}
for sim in sims:
    #core_list=[381]
    core_list={'u501':[323], 'u602a':[381], 'u502':[112, 113]}[sim]
    #core_list=[31,32]
    #core_list=None
    caul_tool[sim]=energy_caul(TL.loops[sim])
    frames='movie'
    frames='some'
    #frames=[0,10,50,100]
    frrt=caul_tool[sim].run(  core_list=core_list, frames=frames)#, mass=mt[sim].unique_mass, dof=mt[sim].dof, volume=mt[sim].volume)

