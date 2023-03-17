from starter2 import *
import xtra_energy

import core_proj_three
reload(core_proj_three)

#import three_loopers_six as TL
import camera_path
import three_loopers_u500 as TL
import pcolormesh_helper as pch 
sim_list=['u501','u502','u503']
sim_list=['u502']

if 0:
    for sim in sim_list:
        loop = TL.loops[sim]
        camera = camera_path.camera_1( loop, 'smooth_zoom_2')
        core_proj_three.core_proj_multiple(loop,axis_list=[0],core_list=[74],frame_list=[0,10],camera=camera, main_core=74)


class vel_plot():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
    def run(self, core_list=None, frame_list=None, tsing=None):
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)
        def get_time_index(time):
            index=np.argmin( np.abs( thtr.times/colors.tff-time))
            return index
        for core_id in core_list:
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            self.cores_used.append(core_id)

            frame_mask = np.zeros_like(thtr.times, dtype='bool')
            frame_mask[get_time_index(tsing.tend_core[core_id])]=True
            frame_mask[get_time_index(tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(0.9*tsing.tsing_core[core_id])]=True
            frame_mask[get_time_index(0.75*tsing.tsing_core[core_id])]=True
            #frame_mask[0]=True
            frame_mask[-1]=True
            frame_list=thtr.frames[frame_mask]
            rm = rainbow_map(len(frame_list))

            fig,axes = plt.subplots(3,len(frame_list), figsize=(12,12))
            fig2,ax2=plt.subplots(1,2)
            img_collector=[]
            ext_rho=extents()
            ext_rho2=extents()
            ext_r=extents()
            for nframe,frame in enumerate(frame_list):
                print("profile on %s c%04d n%04d"%(this_looper.sim_name, core_id, frame))
                ds = this_looper.load(frame)
                #xtra_energy.add_energies(ds)
                #xtra_energy.add_gdotgradrho(ds)
                nf = np.where( this_looper.tr.frames == frame)[0][0]

                center = nar([ms.mean_xc[nf], ms.mean_yc[nf],ms.mean_zc[nf]])
                #center = nar([ms.mean_x[nf], ms.mean_y[nf],ms.mean_z[nf]])
                msR = ms.rc
                msR[ msR<1/2048]=1/2048
                
                MaxRadius=msR[:,nf].max()
                Radius = max([4/2048, MaxRadius])
                #Radius = 32/2048
                print('Max Radius',MaxRadius, Radius)
                rsph = ds.arr(Radius,'code_length')
                sph = ds.sphere(center,rsph)
                left = center - Radius
                right = center + Radius
                dx=1/2048
                nzones = np.floor((right-left)/dx).astype('int')
                print(nzones, left, right)
                print(right-left)
                cg = ds.covering_grid(4,left,nzones, num_ghost_zones=1)
                rho = cg[YT_density]
                #coordinates
                LOS_velocity = YT_velocity_x
                proj_axis=0
                h_axis = ds.coordinates.x_axis[proj_axis]
                v_axis = ds.coordinates.y_axis[proj_axis]
                vx = cg[LOS_velocity]
                dv = cg[YT_cell_volume]
                rho_proj = cg[YT_density].mean(axis=proj_axis)
                proj=axes[0][nframe].imshow( rho_proj,norm=mpl.colors.LogNorm(), extent=[left[v_axis],right[v_axis],left[h_axis],right[h_axis]])


                if 1:
                    cg.set_field_parameter( 'center',ds.arr(center, 'code_length'))
                    R = cg['radius'].flatten()
                    argsort = np.argsort(R)
                    sl = slice(0,50)
                    #mean_vx = (vx*rho*dv).flatten()[argsort][sl].sum()/(rho*dv).flatten()[sl].sum()
                    subset = rho > 100
                    mean_vx = (vx*rho*dv)[subset].sum()/(rho*dv)[subset].sum()
                    vx_nomean = vx-mean_vx
                if 0:
                    vx_proj=( rho*dv*(vx_nomean)).sum(axis=proj_axis)/(rho*dv).sum(axis=proj_axis)
                    vx_proj=np.abs(vx_proj)
                    proj=axes[1][nframe].imshow( vx_proj,norm=mpl.colors.Normalize(vmin=0,vmax=+1))
                    fig.colorbar(proj,ax=axes[1][nframe])

                if 0:
                    b=axes[1][nframe].contour(rho.mean(axis=proj_axis), 1, colors='r', alpha=0.1)
                    #pdb.set_trace()

                if 0:
                    g1 = rho
                    g2 = np.abs(vx_nomean)
                    ext=extents()
                    #ext( g1[g1>0])
                    #ext( g2[g2>0])
                    ext(nar([1e-2,5e4]))
                    #pdb.set_trace()
                    norm = mpl.colors.LogNorm( vmin=ext.minmax[0], vmax=ext.minmax[1])
                    if 0:
                        binsx = np.geomspace( ext.minmax[0], ext.minmax[1], 64)
                        binsy = np.geomspace( ext.minmax[0], ext.minmax[1], 64)
                    if 1:
                        binsx = np.geomspace( 1e-2,5e4,64)
                        #binsy = np.geomspace( 0.1,
                        #binsy = np.linspace( 1e-3,10,64)
                        binsy = np.geomspace( 1e-3,10,64)

                    #do_plot( p1, axlist[0], norm)
                    #do_plot( p2, axlist[1], norm)

                    hist, xbins, ybins = np.histogram2d( g1.flatten(), g2.flatten(), bins=[binsx,binsy])

                    aaax=axes[2][nframe]
                    pch.helper( hist, xbins, ybins, ax=aaax)
                    #axbonk(aaax,xscale='log',yscale='linear', xlabel='rho',ylabel='v')
                    axbonk(aaax,xscale='log',yscale='log', xlabel='rho',ylabel='v')
                if 1:
                    binsx = np.geomspace( 1e-2,5e4,64)
                    binsy = np.geomspace( 1e-3,10,64)
                    ax2[0].hist( rho.flatten().v, histtype='step',bins=binsx,density=True,cumulative=False)
                    ax2[1].hist( np.abs(vx_nomean).flatten().v, histtype='step',bins=binsy,density=True,cumulative=False)
                    ax2[0].set(xscale='log',yscale='log',xlabel='rho')
                    ax2[1].set(xscale='log',yscale='log',xlabel='|V|')



                if 0:
                    this_x=[ms.this_x[:,nf],ms.this_y[:,nf],ms.this_z[:,nf]][v_axis]#note transpose
                    this_y=[ms.this_x[:,nf],ms.this_y[:,nf],ms.this_z[:,nf]][h_axis]#note transpose
                    axes[0][nframe].scatter( this_x,this_y,s=0.1,alpha=0.5)
                if 0:
                    axes[2][nframe].hist( np.log10((rho).v.flatten()), histtype='step')
                if 0:
                    axes[2][nframe].hist( (vx).v.flatten(), histtype='step')
                    axes[2][nframe].hist( (vx-mean_vx).v.flatten(), histtype='step')

                if 0:
                    #just to check 
                    pw = yt.ProjectionPlot(ds, 0, YT_density, data_source=sph, center=center, origin='window')
                    frb =  pw.data_source.to_frb( 2*sph.radius,512)
                    axes[1][nframe].imshow( frb[YT_density], norm=mpl.colors.LogNorm())

                if 0:
                    cg.set_field_parameter( 'center',ds.arr(center, 'code_length'))

                    proj=axes[1][nframe].imshow( cg['radius'].mean(axis=0).v)#,norm=mpl.colors.Normalize(vmin=-2,vmax=+2))
                    fig.colorbar(proj,ax=axes[1][nframe])

                if 0:
                    import other_scrubber
                    reload(other_scrubber)
                    print('wut')
                    R0 = cg['radius'].v
                    print('fl')
                    R1=R0.flatten()
                    print('argsrt', R1.size, R0.shape)
                    order = np.argsort(R1)
                    print('hey')
                    vel = []
                    for axis in 'xyz':
                        #vel.append( cg['velocity_%s'%axis].flatten()[order][:10].mean())
                        v = cg['velocity_%s'%axis].flatten()
                        D = cg[YT_density].flatten()
                        dv= cg[YT_cell_volume].flatten()
                        mean_v = (v*D*dv).sum()/(D*dv).sum()
                        vel.append( mean_v)
                    print('scrub')
                    scrub = other_scrubber.scrubber(cg, reference_velocity = vel)
                    print('proj')


                    proj=axes[1][nframe].imshow( scrub.vr_rel.sum(axis=0) ,norm=mpl.colors.Normalize(vmin=-2,vmax=+2))
                    fig.colorbar(proj,ax=axes[1][nframe])



            fig.savefig( 'plots_to_sort/subsonic_%s_c%04d'%(this_looper.sim_name,core_id))
            fig2.savefig( 'plots_to_sort/rho_v_hist_%s_c%04d'%(this_looper.sim_name,core_id))




import tsing
reload(tsing)
if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

#import anne
#reload(anne)
##anne.make_inflection()
if 'multi_proj' not in dir() or True :
    for sim in sim_list:
        all_cores=np.unique( TL.loops[sim].tr.core_ids)
        core_list=list(all_cores)
        core_list=None#[323]
        core_list=[323]
        core_list=[25]
        core_list=[114]
        core_list=[195]
        core_list=[74]
        #core_list = TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[10:11]
        #core_list=[114]


        #core_list = [114]
        #core_list=[361]
        #core_list=[8]
        #core_list=[381]
        #core_list=[323]
        TF=TL.loops[sim].target_frame
        all_frames = TL.loops[sim].tr.frames
        nframes = all_frames.size
        frame_list = all_frames #[0, TF, all_frames[int(nframes/2)]]
        #frame_list = [100]

        for core_id in core_list:
            frame_list = [tsing_tool[sim].tend_frame[core_id]]
            vp=vel_plot(TL.loops[sim])
            vp.run(core_list=[core_id],tsing=tsing_tool[sim], frame_list=frame_list)#, r_inflection=anne.inflection[sim])
