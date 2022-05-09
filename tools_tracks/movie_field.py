from starter2 import *
import pcolormesh_helper as pch
reload(pch)
import annotate_hair 
if 1:
    import three_loopers_u500 as TL5
    this_looper = TL5.loops['u501']

YT_mag_vel = "gas","cos_b_v"
def add_stuff(obj):
    def magvel(field,data):
        bx= data['magnetic_field_x']
        by= data['magnetic_field_y']
        bz= data['magnetic_field_z']
        vx= data['velocity_x']
        vy= data['velocity_y']
        vz= data['velocity_z']
        bmag = data['magnetic_field_strength']
        vmag = data['velocity_magnitude']
        return (bx*vx+by*vy+bz*vz)/(bmag*vmag)
    obj.add_field( YT_mag_vel, magvel, units='dimensionless', sampling_type='cell',take_log=False)

if 1:
    thtr=this_looper.tr
    for nframe,frame in enumerate(thtr.frames):
        for core_id in [323]:
            ds = this_looper.load(frame)
            add_stuff(ds)
            ms = trackage.mini_scrubber(thtr,core_id)
            ms.particle_pos(core_id)

            all_x,all_y,all_z=ms.particle_x,ms.particle_y, ms.particle_z
            all_p = [all_x,all_y,all_z]
            all_p_s = np.stack(all_p)
            max_max = all_p_s.max(axis=1).max(axis=1)
            min_min = all_p_s.min(axis=1).min(axis=1)
            #Bump to zones
            quant_res=128
            min_min = np.round( min_min*quant_res)/quant_res
            max_max = np.round( max_max*quant_res)/quant_res
            cen = 0.5*(min_min+max_max)
            scale = (max_max-min_min).max()


            proj_axis=0
            rect  = ds.region(cen, min_min,max_max)

            xax = ds.coordinates.x_axis[proj_axis]
            yax = ds.coordinates.y_axis[proj_axis]
            xaf = 'xyz'[xax]
            yaf = 'xyz'[yax]

            #px = ds.coordinates.pixelize( 
            field=YT_mag_vel
            proj = ds.proj(field,proj_axis,center=cen,data_source=rect)
            pw=proj.to_pw(center=cen,origin='domain')
            pw.annotate_hair('magnetic_field_%s'%xaf,'magnetic_field_%s'%yaf, this_looper=this_looper,frame=frame, core_id=core_id)
            pw.zoom(1./scale)
            pw.set_axes_unit('code_length')
            pw.set_zlim(field,-1,1)
            outname='plots_to_sort/derp_%s_c%04d_i%04d'%(this_looper.sim_name,core_id,nframe)
            pw.save(outname)
            print(outname)

