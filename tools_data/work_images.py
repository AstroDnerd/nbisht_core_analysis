
from starter2 import *
import xtra_energy
reload(xtra_energy)
core_list=[79]
frame_list = range(0,125,5)
#frame_list=[10]
if 1:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_force_terms],
                                     sim_name = 'u05',
                                     out_prefix = 'work1',
                                     target_frame = 125,
                                     frame_list = frame_list, #[0,120],#list(range(11,120))+[125],
                                     core_list = core_list, #[THIS_CORE],
                                     fields_from_grid=['x','y','z','mag_work'],
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='datasets_small/bad_particles.h5')
@looper.core_loop
def core_proj_follow(looper,snapshot, field='density', axis_list=[0,1,2], color='r'):
    for ax in axis_list:
        scale_min = snapshot.ds.arr(0.05,'code_length')
        scale = max([snapshot.R_mag.max(), scale_min])
        sph = snapshot.ds.sphere(snapshot.R_centroid,scale)
        proj = snapshot.ds.proj(field,ax,center=snapshot.R_centroid, data_source = sph)
        pw = proj.to_pw(center = snapshot.R_centroid,width=(1.0,'code_length'), origin='domain')
        pw.zoom(1./(2*scale.v))
        #pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(field,'gray')
        #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max(), circle_args={'color':color} ) #R_mag.max())
        #pw.annotate_text(snapshot.R_centroid,
        #                 "%d"%snapshot.core_id,text_args={'color':color}, 
        #                 inset_box_args={'visible':False},
        #                 coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
    return pw
core_proj_follow( this_looper,axis_list=[0],field='density')
