
from starter2 import *

import pyximport; pyximport.install()
import particle_ops
all_cores  =looper.get_all_nonzero()
all_cores.sort()
#ni = int(sys.argv[1])
core_list=[79]# all_cores[ni:ni+10]
frame_list=[100]
print("CORE_LIST", core_list)
print(len(core_list))
print(np.where(core_list==154))
rm = rainbow_map(len(all_cores))
#make an instance
#for THIS_CORE in core_list
if 'this_looper' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
    this_looper = looper.core_looper(directory= directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list, #[THIS_CORE],
                                     fields_from_grid=['x'] #,'y','z']
                                  )
    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    #this_looper.get_tracks()

dbg = 0
import pdb

def add_tracer_density(obj):
    def tracer_number_density(field,data):
        norm = data.get_field_parameter('particle_norm')
        cic = data[('deposit','all_cic')].in_units('code_density')
        number_density = cic/norm
        return number_density
    these_units='1' #/code_length**3'
    obj.add_field('tracer_number_density',function=tracer_number_density, #units=these_units,
                   validators=[yt.ValidateParameter("particle_norm")])

    def deposit_target_particles(field, data):
        #pdb.set_trace()
        target_indices= data.get_field_parameter('target_indices')
        blank = np.zeros(data.ActiveDimensions, dtype='float64')
        if not hasattr(target_indices,'sum'):
            return blank
        if data.NumberOfParticles == 0: return blank
        if dbg > 0:
            if hasattr(target_indices,'sum'):
                print("indices inside", target_indices.sum())
            else:
                print("indices inside", target_indices)
                return blank
        #int 32 because it's just a flag.
        #mask_to_get = np.zeros(target_indices.shape, dtype='int32')
        mask_to_get = data.get_field_parameter('mask_to_get')
        my_indices = data['particle_index'].astype('int64')
        found_any, mask = particle_ops.mask_particles(
            target_indices, my_indices, mask_to_get)
        data.set_field_parameter('mask_to_get',mask_to_get) #keeping the mask is faster
        pos_to_get = data['particle_position'][mask == 1]
        d = data.deposit(pos_to_get, method = "count")
        d = data.ds.arr(d, input_units = "cm**-3")
        #print("ahoy.", d.sum())
        return data.apply_units(d, field.units)

    obj.add_field(      ("deposit","deposit_target_particles"),
             function = deposit_target_particles,
             validators = [yt.ValidateSpatial(), 
                           yt.ValidateParameter('target_indices'), 
                           yt.ValidateParameter('mask_to_get'), 
                           yt.ValidateGridType()],
             display_name = "target_particles")

#    ad2 = ds2.all_data()
#    my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
#    ad.set_field_parameter('particle_norm',my_norm)


mask_stash = None
@looper.particle_loop
def mask_test(looper,snapshot,axis_list=[0,1,2]):
    global mask_stash
    if snapshot.frame not in looper.all_data:
        ds = looper.load(frame=snapshot.frame,derived=[add_tracer_density])
        looper.all_data[snapshot.frame] = ds.all_data()
    ad = looper.all_data[snapshot.frame]
    ds = looper.ds_list[snapshot.frame]
    seconds = [b for a,b in ds.field_info.keys()]
    if 'tracer_number_density' not in seconds:
        add_tracer_density(ds)
    if mask_stash is None:
        mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
    ad.set_field_parameter('target_indices',snapshot.target_indices+0)
    mask_stash[:] *= 0 #this should happen unless you're looking for all the cores.
    ad.set_field_parameter('mask_to_get',mask_stash)
    my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
    ad.set_field_parameter('particle_norm',my_norm)
    axis=0
    deposit_tuple=("deposit","deposit_target_particles")
    proj = snapshot.ds.proj(deposit_tuple,axis,center='c',data_source=ad)
    if 1:
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(deposit_tuple,'gray')
        radius_from_core = []
        core_label = ""
        #pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
        #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
        outname = '%s_cic_test_%sn%04d'%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(outname)
        print( pw.save(outname))
looper.core_looper.mask_test=mask_test
this_looper.mask_test()

