
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
#if 'this_looper' not in dir():
#    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
#    this_looper = looper.core_looper(directory= directory,
#                                     sim_name = 'u05',
#                                     out_prefix = 'test',
#                                     target_frame = 125,
#                                     frame_list = frame_list,
#                                     core_list = core_list, #[THIS_CORE],
#                                     fields_from_grid=['x'] #,'y','z']
#                                  )
#    this_looper.get_target_indices(h5_name='u05_0125_peaklist.h5',
#                                     bad_particle_list='bad_particles.h5')
#    #this_looper.get_tracks()
#
#dbg = 0
#import pdb
def deposit_target_particles_tmp(field, data):
    d=data['density']
    return data.apply_units(d, field.units)
def deposit_target_particles(field, data):
    #pdb.set_trace()

    target_indices= data.get_field_parameter('target_indices')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if not hasattr(target_indices,'sum'):
        return blank
    if data.NumberOfParticles == 0: return blank
#if dbg > 0:
#    if hasattr(target_indices,'sum'):
#        print("indices inside", target_indices.sum())
#    else:
#        print("indices inside", target_indices)
#        return blank
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
    #d=data['density']
    return  data.apply_units(d, field.units)

def deposit_target_particles_boolean(field,data):
    target_indices= data.get_field_parameter('target_indices')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if not hasattr(target_indices,'sum'):
        return blank
    if data.NumberOfParticles == 0: return blank
    deposit_tuple=("deposit","deposit_target_particles")
    output = data[deposit_tuple[1]]
    return output


def add_tracer_density(obj):
    def tracer_number_density(field,data):
        norm = data.get_field_parameter('particle_norm')
        if norm == 0:
            norm = 1
        cic = data[('deposit','all_cic')].in_units('code_density')
        number_density = cic/norm
        return number_density
    these_units='1' #/code_length**3'
    obj.add_field('tracer_number_density',function=tracer_number_density, #units=these_units,
                   validators=[yt.ValidateParameter("particle_norm")],sampling_type='cell')

    obj.add_field(      ("deposit","deposit_target_particles"),
             function = deposit_target_particles,
             validators = [yt.ValidateSpatial(), 
                           yt.ValidateParameter('target_indices'), 
                           yt.ValidateParameter('mask_to_get'), 
                           yt.ValidateGridType()],
             display_name = "target_particles",sampling_type='cell')
    obj.add_field(      ("deposit","deposit_target_particles_boolean"),
             function = deposit_target_particles_boolean,
             validators = [yt.ValidateSpatial(), 
                           yt.ValidateParameter('target_indices'), 
                           yt.ValidateParameter('mask_to_get'), 
                           yt.ValidateGridType()],
             display_name = "target_particles",sampling_type='cell')

#    ad2 = ds2.all_data()
#    my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
#    ad.set_field_parameter('particle_norm',my_norm)


mask_stash = None
@looper.core_loop
def mask_test(looper,snapshot,axis_list=[0,1,2]):
    print("SHIT GUYS")
    return
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
#looper.core_looper.mask_test=mask_test
#this_looper.mask_test()
mask_stash = None
@looper.frame_loop
def mask_test2(myloop,axis_list=[0,1,2],core_list=None):
    if core_list is None:
        core_list = myloop.core_list
    global mask_stash
    axis=0
    ds = myloop.load(frame=myloop.current_frame,derived=[add_tracer_density])
    ad=ds.all_data()
    add_tracer_density(ds)
    if mask_stash is None:
        mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
    all_target_indices = np.concatenate( [myloop.target_indices[core_id] for core_id in core_list])
    mask_stash[:] *= 0 
    my_norm = (8*ds.index.grids[0]['particle_mass'][0]/ds.index.grids[0].dds.prod()).in_units('code_density')
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',mask_stash)
    ad.set_field_parameter('particle_norm',my_norm)
    myloop.ad = ad
    deposit_tuple=("deposit","deposit_target_particles")
    proj = ds.proj(deposit_tuple,axis,center='c',data_source=ad)
    if 1:
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(deposit_tuple,'gray')
        radius_from_core = []
        core_label = ""
        #pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
        #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
        outname = '%s_cic_test_n%04d'%(myloop.out_prefix,myloop.current_frame)
        print(outname)
        print( pw.save(outname))
#looper.core_looper.mask_test=mask_test
#this_looper.mask_test()

def get_deposit_field(myloop,frame=None,core_list=None, mask_stash=None):
    if frame is None: frame = myloop.current_frame
    if core_list is None: core_list = myloop.core_list
    ds = myloop.load(frame=frame,derived=[add_tracer_density])
    ad = ds.all_data()
    all_target_indices = np.concatenate( [myloop.target_indices[core_id] for core_id in core_list])
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',mask_stash)
    return ad


@looper.frame_loop
def mask_test3(myloop,axis_list=[0,1,2],core_list=None):
    global mask_stash
    axis=0
    ds = myloop.load(frame=myloop.current_frame,derived=[add_tracer_density])
    add_tracer_density(ds)
    if mask_stash is None:
        mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
    ad=get_deposit_field(myloop,frame=myloop.current_frame, core_list=core_list, mask_stash=mask_stash)
    mask_stash[:] *= 0 
    deposit_tuple=("deposit","deposit_target_particles")
    proj = ds.proj(deposit_tuple,axis,center='c',data_source=ad)
    pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
    pw.set_cmap(deposit_tuple,'gray')
    radius_from_core = []
    core_label = ""
    #pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
    #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
    outname = '%s_cic_test_3_n%04d'%(myloop.out_prefix,myloop.current_frame)
    print( pw.save(outname))


@looper.frame_loop
def mask_test4(myloop,axis_list=[0,1,2],core_list=None):
    global mask_stash
    axis=0
    ds = myloop.load(frame=myloop.current_frame,derived=[add_tracer_density])
    add_tracer_density(ds)
    if mask_stash is None:
        mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
    ad=get_deposit_field(myloop,frame=myloop.current_frame, core_list=core_list, mask_stash=mask_stash)
    mask_stash[:] *= 0 
    deposit_tuple=("deposit","deposit_target_particles_boolean")
    proj = ds.proj(deposit_tuple,axis,center='c',data_source=ad)
    pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
    pw.set_cmap(deposit_tuple,'gray')
    radius_from_core = []
    core_label = ""
    #pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
    #pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
    outname = '%s_cic_test_4_n%04d'%(myloop.out_prefix,myloop.current_frame)
    print( pw.save(outname))


