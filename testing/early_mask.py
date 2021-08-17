
from starter2 import *

import pyximport; pyximport.install()
import particle_ops
all_cores  =looper.get_all_nonzero()
def target_particle_volume_fast(field, data):
    #!!! not tested.
    #pdb.set_trace()
    target_indices= data.get_field_parameter('target_indices')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if type(target_indices) == float:
        return blank
    if data.NumberOfParticles == 0: return blank
    mask_to_get = data.get_field_parameter('mask_to_get')
    my_indices = data['particle_index'].astype('int64')
    print("target_particle_volume nparticles %d"%my_indices.size)
    #This search needs both lists to be sorted.
    #They're probably not sorted, so we sort them first. 
    #Since they're not sorted, we need to return to the original 
    #order.  Thus the two argsorts.  We do not need to acces
    #the targets, so we just sort those directly
    args_sort = np.argsort(my_indices) 
    args_return =np.argsort(args_sort) 
    my_ind_sorted = my_indices[args_sort]
    targets_sorted = np.sort(target_indices)
    found_any, mask_sorted = particle_ops.mask_particles_sorted_t7(
        targets_sorted, my_ind_sorted, mask_to_get)
    mask = mask_sorted[args_return]
    data.set_field_parameter('mask_to_get',mask_to_get) #keeping the mask is faster
    pos_to_get = data['particle_position'][mask == 1]
    d = data.deposit(pos_to_get, method = "count")
    d = data.ds.arr(d, input_units = "cm**-3")
    mask = d>0
    d[mask]  = data['cell_volume'][mask]
    return d
    #return  data.apply_units(d, field.units)

def target_particle_volume(field, data):
    #pdb.set_trace()
    target_indices= data.get_field_parameter('target_indices')
    blank = np.zeros(data.ActiveDimensions, dtype='float64')
    if type(target_indices) == float:
        return blank
    if data.NumberOfParticles == 0: return blank
    mask_to_get = data.get_field_parameter('mask_to_get')
    my_indices = data['particle_index'].astype('int64')
    print("target_particle_volume nparticles %d"%my_indices.size)
    found_any, mask = particle_ops.mask_particles(
        target_indices, my_indices, mask_to_get)
    data.set_field_parameter('mask_to_get',mask_to_get) #keeping the mask is faster
    pos_to_get = data['particle_position'][mask == 1]
    d = data.deposit(pos_to_get, method = "count")
    d = data.ds.arr(d, input_units = "cm**-3")
    mask = d>0
    d[mask]  = data['cell_volume'][mask]
    return d
    #return  data.apply_units(d, field.units)

def add_tracer_density(obj):
    obj.add_field(      ("deposit","target_particle_volume"),
             function = target_particle_volume,
             validators = [yt.ValidateSpatial(), 
                           yt.ValidateParameter('target_indices'), 
                           yt.ValidateParameter('mask_to_get'), 
                           yt.ValidateGridType()],
                  units = '1/cm**3',
             display_name = "target_particle_volume",sampling_type='cell')
    obj.add_field(      ("deposit","target_particle_volume_fast"),
             function = target_particle_volume_fast,
             validators = [yt.ValidateSpatial(), 
                           yt.ValidateParameter('target_indices'), 
                           yt.ValidateParameter('mask_to_get'), 
                           yt.ValidateGridType()],
                  units = '1/cm**3',
             display_name = "target_particle_volume",sampling_type='cell')

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
def project_particle_mask(myloop,axis_list=[0,1,2],core_list=None,mask_stash=None):
    for axis in axis_list:
        ds = myloop.load(frame=myloop.current_frame,derived=[add_tracer_density])
        add_tracer_density(ds)
        if mask_stash is None:
            mask_stash = np.zeros(ds['NumberOfParticles'], dtype='int32')
        if core_list is None:
            core_list = myloop.target_indices.keys()
        mask_stash[:] *= 0 
        deposit_tuple=("deposit","target_particle_volume")
        all_target_indices = np.concatenate( [myloop.target_indices[core_id] for core_id in core_list])
        field_parameters={}
        field_parameters['target_indices']=all_target_indices
        field_parameters['mask_to_get']=mask_stash
        proj = ds.proj(deposit_tuple,axis,center='c',field_parameters=field_parameters)
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(deposit_tuple,'gray')
        outname = '%s_cic_test_4_n%04d'%(myloop.out_prefix,myloop.current_frame)
        print( pw.save(outname))

