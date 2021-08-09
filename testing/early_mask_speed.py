from starter2 import *

import three_loopers_mountain_top as TLM
import testing.early_mask as em
reload(em)
this_looper=TLM.loops['u301']
ds = this_looper.load(0)
em.add_tracer_density(ds)
ad = ds.all_data()
all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in this_looper.core_list])
all_target_indices = all_target_indices.astype('int64')
ad.set_field_parameter('target_indices',all_target_indices)
ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
t1=time.time()
tpvf=ad[('deposit','target_particle_volume_fast')]
t2=time.time()
tpv=ad[('deposit','target_particle_volume')]
t3=time.time()
print('did it work? %0.2e'%( np.abs( tpv-tpvf).sum()))
print('did it work? %0.2e %0.2e'%(t3-t2,t2-t2))
