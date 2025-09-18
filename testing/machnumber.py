

#ad=mp.ad
vx = ad['x-velocity']
vy = ad['y-velocity']
vz = ad['z-velocity']
vv = ad['velocity_magnitude']
dv_real = ad['cell_volume']
if 1:
    all_target_indices = np.concatenate( [this_looper.target_indices[core_id] for core_id in core_list])
    ad.set_field_parameter('target_indices',all_target_indices)
    ad.set_field_parameter('mask_to_get',np.zeros_like(all_target_indices,dtype='int32'))
    deposit_tuple = ("deposit","target_particle_volume")
    dv_particles = ad[deposit_tuple]
dv = dv_particles
vx0=np.mean(vx)
vy0=np.mean(vy)
vz0=np.mean(vz)
density_mean = ((dv*ad['density']).sum()/dv.sum()).v
velocity_variance = np.sum(dv*((vx-vx0)**2+(vy-vy0)**2+(vz-vz0)**2))**0.5
density_variance=(dv*(ad['density'].v-density_mean)**2).sum()
print("sigma_v %0.2f sigma_rho %0.2f sigma_rho/sigma_v %0.2f"%(velocity_variance, density_variance,density_variance.v/velocity_variance.v))

