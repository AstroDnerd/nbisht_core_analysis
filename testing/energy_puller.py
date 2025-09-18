from starter2 import *
import xtra_energy
reload(xtra_energy)
import data_locations as dl
all_nonzero = looper.get_all_nonzero()
#many1: frame_list = list(range(10,130,10)) + [125]
#many1: core_list =  [ 120, 10, 308, 41, 44, 110, 1, 65], 
#frame 120 core 96
core_list = all_nonzero.astype('int')[-3:]
core_list=[79]
frame_list = list(range(10,130,10)) + [125]

fields = ['x','y','z','velocity_magnitude','magnetic_field_strength', 'density']
fields += ['kinetic_energy','magnetic_energy','grav_energy','therm_energy']
fields += ['momentum_flux','gravity_work','mag_work','pressure_work']
#fields = ['grav_energy']#,'therm_energy']
if 0:
#h5ptr = h5py.File('u05_0125_peaklist.h5','r')
#cen = h5ptr['peaks'][79]
    cen=np.array([0.40258789, 0.6862793 , 0.05249023])

    ds = yt.load(dl.enzo_directory+"/DD0125/data0125")
#xtra_energy.add_energies(ds)
#xtra_energy.add_test_energies(ds)
    xtra_energy.add_force_terms(ds)
    xtra_energy.add_energies(ds)

    radius=0.1
    sphere=ds.sphere(cen,radius)
    for field in fields:
        proj = ds.proj(field,0,center=cen,data_source=sphere)
        pw = proj.to_pw(center=cen,width=2.*radius)
        #pw.set_zlim('pressure_work',-1e6,1e6)
        pw.set_log(field,True,linthresh=0.1)
        pw.save('/home/dcollins4096/PigPen/test_c0079')

#print(gx['momentum_flux'])
#print(gx['grav_x'])
#for frame in this_looper.ds_list:
#    xtra_energy.add_test_energies(this_looper.ds_list[frame])
#sys.exit(0)
output_base = "energy_test_full"
core_list=all_nonzero
for core in core_list:
    output_name = '%s_c%04d.h5'%(output_base,core)
    if os.path.exists(output_name):
        print("Skipping "+output_name)
        continue
    this_looper = looper.core_looper(directory= dl.enzo_directory,
                                     sim_name = 'u05',
                                     out_prefix = 'test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list =  [core],# core_list,
                                     derived=[xtra_energy.add_energies,xtra_energy.add_force_terms],
                                     fields_from_grid=fields,
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='datasets_small/bad_particles.h5')
    this_looper.get_tracks()
    this_looper.save(output_name)
