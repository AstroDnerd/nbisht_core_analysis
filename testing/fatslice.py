

from starter2 import *
import loop_tools
reload(loop_tools)
import pyximport; pyximport.install()
import particle_ops
all_cores  =looper.get_all_nonzero()
core_list=all_cores
save_field = 'all_cores_n0000.h5'
if 'this_looper' not in dir() and os.path.exists(save_field) and False:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,savefile=save_field)
    #this_looper.derived = [xtra_energy.add_force_terms] #for some reason derived quantities aren't saved.
if 'this_looper' not in dir() and True:
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_force_terms],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()

grid_quan_temp={}
grid_quan_temp["DomainLeft"] =np.zeros(3)
grid_quan_temp["DomainRight"] =np.ones(3)
grid_quan_temp["DomainWidth"] =np.ones(3)
grid_quan_temp["max_dx"]      =1./128
grid_quan_temp["min_dx"]      =1./128/2**4


@looper.core_loop
def fat_core_slice(myloop,snapshot,field='density',axis_list=[0],core_list=None,mask_stash=None,):
    MIN_DX = snapshot.ds.quan(1./2048,'code_length')

    for axis in axis_list:
        grid_quan=None
        #This grid_quan is not a good workaround for an oversight in my datastructure.
        if snapshot.ds is None:
            grid_quan=grid_quan_temp
        shifted = loop_tools.shift_particles(snapshot.ds,snapshot.pos,shiftRight=False,grid_quan=grid_quan)
        these_coords = shifted[:,axis]
        regleft = np.zeros(3)
        regright= np.ones(3)
        regleft[axis] = these_coords.min()
        regright[axis] = max([these_coords.max(), these_coords.min()+MIN_DX])
        center = 0.5*(regleft+regright)
        reg = snapshot.ds.region(center,regleft,regright)
        proj = snapshot.ds.proj(field,axis,center=center,data_source=reg)
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap(field,'gray')
        pw.annotate_select_particles(1.0, col='r', indices=snapshot.target_indices)
        outname = '%s_fatslice_c%04d_n%04d'%(myloop.out_prefix,snapshot.core_id,myloop.current_frame)
        print( pw.save(outname))

#this_looper.core_list=[79,31]
fat_core_slice(this_looper,axis_list=[0,1,2],field='PotentialField')
