"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import xtra_energy
reload(loop_apps)
core_list=[79,83]
frame_list=[40,50] #range(0,130,10)
fields=['density']  

#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
if os.path.exists('app_test.h5'):
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,savefile='app_test.h5')
    this_looper.derived = [xtra_energy.add_force_terms] #for some reason derived quantities aren't saved.
if 'this_looper' not in dir():

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
    this_looper.save('app_test.h5')

if 0:
    #print the centroid of each core.
    loop_apps.print_centroid(this_looper)

if 0:
    #project the entire domain, with particles form cores in the list plotted
    #All cores together.
    loop_apps.proj_cores(this_looper,axis_list=[0],core_list=core_list,field='density')

if 0:
    #this takes a list of particle indices.  
    #mostly for debugging purposes.
    loop_apps.select_particles(this_looper,axis_list=[0],
                               these_particles=this_looper.target_indices[ core_list[0] ])

if 1:
    #Plot each core individually, zoomed in to the core itself.
    loop_apps.core_proj_follow(this_looper,field='mag_work',axis_list=[0],force_log=True)


if 0:
    """Draws circles around each core."""
    loop_apps.core_circle( this_looper, axis_list=[0])

if 0:
    loop_apps.proj_cores_with_annotations( this_looper, axis_list=[0],color_dict={79:'g',83:'b'})
