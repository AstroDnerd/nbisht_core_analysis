
"""
A container for code portability.
Put
export machine = my_machine_name
in your .bashrc
Shipsterns and Mullaghmore are two of my computers, copy that.
"""
from starter2 import *


sim_u05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128' 
sim_u10 = '/archive2/luzlourdes/u10'
sim_u11 = '/archive2/luzlourdes/u11' 
sims= {'u05':sim_u05,'u10':sim_u10,'u11':sim_u11}


peaks_u05 = 'datasets_small/u05_125_peaklist.h5'
peaks_u10 = 'datasets_small/u10_082_peaklist.h5'
peaks_u11 = 'datasets_small/u11_088_peaklist.h5'
peak_list = {'u05':peaks_u05,'u10':peaks_u10,'u11':peaks_u11}


target_frames = {'u05':125,'u10':82,'u11':88}


u05_sixteen_frame ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe/*h5'
u10_every_ten = "/home/luzlourdes/scripts/p19_newscripts/p19_newscripts/core_storage/u10_tens/all_primitives_c*_nXXX0.h5"  #TEST
u11_every_ten = "/home/luzlourdes/scripts/p19_newscripts/p19_newscripts/core_storage/u11_tens/all_primitives_c*_nXXX0.h5"  #TEST
every_ten = {'u05':u05_sixteen_frame,'u10':u10_every_ten,'u11':u11_every_ten}


bad_particles_u05='datasets_small/u05_bad_particles.h5'
bad_particles={'u05':bad_particles_u05,'u10':None,'u11':None}


n_particles={'u05':'datasets_small/u05c_particles.txt',
             'u10':'datasets_small/u10c_particles.txt',
             'u11':'datasets_small/u11c_particles.txt'}
