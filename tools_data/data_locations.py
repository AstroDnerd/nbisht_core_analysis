
"""
A container for code portability.
Put
export machine = my_machine_name
in your .bashrc
Shipsterns and Mullaghmore are two of my computers, copy that.
"""
from starter2 import *

machine = None
if os.path.exists("/scratch1"):
    machine = 'Nazare'
elif os.path.exists("/data/cb1"):
    machine = 'Cloudbreak'
else:
    print("Bad error: cannot detrmine machine")


if machine == 'Cloudbreak':
    sim_u05 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u10 = '/data/cb1/Projects/P19_CoreSimulations/u10_r4_l4_128-Beta2/GravPotential'
    sim_u11 = '/data/cb1/Projects/P19_CoreSimulations/u11_r4_l4_128-Beta20/GravPotential'

    u05_every_ten = '/data/cb1/Projects/P19_CoreSimulations/CoreSets/u05_every_ten/*h5'
    u10_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u10_every_ten/u10_all_primitives_primitives_c*_nXXX0.h5"
    u11_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u11_every_ten/u11_all_primitives_primitives_c*_nXXX0.h5"

    n_particles={'u05':'datasets_small/u05_n_particles.txt',
                 'u10':'datasets_small/u10_n_particles.txt',
                 'u11':'datasets_small/u11_n_particles.txt'}
    
    peaks_u05 = 'datasets_small/u05_0125_peaklist.h5'
    peaks_u10 = 'datasets_small/u10_0082_peaklist.h5'
    peaks_u11 = 'datasets_small/u11_0088_peaklist.h5'

elif machine == 'Nazare':
    sim_u05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128' 
    sim_u10 = '/archive2/luzlourdes/u10'
    sim_u11 = '/archive2/luzlourdes/u11' 
    
    u05_every_ten ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe/*h5'
    u10_every_ten = "/home/luzlourdes/scripts/p19_newscripts/p19_newscripts/core_storage/u10_tens/all_primitives_c*_nXXX0.h5"  #TEST
    u11_every_ten = "/home/luzlourdes/scripts/p19_newscripts/p19_newscripts/core_storage/u11_tens/all_primitives_c*_nXXX0.h5"  #TEST

    n_particles={'u05':'datasets_small/u05c_particles.txt',
                 'u10':'datasets_small/u10c_particles.txt',
                 'u11':'datasets_small/u11c_particles.txt'}

    peaks_u05 = 'datasets_small/u05_125_peaklist.h5'
    peaks_u10 = 'datasets_small/u10_082_peaklist.h5'
    peaks_u11 = 'datasets_small/u11_088_peaklist.h5'


peak_list = {'u05':peaks_u05,'u10':peaks_u10,'u11':peaks_u11}
sims= {'u05':sim_u05,'u10':sim_u10,'u11':sim_u11}
every_ten = {'u05':u05_every_ten,'u10':u10_every_ten,'u11':u11_every_ten}
target_frames = {'u05':125,'u10':82,'u11':88}

bad_particles_u05='datasets_small/u05_bad_particles.h5'
bad_particles={'u05':bad_particles_u05,'u10':None,'u11':None}


