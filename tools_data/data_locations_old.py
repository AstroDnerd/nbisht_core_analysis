"""
A container for code portability.
Put
export machine = my_machine_name
in your .bashrc
Shipsterns and Mullaghmore are two of my computers, copy that.
"""
from starter2 import *

#enzo_directory = '/archive2/dcollins4096/Paper19/u10_r4_l4_128'  #TO BEGIN WORK WITH NEW SIMS
enzo_directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'

sixteen_frame ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe'
snapshot_base ='/scratch1/dcollins/Paper19/Datasets/'
output_directory = "./plots_to_sort/"

if 'machine' in os.environ:
    if os.environ['machine'] in ['shipsterns']:
        enzo_directory = '/home/dcollins/scratch/u05-r4-l4-128/GraviPotential'  #CHANGE TO NEW SIMS
        snapshot_location = '/home/dcollins/scratch/Paper19/track_index_fix'
    if os.environ['machine'] in ['mullaghmore']:
        enzo_directory = '/Users/dcollins/scratch/P19/u05-r4-l4-128/GraviPotential'
        snapshot_location = '/home/dcollins/scratch/Paper19/track_index_fix'

