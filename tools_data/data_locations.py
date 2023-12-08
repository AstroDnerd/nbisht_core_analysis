from starter2 import *

#
# guess the computer
#

machine = os.environ.get('machine',None)
if machine is None:
    if 'machine' in os.environ:
        machine = os.environ['machine']
    elif os.path.exists("/scratch1"):
        machine = 'Nazare'
    elif  os.path.exists("/data/cb1"):
        machine = 'cloudbreak'
    else:
        print("data_locations.py Bad error: cannot detrmine machine")

output_directory = "./plots_to_sort"
sim_u14=None        
if machine == 'Nazare':
    sf_path = "/scratch3/dcollins/Paper19/Datasets/VelocitySF/"
    coresets = {'ourset':'/scratch3/dcollins/Paper19/Datasets/'}  #ADDED here even though it should work in the place of the last 'else'
    sim_u05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128-Beta0.2'
    sim_u10 = '/archive2/dcollins4096/Paper19/u10_r4_l4_128-Beta2/GravPotential'
    sim_u11 = '/archive2/dcollins4096/Paper19/u11_r4_l4_128-Beta20/GravPotential'
    
    sim_u101 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128-Beta0.2'
    sim_u102 = '/archive2/dcollins4096/Paper19/u10_r4_l4_128-Beta2/GravPotential'
    sim_u103 = '/archive2/dcollins4096/Paper19/u11_r4_l4_128-Beta20/GravPotential'

    sim_u201 =  '/scratch3/dcollins/Paper19/Simulations_Good/u05-beta0.2'
    sim_u202 =  '/scratch3/dcollins/Paper19/Simulations_Good/u202-beta2'
    sim_u203 =  '/scratch3/dcollins/Paper19/Simulations_Good/u203-beta20'

    sim_u301 = sim_u201
    sim_u302 = sim_u202
    sim_u303 = sim_u203

    sim_u501 = sim_u201
    sim_u502 = sim_u202
    sim_u503 = sim_u203

    mountain_top = {'u301':"/scratch3/dcollins/Paper19/Datasets/mountain_tops/u301_new_tracks_take_9c.h5",
                    "u302":"/scratch3/dcollins/Paper19/Datasets/mountain_tops/u302_new_tracks_take_9b.h5",
                    "u303":"/scratch3/dcollins/Paper19/Datasets/mountain_tops/u303_new_tracks_take_9b.h5",
                    "t02":"/scratch3/dcollins/Paper19/Datasets/mountain_tops/t02_mountain_top.h5"}

elif machine == 'cloudbreak':
    sim_u05 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u10 = '/data/cb1/Projects/P19_CoreSimulations/u10_r4_l4_128-Beta2/GravPotential'
    sim_u11 = '/data/cb1/Projects/P19_CoreSimulations/u11_r4_l4_128-Beta20/GravPotential'
    sim_u14 = '/data/cb1/Projects/P19_CoreSimulations/u14_sphere'


    sim_u101 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u102 = '/data/cb1/Projects/P19_CoreSimulations/u10_r4_l4_128-Beta2/GravPotential'  
    sim_u103 = '/data/cb1/Projects/P19_CoreSimulations/u11_r4_l4_128-Beta20/GravPotential' 

    #The 200 series continues u05, u10, and u11 to 1tff, 
    #and ties the analysis to cores found at 1 tff.
    sim_u201 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u202 = '/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential'  
    sim_u203 = '/data/cb1/Projects/P19_CoreSimulations/u203-Beta20/GravPotential' 
    sim_t02 = sim_u202

    #The 300 series uses the 200 series data, but extracts cores in a different manner.
    sim_u301 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u302 = '/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential'  
    sim_u303 = '/data/cb1/Projects/P19_CoreSimulations/u203-Beta20/GravPotential' 

    mountain_top = {'u301':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u301_new_tracks_take_9c.h5",
                    "u302":"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u302_new_tracks_take_9b.h5",
                    "u303":"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/u303_new_tracks_take_9b.h5",
                    "t02":"/data/cb1/Projects/P19_CoreSimulations/CoreSets/mountain_tops/t02_mountain_top.h5"}
    sim_u501 = sim_u201
    sim_u502 = sim_u202
    sim_u503 = sim_u203

    #path for structure functions
    sf_path = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/VelocitySF"
    coresets = {'ourset':'/data/cb1/Projects/P19_CoreSimulations/CoreSets/'}
else:
    sf_path = "/scratch3/dcollins/Paper19/Datasets/VelocitySF/"
    coresets = {'ourset':'/scratch3/dcollins/Paper19/Datasets/'}
    sim_u05 = None
    sim_u10 = None
    sim_u11 = None


    sim_u101 = None
    sim_u102 = None
    sim_u103 = None
    sim_u201 = None
    sim_u202 = None
    sim_u203 = None
    sim_u301 = None
    sim_u302 = None
    sim_u303 = None


sims={'u05':sim_u05,'u10':sim_u10,'u11':sim_u11,'u101':sim_u101,'u102':sim_u102,'u103':sim_u103}
sims.update({'u201':sim_u201,'u202':sim_u202,'u203':sim_u203})
sims.update({'u301':sim_u201,'u302':sim_u202,'u303':sim_u203})
sims.update({'u501':sim_u501,'u502':sim_u502,'u503':sim_u503})
sims.update({'u601':sim_u501,'u602':sim_u502,'u603':sim_u503})
sims.update({'u701':sim_u501,'u702':sim_u502,'u703':sim_u503})
sims.update({'u901':sim_u501,'u902':sim_u502,'u903':sim_u503})
sims['u14']=sim_u14
sims['t02']=sim_t02

new_tracks = coresets['ourset']+"/NewTracks"


#peaks_u05 = 'datasets_small/u05_0125_peaklist.h5'
#peaks_u10 = 'datasets_small/u10_0082_peaklist.h5'
#peaks_u11 = 'datasets_small/u11_0088_peaklist.h5'
#peaks_u101 = 'datasets_small/u101_0080_peaklist.h5'
#peaks_u102 = 'datasets_small/u102_0080_peaklist.h5'
#peaks_u103 = 'datasets_small/u103_0080_peaklist.h5'
#peaks_u202 = 'datasets_small/u202_0118_peaklist.h5'
#peaks_u203 = 'datasets_small/u203_0107_peaklist.h5'
#peaks_u14  = 'datasets_small/u14_0025_peaklist.h5'
#peaks_u301 = 'datasets_small/u301_0125_peaklist.h5'
#peaks_u302 = 'datasets_small/u302_0118_peaklist.h5'
#peaks_u303 = 'datasets_small/u303_0107_peaklist.h5'
#peak_list = {'u05':peaks_u05,'u10':peaks_u10,'u11':peaks_u11, 'u101':peaks_u101,'u102':peaks_u102,'u103':peaks_u103}
#peak_list.update( {'u201':peaks_u05,'u202':peaks_u202, 'u203':peaks_u203})
#peak_list.update( {'u301':peaks_u301,'u302':peaks_u302, 'u303':peaks_u303})
#peak_list.update( {'u701':peaks_u301,'u702':peaks_u302, 'u703':peaks_u303})
#peak_list['u14']=peaks_u14
#peak_list['t02']="datasets_small/t02_0060_peaklist.h5"
#
#
#target_frames={'u05':125,'u10':82,'u11':88,'u101':80,'u102':80,'u103':80,'u14':25}
#target_frames.update({'u201':125,'u202':118,'u203':107})
#target_frames.update({'u301':125,'u302':118,'u303':107})
#target_frames.update({'u401':125,'u402':118,'u403':107})
#target_frames.update({'u501':125,'u502':118,'u503':107})
#target_frames.update({'u601':125,'u602':118,'u603':107})
#target_frames.update({'u701':125,'u702':118,'u703':107})
#target_frames['t02']=60 #half a free fall time
#
#bad_particles_u05='datasets_small/u05_bad_particles.h5'
#bad_particles={'u05':bad_particles_u05,'u10':None,'u11':None}
#
#n_particles={'u05':'datasets_small/u05_n_particles.txt',
#             'u10':'datasets_small/u10_n_particles.txt',
#             'u11':'datasets_small/u11_n_particles.txt',
#             'u101':'datasets_small/u101_n_particles.txt',
#             'u102':'datasets_small/u102_n_particles.txt',
#             'u103':'datasets_small/u103_n_particles.txt',
#             'u301':'datasets_small/u301_n_particles.txt',
#             'u302':'datasets_small/u302_n_particles.txt',
#             'u303':'datasets_small/u303_n_particles.txt',
#             'u14':'datasets_small/u14_n_particles.txt'}
##
