
from starter2 import *
import xtra_energy
import pcolormesh_helper as pch 
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import pcolormesh_helper as pch
import colors
import movie_frames 
G = colors.G
plt.close('all')

from collections import defaultdict


class theparticles():
    def __init__(self,the_loop):
        print('init complete')
        self.this_looper = the_loop


# exploring this tool, ignore
class breakallloops(Exception):
    pass
#if nd-1 == numOne:
#    print('breaking')
#    raise breakallloops

#except breakallloops:
#    pass



import track_loader as TL
#sim_list=['m0230', 'm0231', 'm0232', 'm0233', 'm0234', 'm0235', 'm0236', 'm0237', 'm0238', 'm0239', 'm0240', 'm0241', 'm0242',\
#      'm0243', 'm0244', 'm0245', 'm0246', 'm0247', 'm0250', 'm0260', 'm0270', 'm0280', 'm02100', 'm02110']  
sim_list = ['m0230', 'm0231', 'm0232']
tframe_230 = defaultdict(set)
tframe_231 = defaultdict(set)
tframe_232 = defaultdict(set)
listofdicts = [tframe_230, tframe_231, tframe_232]

if 'targets' not in dir():
    targets=True
TL.load_tracks(sim_list)

if 'targets2' not in dir() or targets:
    for num, sim in enumerate(sim_list):
        all_cores=np.unique(TL.loops[sim].tr.core_ids)
        core_list=list(all_cores) 
        #core_list = [123]
        if 0:  #announce!
            print("The core count of target frame %s is %d"%(sim, len(core_list)))
            print("The newborn tags are: ") 
            print(core_list) 
        
        # in a definition of the class
        for core_id in core_list:
            this_looper = TL.loops[sim]
            ms = trackage.mini_scrubber(this_looper.tr,core_id,do_velocity=False)
            the_particles = ms.particle_id
            
            listofdicts[num][core_id].update(the_particles)
            #print('the particles for %d are: '%core_id,the_particles)
    
    def rename_keys_by_numbers(d):
        return {i: value for i, (key, value) in enumerate(d.items())}
    print('keys by core_id', listofdicts)

    numbereddicts = []
    for numof, dicts in enumerate(listofdicts):
        rename = rename_keys_by_numbers(dicts)
        numbereddicts.append(rename)
    print('keys by number', numbereddicts)




    # PLAY/FIGURE THINGS OUT
    sharingsall = defaultdict(list)
    sharingscores = defaultdict(list)
    
    nd = len(numbereddicts)  #number of dictionaries 
    #ncd = len(numbereddicts[nd].keys())  #number of cores per dictionaries
    #nppc = numbereddicts[nd][ndc]  #len for number of particles per cores

    for numOne, dictio in enumerate(numbereddicts):
        ncd = numbereddicts[numOne].keys()
        print('dictionary number ',numOne)  
        
        for numTwo, dictioCore in enumerate(ncd):
            print('core number ',numTwo)  
            partsInCoreA = numbereddicts[numOne][numTwo] 
            
            if numOne+1 < nd:
                partsInCoreB = numbereddicts[numOne+1][numTwo] 
                
                N = len(partsInCoreA.intersection(partsInCoreB))
                print('firstframe sharing',N)
                sharingscores[numTwo].append(N)

            # test! edit, need to account for when numOne and numTwo reaches limit
            if numOne+1 == nd: 
                partsInCoreC = numbereddicts[numOne][numTwo] 
                N = len(partsInCoreB.intersection(partsInCoreC))
                print('final sharing',N)
                sharingscores[numTwo].append(N)
               
            if numOne+1 < nd:
                ncdn = numbereddicts[numOne+1].keys()
                if numTwo+1 < len(ncdn): 
                    partsInCoreB = numbereddicts[numOne+1][numTwo+1] 
                    N = len(partsInCoreA.intersection(partsInCoreB))
                    print('next frame sharing',N)
                    sharingscores[numTwo].append(N)
    
            # test, edit, need to account for when numOne and numTwo reaches limit 
            if numOne+1 == nd: 
                if numTwo+1 == len(ncdn): 
                    partsInCoreC = numbereddicts[numOne][numTwo] 
                    N = len(partsInCoreB.intersection(partsInCoreC))
                    print('next frame sharing',N)
                    sharingscores[numTwo].append(N)

        sharingsall[numOne].append(sharingscores[numTwo])
    pdb.set_trace()
    




# PURPOSE
        # want to track: targetframe, cores, particles in cores
        # in one tframe, if two cores, are the cores particles all different? they should.
        # in the next tframe, what happened to the previous tframes's cores?
        #     do they have the same number of particles?
        #      did they loose particles?
        #      did they gain new particles?
        # in the next tframe, are there new cores?
        #     are all of the new cores made of new particles? or do they come from other cores (EDIT)

'''
        if 0:  # for all particles in a target frame!
            all_particles = np.unique(TL.loops[sim].tr.particle_ids) 
            particle_list = list(all_particles) 
            print("The particle count of target frame %s is %d"%(sim, len(particle_list)))
            print("This does not give us information about which particles belong to which cores")
            print("The particles are: ") 
            print(particle_list)
'''


