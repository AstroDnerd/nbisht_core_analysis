import warnings
warnings.filterwarnings("ignore")

import matplotlib.ticker as tck
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=False) # Use LaTeX font
import h5py
import numpy as np
import random 
from starter2 import *
import track_loader as TL
import json
import ucore
import monster
import pandas as pd


path_to_output_plots = "/home/nbisht/cdbreak_desktop/nikhilb_home/results/plots_to_sort"
path_to_sink_clump = '/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/sink_clump_position_with_sink_totaldensity.json'
nonsink_trackname = 'nb101'
sink_trackname = 'nb102'

sys.path.insert(1, 'analysis_pipe/')
import rho_time

#get core data only for valid tsing cores
def export_data_for_ML(trackname, plot_rho = False):
    this_track = track_info.tracks[trackname]
    with open(path_to_sink_clump, 'r') as fp:
        sink_core_dict = json.load(fp)
    monster.load([trackname])
    mon = monster.closet[trackname]
    core_list = list(mon.tsing.tsing_core.keys())
    if plot_rho == True:
        rho_time.run(trackname, core_list, tsing_dict = mon.tsing.tsing_core, tend_dict = mon.tsing.tend_core)

    #Making the dataset, Clump_id is -1 if not a core_tracer and other clump values are individual values
    data = {'Clump_id':[], 'Particle_id':[], 'X':[], 'Y':[], 'Z':[], 'Density': [], 'Vx':[], 'Vy':[], 'Vz':[], 't_hard': [],
            'O_Clump_X':[], 'O_Clump_Y':[], 'O_Clump_Z':[], 'O_Clump_Vx':[], 'O_Clump_Vy':[], 'O_Clump_Vz':[],  'O_Clump_density':[], 'O_t_end':[]}
    t_sing_actual_frame = []
    t_end_actual_frame = []
    for core_id in core_list:
        ms = mon.get_ms(core_id)
        print('data_exporter',core_id)
        ms.particle_pos(core_id)
        t_sing_actual_frame.append(mon.tsing.tsing_frame[core_id])
        tsing_fr = mon.get_frame_index(mon.tsing.tsing_frame[core_id])
        t_end_actual_frame.append(mon.tsing.tend_frame[core_id])
        tend_fr = mon.get_frame_index(mon.tsing.tend_frame[core_id])
        tsing_val = mon.tsing.tsing_core[core_id]
        tend_val = mon.tsing.tend_core[core_id]
        pid = ms.particle_id
        number_of_particles = int(ms.nparticles)
        xpos = ms.particle_x[:,tsing_fr]
        ypos = ms.particle_y[:,tsing_fr]
        zpos = ms.particle_z[:,tsing_fr]
        rho = ms.density[:,tsing_fr]
        vx = ms.raw_vx[:,tsing_fr]
        vy = ms.raw_vy[:,tsing_fr]
        vz = ms.raw_vz[:,tsing_fr]
        clump_pos = ms.mean_center[:,tend_fr]
        clump_vx = ms.mean_vx[tend_fr]
        clump_vy = ms.mean_vx[tend_fr]
        clump_vz = ms.mean_vz[tend_fr]
        clump_dens = np.mean(ms.density[:,tend_fr])
        #insertion
        data['Clump_id'].extend([core_id]*number_of_particles)
        data['Particle_id'].extend(list(pid))
        data['X'].extend(list(xpos))
        data['Y'].extend(list(ypos))
        data['Z'].extend(list(zpos))
        data['Density'].extend(list(rho))
        data['Vx'].extend(list(vx))
        data['Vy'].extend(list(vy))
        data['Vz'].extend(list(vz))
        data['t_hard'].extend([tsing_val]*number_of_particles)
        data['O_Clump_X'].extend([clump_pos[0]]*number_of_particles)
        data['O_Clump_Y'].extend([clump_pos[1]]*number_of_particles)
        data['O_Clump_Z'].extend([clump_pos[2]]*number_of_particles)
        data['O_Clump_Vx'].extend([clump_vx]*number_of_particles)
        data['O_Clump_Vy'].extend([clump_vy]*number_of_particles)
        data['O_Clump_Vz'].extend([clump_vz]*number_of_particles)
        data['O_Clump_density'].extend([clump_dens]*number_of_particles)
        data['O_t_end'].extend([tend_val]*number_of_particles)

    #Now for non core particles
    import time
    start_time = time.time()
    t_begin = int(np.around(np.mean(t_sing_actual_frame))) #Choose this as the initial frame for non_core particles
    t_final = int(np.around(np.max(t_end_actual_frame))) #Choose this as the final frame for non_core particles
    #Open these frames and get particle data
    ds_begin = yt.load("%s/DD%04d/data%04d"%(this_track.sim_directory,t_begin,t_begin))
    all_data_begin = ds_begin.all_data()
    core_pids = data['Particle_id']
    all_pids = np.array(all_data_begin[ds_begin.fields.all.particle_index])
    indices_core_pids = [np.where(all_pids == i)[0][0] for i in core_pids]
    non_core_pids = np.delete(all_pids,indices_core_pids)
    indices_non_core_pids = [np.where(all_pids == i)[0][0] for i in non_core_pids]
    

    data['Clump_id'].extend([-1]*len(non_core_pids))
    data['Particle_id'].extend(list(non_core_pids))
    data['X'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_position_x])[indices_non_core_pids].tolist())
    data['Y'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_position_y])[indices_non_core_pids].tolist())
    data['Z'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_position_z])[indices_non_core_pids].tolist())

    non_pos = np.array(all_data_begin[ds_begin.fields.all.particle_position])[indices_non_core_pids]
    peak_density = ds_begin.find_field_values_at_points([("gas", "density")],non_pos)
    data['Density'].extend(peak_density.tolist())

    data['Vx'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_velocity_x])[indices_non_core_pids].tolist())
    data['Vy'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_velocity_y])[indices_non_core_pids].tolist())
    data['Vz'].extend(np.array(all_data_begin[ds_begin.fields.all.particle_velocity_z])[indices_non_core_pids].tolist())
    data['t_hard'].extend([np.mean(data['t_hard'])]*len(non_core_pids))

    del ds_begin
    del all_data_begin

    ds_final = yt.load("%s/DD%04d/data%04d"%(this_track.sim_directory,t_final,t_final))
    all_data_final = ds_final.all_data()
    all_pids = np.array(all_data_final[ds_final.fields.all.particle_index])
    indices_non_core_pids = [np.where(all_pids == i)[0][0] for i in non_core_pids]
    
    data['O_Clump_X'].extend(np.array(all_data_final[ds_final.fields.all.particle_position_x])[indices_non_core_pids].tolist())
    data['O_Clump_Y'].extend(np.array(all_data_final[ds_final.fields.all.particle_position_y])[indices_non_core_pids].tolist())
    data['O_Clump_Z'].extend(np.array(all_data_final[ds_final.fields.all.particle_position_z])[indices_non_core_pids].tolist())
    data['O_Clump_Vx'].extend(np.array(all_data_final[ds_final.fields.all.particle_velocity_x])[indices_non_core_pids].tolist())
    data['O_Clump_Vy'].extend(np.array(all_data_final[ds_final.fields.all.particle_velocity_y])[indices_non_core_pids].tolist())
    data['O_Clump_Vz'].extend(np.array(all_data_final[ds_final.fields.all.particle_velocity_z])[indices_non_core_pids].tolist())

    non_pos = np.array(all_data_final[ds_final.fields.all.particle_position])[indices_non_core_pids]
    peak_density = ds_final.find_field_values_at_points([("gas", "density")],non_pos)
    data['O_Clump_density'].extend(peak_density.tolist())
    
    data['O_t_end'].extend([np.max(data['O_t_end'])]*len(non_core_pids))

    del ds_final
    del all_data_final

    #Open dataframe
    df_name = this_track.export_to_ML_fname
    df = pd.DataFrame(data)
    df.to_csv(df_name)
    
    diff = time.time() - start_time
    print('Time taken by looper:',diff)

    pdb.set_trace()
    del df

    
    
export_data_for_ML(nonsink_trackname, plot_rho = False)



