import warnings
warnings.filterwarnings("ignore")
import h5py
import numpy as np
import random 
from starter2 import *
import json

def SinkVsNonSink_BRhoPlotter(trackname, sink_trackname):
    this_track = track_info.tracks[trackname]
    this_simname = this_track.sim_directory

    this_frame = this_track.target_frame
    ds = yt.load("%s/DD%04d/data%04d"%(this_simname,this_frame,this_frame))
    all_data = ds.all_data()

    #load sink simulation
    sink_track = track_info.tracks[sink_trackname]
    sink_simname = sink_track.sim_directory
    sink_frame = sink_track.target_frame
    got_right_frame = False
    #make sure sink frame chosen is same or close to non sink chosen frame
    while got_right_frame==False:
        ds_sink = yt.load("%s/DD%04d/data%04d"%(sink_simname,sink_frame,sink_frame))
        if ds_sink.current_time-ds.current_time >1e-6:
            sink_frame-=1
        elif ds_sink.current_time-ds.current_time <-1e-6:
            sink_frame+=1
        else:
            got_right_frame = True
    all_data_sink = ds_sink.all_data()

    solar_mass = 5900*1.91e33
    sink_particle_mass = []
    sink_particle_xpos = []
    sink_particle_ypos = []
    sink_particle_zpos = []
    number_of_sink_particles = []
    sink_particle_index = []
    time_arr = []

    type_4_indices = np.where(np.array(all_data_sink[ds_sink.fields.all.particle_type])==4)[0]
    sink_particle_mass.append(np.array(all_data_sink[ds_sink.fields.all.particle_mass])[type_4_indices].tolist())
    sink_particle_xpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_x])[type_4_indices].tolist())
    sink_particle_ypos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_y])[type_4_indices].tolist())
    sink_particle_zpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_z])[type_4_indices].tolist())
    number_of_sink_particles.append(len(type_4_indices))
    sink_particle_index.append(np.array(all_data_sink[ds_sink.fields.all.particle_index])[type_4_indices].tolist())
    time_arr.append(ds.current_time)

    #setup plot
    plot = yt.PhasePlot(all_data, ("gas", "density"), ("gas", "magnetic_field_strength"),  ("gas", "velocity_magnitude"), fontsize = 10)
    suppress_annotation = False
    
    '''
    for sink_i in range(len(sink_particle_mass[-1])):
        plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i]),sink_particle_zpos[-1][sink_i]),
            '.', coord_system="data", text_args={"fontsize": 100, "color":'black'})
    if suppress_annotation ==False:
        offset = 0.01
        for sink_i in range(len(sink_particle_mass[-1])):
            plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i])+random.choice([-offset,0,offset]), sink_particle_zpos[-1][sink_i]),
                '{:0.1e}'.format(sink_particle_mass[-1][sink_i]*solar_mass), coord_system="data", text_args={"fontsize": 50, "color":'blue'})
    '''

    #save plot
    plot.save("plots_to_sort/B_Rho_phasediag_DD"+str(this_frame)+".png")


SinkVsNonSink_BRhoPlotter('nb101', 'nb102')