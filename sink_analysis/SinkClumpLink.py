import warnings
warnings.filterwarnings("ignore")
import h5py
import numpy as np
import random 
from starter2 import *
import json

def make_dir(dir_path):
        if not os.path.exists(dir_path):
                print("making directory:",dir_path)
                os.makedirs(dir_path)

def SinkClumpLinker(trackname, sink_trackname, do_projections = False):
    this_track = track_info.tracks[trackname]
    this_simname = this_track.sim_directory
    if os.path.exists(this_track.SinkClumpLink_fname):
        print("Sink Clump File exists!")
        return True
    #let's get all peaks
    peaklist_name = this_track.peak_fname
    f = h5py.File(peaklist_name,"r")
    peak = np.array(f['peaks'])
    f.close()

    this_frame = this_track.target_frame
    ds = yt.load("%s/DD%04d/data%04d"%(this_simname,this_frame,this_frame))
    all_data = ds.all_data()

    #load sink simulation
    sink_track = track_info.tracks[sink_trackname]
    sink_simname = sink_track.sim_directory
    sink_frame = sink_track.target_frame

    solar_mass = 5900*1.91e33
    sink_particle_mass = []
    sink_particle_xpos = []
    sink_particle_ypos = []
    sink_particle_zpos = []
    number_of_sink_particles = []
    sink_particle_index = []
    time_arr = []
    
    ds_sink = yt.load("%s/DD%04d/data%04d"%(sink_simname,sink_frame,sink_frame))
    all_data_sink = ds_sink.all_data()
    type_4_indices = np.where(np.array(all_data_sink[ds_sink.fields.all.particle_type])==4)[0]
    sink_particle_mass.append(np.array(all_data_sink[ds_sink.fields.all.particle_mass])[type_4_indices].tolist())
    sink_particle_xpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_x])[type_4_indices].tolist())
    sink_particle_ypos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_y])[type_4_indices].tolist())
    sink_particle_zpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_z])[type_4_indices].tolist())
    number_of_sink_particles.append(len(type_4_indices))
    sink_particle_index.append(np.array(all_data_sink[ds_sink.fields.all.particle_index])[type_4_indices].tolist())
    time_arr.append(ds.current_time)

    if do_projections==True:
        make_dir("plots_to_sort/SinkClumpPlots")
        #setup plot
        plot = yt.ProjectionPlot(ds, 'z', ("enzo","Density"), fontsize = 10)
        plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
        plot.annotate_scale(corner="upper_right")
        suppress_annotation = False

        #plot peak locations
        peak_densities = []
        peak_densities_abovecrit = []
        for peak_index in range(len(peak)):
            peak_pos = peak[peak_index]
            point_obj =  ds.point([peak_pos[0], peak_pos[1],peak_pos[2]])
            peak_density = point_obj['enzo', 'Density']
            peak_densities.append(float(peak_density))
            if float(peak_density)>10000:
                peak_densities_abovecrit.append(float(peak_density))
                plot.annotate_sphere((peak_pos[0], peak_pos[1],peak_pos[2]), radius=(10*np.log10(float(peak_density)), "um"), coord_system="data", circle_args={"color": "black"})

        print("Considering All peaks, (min,max,avg,number) of densities: ",min(peak_densities),max(peak_densities), sum(peak_densities)/len(peak_densities),len(peak_densities))
        print("Considering peaks above critical, (min,max,avg,number) of densities: ",min(peak_densities_abovecrit),max(peak_densities_abovecrit), sum(peak_densities_abovecrit)/len(peak_densities_abovecrit),len(peak_densities_abovecrit))

        #plot sink particles

        if len(type_4_indices)>0:
            print(len(sink_particle_mass[-1]), min(sink_particle_mass[-1]),max(sink_particle_mass[-1]))
        for sink_i in range(len(sink_particle_mass[-1])):
            plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i]),sink_particle_zpos[-1][sink_i]),
                '.', coord_system="data", text_args={"fontsize": 100, "color":'black'})
        if suppress_annotation ==False:
            offset = 0.01
            for sink_i in range(len(sink_particle_mass[-1])):
                plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i])+random.choice([-offset,0,offset]), sink_particle_zpos[-1][sink_i]),
                    '{:0.1e}'.format(sink_particle_mass[-1][sink_i]*solar_mass), coord_system="data", text_args={"fontsize": 50, "color":'blue'})

        #save plot
        plot.save("plots_to_sort/SinkClumpPlots/CorePositionsAll_DD"+str(this_frame)+".png")

        fig = plt.figure()
        plt.hist(peak_densities)
        plt.savefig("plots_to_sort/SinkClumpPlots/CorePeakDensityDistbn_DD"+str(this_frame)+".png")

        fig = plt.figure()
        plt.hist(peak_densities_abovecrit)
        plt.savefig("plots_to_sort/SinkClumpPlots/CorePeakHighDensityDistbn_DD"+str(this_frame)+".png")

        fig = plt.figure()
        plt.hist(sink_particle_mass[-1])
        plt.savefig("plots_to_sort/SinkClumpPlots/SinkMassDistbn_DD"+str(this_frame)+".png")


    #now save core ids and respective sink indices by location in dictionary
    sink_core_dict = {}

    #used_sinks = []
    #This is for last frame
    for peak_index in range(len(peak)):
        sink_core_dict[peak_index] = -1
        peak_pos = peak[peak_index]
        point_obj = ds.point([peak_pos[0], peak_pos[1],peak_pos[2]])
        peak_density = float(point_obj['enzo', 'Density'])
        min_dist = 1e20
        min_index = -1
        for sink_i in range(len(sink_particle_mass[-1])):
            nearest = np.sqrt((sink_particle_xpos[-1][sink_i]-peak_pos[0])**2 + (sink_particle_ypos[-1][sink_i]-peak_pos[1])**2 +(sink_particle_zpos[-1][sink_i]-peak_pos[2])**2)
            if nearest<min_dist:
                min_dist = nearest
                min_index = sink_i
        if min_dist<0.015:
            sink_core_dict[peak_index] = {'sink_mass':[sink_particle_mass[-1][min_index]], 'sink_particle_id':[int(sink_particle_index[-1][min_index])], 
                                        'sink_xpos': [sink_particle_xpos[-1][min_index]], 'sink_ypos': [sink_particle_ypos[-1][min_index]], 'sink_zpos': [sink_particle_zpos[-1][min_index]]}
            #used_sinks.append((int(sink_particle_index[-1][min_index]),sink_particle_xpos[-1][min_index],sink_particle_ypos[-1][min_index],sink_particle_zpos[-1][min_index],sink_particle_mass[-1][min_index]))

    all_sinks = []
    for i in range(len(sink_particle_mass[-1])): 
        all_sinks.append(((int(sink_particle_index[-1][i]),sink_particle_xpos[-1][i],sink_particle_ypos[-1][i],sink_particle_zpos[-1][i],sink_particle_mass[-1][i])))

    #print(used_sinks)
    #print(all_sinks)
    for frame_num in range(sink_frame-1,-1,-1):
        sink_particle_mass = []
        sink_particle_xpos = []
        sink_particle_ypos = []
        sink_particle_zpos = []
        number_of_sink_particles = []
        sink_particle_index = []
        ds_sink = yt.load("%s/DD%04d/data%04d"%(sink_simname,frame_num,frame_num))
        all_data_sink = ds_sink.all_data()
        type_4_indices = np.where(np.array(all_data_sink[ds_sink.fields.all.particle_type])==4)[0]
        sink_particle_mass.append(np.array(all_data_sink[ds_sink.fields.all.particle_mass])[type_4_indices].tolist())
        sink_particle_xpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_x])[type_4_indices].tolist())
        sink_particle_ypos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_y])[type_4_indices].tolist())
        sink_particle_zpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_z])[type_4_indices].tolist())
        number_of_sink_particles.append(len(type_4_indices))
        sink_particle_index.append(np.array(all_data_sink[ds_sink.fields.all.particle_index], dtype = int)[type_4_indices].tolist())
        for sp in sink_core_dict.keys():
            if sink_core_dict[sp]!=-1:
                last_x = sink_core_dict[sp]['sink_xpos'][-1]
                last_y = sink_core_dict[sp]['sink_ypos'][-1]
                last_z = sink_core_dict[sp]['sink_zpos'][-1]
                min_dist = 1e20
                min_index = -1
                for sink_index in range(len(type_4_indices)):
                    nearest = np.sqrt((sink_particle_xpos[-1][sink_index]-last_x)**2 + (sink_particle_ypos[-1][sink_index]-last_y)**2 +(sink_particle_zpos[-1][sink_index]-last_z)**2)
                    if nearest<min_dist:
                        min_dist = nearest
                        min_index = sink_index
                if min_dist<0.015:      
                    sink_core_dict[sp]['sink_mass'].append(sink_particle_mass[-1][min_index])
                    sink_core_dict[sp]['sink_particle_id'].append(int(sink_particle_index[-1][min_index]))
                    sink_core_dict[sp]['sink_xpos'].append(sink_particle_xpos[-1][min_index])
                    sink_core_dict[sp]['sink_ypos'].append(sink_particle_ypos[-1][min_index])
                    sink_core_dict[sp]['sink_zpos'].append(sink_particle_zpos[-1][min_index])

    for sp in sink_core_dict.keys():
        if sink_core_dict[sp]!=-1:
            sink_core_dict[sp]['sink_mass'].reverse()
            sink_core_dict[sp]['sink_particle_id'].reverse()
            sink_core_dict[sp]['sink_xpos'].reverse()
            sink_core_dict[sp]['sink_ypos'].reverse()
            sink_core_dict[sp]['sink_zpos'].reverse()
    with open(this_track.SinkClumpLink_fname, 'w') as fp:
        json.dump(sink_core_dict, fp)

