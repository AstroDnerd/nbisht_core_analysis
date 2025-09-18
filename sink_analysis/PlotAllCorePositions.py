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

def make_dir(dir_path):
        if not os.path.exists(dir_path):
                print("making directory:",dir_path)
                os.makedirs(dir_path)

path_to_output_plots = "/home/nbisht/cdbreak_desktop/nikhilb_home/results/plots_to_sort"
#get core positions
f = h5py.File("datasets_small/nb101_0088_peaklist.h5","r")
peak = np.array(f['peaks'])
f.close()

#get number of tracers in core
f = open("datasets_small/nb101_n_particles.txt",'r')
lines = [line.rstrip() for line in f]
num = [int(l.split(' ')[1]) for l in lines]
f.close()

if 'this_simname' not in dir():
    this_simname='nb101'

sink_simname = '/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare'

frame = dl.target_frames[this_simname]
ds = yt.load("%s/DD%04d/data%04d"%(dl.sims[this_simname],frame,frame))
all_data = ds.all_data()
print(ds.field_list)

#setup plot
plot = yt.ProjectionPlot(ds, 'z', ("enzo","Density"), fontsize = 10)
plot.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
plot.annotate_scale(corner="upper_right")

#plot peak locations
peak_densities = []
peak_densities_above10000 = []
for peak_index in range(len(num)):
    peak_pos = peak[peak_index]
    point_obj = ds.point([peak_pos[0], peak_pos[1],peak_pos[2]])
    peak_density = point_obj['enzo', 'Density']
    peak_densities.append(float(peak_density))
    if float(peak_density)>10000:
        peak_densities_above10000.append(float(peak_density))
        plot.annotate_sphere((peak_pos[0], peak_pos[1],peak_pos[2]), radius=(10*np.log10(float(peak_density)), "um"), coord_system="data", circle_args={"color": "black"})

print(min(peak_densities),max(peak_densities), sum(peak_densities)/len(peak_densities),len(peak_densities))
print(min(peak_densities_above10000),max(peak_densities_above10000), sum(peak_densities_above10000)/len(peak_densities_above10000),len(peak_densities_above10000))

#plot sink particles
solar_mass = 5900*1.91e33
sink_particle_mass = []
sink_particle_xpos = []
sink_particle_ypos = []
sink_particle_zpos = []
number_of_sink_particles = []
sink_particle_index = []
time_arr = []
suppress_annotation = False
ds_sink = yt.load("%s/DD%04d/data%04d"%(sink_simname,frame,frame))
all_data_sink = ds_sink.all_data()
type_4_indices = np.where(np.array(all_data_sink[ds_sink.fields.all.particle_type])==4)[0]
sink_particle_mass.append(np.array(all_data_sink[ds_sink.fields.all.particle_mass])[type_4_indices].tolist())
sink_particle_xpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_x])[type_4_indices].tolist())
sink_particle_ypos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_y])[type_4_indices].tolist())
sink_particle_zpos.append(np.array(all_data_sink[ds_sink.fields.all.particle_position_z])[type_4_indices].tolist())
number_of_sink_particles.append(len(type_4_indices))
sink_particle_index.append(np.array(all_data_sink[ds_sink.fields.all.particle_index])[type_4_indices].tolist())
if len(type_4_indices)>0:
    print(len(sink_particle_mass[-1]), min(sink_particle_mass[-1]),max(sink_particle_mass[-1]))
time_arr.append(ds.current_time)
for sink_i in range(len(sink_particle_mass[-1])):
    plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i]),sink_particle_zpos[-1][sink_i]),
        '.', coord_system="data", text_args={"fontsize": 100, "color":'black'})
if suppress_annotation ==False:
    offset = 0.01
    for sink_i in range(len(sink_particle_mass[-1])):
        plot.annotate_text((sink_particle_xpos[-1][sink_i], float(sink_particle_ypos[-1][sink_i])+random.choice([-offset,0,offset]), sink_particle_zpos[-1][sink_i]),
            '{:0.1e}'.format(sink_particle_mass[-1][sink_i]*solar_mass), coord_system="data", text_args={"fontsize": 50, "color":'blue'})

#plot.zoom(zoom_level)
#save plot
plot.save(path_to_output_plots+"/CorePositionsAll_DD"+str(frame)+".png")

fig = plt.figure()
plt.hist(peak_densities)
plt.savefig(path_to_output_plots+"/CorePeakDensityDistbn_DD"+str(frame)+".png")

fig = plt.figure()
plt.hist(peak_densities_above10000)
plt.savefig(path_to_output_plots+"/CorePeakHighDensityDistbn_DD"+str(frame)+".png")

fig = plt.figure()
plt.hist(sink_particle_mass[-1])
plt.savefig(path_to_output_plots+"/SinkMassDistbn_DD"+str(frame)+".png")

sink_core_dict = {}

for peak_index in range(len(num)):
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

for frame_num in range(frame-1,-1,-1):
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
import json
with open('/data/cb1/nbisht/anvil_scratch/projects/128/B2_sink_nazare/datasets/sink_clump_position.json', 'w') as fp:
    json.dump(sink_core_dict, fp)

