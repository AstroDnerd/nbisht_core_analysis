from starter2 import *

#import three_loopers_mountain_top as TLM

#this_looper=new_looper
if 1:
    frame_ind=0
    #this_looper = TLM.loops['u301']
    x = this_looper.tr.track_dict['x'][:,frame_ind]
    y = this_looper.tr.track_dict['y'][:,frame_ind]
    z = this_looper.tr.track_dict['z'][:,frame_ind]
    #x = this_looper.tr.track_dict['particle_pos_x'][:,frame_ind]
    #y = this_looper.tr.track_dict['particle_pos_y'][:,frame_ind]
    #z = this_looper.tr.track_dict['particle_pos_z'][:,frame_ind]
    dx=dy=dz=1./128
    i = np.floor( x/dx).astype('int')
    j = np.floor( y/dx).astype('int')
    k = np.floor( z/dx).astype('int')
    nx=128
    ny=128
    index=i + nx*(j+ny*k)
    print("Number of repeated particles:", index.size-np.unique(index).size)
    args_index = np.argsort(index)
    isrt = index[args_index]
    idif = isrt[1:]-isrt[:-1]
    repeat = idif == 0
    vals = isrt[1:][repeat]
    also = isrt[:-1][repeat]
    print("Check that we picked up repeaters: ", (vals==also).all())
    cores1 = this_looper.tr.core_ids[args_index][1:][repeat]
    cores2 = this_looper.tr.core_ids[args_index][:-1][repeat]
    pairs=defaultdict(list)
    for c1,c2 in zip(cores1,cores2):
        pairs[c1].append(c2)
    particles1 = this_looper.tr.particle_ids[args_index][1:][repeat]
    particles2 = this_looper.tr.particle_ids[args_index][:-1][repeat]

    tr = this_looper.tr
    td = this_looper.tr.track_dict
if 1:
    delta_particle = td['particle_pos_y']-td['y']
    #delta_index = td['bucket'][:,0]-i
    bad = delta_index != 0
    fig,ax=plt.subplots(1,1)
    ax.scatter(td['particle_pos_x'][bad], td['particle_pos_y'][bad])
    fig.savefig('plots_to_sort/bad_positions.png')
    


