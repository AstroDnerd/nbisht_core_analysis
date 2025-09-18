from go import *
# for new ds.all_data() profile plot - 
import core_dump
def rs(arr):
    return arr.reshape(arr.shape + (1,)) 
def rs2(arr):
    return arr.reshape(arr.shape + (1,) + (1,)) 

import trackage
reload(trackage)
track = trackage.track_manager(None)
track.read('all_large_0000_0125.h5')
den0 = rs2(track['density'][:,0]) #track['density'].reshape(track['density'].shape + (1,)) 
den1 = rs2(track['density'][:,1]) #track['density'].reshape(track['density'].shape + (1,)) 
ke0 = rs2(track['kinetic_energy'][:,0])
ke1 = rs2(track['kinetic_energy'][:,1])
cell_v0= rs2(track['cell_volume'][:,0] )
cell_v1 = rs2(track['cell_volume'][:,1] )
# SWITCHED TO CELL_BACON
#data = {'density':den, 'my_vol':cell_v, 'kinetic_energy':ke} 
data = {'d1':den1,'d0':den0, 'my_vol':cell_v, 'kinetic_energy':ke} 
bbox = np.array([[0.,1.]]*3) 
ds = yt.load_uniform_grid(data, den.shape, length_unit="cm", bbox=bbox)
ad = ds.all_data() 
if 1:
    p=yt.PhasePlot(ad,'d0','d1',['my_vol'], x_bins=32,y_bins=32,weight_field=None,fractional=False)
    print(p.save('test_top81'))
if 0:
    # profile of particles that end in a specific core
    pdf = yt.create_profile(ad,field[0],'my_vol',weight_field=None,fractional=False)
    the_x = pdf.x
    the_y = pdf['my_vol']
    plt.clf()
    plt.plot(the_x,the_y)
    plt.xscale('log'); plt.yscale('log')
    plt.savefig('test_c14_c17_writeplot.png')
if 0:
    bbox = np.array([[0.,1.]]*3) 
    den1 = new_track['density'][:,0]
    den2 = new_track['density'][:,1]
    cell_v = new_track['cell_volume'][:,1]
    den1=rs(den1)
    den2=rs(den2)
    cell_v=rs(cell_v)
# SWITCHED TO CELL_BACON
    data = {'den1':den1,'den2':den2, 'my_vol':cell_v}  
    ds = yt.load_uniform_grid(data, den.shape, length_unit="cm", bbox=bbox)
    ad = ds.all_data() 

    pdf = yt.create_profile(ad,'den1',"my_vol",weight_field=None,fractional=False)
    plt.clf()
    plt.plot(pdf.x,pdf['my_vol'])
    plt.savefig('test2.png')
#p=yt.PhasePlot(ad,'den1','den2',['my_vol'],weight_field=None,fractional=False)
#outname = 'bayes_den_0000_0125_top81'
#p.save(outname)
#print(outname)
