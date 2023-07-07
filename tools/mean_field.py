from starter2 import *

#import three_loopers_u500 as TL
import track_loader as TL
sim_list=['u501','u502','u503']
if 'B' not in dir():
    B = {}
def compute_mean_field(trackname):
    print("Compute Mean Field.")
    TL.load_tracks(trackname)
    loop = TL.loops[trackname]
    ds = loop.load(0)
    ad = ds.all_data()
    Bx = ad['magnetic_field_x'].mean()
    By = ad['magnetic_field_y'].mean()
    Bz = ad['magnetic_field_z'].mean()
    B[trackname] = np.sqrt(Bx*Bx+By*By+Bz*Bz)

B['u501']=11.20998243
B['u502']=3.5449077
B['u503']=1.12099824

B['b001']=B['u501']
B['b002']=B['u502']
B['b003']=B['u503']
