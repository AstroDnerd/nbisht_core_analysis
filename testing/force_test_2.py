import xtra_energy
reload(xtra_energy)
frame = 0
pig = '/home/dcollins4096/PigPen'
directory = "/scratch1/dcollins/Paper19/u13_OrszagTang"

ke=[]
be=[]
te=[]
bwork=[]
twork=[]

for frame in range(5):
    fname = "%s/DD%04d/data%04d"%(directory,frame,frame)
    ds = yt.load(fname)
    xtra_energy.add_force_terms(ds)
    ad = ds.all_data()
    ke.append( (ad['cell_volume']*ad['kinetic_energy']).sum())
    te.append( (ad['cell_volume']*ad['thermal_energy']*ad['density']).sum())
    be.append( (ad['cell_volume']*ad['magnetic_energy']).sum())
    twork.append( (ad['cell_volume']*ad['pressure_work']).sum())
    bwork.append( (ad['cell_volume']*ad['mag_work']).sum())
plt.clf()
plt.plot(ke,c='r')
plt.plot(be,c='b')
plt.plot(te,c='y')
plt.plot(bwork,'b:')
plt.plot(twork,'y:')
plt.savefig('%s/energy.png'%pig)
