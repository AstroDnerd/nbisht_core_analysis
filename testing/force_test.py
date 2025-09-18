
import xtra_energy
reload(xtra_energy)
frame = 0
pig = '/home/dcollins4096/PigPen'
directory = "/scratch1/dcollins/Paper19/u13_OrszagTang"
fname = "%s/DD%04d/data%04d"%(directory,frame,frame)
ds = yt.load(fname)
xtra_energy.add_force_terms(ds)
density = ds.all_data()['density']
g=ds.index.grids[0]
dx = 1./16
X,Y = np.mgrid[0.5*dx:(1):dx, 0.5*dx:(1):dx]
#B0 = 1./np.sqrt(4*np.pi)
B0 = 1.
fourpi = np.pi*4
twopi = np.pi*2
V0 = 1.0
vmag = g['velocity_magnitude'][:,:,0].v
#plave(vmag,'%s/ot_vel_%04d.png'%(pig,frame))
vx0 = -V0*np.sin(2*np.pi*Y)
vy0 =  V0*np.sin(2*np.pi*X)
vm0 = np.sqrt(vx0**2+vy0**2)
bx = -B0*np.sin(twopi*Y)
by =  B0*np.sin(fourpi*X)
bm = np.sqrt(bx**2+by**2)
jztest = B0*fourpi*np.cos(fourpi*X)+B0*twopi*np.cos(twopi*Y)
work =( vx0*(-1*jztest*by) + vy0*(jztest*bx) ) /fourpi
#plave(vm0,'%s/ot_vel_an_%04d.png'%(pig,frame))
wrk=g['mag_work'][:,:,0].v
plave(wrk,"%s/ot_mag_work_%04d.png"%(pig,frame))
plave(work, "%s/ot_mag_work_test_%04d.png"%(pig,frame))
#plave(g['current_z'][:,:,0].v,"%s/ot_current_z_%04d.png"%(pig,frame))
#plave(jztest,"%s/ot_current_z_test_%04d.png"%(pig,frame))
#plave((jztest-g['current_z'][:,:,0].v)/jztest.max(),"%s/ot_current_test2_%04d.png"%(pig,frame))
#plave(jztest,"%s/ot_current_test2_%04d.png"%(pig,frame))
#bmag = g['magnetic_field_strength'][:,:,0].v
#plave(bmag,'%s/ot_bmag_%04d.png'%(pig,frame))
#plave(bm,'%s/ot_bmag_2_%04d.png'%(pig,frame))

