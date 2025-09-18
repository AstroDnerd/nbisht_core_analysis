from dtools.starter1 import *
import yt
import read_ak as rak
reload(rak)
import build_dataset as bs
reload(bs)

import add_tracer_tool as AT
reload(AT)

fields = [('enzo','Density'), ('enzo','x-velocity'),('enzo','y-velocity'),('enzo','z-velocity'),('enzo','Bx'),('enzo','By'),('enzo','Bz')]
fields = [('enzo','Density'), ('enzo','x-velocity'),('enzo','y-velocity'),('enzo','z-velocity')]
if 1:
    suite = "mach4"
    size = 1024
    setname = "/anvil/scratch/x-ux454321/p19b_high_res/IC_1024/mach4_xi0.5_b0/DD0090/data0090"
    ds= yt.load(setname)
    cg = ds.covering_grid(0,[0.0]*3,[size]*3)
    print('make cg')
    stuff1024=[]
    for field in fields:
        print('make',field)
        stuff1024.append( cg[field])

if 0:
    #read the cubes
    directory = "/anvil/scratch/x-ux454321/p19b_high_res/IC_1024/B02"
    base = "%s/cube080"%directory
    if 'stuff1024' not in dir():
        print('read the cubes.')
        stuff1024 = rak.read_alexei(base,Isothermal=True, dtype='<f8')
    output_directory='/anvil/scratch/x-ux454321/p19b_high_res/assembler_1024/B02'
    for n in range(len(stuff1024)):
        stuff1024[n].shape=[1024]*3

if 1:
    output_directory = "/anvil/scratch/x-ux454321/p19b_high_res/take1/mach4_1024"

    sl = tuple([slice(0,512)]*3)
    sl = tuple([slice(None)]*3)
    bs.write_enzo_set(stuff1024[0][sl],output_directory,'p19b_%s_density.1024'%suite)
    bs.write_enzo_set(stuff1024[1][sl],output_directory,'p19b_%s_vx.1024'%suite)
    bs.write_enzo_set(stuff1024[2][sl],output_directory,'p19b_%s_vy.1024'%suite)
    bs.write_enzo_set(stuff1024[3][sl],output_directory,'p19b_%s_vz.1024'%suite)
    bs.write_enzo_set(stuff1024[4][sl],output_directory,'p19b_%s_bx.1024'%suite,extend=0 )
    bs.write_enzo_set(stuff1024[5][sl],output_directory,'p19b_%s_by.1024'%suite,extend=1 )
    bs.write_enzo_set(stuff1024[6][sl],output_directory,'p19b_%s_bz.1024'%suite,extend=2 )

if 0:
    #DO ONE STEP.
    #ADD TRACERS
    output_directory='/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Run1024_4096'
    particle_base="/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Particles1024_4096"
    setname = "%s/DD0000/data0000"%output_directory
    frame=0
    ds = yt.load(setname)
    AT.add_particles(ds,setname,particle_base, outnumber = frame)


if 0:
    #read the dataset, make sure its ok.
    names = ['density','vx','vy','vz','Bx','By','Bz']
    fig,axes=plt.subplots(3,3)
    axlist=axes.flatten()
    for nc,cube in enumerate(stuff1024):
        print(cube.size/1024**3)
        print(names[nc])
        cube.shape = 1024,1024,1024
        axlist[nc].imshow(cube.sum(axis=0))
    fig.savefig('%s/cubes'%(plot_dir))

if 0:
    #compare to B02 along z, which works
    #ds128  = yt.load('/anvil/scratch/x-ux454321/p78c_high_res/128/B02/DD0000/data0000')
    #sim = 'b2'
    #ds128  = yt.load('/anvil/scratch/x-ux454321/p78c_high_res/128/B2/DD0000/data0000')
    sim = 'b02'
    ds128  = yt.load('/anvil/scratch/x-ux454321/p78c_high_res/128/B02/DD0000/data0000')

    axdir=2
    proj = ds128.proj('density',axdir)
    fields = [('enzo','Density'), ('enzo','x-velocity'),('enzo','y-velocity'),('enzo','z-velocity'),('enzo','Bx'),('enzo','By'),('enzo','Bz')]
    fig,axes = plt.subplots(3,3)
    axlist=axes.flatten()
    frb = proj.to_frb(1,[128,128],[0.5]*3)
    for nf, field in enumerate(fields):
        print('field',field)
        axlist[nf].imshow(frb[field])
    fig.savefig('%s/cubes_enzo_%s_ax%d'%(plot_dir,sim,axdir))



if 0:
    #make a test cube.
    #dx=1/64
    dx=1/1024
    directory='/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Run1024'
    x,y,z=np.mgrid[0:1:dx, 0:1:dx,0:1:dx]
    rho = 1+0.1*np.sin(2*np.pi*(np.sqrt(x**2+y**2+z**2)))

    print('vx')
    vx = 0.1*np.sin(2*np.pi*(np.sqrt(x**2+y**2+z**2)))
    print('vy')
    vy = 0.1*np.sin(2*np.pi*(np.sqrt(x**2+y**2+z**2)))
    print('vz')
    vz = 0.1*np.sin(2*np.pi*(np.sqrt(x**2+y**2+z**2)))
    print('bx')
    bx = np.zeros( nar(rho.shape)+nar([1,0,0])) + 0.1
    print('by')
    by = np.zeros( nar(rho.shape)+nar([0,1,0])) + 0.2
    print('bz')
    bz = np.zeros( nar(rho.shape)+nar([0,0,1])) + 0.3



    bs.write_enzo_set(rho,directory,'density')
    bs.write_enzo_set(vx,directory,'vx')
    bs.write_enzo_set(vy,directory,'vy')
    bs.write_enzo_set(vz,directory,'vz')
    bs.write_enzo_set(bx,directory,'bx',dims = nar(rho.shape) + nar([1,0,0]) )
    bs.write_enzo_set(by,directory,'by',dims = nar(rho.shape) + nar([0,1,0]) )
    bs.write_enzo_set(bz,directory,'bz',dims = nar(rho.shape) + nar([0,0,1]) )






