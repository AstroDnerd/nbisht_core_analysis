
from dtools.starter1 import *
import yt
import read_ak as rak
reload(rak)
import build_dataset as bs
reload(bs)

fields = [('enzo','Density'), ('enzo','x-velocity'),('enzo','y-velocity'),('enzo','z-velocity'),('enzo','Bx'),('enzo','By'),('enzo','Bz')]

if 1:
    suite = "B2"
    outdir = "P19B_B2"
    setname = "/anvil/scratch/x-ux454321/p19b_high_res/IC_512/%s/512/RS0000/restart0000"%suite
    ds512  = yt.load(setname)
    cg = ds512.covering_grid(0,[0.0]*3,[512]*3)
    print('make cg')
    stuff=[]
    for field in fields:
        print('make',field)
        stuff.append( cg[field])
    #output_directory='/anvil/scratch/x-ux454321/p19b_high_res/take1/Particles_%s'%suite
    output_directory = "/anvil/scratch/x-ux454321/p19b_high_res/take1/%s"%outdir


if 1:
    #build the set
    bs.write_enzo_set(stuff[0],output_directory,field='p19b_%s_density.512'%suite)
    bs.write_enzo_set(stuff[1],output_directory,field='p19b_%s_vx.512'%suite)
    bs.write_enzo_set(stuff[2],output_directory,field='p19b_%s_vy.512'%suite)
    bs.write_enzo_set(stuff[3],output_directory,field='p19b_%s_vz.512'%suite)
    bs.write_enzo_set(stuff[4],output_directory,field='p19b_%s_bx.512'%suite,extend=0 )
    bs.write_enzo_set(stuff[5],output_directory,field='p19b_%s_by.512'%suite,extend=1 )
    bs.write_enzo_set(stuff[6],output_directory,field='p19b_%s_bz.512'%suite,extend=2 )

if 0:
    #DO ONE STEP.
    #ADD TRACERS
    import p78_assemble.add_tracer_tool as AT
    reload(AT)
    particle_base="/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Particles512"
    setname = "%s/DD0000/data0000"%output_directory
    frame=0
    ds = yt.load(setname)
    AT.add_particles(ds,setname,particle_base, outnumber = frame)



if 0:
    #check that the hdf5 files are correct

    output_directory='/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Run128'
    fieldname="p19b_b02_bx.128"
    fname = "%s/%s"%(output_directory,fieldname)
    fptr=h5py.File(fname,'r')
    b = fptr[fieldname][()]
    fptr.close()
    print('B',b.shape)

    #dsnew = yt.load('/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Run128/DD0000/data0000')
    #cgnew = dsnew.covering_grid(0,[0.0]*3,[128]*3)
    fname2 = "%s/DD%04d/data%04d.cpu0000"%(output_directory,0,0)
    fptr = h5py.File(fname2,'r')
    c = fptr['Grid00000001']['BxF'][()]
    fptr.close()

    print('C',c.shape)


if 0:
    #compare new vs old cubes.
    dsnew = yt.load('/anvil/scratch/x-ux454321/p78c_high_res/IC_assembler/Run128/DD0000/data0000')
    cgnew = dsnew.covering_grid(0,[0.0]*3,[128]*3)

    nx = 2
    ny = len(fields)
    s = 3
    plt.close('all')
    fig,axes=plt.subplots(nx, ny, figsize=(s*ny, s*nx))
    plotax=2
    for nf,field in enumerate(fields):
        d1 = cg[field]
        d2 = cgnew[field]
        p1 = d1.sum(axis=plotax)
        #p2 = d2.sum(axis=2).transpose()
        p2 = d2.sum(axis=plotax)
        ax = axes[0][nf]
        p=ax.imshow(p1)
        fig.colorbar(p,ax=ax)
        ax = axes[1][nf]
        p=ax.imshow(p2)
        fig.colorbar(p,ax=ax)
    fig.savefig('%s/compare'%plot_dir)

