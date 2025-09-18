from starter2 import *
import p78_assemble.add_tracer_tool as AT




#ds1 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/u05/DD0000_fresh/data0000')
#ds2 = yt.load('/scratch1/dcollins/Paper48_DrivingFiller/Bred8mhd/DD0001/data0001')
#boundary(ds1,ds2)
#copy_all_fields(ds1,ds2,fields = ['BxF', 'ByF', 'BzF','Density', 'x-velocity','y-velocity','z-velocity', 'Bx','By','Bz'])


if 1:
    """Usage."""
    #dirname = '/scratch1/dcollins/Paper19/SphereTest/s05b_uni_repeat'
    #outdir = '/scratch1/dcollins/Paper19/SphereTest/s05b_uni_repeat/tracers'
    #dirname = '/scratch1/dcollins/Paper06_Multiphase/fd03_add_tracer_amr'
    #outdir  = '/scratch1/dcollins/Paper06_Multiphase/fd03_add_tracer_amr/DD0001t'
    #dirname = '/scratch/00369/tg456484/Paper42_runs/eq44_Actually9_d0.5_p59_L4_J8/transfer/'
    #outdir  = '/scratch/00369/tg456484/Paper42_runs/eq44_Actually9_d0.5_p59_L4_J8/transfer/D1'
    #frame = 0;setname = '%s/DD%04d/data%04d'%(dirname,frame,frame)

    #p78c take 1
    #dirname = '/anvil/scratch/x-ux454321/p78c_high_res/512/B02/512/'
    #outdir  = '/anvil/scratch/x-ux454321/p78c_high_res/512/B02/512/WithTracers'
    #dirname = '/anvil/scratch/x-ux454321/p78c_high_res/512/B2'
    #outdir  = '/anvil/scratch/x-ux454321/p78c_high_res/512/B2'
    #dirname = '/anvil/scratch/x-ux454321/p78c_high_res/512/B20'
    #outdir  = '/anvil/scratch/x-ux454321/p78c_high_res/512/B20'
    #frame=0
    #setname = "%s/RS0000_source/restart0000"%dirname
    #p78c take 2
    frame=0
    dirname = "/anvil/scratch/x-ux454321/p78c_high_res/512_take2"
    outdir  = "/anvil/scratch/x-ux454321/p78c_high_res/512_take2/DD0000_blank_particles"
    setname = "%s/DD0000_512_blank/data0000"%dirname
    ds = yt.load(setname)
    AT.add_particles(ds,setname,outdir, outnumber = frame)
