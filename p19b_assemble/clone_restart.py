

from dtools.starter1 import *
import yt

def clone(old_ds_fname,blank_ds_fname):
    ds_old = yt.load(old_ds_fname)
    ds_blank = yt.load(blank_ds_fname)

    for ng, grid in enumerate(ds_old.index.grids):
        left_old = grid.LeftEdge
        right_old = grid.RightEdge
        grid_blank = ds_blank.index.grids[ng]
        left_blank =  grid_blank.LeftEdge
        right_blank = grid_blank.RightEdge
        L1 = ((left_old-left_blank)**2).sum()
        R1 = ((right_old-right_blank)**2).sum()
        if (np.abs(L1)+np.abs(R1)) >0:
            print("OH NO")
            pdb.set_trace()
        fname_old = grid.filename
        fname_blank = grid_blank.filename
        h5ptr_old = h5py.File(fname_old,'r')
        h5ptr_blank=h5py.File(fname_blank,'r+')
        try:
            grid_old = h5ptr_old['Grid%08d'%grid.id]
            grid_blank = h5ptr_blank['Grid%08d'%grid.id]
            for field in grid_old:
                if field in grid_blank:
                    del grid_blank[field]
                grid_blank[field]=grid_old[field][()]
            print( 'Fixed grid',grid.id)
        except:
            raise
        finally:
            h5ptr_old.close()
            h5ptr_blank.close()



old_set = "/anvil/scratch/x-ux454321/p19b_high_res/assembler/OLD_B02/RS0000/restart0000"
blank_set = "/anvil/scratch/x-ux454321/p19b_high_res/assembler/New_B02/DD0000/data0000"
particle_set = "/anvil/scratch/x-ux454321/p19b_high_res/assembler/Particles_B02"

if 0:
    clone(old_set,blank_set)
if 1:
    #DO ONE STEP.
    #ADD TRACERS
    import add_tracer_tool as AT
    reload(AT)
    setname = blank_set
    frame=0
    ds = yt.load(setname)
    AT.add_particles(ds,setname,particle_set, outnumber = frame)
