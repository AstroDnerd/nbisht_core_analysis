

from dtools.starter1 import *
import yt



def plot_two(old,new,name):
    ds_old = yt.load(old)
    ds_new = yt.load(new)
    for field in ['density','magnetic_field_x','velocity_z']:
        yt.ProjectionPlot(ds_old,0,field).save('%s/%s_old'%(plot_dir,name))
        yt.ProjectionPlot(ds_new,0,field).save('%s/%s_new'%(plot_dir,name))

old_set_b02 = "/anvil/scratch/x-ux454321/p19b_high_res/take_1_restart_512/OLD_B02/RS0000/restart0000"
blank_set_b02 = "/anvil/scratch/x-ux454321/p19b_high_res/take_1_restart_512/New_B02/DD0000/data0000"
plot_two(old_set_b02,blank_set_b02,'b02')
