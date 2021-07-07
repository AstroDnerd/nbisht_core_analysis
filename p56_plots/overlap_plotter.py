from starter2 import *

reload(loop_apps)


fractions,cores=get_overlapping_cores(ht3,185)
catman = np.concatenate
cores = catman([cores,[185]])[::-1]
#ht3b = hull_tool(tl.looper3)
#ht3b.plot_2d(core_list = cores, accumulate=True)


frame_list = list(range(0,88,10))+[88]
loop_apps.proj_select_particles(tl.looper3,axis_list=[2],core_list=cores,field='density',frame_list=frame_list)
#frame_list = list(range(0,125,10))+[125]
#field_list=[0,125]
#loop_apps.proj_select_particles(tl.looper1,axis_list=[2],core_list=[10,11],field='density',frame_list=frame_list)
