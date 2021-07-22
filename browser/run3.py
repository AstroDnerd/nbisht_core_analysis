

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/Users/dcollins/Dropbox/RESEARCH5/Paper19d_u300"
sim_list=['u301','u302','u303']
#sim_list=['u301']
for this_simname in sim_list:

    print('run',this_simname)

    #p_tc_raw = product.product('Tc_raw',fname='%s/data_small/%s_ct.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)
    #p_tc_raw_log = product.product('Log Tc_raw',fname='%s/data_small/%s_mean_collapse_log.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400,number_format="%0.2f")

    #p_distance = product.product('Log Distance to Next',fname='%s/data_small/%s_neighbors_distance_0_log.h5'%(basedir,this_simname),field='distance_0',style='value',width=400,number_format='%0.2f')
    #p_mean_density = product.product('Log Core Density',fname='%s/data_small/%s_neighbors_mean_density_log.h5'%(basedir,this_simname),field='mean_density',style='value',width=400, number_format='%0.2f')
    #u05_neighbors_distance_0.h5
    #u05_neighbors_mean_density.h5

 

    p_peak_density = product.product('Peak Density',fname="datasets_small/%s_mountain_tops_take_9.h5"%this_simname,
                                     field='peak_density',style='value_target_file',width=400)

    r_mountain = r"%s/mountain_tops/%s/%s_peak_p(\d\d\d\d)_34_Projection_x_density.png"%(basedir, this_simname, this_simname)
    p_mountain = product.product("mountain top", regexp=r2, parameters=['core_id'],style='single',width=400)
    p_mountain.get_frames()

    r_density = r"%s/density_time/%s/%s_density_6_c(\d\d\d\d).png"%(basedir, this_simname, this_simname)
    p_density = product.product("density time", regexp=r_density, parameters=['core_id'],style='single',width=400)
    p_density.get_frames()

    r_core_proj_follow   = r"%s/Cores/%s/%s_core_zoom_annotate_c(\d\d\d\d)_n(\d\d\d\d)_Projection_x_density.png"%(basedir,this_simname,this_simname)
    p_core_proj_follow = product.product('core_proj_follow',regexp=r_core_proj_follow,parameters=['core_id','frame'],style='frames',width=400)
#p2.check_glob()
    p_core_proj_follow.get_frames()



    product_list=[p_peak_density,p_mountain,p_density, p_core_proj_follow]
    cl=make_page.make_page(product_list, core_list=None,htmlname='browser/output_%s.html'%(this_simname))


