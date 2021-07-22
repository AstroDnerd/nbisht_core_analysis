

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/Users/dcollins/RESEARCH3/Paper19c_u10u11/0000_more_good_plots/"
basedir = "/Users/dcollins/Dropbox/RESEARCH5/Paper19c_u10u11/0000_more_good_plots/"
#basedir = "/Users/dcollins/RESEARCH4/Paper19c_u10u11/0000_more_good_plots/"
for this_simname in ['u05','u10','u11']:

    glob1 = "%s/density_time/%s/%s_density_6_c????.png"%(basedir,this_simname,this_simname)
    reg1 = re.compile(r'%s/density_time/%s/%s_density_6_c(\d\d\d\d).png'%(basedir,this_simname,this_simname))
    p_rho_t = product.product('rho_t', regexp=reg1, myglob=glob1, parameters=['core_id'],style='single',width=400)
    p_rho_t.get_frames()

    p_tc_raw = product.product('Tc_raw',fname='%s/data_small/%s_ct.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)

    p_tc_global = product.product('Tc_global',fname='%s/data_small/%s_ct_glob.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)
    p_tc_raw = product.product('Tc_raw',fname='%s/data_small/%s_ct.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)
    #p_tc_raw_log = product.product('Log Tc_raw',fname='%s/data_small/%s_ct_log.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400,number_format="%0.2f")
    p_tc_raw_log = product.product('Log Tc_raw',fname='%s/data_small/%s_mean_collapse_log.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400,number_format="%0.2f")

    p_distance = product.product('Log Distance to Next',fname='%s/data_small/%s_neighbors_distance_0_log.h5'%(basedir,this_simname),field='distance_0',style='value',width=400,number_format='%0.2f')
    p_mean_density = product.product('Log Core Density',fname='%s/data_small/%s_neighbors_mean_density_log.h5'%(basedir,this_simname),field='mean_density',style='value',width=400, number_format='%0.2f')
    #u05_neighbors_distance_0.h5
    #u05_neighbors_mean_density.h5


    g2 = "%s/alpha_time/%s/%s_density_radius_c????.png"%(basedir, this_simname, this_simname)
    r2 = re.compile(r"%s/alpha_time/%s/%s_density_radius_c(\d\d\d\d).png"%(basedir, this_simname, this_simname))
    p_alpha_time = product.product("alpha-time", regexp=r2, myglob=g2, parameters=['core_id'],style='single',width=400)
    p_alpha_time.get_frames()

    g_neighbor = r"%s/proj_final/zoom4/%s/%s_core_zoom_annotate_c????_n????_Projection_x_density.png"%(basedir, this_simname, this_simname)
    r_neighbor =  re.compile(r"%s/proj_final/zoom4/%s/%s_core_zoom_annotate_c(\d\d\d\d)_n\d\d\d\d_Projection_x_density.png"%(basedir, this_simname, this_simname))
    p_neighbor = product.product("Proj Neighbors", regexp=r_neighbor,myglob=g_neighbor,parameters=['core_id'],style='single',width=400)
    p_neighbor.get_frames()

    g_hull = r"%s/convex_hull/FULL/%s/%s_hull_3d_t_c????_n0000.png"%(basedir,this_simname,this_simname)
    r_hull = re.compile(r"%s/convex_hull/FULL/%s/%s_hull_3d_t_c(\d\d\d\d)_n0000.png"%(basedir,this_simname,this_simname))
    p_hull = product.product("Hull", regexp=r_hull,myglob=g_hull,parameters=['core_id'],style='single',width=400)
    p_hull.get_frames()

    g_hull = r"%s/convex_hull/FULL/%s/%s_hull_3d_t_c????_n0000.png"%(basedir,this_simname,this_simname)
    r_hull = re.compile(r"%s/convex_hull/FULL/%s/%s_hull_3d_t_c(\d\d\d\d)_n0000.png"%(basedir,this_simname,this_simname))
    p_hull = product.product("Hull", regexp=r_hull,myglob=g_hull,parameters=['core_id'],style='single',width=400)
    p_hull.get_frames()

    p_core = product.product("ID",style='core_id')



    g_vorticity = "%s/vorticity/%s/%s_vorticity_c????.png"%(basedir, this_simname, this_simname)
    r_vorticity = re.compile(r"%s/vorticity/%s/%s_vorticity_c(\d\d\d\d).png"%(basedir, this_simname, this_simname))
    p_vorticity = product.product("vorticity", regexp=r_vorticity, myglob=g_vorticity, parameters=['core_id'],style='single',width=400)
    p_vorticity.get_frames()

    g3  = r"%s/proj_follow/%s/%s_c????_n????_centered_Projection_x_density.png"%(basedir,this_simname,this_simname)
    r3   = re.compile(r"%s/proj_follow/%s/%s_c(\d\d\d\d)_n(\d\d\d\d)_centered_Projection_x_density.png"%(basedir,this_simname,this_simname))
    p_core_proj_follow = product.product('core_proj_follow',regexp=r3,myglob=g3,parameters=['core_id','frame'],style='frames',width=400)
#p2.check_glob()
    p_core_proj_follow.get_frames()

    gx = "%s/velocity_time/%s/%s_vi_t_rel_c????.png"%(basedir, this_simname, this_simname)
    rx = re.compile(r"%s/velocity_time/%s/%s_vi_t_rel_c(\d\d\d\d).png"%(basedir, this_simname, this_simname))
    p_vel_time = product.product("velocity-time", regexp=rx, myglob=gx, parameters=['core_id'],style='single',width=400)
    p_vel_time.get_frames()


    #future self: merge conflict, not sure which was right.
    #product_list=[p_tc_global,p_rho_t,p_alpha_time,p_vel_time,p_core,p_core_proj_follow,p_neighbor,p_hull]
    #product_list=[p_distance, p_mean_density,p_tc_raw,p_tc_raw_log,p_vorticity,p1,p2,p3]
    cl=make_page.make_page(product_list, core_list=None,htmlname='browser/output_%s.html'%(this_simname))


