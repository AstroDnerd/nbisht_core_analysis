
from go import *
import product
reload(product)
import make_page
reload(make_page)
basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
glob1 = "%s/density_radius/density_radius_c????.png"%basedir
reg1  = re.compile(r"%s/density_radius/density_radius_c(\d\d\d\d).png"%basedir)
p1 = product.product('rho_r',regexp=reg1,myglob=glob1,parameters=['core_id'],style='single',width=400)
p1.get_frames()

glob1 = "%s/core_proj_follow/density_radius_c????.png"%basedir
glob2  = r"%s/core_proj_follow/follow_c????_n????_centered_Projection_x_density.png"%(basedir)
reg2   = re.compile(r"%s/core_proj_follow/follow_c(\d\d\d\d)_n(\d\d\d\d)_centered_Projection_x_density.png"%(basedir))
p2 = product.product('core_proj_follow',regexp=reg2,myglob=glob2,parameters=['core_id','frame'],style='frames',width=400)
#p2.check_glob()
p2.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/density_time_all_cores/with_proper_tff/all/test_density_5_c????.png"%basedir
regn  = re.compile(r"%s/density_time_all_cores/with_proper_tff/all/test_density_5_c(\d\d\d\d).png"%basedir)
p3 = product.product('rho_t',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p3.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/density_radius_subfits/density_radius_c????.png"%basedir
regn = re.compile("%s/density_radius_subfits/density_radius_c(\d\d\d\d).png"%basedir)
p4 = product.product('alpha_t',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p4.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/velocity/velocity_radius_components_rel/vi_r_c????.png"%basedir
regn = re.compile("%s/velocity/velocity_radius_components_rel/vi_r_c(\d\d\d\d).png"%basedir)
p5 = product.product('v_r_rel',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p5.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/velocity/velocity_time_components_rel/vi_t_rel_c????.png"%basedir
regn  = re.compile("%s/velocity/velocity_time_components_rel/vi_t_rel_c(\d\d\d\d).png"%basedir)
p6 = product.product('v_t_rel',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p6.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s//velocity/velocity_time_components_avg/vi_t_rel_c????.png"%basedir
regn  = re.compile("%s//velocity/velocity_time_components_avg/vi_t_rel_c(\d\d\d\d).png"%basedir)
p7 = product.product('v_t_rel_sigma',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p7.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s//velocity/velocity_radius_symlog/vi_r_symlog_c????.png"%basedir
regn = re.compile(r"%s//velocity/velocity_radius_symlog/vi_r_symlog_c(\d\d\d\d).png"%basedir)
p8 = product.product('v_r_sym',regexp=regn,myglob=globn,parameters=['core_id'],style='single',width=400)
p8.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/zoom1/x/c????/test_full_particles_c????_n????_Projection_x_density.png"%basedir
regn = re.compile(r"%s/zoom1/x/c(\d\d\d\d)/test_full_particles_c(\d\d\d\d)_n(\d\d\d\d)_Projection_x_density.png"%basedir)
p9 = product.product('v_r_sym',regexp=regn,myglob=globn,parameters=['core_id','core_id','frame'],style='frames',width=400)
p9.get_frames()

basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/zoom2/x/c????/u05_centroids_zoom2_c????_n????_Projection_x_density.png"%basedir
regn = re.compile(r"%s/zoom2/x/c(\d\d\d\d)/u05_centroids_zoom2_c(\d\d\d\d)_n(\d\d\d\d)_Projection_x_density.png"%basedir)
p10 = product.product('v_r_sym',regexp=regn,myglob=globn,parameters=['core_id','core_id','frame'],style='frames',width=400)
p10.get_frames()

#zoom 2, just the last frame                  
basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/zoom2/x/c????/u05_centroids_zoom2_c????_n0125_Projection_x_density.png"%basedir
regn = re.compile(r"%s/zoom2/x/c(\d\d\d\d)/u05_centroids_zoom2_c(\d\d\d\d)_n0125_Projection_x_density.png"%basedir)
p10 = product.product('v_r_sym',regexp=regn,myglob=globn,parameters=['core_id','core_id'],style='single',width=400)
p10.get_frames()

#zoom 4, just the last frame                  
basedir='/Users/dcollins/RESEARCH3_USE_R4/Paper19_47_overlap/0000_main_plots/'
globn = "%s/zoom4/x/c????/u05_centroids_zoom4_c????_n0125_Projection_x_density.png"%basedir
regn = re.compile(r"%s/zoom4/x/c(\d\d\d\d)/u05_centroids_zoom4_c(\d\d\d\d)_n0125_Projection_x_density.png"%basedir)
p11 = product.product('v_r_sym',regexp=regn,myglob=globn,parameters=['core_id','core_id'],style='single',width=400)
p11.get_frames()

large1 = [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 47, 49, 50, 51, 53, 56, 63, 220]
suburbs1 = [41, 42, 43, 44, 45, 46, 54, 55, 57, 58, 64, 108, 110, 111, 113]
core_list = large1 + suburbs1   
core_list=None

cl=make_page.make_page([p3,p4,p7,p8,p5,p2,p9,p11], core_list=core_list)
