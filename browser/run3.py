

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/data/cb1/Public/P19_FlatBrowser"
sim_list=['u601','u602','u603']
#sim_list=['u301']

for isim,this_simname in enumerate(sim_list):
    nsim=isim+1

    RE = []
    NAMES = []
    PROD = []

    print('run',this_simname)

    #p_tc_raw = product.product('Tc_raw',fname='%s/data_small/%s_ct.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)
    #p_tc_raw_log = product.product('Log Tc_raw',fname='%s/data_small/%s_mean_collapse_log.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400,number_format="%0.2f")

    #p_distance = product.product('Log Distance to Next',fname='%s/data_small/%s_neighbors_distance_0_log.h5'%(basedir,this_simname),field='distance_0',style='value',width=400,number_format='%0.2f')
    #p_mean_density = product.product('Log Core Density',fname='%s/data_small/%s_neighbors_mean_density_log.h5'%(basedir,this_simname),field='mean_density',style='value',width=400, number_format='%0.2f')
    #u05_neighbors_distance_0.h5
    #u05_neighbors_mean_density.h5

 

    p_peak_density = product.product('Peak Density',fname="browser_data/%s_mountain_top.h5"%this_simname,
                                     field='peak_density',style='value_target_file',width=400)
    p_nparticles = product.product('Nparticles',fname="browser_data/u50%d_nparticles.h5"%nsim,
                                     field='n_particles',style='value')
    p_neighbor = product.product('Neighborhood',fname="browser_data/neighborhood_by_core_u50%d.h5"%nsim,
                                     field='neighborhood',style='value',number_format='%d')
    PROD.append(p_peak_density)
    PROD.append(p_nparticles)
    PROD.append(p_neighbor)

    NAMES.append("mountain top")
    RE.append(r"/*mountain_tops/%s/mountain_top_%s_c(\d\d\d\d)_Projection_x_density.png"%( this_simname, this_simname))

    NAMES.append("Pathlines")
    RE.append(r"/*blowing_hair_u500/u50%d/u50%d_hair_x_c(\d\d\d\d).png"%(nsim  ,nsim))

    NAMES.append("Density paths")
    RE.append(r"/*density_time_hair/u50%d/u50%d_rho_t_c(\d\d\d\d).png"%(nsim,nsim))

    NAMES.append("Density heatmap")
    RE.append( r"/*density_time_heat/u50%d/u50%d_rho_t_heat_c(\d\d\d\d).png"%(nsim,nsim))
    NAMES.append("Velocity path")
    RE.append(r"/*velocity_hair/u50%d/u50%d_v_t_c(\d\d\d\d).png"%(nsim,nsim))

    NAMES.append("Velocity heatmap")
    RE.append(r"/*velocity_heat/u50%d/u50%d_v_t_heat_c(\d\d\d\d).png"%(nsim,nsim))

    NAMES.append("VR VT")
    RE.append(r'/*vr_vt/u50%d/u50%d_vr_vt_c(\d\d\d\d).png'%(nsim,nsim))

    RE.append(r"/*otherones/u40%d/otherones_u40%d_c(\d\d\d\d)__0000_0125.png"%(nsim,nsim))
    NAMES.append('Otherones')

    NAMES.append('(Grad phi)^2')
    RE.append(r"/*BindingEnergy/u50%d/u50%d_c(\d\d\d\d)_potfit.png"%(nsim,nsim))

    link_dir=''
    for N in range(len(RE)):
        print(RE[N])
        newprod = product.product(NAMES[N], regexp=RE[N], parameters=['core_id'],style='single',width=400,data_dir=basedir,link_dir=link_dir)
        newprod.get_frames()
        PROD.append(newprod)





#   r_density = r"%s/density_time/%s/%s_density_6_c(\d\d\d\d).png"%(basedir, this_simname, this_simname)
#   p_density = product.product("density time", regexp=r_density, parameters=['core_id'],style='single',width=400)
#   p_density.get_frames()

#   r_core_proj_follow   = r"%s/Cores/%s/%s_core_zoom_annotate_c(\d\d\d\d)_n(\d\d\d\d)_Projection_x_density.png"%(basedir,this_simname,this_simname)
#   p_core_proj_follow = product.product('core_proj_follow',regexp=r_core_proj_follow,parameters=['core_id','frame'],style='frames',width=400)
#p2.check_glob()
#   p_core_proj_follow.get_frames()



    #product_list=[p_peak_density,p_mountain,p_density, p_core_proj_follow]
    #product_list=[p_nparticles,p_peak_density,p_mountain,p_hair, p_other,p_rho_hair, p_rho_heat]
    product_list=PROD

    cl=make_page.make_page(product_list, core_list=None,htmlname='%s/output_%s.html'%(basedir,this_simname))


