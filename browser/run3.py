

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/data/cb1/Public/P19_FlatBrowser"
sim_list=['u601','u602','u603']
#sim_list=['u301']

for isim,this_simname in enumerate(sim_list):
    frame = [125,118, 107][isim] 
    nsim=isim+1

    RE = []
    NAMES = []
    PROD = []

    print('run',this_simname)

    p_peak_density = product.product('Peak Density',fname="browser_data/%s_mountain_top.h5"%this_simname,
                                     field='peak_density',style='value_target_file',width=400)
    p_nparticles = product.product('Nparticles',fname="browser_data/u50%d_nparticles.h5"%nsim,
                                     field='n_particles',style='value')
    p_neighbor = product.product('Neighborhood',fname="browser_data/neighborhood_by_core_u50%d.h5"%nsim,
                                     field='neighborhood',style='value',number_format='%d')
    p_mode = product.product('mode',fname="browser_data/core_formation_mode_old_u50%d.h5"%nsim,
                                     field='modes',style='value',number_format='%s')
    p_mode_new = product.product('mode',fname="browser_data/core_formation_mode_new_u60%d.h5"%nsim,
                                     field='modes',style='value',number_format='%s')
    PROD.append(p_peak_density)
    PROD.append(p_nparticles)
    PROD.append(p_neighbor)
    PROD.append(p_mode)
    PROD.append(p_mode_new)

    NAMES.append("mountain top")
    RE.append(r"/*mountain_tops/%s/mountain_top_%s_c(\d\d\d\d)_Projection_x_density.png"%( this_simname, this_simname))

    #context, zoomed to 8 root grid zones
    NAMES.append("friends zoom8")
    RE.append(r"/*friends_zoom8/u50%d/friends_u50%d_c(\d{4})_n%04d_Projection_x_density.png"%(nsim,nsim, frame))
    #context, zoomed to 1 root grid zones
    NAMES.append("friends zoom1")
    RE.append(r"/*friends_zoom1/u50%d/friends_u50%d_c(\d{4})_n%04d_Projection_x_density.png"%(nsim,nsim, frame))


    NAMES.append("Buddy centroids")
    RE.append(r"/*buddy_centroids/u50%d/u50%d_buddies__centroid_xyz_c(\d{4}).png"%(nsim,nsim))

    #hair
    NAMES.append("Pathlines")
    RE.append(r"/*blowing_hair_u500/u50%d/u50%d_hair_x_c(\d\d\d\d).png"%(nsim  ,nsim))

    #3 projections
    NAMES.append('3 proj')
    RE.append(r"/*3proj/u50%d/u50%d_c(\d{4}).mp4"%(nsim,nsim))


    #anatomy.
    NAMES.append("Anatomy")
    RE.append(r"/*Anatomy3/u50%d/u50%d_rho_vel_hist_t_c(\d{4}).png"%(nsim,nsim))
    #NAMES.append("Anatomy")
    #RE.append(r"/*Anatomy/u50%d/u50%d_rho_vel_hist_t_c(\d{4}).png"%(nsim,nsim))

    #with density plots
    #NAMES.append("Anatomy 4")
    #RE.append(r"/*Anatomy_density/u50%d/u50%d_rho_vel_hist_t_c(\d{4}).png"%(nsim,nsim))
    NAMES.append("Anatomy play")
    RE.append(r"/*Anatomy_play/u50%d/u50%d_rho_vel_hist_t_c(\d{4}).png"%(nsim,nsim))
    print(RE[-1])

    NAMES.append("EG to EK")
    RE.append(r"/*Anatomy_play/u50%d/eg_to_ek_vs_radius_u50%d_c(\d{4}).png"%(nsim,nsim))

    ##Why Bound
    ##Not the right version
    #NAMES.append("GE/KE Radius Time")
    #RE.append("/*why2/u50%d/why2_u50%d_c(\d{4}).png"%(nsim,nsim))

    #Why Bound
    #good version, seismic colormap
    #NAMES.append("GE/KE Radius Time")
    #RE.append("/*why2_centroid/u50%d/why2_u50%d_c(\d{4}).png"%(nsim,nsim))
    NAMES.append("GE/KE Radius Time")
    RE.append("/*Anatomy_play/u50%d/why2_u50%d_c(\d{4}).png"%(nsim,nsim))

    NAMES.append("Avg interior mass")
    RE.append("/*MassInterior/u50%d/mass_interior_vs_radius_u50%d_c(\d{4}).png"%(nsim,nsim))

    NAMES.append("Radial/Tangential Velocity")
    RE.append("/*velocity_spheres/u50%d/velocity_color_u50%d_c(\d{4}).png"%(nsim,nsim))

    NAMES.append("Radial/Tangential Velocity")
    RE.append("/*velocity_spheres/u50%d/velocity_vs_radius_u50%d_c(\d{4}).png"%(nsim,nsim))

    RE.append(r"/*otherones/u40%d/otherones_u40%d_c(\d\d\d\d)__0000_%04d.png"%(nsim,nsim, frame))
    NAMES.append('Otherones')

    #NAMES.append("Density Velocity")
    #RE.append(r"/*density_velocity/u50%d/u50%d_rho_vnorm_t_c(\d{4}).png"%(nsim,nsim))

    #Clearly shows that mass conservation is only as good as the resolution is steady.
    #NAMES.append("Density mass")
    #RE.append(r"/*density_mass_dof/u50%d/u50%d_mass_volume_t_c(\d{4}).png"%(nsim,nsim))
    #Shows that mass is not conserved.
    #NAMES.append("Density Mass")
    #RE.append(r"/*density_mass/u50%d/u50%d_rho_mass_t_c(\d{4}).png"%(nsim,nsim))

    #proto anatomy
    #NAMES.append("Density vel hist")
    #RE.append(r"/*density_vel_hist/u50%d/u50%d_rho_vel_hist_t_c(\d{4}).png"%(nsim,nsim))

    #good velocity plots.
    NAMES.append("V extras")
    RE.append(r"/*velocity_mean_with_pearson/u50%d/u50%d_v_t_c(\d\d\d\d).png"%(nsim,nsim))
    #just velocity
    #NAMES.append("V-Vmean")
    #RE.append(r"/*velocity_mean/u50%d/u50%d_v_t_c(\d\d\d\d).png"%(nsim,nsim))



    if 1:
        #good movies
        #bulky
        NAMES.append("GE/KE")
        RE.append(r"/*GE_KE_zoom4/u50%d/u50%d_c(\d{4}).mp4"%(nsim,nsim))


    #nice mountain top image with the acceleration field
    NAMES.append("Gravity Density Ring")
    RE.append(r"/*proj_grav/u60%d/mountain_top_u60%d_c(\d\d\d\d)_acceleration_Projection_x_density.png"%(nsim,nsim))

    #energies.
    NAMES.append("GE/KE 2/128")
    RE.append(r"/*GEKEmass_larger/u60%d/eng_x_u60%d_c(\d\d\d\d)_n%04d_ge_ke_t2.png"%(nsim,nsim,frame))

    #velocity streamlines, v=v-v_mean_particles
    NAMES.append("Velocity Streamlines M1")
    RE.append(r"/*swirl_1/u60%d/proj_x_u60%d_c(\d\d\d\d)_n%04d_density_velocity.png"%(nsim,nsim, frame))


    NAMES.append('(Grad phi)^2')
    RE.append(r"/*BindingEnergy/u50%d/u50%d_c(\d\d\d\d)_potfit.png"%(nsim,nsim))

    if 0:
        NAMES.append('B and rho') #u503_b_and_rho_t_c0248
        RE.append(r"/*B_and_rho/u50%d/u50%d_b_and_rho_t_c(\d\d\d\d).png"%(nsim,nsim))



    #NAMES.append("mountain rings")
    #RE.append(r"/*mountain_rings/u60%d/mountain_top_u60%d_c(\d\d\d\d)_Projection_x_density.png"%(nsim,nsim))
    #NAMES.append('GE mass')
    #RE.append(r"/*ge_cuml/u50%d/u50%d_cuml_c(\d\d\d\d).png"%(nsim,nsim))
    #NAMES.append('GE cuml')
    #RE.append(r"/*GEcuml/u60%d/u60%d_cuml_c(\d\d\d\d).png"%(nsim,nsim))

    NAMES.append("GE/KE")
    RE.append(r"/*GEKEmass/u60%d/eng_x_u60%d_c(\d\d\d\d)_n%04d_ge_ke_t1.png"%(nsim,nsim,frame))
    #NAMES.append("Swirl M2")
    #RE.append(r"/*swirl_2/u60%d/proj_x_u60%d_c(\d\d\d\d)_n%04d_density_velocity_whole_mean.png"%(nsim,nsim,frame))
    #NAMES.append("Swirl M3")
    #RE.append(r"/*swirl_3/u60%d/proj_x_u60%d_c(\d\d\d\d)_n%04d_density_velocity_vcentral.png"%(nsim,nsim,frame))


    #NAMES.append("Density paths")
    #RE.append(r"/*density_time_hair/u50%d/u50%d_rho_t_c(\d\d\d\d).png"%(nsim,nsim))

    #NAMES.append("Density heatmap")
    #RE.append( r"/*density_time_heat/u50%d/u50%d_rho_t_heat_c(\d\d\d\d).png"%(nsim,nsim))
    #NAMES.append("Velocity path")
    #RE.append(r"/*velocity_hair/u50%d/u50%d_v_t_c(\d\d\d\d).png"%(nsim,nsim))

    #NAMES.append("Velocity heatmap")
    #RE.append(r"/*velocity_heat/u50%d/u50%d_v_t_heat_c(\d\d\d\d).png"%(nsim,nsim))

    #NAMES.append("VR VT")
    #RE.append(r'/*vr_vt/u50%d/u50%d_vr_vt_c(\d\d\d\d).png'%(nsim,nsim))
    link_dir=''
    for N in range(len(RE)):
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


