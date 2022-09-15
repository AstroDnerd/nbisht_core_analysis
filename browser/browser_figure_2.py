

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/data/cb1/Public/P19_FlatBrowser"
sim_list=['u601','u602','u603']
#sim_list=['u301']

all_RE={}
all_NAMES={}
all_PROD={}
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
    p_mode = product.product('mode',fname="browser_data/core_formation_mode_u50%d.h5"%nsim,
                                     field='modes',style='string', number_format='%s')
    PROD.append(p_peak_density)
    PROD.append(p_nparticles)
    #PROD.append(p_neighbor)
    PROD.append(p_mode)

    NAMES.append("mountain top")
    RE.append(r"/*mountain_tops/%s/mountain_top_%s_c(\d\d\d\d)_Projection_x_density.png"%( this_simname, this_simname))

    #context, zoomed to 8 root grid zones
    NAMES.append("friends zoom 8x")
    RE.append(r"/*friends_zoom8/u50%d/friends_u50%d_c(\d{4})_n%04d_Projection_x_density.png"%(nsim,nsim, frame))

    #hair
    NAMES.append("Pathlines")
    RE.append(r"/*blowing_hair_u500/u50%d/u50%d_hair_x_c(\d\d\d\d).png"%(nsim  ,nsim))

    #3 projections
    NAMES.append('Collapse Movies')
    RE.append(r"/*3proj/u50%d/u50%d_c(\d{4}).mp4"%(nsim,nsim))


    RE.append(r"/*otherones/u40%d/otherones_u40%d_c(\d\d\d\d)__0000_%04d.png"%(nsim,nsim, frame))
    NAMES.append('Other ones')

    link_dir=''
    for N in range(len(RE)):
        newprod = product.product(NAMES[N], regexp=RE[N], parameters=['core_id'],style='single',width=400,data_dir=basedir,link_dir=link_dir)
        newprod.get_frames()
        PROD.append(newprod)

    all_RE[this_simname] = RE
    all_NAMES[this_simname] = NAMES
    all_PROD[this_simname] = PROD

core_ids={}
core_ids['u601']=[323, 8, 27, 32, 37, 44, 84, 275]
core_ids['u602']=[32, 378]
core_ids['u603']=[233, 233, 235]
make_page.make_page_multi(all_PROD, core_ids, "%s/PaperI_figure2.html"%basedir)
