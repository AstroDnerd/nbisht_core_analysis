

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/data/cb1/Public/P19_FlatBrowser"
sim_list=['u601','u602','u603']
#sim_list=['u603']
#sim_list=['u301']

numbers={'u601':1,'u602':2,'u603':3}
for isim,this_simname in enumerate(sim_list):
    frame = [125,118, 107][isim] 
    nsim=numbers[this_simname]

    RE = []
    NAMES = []
    PROD = []
    CAPTIONS=[]

    print('run',this_simname)

    p_peak_density = product.product('Peak Density',fname="browser_data/%s_mountain_top.h5"%this_simname,
                                     field='peak_density',style='value_target_file',width=400)
    p_nparticles = product.product('Nparticles',fname="browser_data/u50%d_nparticles.h5"%nsim,
                                     field='n_particles',style='value')
    p_neighbor = product.product('Neighborhood',fname="browser_data/neighborhood_by_core_u50%d.h5"%nsim,
                                     field='neighborhood',style='value',number_format='%d')
    p_mode = product.product('mode',fname="browser_data/core_formation_mode_u50%d.h5"%nsim,
                                     field='modes',style='string',number_format='%s')
    PROD.append(p_peak_density)
    CAPTIONS.append("The Peak Density")
    PROD.append(p_nparticles)
    CAPTIONS.append("Number of Particles.")
    PROD.append(p_neighbor)
    CAPTIONS.append("<i>Neighborhoods </i> are defined by cores with overlapping preimages.")
    PROD.append(p_mode)
    CAPTIONS.append("Cores form Alone, in Binaries (Bn) or in clusters (Cn)")

    NAMES.append("mountain top")
    RE.append(r"/*mountain_tops/%s/mountain_top_%s_c(\d\d\d\d)_Projection_x_density.png"%( this_simname, this_simname))
    CAPTIONS.append(  "The core and its tracer particles (red).  The yellow contour only serves to locate the particles. 32x zoom.")


    #context, zoomed to 8 root grid zones
    NAMES.append("friends zoom 8x")
    RE.append(r"/*friends_zoom8/u50%d/friends_u50%d_c(\d{4})_n%04d_Projection_x_density.png"%(nsim,nsim, frame))
    CAPTIONS.append("The core and its friends.  8x zoom.")

    #hair
    NAMES.append("Pathlines")
    RE.append(r"/*blowing_hair_u500/u50%d/u50%d_hair_x_c(\d\d\d\d).png"%(nsim  ,nsim))
    CAPTIONS.append("Path lines.  The particles begin at the black points and end on the red ones.  The <i>preimage</i> is defined by positions at the start.")

    #3 projections
    NAMES.append('Collapse Movies')
    RE.append(r"/*3proj/u50%d/u50%d_c(\d{4}).mp4"%(nsim,nsim))
    CAPTIONS.append("Movies of projections of density along 3 axes.  This core is green, other nearby cores are in red and orange.  Projections along X, Y, and Z are seen.")


    RE.append(r"/*otherones/u40%d/otherones_u40%d_c(\d\d\d\d)__0000_%04d.png"%(nsim,nsim, frame))
    NAMES.append('Other ones')
    CAPTIONS.append("The Other Ones.  Particles that do not end in dense cores, but begin near particles that do.  The top row is t=0, the bottom row is t=tfinal.  The left column is the spatial extent, and the right column density vs. radius.  The blue (green) points initially fill the <i>convex hull</i> defined by the initial particle position.  ")

    link_dir=''
    for N in range(len(RE)):
        newprod = product.product(NAMES[N], regexp=RE[N], parameters=['core_id'],style='single',width=400,data_dir=basedir,link_dir=link_dir)
        newprod.get_frames()
        PROD.append(newprod)


    core_list=[57]
    cl=make_page.make_page(PROD, core_list=core_list,htmlname='%s/captions.html'%basedir, captions=CAPTIONS)


