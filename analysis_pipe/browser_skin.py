from starter2 import *

import track_loader as TL

import product 
reload(product)
import make_page 
reload(make_page)

basedir = "./"
plotdir="./plots_to_sort"

def make_browser(trackname):
    track=track_info.tracks[trackname]
    mountain_top=track.mountain_top

    RE = []
    NAMES = []
    PROD = []

    p_peak_density = product.product('Peak Density',fname=mountain_top,
                                     field='peak_density',style='value_target_file',width=400)
    p_mode = product.product('mode',fname=track.mode_fname,
                                     field='modes',style='value',number_format='%s')
    PROD.append(p_peak_density)
    PROD.append(p_mode)

    NAMES.append("mountain top")
    RE.append(r"plots_to_sort/mountain_top_%s_c(\d\d\d\d)_Projection_x_density.png"%( trackname))

    link_dir='./'
    for N in range(len(RE)):
        print('wut')
        newprod = product.product(NAMES[N], regexp=RE[N], parameters=['core_id'],style='single',width=400,data_dir=basedir,link_dir=link_dir)
        newprod.get_frames()
        PROD.append(newprod)


    product_list=PROD

    cl=make_page.make_page(product_list, core_list=None,htmlname='%s/output_%s.html'%(plotdir,trackname))




