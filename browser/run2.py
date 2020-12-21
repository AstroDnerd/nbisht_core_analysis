

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/Users/dcollins/RESEARCH3/Paper19c_u10u11/0000_more_good_plots/"
basedir = "/Users/dcollins/Dropbox/RESEARCH5/Paper19c_u10u11/0000_more_good_plots/"
for this_simname in ['u05','u10','u11']:
    glob1 = "%s/density_time/%s/%s_density_6_c????.png"%(basedir,this_simname,this_simname)
    reg1 = re.compile(r'%s/density_time/%s/%s_density_6_c(\d\d\d\d).png'%(basedir,this_simname,this_simname))
    p1 = product.product('rho_t', regexp=reg1, myglob=glob1, parameters=['core_id'],style='single',width=400)
    p1.get_frames()

    p13 = product.product('Tc_raw',fname='%s/data_small/%s_ct.h5'%(basedir,this_simname),field='collapse_times',style='value',width=400)

    g2 = "%s/alpha_time/%s/%s_density_radius_c????.png"%(basedir, this_simname, this_simname)
    r2 = re.compile(r"%s/alpha_time/%s/%s_density_radius_c(\d\d\d\d).png"%(basedir, this_simname, this_simname))
    p2 = product.product("alpha-time", regexp=r2, myglob=g2, parameters=['core_id'],style='single',width=400)
    p2.get_frames()

    g3  = r"%s/proj_follow/%s/%s_c????_n????_centered_Projection_x_density.png"%(basedir,this_simname,this_simname)
    r3   = re.compile(r"%s/proj_follow/%s/%s_c(\d\d\d\d)_n(\d\d\d\d)_centered_Projection_x_density.png"%(basedir,this_simname,this_simname))
    p3 = product.product('core_proj_follow',regexp=r3,myglob=g3,parameters=['core_id','frame'],style='frames',width=400)
#p2.check_glob()
    p3.get_frames()


    product_list=[p13,p1,p2,p3]
    cl=make_page.make_page(product_list, core_list=None,htmlname='browser/output_%s.html'%(this_simname))


