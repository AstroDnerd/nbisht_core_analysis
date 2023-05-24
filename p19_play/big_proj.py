from starter2 import *
reload(loop_apps)

#TODO
#Why do the figure labels come out in not latex?
#Needs to not be in g/cm^2
if 0:
    ds = yt.load(tll.directory+"/DD0125/data0125")
    peaks = h5py.File('datasets_small/u05_0125_peaklist.h5','r')
    peak_list = peaks['peaks'][:]
    peaks.close()
    proj = ds.proj('density',0)
if 1:
    pw=proj.to_pw()
    pw.set_cmap('density','Greys')
if 0:
    peak, center = ds.h.find_max('density')
if 1:
    pw.set_center([center[1],center[2]])
if 0:
    pd = looper.count_particles()
    for n,p in enumerate(peak_list):
        if pd[n] > 1:
            point_x = p[1]
            point_y = p[2]
            if point_x < pw.bounds[0]: point_x += 1
            if point_x > pw.bounds[1]: point_x -= 1
            if point_y < pw.bounds[2]: point_y += 1
            if point_y > pw.bounds[3]: point_y -= 1
            pw.annotate_text(p,"%d"%n,coord_system='data')
if 0:
    #this doesn't quite work.
    #Still unhappy with plot labels.
    print(pw.save('plots_to_sort/big_proj_n%04d.png'%125))
    pw.plots['density'].axes.set_xlabel(r'$\rm{yy\ \ (\rm{cm})}$')
    pw.plots['density'].axes.set_ylabel(r'$\rm{z}$')
    pw.plots['density'].cb.set_label(r'$\Sigma$')
    pw.plots['density'].figure.savefig('plots_to_sort/test2.png')
if 1:
    print(pw.save('plots_to_sort/big_proj_n%04d.pdf'%125))
