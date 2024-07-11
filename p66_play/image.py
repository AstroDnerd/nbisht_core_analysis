from starter2 import *
import track_loader as TL
#import ucore
#reload(ucore)


import dtools.davetools as dt


sim_list=['u501','u502','u503']
import track_loader as TL
TL.load_tracks(sim_list)
import monster
monster.load(sim_list)


def image(mon,core_list,frames):

    if type(frames) == str:
        if frames == 'movie':
            frame_slice=slice(None)
        elif frames == 'short':
            frame_slice = np.zeros_like(mon.frames,dtype='bool')
            frame_slice[::10]=True
            frame_slice[-1]=True
        elif frames == 'realshort':
            frame_slice=slice(-1,None,10)
            frame_slice = np.zeros_like(mon.frames,dtype='bool')
            frame_slice[0]=True;frame_slice[-1]=True
        frame_list=mon.frames[frame_slice]


    for frame in frame_list:
        for core_id in core_list:
            print('image',mon.name,core_id,frame)
            sph1 = mon.get_sphere(core_id,frame,'r1')
            sphmax = mon.get_sphere(core_id,frame,'rmax')
            sphinf = mon.get_sphere(core_id,frame,'rinf')
            sphsmrt = mon.get_sphere(core_id,frame,'rsmart')
            rad = [sph1.radius, sphmax.radius,sphinf.radius]
            sss = [sph1, sphmax,sphinf]
            which = np.argmax(rad)
            sph = sss[which]
            ds = mon.get_ds(frame)
            proj = ds.proj('density',0,data_source=sph, center=sph.center)
            pw = proj.to_pw()
            pw.annotate_sphere(sph1.center,sph1.radius,circle_args={'color':'r'})
            pw.annotate_sphere(sphmax.center,sphmax.radius,circle_args={'color':'g'})
            pw.annotate_sphere(sphinf.center,sphinf.radius,circle_args={'color':'orange'})
            pw.annotate_sphere(sphsmrt.center,sphsmrt.radius,circle_args={'color':'cyan'})

            pw.set_width(2*sph.radius)
            pw.save('%s/proj_%s_n%04d_c%04d'%(plot_dir,mon.name,frame,core_id))
for sim in sim_list[:1]:
    mon = monster.closet[sim]
    core_list =  mon.this_looper.core_by_mode['A'][3:4]
    pw=image(mon,core_list,frames='short')
