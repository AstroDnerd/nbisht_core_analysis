from starter2 import *
from collections import defaultdict
import scipy
import colors

if 'set_looper' not in dir():
    savefile='u301_long_pos_only.h5'
    set_looper=looper.core_looper(directory= dl.sims['u301'],savefile_only_trackage=savefile)
    thtr = set_looper.tr
    if 1:
        bad = np.where(thtr.track_dict['density'] <= 0)
        bad_pids = thtr.particle_ids[bad[0]]
        for bad_id in bad_pids:
            print("go ", bad_id)
            n_bad_densities= (thtr.p([bad_id],'density')  <= 0).sum()
            if n_bad_densities == 0:
                print("No bad densities.")
                raise
            particle_index = np.where( thtr.particle_ids == bad_id)
            for field in thtr.track_dict.keys():
                arr = thtr.track_dict[field]
                smaller = np.delete( arr, particle_index,axis=0)
                thtr.track_dict[field]=smaller
            thtr.particle_ids =  np.delete(thtr.particle_ids, particle_index)
            thtr.core_ids =  np.delete(thtr.core_ids, particle_index)
    stillbad = np.where(thtr.track_dict['density'] <= 0)
    print("STILL BAD", stillbad)
    thtr = set_looper.tr
    set_looper.out_prefix='core_13etal'
    thtr.sort_time()
    core_id=24

if 1:
    thtr = set_looper.tr
    core_id = 24

if 0:
    thtr = TLM.loops['u303'].tr 
    core_id = 186
if 0:
    thtr = TLM.loops['u301'].tr 
    core_id = 24


if 1:
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    reload(trackage)
    ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
    ms.particle_pos(core_id)
    rx = ms.particle_x - ms.meanx2
    ry = ms.particle_y - ms.meany2
    rz = ms.particle_z - ms.meanz2
    #rx = ms.this_x - ms.meanx2
    #ry = ms.this_y - ms.meany2
    #rz = ms.this_z - ms.meanz2
    r = np.sqrt( rx**2+ry**2 +rz**2)[:,1:] #because the first snap is repeated.

    norm = mpl.colors.Normalize()
    norm.autoscale( r[:,0])
    ax[1][0].hist( r[:,0], histtype='step')
    colormap = mpl.cm.ScalarMappable(norm=norm,cmap='jet')

    #mask = r[:,0] > 0.15
    #r = r[mask,:]




    t = thtr.times[1:]
    dt = (t[1:]-t[:-1])
    mt = 0.5*(t[1:]+t[:-1])
    drdt = ((r[:,1:]-r[:,:-1])/dt)
    dr0 = np.abs(drdt[:,:3].mean(axis=1))
    dr0.shape = (dr0.shape[0],1) #reshape to divide
    part = np.log(np.abs(drdt/dr0))
    exponent = np.log(np.abs(drdt/dr0)).mean(axis=0)
    mask= np.where(drdt[:,0] >-1e7)

    p4 = np.log(np.abs(drdt/dr0))
    for ip in mask[0]:
        c = colormap.to_rgba( r[ip,0])
        ax[0][0].plot( t, r[ip,:],c=c, linewidth=0.1)
        ax[0][1].plot( mt, drdt[ip,:], c=c)
        #ax[1][1].plot( mt, part[ip,:])
        ax[1][1].plot( mt, p4[ip,:],c=c)
        #ax.plot( ms.raw_y[ip,:])
    mp4=p4.mean(axis=0)
    ax[1][1].plot(mt, mp4,c='k')
    ax[0][1].plot(mt,drdt.mean(axis=0), c='k')
    #ax[0][1].plot(mt,(drdt/dr0).mean(axis=0))

    axbonk(ax[0][1],ylabel='<v>')

    #for nt, time in enumerate(dt):
    #    if nt%10 > 0:
    #        continue
    #    color=rainbow_map(len(t))(nt)
    #    ax[1][0].hist( drdt[:,nt], histtype='step',color=color)
        
    fig.savefig('plots_to_sort/shift_test.png')
    


