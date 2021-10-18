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
    name = 'u301_long'
    core_id = 24

if 0:
    thtr = TLM.loops['u303'].tr 
    name='u303'
    core_id = 186
if 0:
    thtr = TLM.loops['u301'].tr 
    name='u301'
    core_id = 24


if 1:
    fig,ax=plt.subplots(2,2,figsize=(8,8))
    reload(trackage)
    ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
    ms.particle_pos(core_id)
    rx = ms.particle_x# - ms.meanx2
    ry = ms.particle_y# - ms.meany2
    rz = ms.particle_z# - ms.meanz2
    #rx = ms.this_x - ms.meanx2
    #ry = ms.this_y - ms.meany2
    #rz = ms.this_z - ms.meanz2
    r = np.sqrt( rx**2+ry**2 +rz**2)[:,1:] #because the first snap is repeated.

    norm = mpl.colors.Normalize()
    norm.autoscale( r[:,0])
    ax[1][0].hist( r[:,0], histtype='step')
    colormap = mpl.cm.ScalarMappable(norm=norm,cmap='jet')

    t = thtr.times[1:]
    dt = (t[1:]-t[:-1])
    mt = 0.5*(t[1:]+t[:-1])
    drdt = ((r[:,1:]-r[:,:-1])/dt)
    dr0 = np.abs(drdt[:,:3].mean(axis=1))
    dr0.shape = (dr0.shape[0],1) #reshape to divide
    part = np.log(np.abs(drdt/dr0))
    exponent = np.log(np.abs(drdt/dr0)).mean(axis=0)
    mask= np.where(drdt[:,0] >-1e7)
    drdt_to_plot = drdt

    if 0:
        #odd slices because we skip the first frame
        rxb = ms.particle_x[:,2:]-ms.particle_x[:,1:-1]
        ryb = ms.particle_y[:,2:]-ms.particle_y[:,1:-1]
        rzb = ms.particle_z[:,2:]-ms.particle_z[:,1:-1]
        drb = rxb #np.sqrt(rxb**2+ryb**2+rzb**2)
        drdt_to_plot=drb/dt
        p4 = np.log( np.abs(rxb/dt)) + np.log( np.abs(ryb/dt))+np.log(np.abs(rzb/dt))
        morename = 't1'
    if 0:
        rxc= ms.particle_x
        ryc= ms.particle_y
        rzc= ms.particle_z
        rc = np.sqrt( rxc**2+ryc**2 +rzc**2)[:,1:] #because the first snap is repeated.
        drdtc = ((rc[:,1:]-rc[:,:-1])/dt)
        dr0 = np.abs(drdtc[:,:3].mean(axis=1))
        dr0.shape = (dr0.shape[0],1) #reshape to divide
        p4 = np.log( np.abs( drdtc))
        morename = 't2'
    if 0:
        rxc= ms.particle_x - ms.meanx2
        ryc= ms.particle_y - ms.meany2
        rzc= ms.particle_z - ms.meanz2
        rc = np.sqrt( rxc**2+ryc**2 +rzc**2)[:,1:] #because the first snap is repeated.
        drdtc = ((rc[:,1:]-rc[:,:-1])/dt)
        drdt_to_plot=drdtc
        dr0 = np.abs(drdtc[:,:3].mean(axis=1))
        dr0.shape = (dr0.shape[0],1) #reshape to divide
        p4 = np.log( np.abs( drdtc))
        morename='t3'
    if 1:
        #Here you can see particles coalessing.
        P0 = np.column_stack([ ms.particle_x[:,0], ms.particle_y[:,0], ms.particle_z[:,0]])
        R3 = np.zeros_like( ms.particle_x[:,1:])
        EPSILON = 5e-4
        morename = 't4'
        for ip,p in enumerate(P0):
            Ptemp = P0 - p
            dp = (Ptemp**2).sum(axis=1)
            args = np.argsort(dp)
            take = args[:5]
            #take = np.where(dp < EPSILON)[0]
            delta = 0
            for i in take:
                dx = ms.particle_x[ip,:]-ms.particle_x[i,:]
                dy = ms.particle_y[ip,:]-ms.particle_y[i,:]
                dz = ms.particle_z[ip,:]-ms.particle_z[i,:]
                dr = np.sqrt( dx**2+dy**2+dz**2)[1:]
                delta = dr + delta
            delta /= len(take)
            R3[ip,:] = delta
        drdt_to_plot = R3[:,1:]
        p4 = np.log( np.abs( drdt_to_plot))

    for ip in mask[0]:
        c = colormap.to_rgba( r[ip,0])
        ax[0][0].plot( t, r[ip,:],c=c, linewidth=0.1)
        ax[0][1].plot( mt, drdt_to_plot[ip,:], c=c)
        #ax[1][1].plot( mt, part[ip,:])
        ax[1][1].plot( mt, p4[ip,:],c=c)
        #ax.plot( ms.raw_y[ip,:])
    mp4=p4.mean(axis=0)
    ax[1][1].plot(mt, mp4,c='k')
    ax[0][1].plot(mt,drdt_to_plot.mean(axis=0), c='k')
    #ax[0][1].plot(mt,(drdt/dr0).mean(axis=0))

    axbonk(ax[0][1],ylabel='<v>')

    #for nt, time in enumerate(dt):
    #    if nt%10 > 0:
    #        continue
    #    color=rainbow_map(len(t))(nt)
    #    ax[1][0].hist( drdt[:,nt], histtype='step',color=color)
        
    fig.savefig('plots_to_sort/%s_lyapunov_test_%s_c%04d.png'%(name,morename,core_id))
    


