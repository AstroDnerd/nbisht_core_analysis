from starter2 import *
#import pdb

import data_locations as dl
reload(dl)

file_list=glob.glob('%s/track_indfix_sixteenframe_core_*.h5'%dl.snapshot_location)

plt.close('all')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/particle_error_test_c0031_threeframes.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/particle_error/track_indfix_sixteenframe_core_0031.h5')
#file_list=glob.glob('/scratch1/dcollins/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
#file_list=glob.glob('/home/dcollins/scratch/Paper19/track_index_fix/track_indfix_sixteenframe_core_*.h5')
#print("my name is Dan")
#print(file_list)
#pdb.set_trace()
for nfile,fname in enumerate(file_list):#[:3]):
    #0164.h5
    t1 = fname.split("/")[-1]
    #l = len("track_three_to_test_core_")
    #l = len("track_sixteen_frames_core_")
    l = len("track_indfix_sixteenframe_core_")

    this_cor = int(t1[l:l+4]) #[fname.index('_'):]
    #this_cor=31
    #if this_cor not in  [12]:#, 31]:
    #    continue
    print(this_cor)
    this_looper=looper.core_looper(directory=dl.enzo_directory)
#directory)
    trw.load_loop(this_looper,fname)
    thtr = this_looper.tr
    all_cores = np.unique(thtr.core_ids)
    core_list=all_cores
    rm = rainbow_map(len(all_cores))

    if 1:
        #time plots
        asort =  np.argsort(thtr.times)
        n0=asort[0]
        tsorted = thtr.times[asort]
        print("my name is dan")
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            density = thtr.c([core_id],'density')
            mag = thtr.c([core_id],'magnetic_field_strength')
            #field_name='magnetic_field_strength'
            #rho = b^x
            #ln rho = x ln b
            #ln rho/ln b = x
            field = np.zeros_like(density)
            print(np.shape(field))
            ok = mag != 0
            field[ok] = np.log(density[ok])/np.log(mag[ok])
            field_name = 'dbx'
            tmap=rainbow_map(ms.ntimes)
            plt.clf()
            norm = mpl.colors.Normalize()
            if 0:
                norm.autoscale( np.log10(field[:,n0]))
            if 1:
                norm.autoscale( field[:,n0])
            cmap = mpl.cm.jet
            color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
            fighist,axhist=plt.subplots(1,1)
            figlines,axlines=plt.subplots(1,1)

            for npart in range(ms.nparticles):
                if 0:
                    c = color_map.to_rgba(np.log10(field[npart,n0]))
                if 1:
                    c = color_map.to_rgba(field[npart,n0])
                #plt.plot( tsorted, field[npart,asort],c='k',linestyle=':',marker='*')
                axlines.plot( tsorted, field[npart,asort],c=c,linewidth=.1)#linestyle=':') #put back
            #err= np.exp(np.log(field).std(axis=0)[asort])
            #err= np.exp(np.log(field).std(axis=0)[asort])
            #plt.plot(tsorted, field.mean(axis=0)[asort],c='k')
            #axlines.errorbar(tsorted, field.mean(axis=0)[asort],c='k',yerr=err)
            #plt.plot(tsorted, field.mean(axis=0),c='k')
            nbins=30
            if 0:
                bins = np.logspace(np.log10(field.min()),np.log10(field.max()),nbins+1)
            if 1:
                bins = np.logspace(field.min(),field.max(),nbins+1)
            binc=0.5*(bins[:-1]+bins[1:])
            nsteps=field.shape[1]
            wire = np.zeros([nbins,nsteps])
            tm = rainbow_map(len(asort))
            field_mean=np.zeros(asort.size)
            for i,n in enumerate(asort):
                print(n)
                hist ,bin_edge = np.histogram(field[:,n],bins=bins)
                binc_2 = 0.5*(bin_edge[:-1]+bin_edge[1:])
                wire[:,i]=hist
                field_mean[i] = (wire[:,i]*binc_2).sum()/wire[:,i].sum()
                axhist.plot(binc_2, wire[:,i],c=tm(i))
                axhist.plot([field_mean[i]]*2,[1,100],c=tm(i))
            thex = np.tile(tsorted,(nbins,1))
            they = np.tile(binc,(len(tsorted),1)).transpose()
            they=np.log10(they)
            axhist.set_xscale('log')
            axhist.set_yscale('symlog',linthreshy=5)
            axhist.set_xlabel(r'$\rho$'); axhist.set_ylabel(r'$N$')
            fighist.savefig('image_tracks/multihist_%s_c%04d.png'%(field_name,core_id))
            if 1:
                figd,axd=plt.subplots(1,1)
                norm = mpl.colors.LogNorm(vmin=1,vmax=0.1*wire.max())
                #norm.autoscale( sorted(wire.flatten()))
                cmap = mpl.cm.jet
                cmap.set_under('w')
                color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
                p=axd.pcolormesh(thex,they,wire,norm=norm)
                if 0:
                    axd.scatter(tsorted, np.log10(field_mean),c='k',label='mean')
                if 1:
                    axd.scatter(tsorted, field_mean,c='k',label='mean')
                #axd.scatter(tsorted, np.log10(density_mean),c='k',label='mean')
                cbar = figd.colorbar(p)
                #figd.savefig('hist.png')
                figd.savefig('image_tracks/%s_hist_c%04d.png'%(field_name,core_id))
            if 0:
                fig3d=plt.figure()
                ax3d=fig3d.add_subplot(111,projection='3d')
                wire[wire>0]=np.log10(wire[wire>0])

                ax3d.plot_wireframe(thex,they,wire)
                ax3d.set_xlabel('time')
                ax3d.set_ylabel('log field')
                fig3d.savefig('fig3d.png')
            

            if 0:
                t0 = thtr.times[asort][0]
                t1 = thtr.times[asort][-1]
                rho0 = np.mean(density[:,asort[0]])
                rho1 = np.mean(density[:,asort[-1]])
                alpha = 1.8
                tc = t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
                G=1620./(4*np.pi)
                tff_global = np.sqrt(3*np.pi/(32*G*1))
                tff_local = np.sqrt(3*np.pi/(32*G*rho0))
                rhot = rho0*(1-(tsorted/tc)**2)**-alpha
                rho_c = 3*np.pi/(32*G*tc**2)
                freefall[core_id].set(rho0=rho0,rho1=rho1,tc=tc, tff_global=tff_global,
                                      tff_local=tff_local,rho_c=rho_c)
                axlines.plot( tsorted, rhot, c='r')

                axd.scatter( tsorted, np.log10(rhot), c='r',label='free fall')
                figd.savefig('image_tracks/density_hist_c%04d.png'%core_id)
                print('wut')
                continue
                axlines.text( tsorted[0], rho1, r'$tc = %0.2e \rho_c = %0.2e$'%(tc,rho_c))
                axlines.text( tsorted[0], 0.5*rho1, r'$tc/tff = %0.2e$'%(tc/tff_global))
                for i,n in enumerate(asort[0:1]):
                    timeline=plt.plot( [tsorted[i]]*2,[1,1e8],c=[0.5]*3,linewidth=0.1)
                    timetext=plt.text( tsorted[i], 1e8, 'n=%d'%thtr.frames[n])
                    outname = 'image_tracks/rho_t_fit2_c%04d_s%04d.png'%(core_id,i)
                    axlines.set_yscale('log')
                    figlines.savefig(outname)
