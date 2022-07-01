from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

plt.close('all')

class dr_thing():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.rho_extents=None
    def run(self,core_list=None):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        if self.rho_extents is None: 
            self.rho_extents=davetools.extents()
            self.r_extents=davetools.extents()
            for nc,core_id in enumerate(all_cores):
                ms = trackage.mini_scrubber(thtr,core_id)
                if ms.nparticles < 20:
                    continue
                density = thtr.c([core_id],'density')
                self.rho_extents(density)
                self.r_extents(ms.r)

        tsorted = thtr.times
        self.core_list=core_list

        ds= self.this_looper.load(0)
        ad = ds.all_data()
        rmcore=rainbow_map(len(core_list))
        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id)
            self.ms = ms
            if ms.nparticles < 10:
                continue
            self.cores_used.append(core_id)
            fig, axd=plt.subplots(2,2)
            axd1 = axd[0][0]; axd2=axd[0][1]#sorry
            axd3 = axd[1][0]; axd4=axd[1][1]

            tmap=rainbow_map(ms.ntimes)
            asort =  np.argsort(thtr.times)
            density = thtr.c([core_id],'density')
            n0=asort[0]
            tsorted = thtr.times[asort]

            gx = thtr.c([core_id],'grav_x')
            gy = thtr.c([core_id],'grav_y')
            gz = thtr.c([core_id],'grav_z')
            grav = 1/(8*np.pi)*(gx*gx+gy*gy+gz*gz)
            grav = np.abs(thtr.c([core_id],'PotentialField'))

            #fig, axd1=plt.subplots(1,1)

            for n_count,n_time in enumerate(asort):
                mask2 = ms.compute_unique_mask(core_id, dx=1./2048,frame=n_count)
                mask2[:]=True
                time=thtr.times[n_time]
                c=tmap(n_count,mask2.sum())
                #c=rmcore(nc, mask2.sum())
                this_r=ms.r[:,n_time]+0
                this_r[ this_r<1/2048]=1/2048
                r_un = nar(sorted(np.unique(this_r)))

                mask2 = mask2 * (this_r > 0)

                axd1.scatter(this_r[mask2],density[mask2,n_time],c=c,label=thtr.times[n_time],s=0.1)
                axd2.hist(np.log10(density[mask2,n_time]), histtype='step',color=c[0])

                axd3.scatter( this_r[mask2], grav[mask2, n_time], c=c,s=0.1)

                the_x = np.log(this_r[mask2]); the_y=np.log(density[mask2,n_time])
                r,p=scipy.stats.pearsonr(the_x,the_y)
                axd4.scatter( thtr.times[n_time], r, c=tmap(n_count,1))
                
                the_x = np.log(this_r[mask2]); the_y=np.log(grav[mask2,n_time])
                r,p=scipy.stats.pearsonr(the_x,the_y)
                axd4.scatter( thtr.times[n_time], r, c=tmap(n_count,1),marker='*')

            davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                             xlim=[1/2048,self.r_extents.minmax[1]], ylim=self.rho_extents.minmax)
            davetools.axbonk(axd2,xscale='linear',yscale='log',xlabel='r',ylabel=r'$\rho$',
                             xlim=np.log10(self.rho_extents.minmax))
            axbonk(axd3,xscale='log',yscale='log',xlabel='r',ylabel='GE')
            outname = '%s/density_radius_c%04d.png'%(dl.output_directory,core_id)
            #outname = '%s/density_radius_n0000.png'%(dl.output_directory)
            fig.savefig(outname)
            print("saved "+outname)
            plt.close(fig)

#import three_loopers_1tff as TL
#
#r1 = dr_thing(TL.looper1)
#r1.run()
import three_loopers_six as TLM
r1 = dr_thing(TLM.loops['u601'])
r1.run()
#r1.run(core_list=[323])
