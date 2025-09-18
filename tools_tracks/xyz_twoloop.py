
from starter2 import *
import data_locations as dl
import davetools
import colors
reload(davetools)
reload(trackage)


class tracks():
    def __init__(self,looper1, looper2):
        self.this_looper=looper1
        self.other_looper=looper2
        self.cores_used=[]


    def run(self,core_list=None, output_suffix=''):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        otr = self.other_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        tsorted = thtr.times/colors.tff
        self.core_list=core_list

        for nc,core_id in enumerate(core_list):
            ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
            ms2 = trackage.mini_scrubber(otr,core_id,do_velocity=True)

            self.ms = ms
            if ms.nparticles < 10:
                continue
            print('tracks %s %d'%(self.this_looper.sim_name,core_id))
            self.cores_used.append(core_id)

            for through in [1]:
                fig, axd=plt.subplots(2,2)
                label='no'
                axl = axd.flatten()
                axd1 = axd[0][0]; axd2=axd[0][1]#sorry
                axd3 = axd[1][0]; axd4=axd[1][1]

                if through:
                    this_p = np.stack([ms.this_x,ms.this_y,ms.this_z])
                    o_p = np.stack([ms2.this_x,ms2.this_y,ms2.this_z])
                    label = 'track_shifted_c%04d'%core_id
                else:
                    this_p = np.stack([ms.raw_x,ms.raw_y,ms.raw_z])
                    label = 'track_raw_c%04d'%core_id

                for nnn in range(3):
                    sl = slice(None)
                    axl[nnn].plot(tsorted,this_p[nnn][sl,:].transpose(),linewidth=0.1,c=[0.5]*4)
                    axl[nnn].plot(tsorted,o_p[nnn][sl,:].transpose(),linewidth=0.1,c=[0.5,1.0,0.5,0.5])
                axl[3].plot( tsorted, ms.rel_vmag[sl,:].transpose(), linewidth=0.1, c=[0.5]*4)
                axl[3].plot( tsorted, ms2.rel_vmag[sl,:].transpose(), linewidth=0.1, c=[0.5,1.0,0.5,0.5])
                axl[3].set(yscale='log')
                outname='plots_to_sort/%s_%s%s.pdf'%(self.this_looper.sim_name,label, output_suffix)
                fig.savefig(outname)
                print(outname)
                plt.close(fig)

#import three_loopers_tenfour as TL4


#for simname in TL4.loops:
#    tr = tracks( TL4.loops[simname])
#    tr.run()


#tr = tracks( new_looper)
#tr.run(core_list=[3])

