
from starter2 import *
import three_loopers_u500 as TL



class radial():
    def __init__(self,this_looper):
        self.this_looper=this_looper

    def run(self, core_list=None,frame_list=None):
        this_looper=self.this_looper
        thtr=this_looper.tr


        if core_list is None:
            core_list=np.unique(this_looper.tr.core_ids)
        if frame_list is None:
            frame_list=this_looper.tr.frames

        for iframe, frame in enumerate(frame_list):
            ds = this_looper.load(frame)
            radius=1e-2
            radius = ds.arr(radius,'code_length')
            nframe=np.where(thtr.frames==frame)[0][0]

            los = 0
            XYZ = 'xyz'[los]
            BI = 'magnetic_field_'+XYZ
            BI = 'magnetic_field_strength'

            for core_id in core_list:
                print('yay',core_id)
                ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True)
                p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                c = ds.arr(p,'code_length')
                r=2/128.

                sph = ds.sphere(c,r)

                rrr = sph['radius']
                dv = sph['cell_volume']
                B  = sph[BI]
                args=np.argsort(rrr)
                rhos = sph['density'][args]
                rrrs = rrr[args]
                dvs  = dv[args]
                Bs   = B[args]
                #Bs = Bs/B[args][0]
                #rhos /= rhos[0]
                mean_rho = np.cumsum( rhos*dvs)/np.cumsum(dvs)
                #mean_B = np.cumsum( rhos*Bs*dvs)/np.cumsum(rhos*dvs)
                mean_B = np.cumsum( Bs*dvs)/np.cumsum(dvs)

                for end in [int(args.size/4), int(args.size/2), int(3/4*args.size)]:
                    fit = np.polyfit( np.log(mean_rho[:end]), np.log(mean_B[:end]),1)
                    print(fit)



                fig,ax=plt.subplots(1,2)
                ax0=ax[0]; ax1=ax[1]
                ax0.plot( rrrs, mean_rho, c='r')
                ax0.plot( rrrs, mean_B, c='b')
                ax1.plot( np.log( mean_rho), np.log( mean_B))
                ax0.set(xscale='log')
                #ax0.set_yscale('symlog',linthresh=1e-1)
                ax0.set_yscale('log')
                fig.tight_layout()
                fig.savefig('plots_to_sort/alpha_r_c%04d_n%04d.png'%(core_id,frame))



sim_list = ['u501']

for sim in sim_list:
    rr = radial( TL.loops[sim])
    this_looper=TL.loops[sim]
    core_list=[323]
    core_list=np.unique( this_looper.tr.core_ids)[20:25]

    frame_list = TL.loops[sim].tr.frames[-1:]
    rrr = radial( TL.loops[sim])
    rrr.run( core_list=core_list,frame_list=frame_list)
