
from starter2 import *
import three_loopers_u500 as TL



class radial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.B = []
        self.rho=[]

        self.pB = []
        self.prho=[]

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
                ms = trackage.mini_scrubber(thtr,core_id, do_magnetic=True)
                p = np.array([ms.mean_x[nframe],ms.mean_y[nframe],ms.mean_z[nframe]])
                c = ds.arr(p,'code_length')
                r=0.5/128.

                sph = ds.sphere(c,r)

                rrr = sph['radius']
                dv = sph['cell_volume']
                B  = sph[BI]
                rho = sph['density']

                self.B.append( (B*dv*rho).sum()/(dv*rho).sum())
                self.rho.append( (rho*dv).sum()/(dv).sum())

                mask = ms.compute_unique_mask(core_id, 1/2048,nframe)
                pB = np.sqrt(ms.b2[mask,nframe])
                prho= ms.density[mask,nframe]
                dv  = ms.cell_volume[mask,nframe]
                pBmean = (pB*prho*dv).sum()/(prho*dv).sum()
                prhomean=(prho*dv).sum()/dv.sum()


                self.pB.append(pBmean)
                self.prho.append(prhomean)

                fig,ax=plt.subplots(1,2)
                ax0=ax[0];ax1=ax[1]

                bins = np.geomspace( rho.min(), rho.max())
                ax0.hist( rho.v, color='g', histtype='step', bins=bins, density=True)
                ax0.hist( prho, color='k', histtype='step', bins=bins, density=True)
                ax0.set_yscale('log')
                fig.savefig('plots_to_sort/hists_c%04d_n%04d'%(core_id,frame))



sim_list = ['u501']

if 'rrr' not in dir() or True:
    for sim in sim_list:
        rr = radial( TL.loops[sim])
        this_looper=TL.loops[sim]
        core_list=[323]
        core_list=np.unique( this_looper.tr.core_ids)[20:25]

        frame_list = TL.loops[sim].tr.frames[-1:]
        rrr = radial( TL.loops[sim])
        rrr.run( core_list=core_list,frame_list=frame_list)
if 1:

    fig,ax=plt.subplots(2,2)
    ax0=ax[0][0];ax1=ax[0][1]; ax2=ax[1][0];ax3=ax[1][1]
    
    the_x,the_y=np.log(rrr.rho), np.log(rrr.B)
    fit = np.polyfit( the_x,the_y,1)
    ax0.plot(the_x,fit[0]*the_x+fit[1])
    ax0.scatter(the_x,the_y)
    ax0.set_title("2d %0.2f"%fit[0])

    the_x,the_y=np.log(rrr.prho), np.log(rrr.pB)
    fit = np.polyfit( the_x,the_y,1)
    ax1.plot(the_x,fit[0]*the_x+fit[1])
    ax1.scatter(the_x,the_y)
    ax1.set_title("3d %0.2f"%fit[0])

    ax2.scatter( rrr.rho, rrr.prho)
    ax3.scatter( rrr.B, rrr.pB)
    ax2.set(xscale='log',yscale='log')
    ax3.set(xscale='log',yscale='log')



    fig.savefig('plots_to_sort/derp.png')
    plt.close('all')
