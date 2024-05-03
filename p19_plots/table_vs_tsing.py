from starter2 import *
#import three_loopers_u500 as TL
import track_loader as TL
from collections import defaultdict
import tsing
reload(tsing)
plt.close('all')
class vs_tsing():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.tsing_list=[]
        self.tsung_list=[]
        self.mode_list=[]
        self.quan_list=defaultdict(list)
        self.cores_used=[]
    def run(self, core_list=None,tsing=None,frame=0):
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        for core_id in core_list:
            print('quan maker',core_id)
            self.cores_used.append(core_id)
            frame_tsing = tsing.tsing_frame[core_id]
            frame_tsung = tsing.tend_frame[core_id]
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            if frame == 'tsing':
                frame = frame_tsing
            if frame == 'tsung':
                frame = frame_tsung

            nframe = np.where( thtr.frames==frame)[0][0]
            ms.compute_ge(core_id)
            ms.compute_ke_rel(core_id)

            D = ms.density[:,nframe]
            DV= ms.cell_volume[:,nframe]
            volume=DV.sum()
            KE = ms.ke_rel[:,nframe]
            GE = ms.ge[:,nframe]
            self.tsing_list.append( tsing.tsing_core[core_id])
            self.tsung_list.append( tsing.tend_core[core_id])
            my_mass= (D*DV).sum()
            self.quan_list['mass'].append( my_mass)
            self.quan_list['volume'].append( DV.sum())
            
            my_ke = (KE*DV).sum()
            my_ge = (GE*DV).sum()

            self.quan_list['ge'].append(my_ge)
            self.quan_list['ke'].append(my_ke)

            my_v3 = (ms.rel_vmag[:,nframe]*DV).sum()/volume
            my_vt = (np.sqrt(ms.vt2_rel[:,nframe])*DV).sum()/volume
            my_vr = (ms.vr_rel[:,nframe]*DV).sum()/volume

            self.quan_list['avg_vr'].append(my_vr)
            self.quan_list['avg_v3'].append(my_v3)
            self.quan_list['avg_vt'].append(my_vt)

            #my_tff = np.sqrt( 3*np.pi/(32*colors.G*D))
            avg_rho= (D*DV).sum()/DV.sum()
            self.quan_list['avg_rho'].append(avg_rho)
            my_tff = np.sqrt( 3*np.pi/(32*colors.G*avg_rho))
            self.quan_list['avg_tff'].append(my_tff.mean())

            my_mode = this_looper.mode_dict[core_id]
            if 'Alone' in my_mode:
                self.mode_list.append('r')
            elif 'Binary' in my_mode:
                self.mode_list.append('g')
            elif 'Cluster' in my_mode:
                self.mode_list.append('b')
            else:
                print("Error in Mode.  Sort it out.")
                pdb.set_trace()

sim_list=['u501','u502','u503']
TL.load_tracks(sim_list)

if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

import convex_hull_tools as CHT
if 'ht' not in dir():
    ht = {}
    for sim in sim_list:
        ht[sim]=CHT.hull_tool( TL.loops[sim])
        ht[sim].make_hulls()

#sim_list=['u502']

if 'vs_0' not in dir():
    vs_0={} #initial frame
    vs_s={} #at tsing
    vs_u={} #at tsung

for sim in sim_list:
    if sim in vs_0:
        continue
    core_list=None
    #core_list=[300]
    #core_list = TL.loops[sim].core_by_mode['Alone']

    if 1:
        vs = vs_tsing(TL.loops[sim])
        vs.run(core_list=core_list,tsing=tsing_tool[sim], frame=0)
        vs_0[sim]=vs

    if 0:
        vs = vs_tsing(TL.loops[sim])
        vs.run(core_list=core_list,tsing=tsing_tool[sim], frame='tsung')
        vs_u[sim]=vs


    if 0:
        vs = vs_tsing(TL.loops[sim])
        vs.run(core_list=core_list,tsing=tsing_tool[sim], frame='tsing')
        vs_s[sim]=vs

if 'the_r' not in dir():
    the_r={}

if 0:
    fig,ax=plt.subplots(1,1)
    ax.plot( sorted(vs.tsung_list))
    fig.savefig('%s/test'%plot_dir)
if 1:
    #plot tffbar
    for frame, vs_tool in [[0,vs_0]]:#,['tsing',vs_s]
        quan_by_sim={}
        ncol=3
        nrow=1
        fig,axes_all=plt.subplots(nrow,ncol, figsize=(8,3))
        axes=[axes_all]
        exty = extents()
        extx = extents()
        for nsim,sim in enumerate(sim_list):
            the_r[sim]={}
            vs = vs_tool[sim]
            htL = ht[sim]

            the_x=nar(vs.tsung_list)
            extx(the_x)
            print('word')
            TheR = {}
            TheLR = {}
            quan_list=list(vs.quan_list.keys()) +['hull_volumes']
            #quan_list = ['density']
            quan_list=['avg_tff']#,'avg_rho']
            for nq,quan in enumerate(quan_list):
                the_r[sim][quan]=defaultdict(list)
                nc=nsim
                nr=nq
                ok_vals = True
                ok_last = True
                the_y = nar(vs.quan_list[quan])
                if 1:
                    rho = nar(vs.quan_list['avg_rho'])
                    rhonorm = (rho-rho.mean())/rho.std()
                    ok_vals = np.abs(rhonorm)<5
                if 1:
                    tlast = vs.this_looper.tr.times[-2]
                    not_last_frame = (the_x*colors.tff)/tlast
                    ok_last = not_last_frame < 1

                for mode in ['r']:#np.unique(vs.mode_list):
                    ok = nar(vs.mode_list)==mode
                    #ok *= ok_vals
                    #ok *= ok_last
                    if 1:
                        okok = nar(vs.quan_list['volume']) < 0.5
                        ok = ok*okok
                    R = scipy.stats.pearsonr( nar(the_x)[ok],nar(the_y)[ok])
                    the_r[sim][quan][mode].append(R)
                    print(R)

                    axes[nr][nc].scatter(nar(the_x)[ok],nar(the_y)[ok],c=nar(vs.mode_list)[ok])
                    axes[nr][nc].text(0.3,0.052,"r=%0.2f"%R[0])
                    axes[nr][nc].set_ylabel(r'$\overline{t_{\rm{ff}}}/t_{\rm{ff}}$')
                    axes[nr][nc].set_xlabel(r'$t_{\rm{sung}}$')
                    exty(nar(the_y)[ok])
        for ax in axes[0]:
            ax.set(xlim=extx.minmax,ylim=exty.minmax)
        fig.tight_layout()
        fig.savefig('%s/tff_tsung'%(plot_dir))
                
if 1:
    #
    # plot all the quantities
    #
    quan_by_sim={}
    for frame, vs_tool in [[0,vs_0]]:#,['tsing',vs_s]
        quan_by_sim={}
        for sim in sim_list:
            quan_by_sim[sim]={}
            the_r[sim]={}
            vs = vs_tool[sim]
            htL = ht[sim]

            the_x=vs.tsung_list
            print('word')
            TheR = {}
            TheLR = {}
            quan_list=list(vs.quan_list.keys()) +['hull_volumes']
            #quan_list = ['density']
            ncol=3
            nrow=len(quan_list)//ncol+1
            fig,axes=plt.subplots(nrow,ncol)
            for nq,quan in enumerate(quan_list):
                the_r[sim][quan]=defaultdict(list)
                nc=nq//ncol
                nr=nq%ncol
                ok_vals = True
                if quan in ['hull_volumes']:
                    the_y = nar(htL.hull_volumes)
                else:
                    the_y = nar(vs.quan_list[quan])
                if quan in ['ge', 'avg_vr']:
                    the_y = np.abs(the_y)
                if quan in ['ge','ke','mass']:
                    the_y = nar(the_y)/nar(vs.quan_list['volume'])
                if 1:
                    rho = nar(vs.quan_list['avg_rho'])
                    rhonorm = (rho-rho.mean())/rho.std()
                    ok_vals = np.abs(rhonorm)<5



                do_log=True
                print(quan)
                for mode in ['r']:#np.unique(vs.mode_list):
                    ok = nar(vs.mode_list)==mode
                    ok *= ok_vals
                    if 1:
                        okok = nar(vs.quan_list['volume']) < 0.5
                        ok = ok*okok
                    R = scipy.stats.pearsonr( nar(the_x)[ok],nar(the_y)[ok])
                    the_r[sim][quan][mode].append(R)
                    if np.abs(R[0]) < 0.01:
                        st = "%0.3f"%R[0]
                    else:
                        st = "%0.2f"%R[0]
                    quan_by_sim[sim][quan]=st
                    axes[nc][nr].scatter(nar(the_x)[ok],nar(the_y)[ok],c=nar(vs.mode_list)[ok])
                    axes[nc][nr].set(title="%0.2f"%R[0])
                axes[nc][nr].set(ylabel=quan)
            fig.tight_layout()
            fig.savefig('%s/pearson_%s'%(plot_dir,sim))
                
if 1:
    #make table
    for frame, vs_tool in [[0,vs_0]]:#,['tsing',vs_s]
        quan_list=['volume','hull_volumes','mass','ge','ke','avg_v3','avg_vr','avg_vt','avg_tff']
        quan_tex={'volume':r'$N_{\rm{particles}}$','mass':r'$\overline{ \rho }$','ge':r'$\overline{ \EG}$', 'ke':r'$\overline{ \EK }$','avg_v3':r'$\sigma_{3d}$','avg_vr':r'$\overline{ v_r }$','avg_vt':r'$\overline{v_t}$','avg_tff':r'$\overline{t_{\rm{ff}}}$','hull_volumes':r'$V_{\rm{hull}}$'}
        import jinja2
        loader=jinja2.FileSystemLoader('.')
        env = jinja2.Environment(loader=loader)
        main_template  = env.get_template('p19_plots/table2_template.tex')
        fptr=open('%s/table2.tex'%plot_dir,'w')
        #quan_tex=dict(zip(quan_list,quan_list))
        fptr.write(main_template.render(quan_list=quan_list,R=quan_by_sim,quan_tex=quan_tex))
        fptr.close()




if 0:
    #hull volume
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            htL = ht[sim]
            #print(nar(htL.cores_used) - nar( vs.cores_used))
            fig,ax=plt.subplots(1,1)
            the_x,the_y=vs.tsing_list, nar(htL.hull_volumes)
            R = scipy.stats.pearsonr( the_x,the_y)
            LR = scipy.stats.pearsonr( the_x,np.log10(the_y))
            ax.scatter( the_x,the_y, color=vs.mode_list)
            ax.set(xlabel='tsing',ylabel='volume',xscale='linear',yscale='log', title="R %0.2f %0.2f"%(R[0],LR[0]))
            fig.savefig('plots_to_sort/hull_volume_vs_tsing_%s_frame_%s'%(sim,frame))
            plt.close(fig)

if 0:
    #mean tff
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,ax=plt.subplots(1,1)
            the_x,the_y=vs.tsing_list, nar(nar(vs.quan_list['avg_tff']))
            R = scipy.stats.pearsonr( the_x,the_y)
            LR = scipy.stats.pearsonr( the_x,np.log10(the_y))
            ax.scatter( the_x,the_y, color=vs.mode_list)
            ax.set(xlabel='tsing',ylabel='volume',xscale='linear',yscale='log', title="R %0.2f %0.2f"%(R[0],LR[0]))
            fig.savefig('plots_to_sort/tff_vs_tsing_%s_frame_%s'%(sim,frame))
            plt.close(fig)

if 0:
    #cell volume
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,ax=plt.subplots(1,1)
            the_x,the_y=vs.tsing_list, nar(nar(vs.quan_list['volume']))
            R = scipy.stats.pearsonr( the_x,the_y)
            LR = scipy.stats.pearsonr( the_x,np.log10(the_y))
            ax.scatter( the_x,the_y, color=vs.mode_list)
            ax.set(xlabel='tsing',ylabel='volume',xscale='linear',yscale='log', title="R %0.2f %0.2f"%(R[0],LR[0]))
            fig.savefig('plots_to_sort/volume_vs_tsing_%s_frame_%s'%(sim,frame))
            plt.close(fig)

if 0:
    #mean density
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,ax=plt.subplots(1,1)
            the_x,the_y=vs.tsing_list, nar(vs.quan_list['mass'])/nar(vs.quan_list['volume'])
            R = scipy.stats.pearsonr( the_x,the_y)
            LR = scipy.stats.pearsonr( the_x,np.log10(the_y))
            ax.scatter( the_x,the_y, color=vs.mode_list)
            ax.set(xlabel='tsing',ylabel='mean density',xscale='linear',yscale='log', title="R %0.2f %0.2f"%(R[0],LR[0]))
            fig.savefig('plots_to_sort/density_vs_tsing_%s_frame_%s'%(sim,frame))
            plt.close(fig)

if 0:
    #mass
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,ax=plt.subplots(1,1)
            the_x,the_y=vs.tsing_list, vs.quan_list['mass']
            R = scipy.stats.pearsonr( the_x,the_y)
            LR = scipy.stats.pearsonr( the_x,np.log10(the_y))
            ax.scatter( the_x,the_y, color=vs.mode_list)
            ax.set(xlabel='tsing',ylabel='InitialMass',xscale='linear',yscale='log', title="R %0.2f %0.2f"%(R[0],LR[0]))
            fig.savefig('plots_to_sort/mass_vs_tsing_%s_frame_%s'%(sim,frame))
            plt.close(fig)
if 0:
    #velocity
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s], ['tsung',vs_u]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,axes=plt.subplots(1,3)
            ax0=axes[0];ax1=axes[1];ax2=axes[2]
            M = vs.quan_list['mass']
            M = 1
            ext=extents()
            the_x,the_y=vs.tsing_list, vs.quan_list['avg_v3']
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,the_y)[0]
            ax0.scatter(the_x,the_y , color=vs.mode_list)
            ax0.set(yscale='linear',xlabel='tsing',title='V3 lin r=%0.2f'%R)
            the_x,the_y=vs.tsing_list, np.abs(vs.quan_list['avg_vr'])
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,the_y)[0]
            #ax1.scatter( vs.tsing_list, np.abs(the_y), color=vs.mode_list)
            ax1.scatter( vs.tsing_list, the_y, color=vs.mode_list)
            ax1.set(yscale='linear',xlabel='tsing',title='VR lin r=%0.2f yv'%R)
            the_x,the_y=vs.tsing_list, vs.quan_list['avg_vt']
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,the_y)[0]
            ax2.scatter( vs.tsing_list, the_y, color=vs.mode_list)
            ax2.set(yscale='linear',xlabel='tsing',title='vt lin r=%0.2f'%R)
            ax0.set(ylim=[0,ext.minmax[1]])
            ax1.set(ylim=[0,ext.minmax[1]])
            ax2.set(ylim=[0,ext.minmax[1]])
            fig.savefig('plots_to_sort/vel_vs_tsing_%s_frame_%s'%(sim,frame))

if 0:
    #energies
    for frame, vs_tool in [[0,vs_0],['tsing',vs_s], ['tsung',vs_u]]:
        for sim in sim_list:
            vs = vs_tool[sim]
            fig,axes=plt.subplots(1,3)
            ax0=axes[0];ax1=axes[1];ax2=axes[2]
            M = vs.quan_list['mass']
            M = 1
            ext=extents()
            the_x,the_y=vs.tsing_list, np.abs(vs.quan_list['ge'])
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,np.log10(the_y))[0]
            ax0.scatter(the_x,the_y , color=vs.mode_list)
            ax0.set(yscale='log',xlabel='tsing',title='EG r=%0.2f'%R)
            the_x,the_y=vs.tsing_list, np.abs(vs.quan_list['ke'])
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,np.log10(the_y))[0]
            #ax1.scatter( vs.tsing_list, np.abs(the_y), color=vs.mode_list)
            ax1.scatter( vs.tsing_list, the_y, color=vs.mode_list)
            ax1.set(yscale='log',xlabel='tsing',title='KE r=%0.2f yv'%R)
            the_x,the_y=vs.tsing_list, np.abs(nar(vs.quan_list['ge'])/nar(vs.quan_list['ke']))
            ext(nar(the_y))
            R = scipy.stats.pearsonr(the_x,the_y)[0]
            ax2.scatter( vs.tsing_list, the_y, color=vs.mode_list)
            ax2.set(yscale='linear',xlabel='tsing',title='EG/EK lin r=%0.2f'%R)
            ax0.set(ylim=[0,ext.minmax[1]])
            ax1.set(ylim=[0,ext.minmax[1]])
            ax2.set(ylim=[0,ext.minmax[1]])
            fig.savefig('plots_to_sort/eng_vs_tsing_%s_frame_%s'%(sim,frame))

if 0:
    #energy vs energy
    for sim in sim_list:
        vs = vs_tool[sim]
        fig,ax=plt.subplots(1,1)
        ax.scatter(np.abs(vs.quan_list['ge']),np.abs(vs.quan_list['ke']), color=vs.mode_list)
        ext = extents()
        ext(np.abs(vs.quan_list['ge']))
        ext(np.abs(vs.quan_list['ke']))
        ax.plot( ext.minmax,ext.minmax,c='k')
        ax.plot( ext.minmax,10*nar(ext.minmax),c='k')
        ax.set(xscale='log',yscale='log',xlabel='GE',ylabel='KE')
        fig.savefig('plots_to_sort/eng_vs_eng_%s'%sim)


if 0:
    #energy vs mass
    for sim in sim_list:
        vs = vs_tool[sim]
        fig,ax=plt.subplots(1,1)
        the_x,the_y=np.abs(vs.quan_list['mass']),np.abs(vs.quan_list['ge'])
        ax.scatter(the_x,the_y, color=vs.mode_list)
        pfit = np.polyfit(np.log10(the_x),np.log10(the_y),1)
        ax.plot( the_x, 10**( pfit[0]*np.log10(the_x)+pfit[1]))
        #ext = extents()
        #ext(np.abs(vs.quan_list['ge']))
        #ext(np.abs(vs.quan_list['ke']))
        #ax.plot( ext.minmax,ext.minmax,c='k')
        #ax.plot( ext.minmax,10*nar(ext.minmax),c='k')
        ax.set(xscale='log',yscale='log',xlabel='M',ylabel='GE')
        fig.savefig('plots_to_sort/eng_vs_mass_%s'%sim)


