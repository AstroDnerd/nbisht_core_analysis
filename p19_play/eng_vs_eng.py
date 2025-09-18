from starter2 import *
import three_loopers_u500 as TL
from collections import defaultdict
import tsing
reload(tsing)
plt.close('all')
class stuff_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.tsing_list=[]
        self.mode_list=[]
        self.quan_list=[defaultdict(list), defaultdict(list)]
    def run(self, core_list=None,tsing=None,frame=0):
        this_looper=self.this_looper
        thtr=this_looper.tr
        if core_list is None:
            core_list = np.unique(this_looper.tr.core_ids)

        for core_id in core_list:
            frame_tsing = tsing.tsing_frame[core_id]
            ms = trackage.mini_scrubber(this_looper.tr,core_id)
            nframe = np.where( thtr.frames==frame)[0][0]
            ntsing = np.where( thtr.frames==frame_tsing)[0][0]
            ms.compute_ge(core_id)
            ms.compute_ke_rel(core_id)
            self.tsing_list.append( tsing.tsing_core[core_id])
            my_mode = this_looper.mode_dict[core_id]
            if 'Alone' in my_mode:
                self.mode_list.append('r')
            elif 'Binary' in my_mode:
                self.mode_list.append('g')
            elif 'Cluster' in my_mode:
                self.mode_list.append('b')

            for nf,FFF in enumerate([0,ntsing]):

                D = ms.density[:,FFF]
                DV= ms.cell_volume[:,FFF]
                KE = ms.ke_rel[:,FFF]
                GE = ms.ge[:,FFF]
                R  = ms.r[:,FFF]
                my_mass= (D*DV).sum()
                self.quan_list[nf]['initial_mass'].append( my_mass)
                self.quan_list[nf]['initial_volume'].append( DV.sum())
                self.quan_list[nf]['Rmax'].append(R.max())
                
                my_ke = (KE*DV).sum()
                my_ge = (GE*DV).sum()

                self.quan_list[nf]['ge'].append(my_ge)
                self.quan_list[nf]['ke'].append(my_ke)




sim_list=['u501','u502','u503']

if 'tsing_tool' not in dir():
    tsing_tool={}
    for ns,sim in enumerate(sim_list):
        obj=tsing.te_tc(TL.loops[sim])
        tsing_tool[sim]=obj
        tsing_tool[sim].run()

if 's_tool' not in dir():
    s_tool={}
for sim in sim_list:
    if sim in s_tool:
        continue
    s = stuff_tool(TL.loops[sim])
    core_list=None
    #core_list=[300]
    #core_list = TL.loops[sim].core_by_mode['Alone']
    s.run(core_list=core_list,tsing=tsing_tool[sim])
    s_tool[sim]=s


if 0:
    #energy vs energy
    for sim in sim_list:
        s = s_tool[sim]
        fig,ax=plt.subplots(1,2)
        ext = extents()
        for nf in [0,1]:
            #weighting by mass is fascinating.
            M = s.quan_list[nf]['initial_mass']
            #M=1
            ax[nf].scatter(np.abs(s.quan_list[nf]['ge'])/M,np.abs(s.quan_list[nf]['ke'])/M, color=s.mode_list)
            ext(np.abs(s.quan_list[nf]['ge'])/M)
            ext(np.abs(s.quan_list[nf]['ke'])/M)
        for nf in [0,1]:
            ax[nf].set(xscale='log',yscale='log',xlabel='GE',ylabel='KE')
            ax[nf].plot( ext.minmax,ext.minmax,c='k')
            ax[nf].plot( ext.minmax,10*nar(ext.minmax),c='k')
        fig.savefig('plots_to_sort/eng_vs_eng_0_tsing_%s'%sim)

if 1:
    #energy vs mass
    for sim in sim_list:
        s = s_tool[sim]
        fig,ax=plt.subplots(1,2)
        ext = extents()
        for nf in [0,1]:
            #weighting by mass is fascinating.
            M = nar(s.quan_list[nf]['initial_mass'])
            R = nar(s.quan_list[nf]['Rmax'])
            GE =np.abs(nar(s.quan_list[nf]['ge']))
            #GE /= colors.G*M*M/R
            GE = GE/M
            #M=1
            ax[nf].scatter(M,GE, color=s.mode_list)
            ext(M)
            ext(GE)
        for nf in [0,1]:
            ax[nf].set(xscale='log',yscale='log',xlabel='mass',ylabel='GE')
            ax[nf].plot( ext.minmax,ext.minmax,c='k')
            ax[nf].plot( ext.minmax,10*nar(ext.minmax),c='k')
        fig.savefig('plots_to_sort/eng_vs_mass_0_tsing_%s'%sim)

