from starter2 import *
import colors
class all_mass_time():
    def __init__(self,this_loop):
        self.this_loop=this_loop
        self.masses=[]
        self.dof=[]
        self.run() #normally I don't do this.
    def run(self):

        tr = self.this_loop.tr
        self.times = tr.times/colors.tff

        x = tr.track_dict['x']
        y = tr.track_dict['y']
        z = tr.track_dict['z']
        dv = tr.track_dict['cell_volume']
        density = tr.track_dict['density']
        Lx = dv**(1./3)
        Lx_max = 1./128
        Lx_min = 1./2048
        Level = np.round(-np.log(Lx/Lx_max)/np.log(2)).astype('int')
        Nx = np.round(1./Lx).astype('int')
        Nx_max = np.round(1./Lx_min).astype('int')
        #Nx = Ny = Nz here. 
        index = x/Lx + Nx*y/Lx + Nx*Nx*(z/Lx)+ Nx_max*Nx_max*Nx_max*Level
        mask = np.ones_like(index, dtype='bool')
        for ntime in range(mask.shape[1]):
            print('nt',ntime)
            this_index=index[:,ntime]
            sort = np.argsort(this_index)
            isorted = this_index[sort]
            unsort = np.argsort(sort)
            mask[1:,ntime] = (isorted[1:]-isorted[:-1] != 0)
            mask[:,ntime] = mask[:,ntime]

        self.masses = (density*dv).sum(axis=0)
        mask = mask.astype('float')
        self.maskedmasses = (density*dv*mask).sum(axis=0)



import three_loopers_otherones as TLO
import three_loopers_tenfour as TL4
simnumb=[2]
if 'tools_core' not in dir():
    tools_core={}
    tools_other={}

    for numb in simnumb:
        tools_core[numb] = all_mass_time(  TL4.loops['u4%02d'%numb])
        tools_other[numb] = all_mass_time( TLO.loops['a%03d'%numb])
if 1:
    fig,ax=plt.subplots(1,1)
    for numb in simnumb:
        tool_o = tools_other[numb]
        tool_c = tools_core[numb]
        ax.plot( tool_o.times, tool_o.masses,'k--',label='others')
        ax.plot( tool_c.times, tool_c.masses,'r--',label='cores')
        ax.plot( tool_o.times, tool_o.maskedmasses,'k-',label='others')
        ax.plot( tool_c.times, tool_c.maskedmasses,'r-',label='cores')
        ax.legend(loc=2)
    axbonk(ax, xlabel=r'$t/t_{\rm{ff}}$', ylabel='M',yscale='log')
    fig.savefig('plots_to_sort/other_masses_%s.png'%'a%03d'%numb)
