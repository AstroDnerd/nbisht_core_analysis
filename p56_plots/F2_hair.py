from starter2 import *
import colors

import hair_dryer_vintage
reload(hair_dryer)


#import three_loopers_mountain_top as TLM
#import three_loopers_tenfour as TL4
#MOD = TLM
#sim_list=['u301','u302','u303']
import three_loopers_six as TL6
MOD = TL6
sim_list=['u601','u602','u603']

import three_loopers_u500 as TL5
MOD = TL5
sim_list=['u501','u502','u503']

hair_tools={}
for ns,this_simname in enumerate(sim_list):
    hair_tools[ns+1]= hair_dryer.hair_tool( MOD.loops[this_simname])
hair_tools[1].sim_name = 'sim 1'
hair_tools[2].sim_name = 'sim 2'
hair_tools[3].sim_name = 'sim 3'

if 0:
    for ns,this_simname in  enumerate(sim_list):
        hair_tools[ns+1].run( )#core_list=[10,11])#core_list=[10])

if 0:
    hair_tools[2].run( core_list = [3,32])
    hair_tools[2].run( core_list = [378])
if 1:
    plt.close('all')
    core_list_1=[0,8,27,32,37,44,84,275, 323]
    core_list_2=[3,32,378]
    #core_list_2=[32]
    core_list_3=[233,235,186]
    fig,ax=plt.subplots(1,1, figsize=(4,3),dpi=400)
    #fig,ax=plt.subplots(1,1,dpi=200)
    #fig,ax=plt.subplots(1,1)
    for sim in [1,2,3]:
        version=''
        core_list = {1:core_list_1, 2:core_list_2,3:core_list_3}[sim]
        for core_id in core_list:
            ax.clear()
            hair_tools[sim].run( core_list = [core_id], external_ax=ax)
            outname = 'plots_to_sort/%s_blowing_hair_%sc%04d%s.png'%('u50%s'%sim,'',core_id,version)
            fig.tight_layout()
            fig.savefig(outname)
            print(outname)
if 0:
    hair_tools[1].run( core_list = [323])
    hair_tools[1].run( core_list = [24])
    hair_tools[2].run( core_list = [3,32])
    hair_tools[2].run( core_list = [378])
    hair_tools[3].run( core_list = [233,235])
    hair_tools[3].run( core_list = [186])
