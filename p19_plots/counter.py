

from starter2 import *
import three_loopers_u500 as TL
from jinja2 import Environment, PackageLoader, select_autoescape
import jinja2
sim_list=['u501','u502','u503']

stuff={}
other_modes=set()
labels={}
for sim in sim_list:
    this_looper=TL.loops[sim]
    stuff[sim]={}


    for mode in ['One','Binary','Cluster']:
        stuff[sim][mode] = len(this_looper.core_by_mode[mode])

    for mode in this_looper.core_by_mode.keys():
        if mode.startswith('S') or mode == 'F1':
            stuff[sim][mode] = len(this_looper.core_by_mode[mode])
            other_modes.add(mode)
            lab = {'S':'Cluster', 'F':'Filament'}[ mode[0]]
            num = mode[1]
            label = "%s %s"%(lab,num)
            labels[mode]=label
            label

                          
if 1:   

    loader=jinja2.FileSystemLoader('.')
    env = jinja2.Environment(loader=loader)
    main_template  = env.get_template('p19_plots/table1_template.tex')
    fptr=open('plots_to_sort/table1.tex','w')
    other_modes = sorted(list(other_modes))
    fptr.write(main_template.render(stuff=stuff, first_modes=['One','Binary','Cluster'], other_modes=other_modes, labels=labels))

    



