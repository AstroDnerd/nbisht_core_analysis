from starter2 import *

import convex_hull_tools as CHT
import matplotlib.colors as colors

reload(CHT)
import hair_dryer
reload(hair_dryer)
import close_tool
#import three_loopers_tenfour as TL4
import three_loopers_six as TL
sim_list=['u601','u602','u603']
#sim_list=['u401','u402','u403']
#sim_list=['u402']
import multiplots
reload(multiplots)
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(TL.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( TL.loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( TL.loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 'overlap_dict' not in dir():
    overlap_dict={}
    overlap_numb={}
    particles = {}
    ratio_matrix={}
    ratio_matrixb={}
    particle_matrix={}

    for ns,this_simname in enumerate(sim_list):
        htool = ht[this_simname]
        overlap_dict[this_simname] = np.zeros( [len(htool.cores_used)]*2) -1
        overlap_numb[this_simname] = np.zeros( [len(htool.cores_used)]*2) -1
        ratio_matrix[this_simname] = np.zeros( [len(htool.cores_used)]*2) 
        ratio_matrixb[this_simname] = []
        particle_matrix[this_simname] = np.zeros( [len(htool.cores_used)]*2)
        particles[this_simname]=[]
        for nc1,core_id_1 in enumerate(htool.cores_used):
            particles[this_simname].append(htool.nparticles[nc1])
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                val = htool.overlaps[core_id_1][nc2]
                num = htool.overlap_number[core_id_1][nc2]
                overlap_dict[this_simname][nc1,nc2] = val
                overlap_numb[this_simname][nc1,nc2] = num

                a,b=htool.overlaps[core_id_1][nc2], htool.overlaps[core_id_2][nc1] 
                ratio = max([a,b])
                rat= sorted( [a,b])
                if rat[1] == 0:
                    ratio = max([a,b])
                else:
                    ratio=rat[0]/rat[1]
                ratio_matrix[this_simname][nc1,nc2]=ratio
                particle_matrix[this_simname][nc1,nc2] = np.sqrt( particles[this_simname][nc1]*particles[this_simname][nc2])




if 1:
    #Probably the good one.
    for ns,this_simname in enumerate(sim_list):
        on = overlap_numb[this_simname]
        ont = on.transpose()
        of = overlap_dict[this_simname]
        oft = overlap_dict[this_simname].transpose()
        #how it works:
        #tool.overlaps[core_id] = list of foreign particles in my contour
        #o[core_id,:] = tool.overlaps[core_id]
        #o[core,:] = other particle in my core contour
        #o[:,core] = other contours with my particles
        both = np.stack([on, ont])
        minmin = both.min(axis=0)

        ofmin = np.stack([of,oft]).min(axis=0)
        of = ofmin


        mask = (on>0.0)*(ont>0.0)
        overlap = ofmin>0

        m1=mask*overlap
        NumberOfOverlap=(m1).sum(axis=1)
        Max = (of*mask).max(axis=1)
        NumberOfOverlap=(m1).sum(axis=1)
        Fraction = (of*mask).sum(axis=1)
        Fraction[NumberOfOverlap>0] /= NumberOfOverlap[ NumberOfOverlap > 0]

        nmax = NumberOfOverlap.max()
        bins_n=np.linspace(-0.5,nmax+0.5)
        bins_f = np.linspace(0,1,20)

        plt.clf()
        print('word')
        figa, axa, axtop,axright = multiplots.three_way_bean(figsize=(4,4), left=0.15, width=0.62, bottom=0.11, height=0.62, histdepth=0.02)
        axa.scatter( NumberOfOverlap, Fraction,c='k')
        axa.scatter( NumberOfOverlap,Max ,c='r')
        axtop.hist( NumberOfOverlap, bins=bins_n, histtype='step',color='k')
        axright.hist( Fraction, bins=bins_f, histtype='step',color='k',orientation='horizontal')
        axright.hist( Max, bins=bins_f, histtype='step',color='r',orientation='horizontal')

        if 1:
            axbonk(axa,xlabel=r'$N_{\rm{overlap}}$', ylabel=r'$\langle O_{i,j}\rangle, \max( O_{i,j})$')
            axa.set_xlim([-0.1,nmax+0.1])
            axbonk(axtop,xlabel='',ylabel=r'$N$')
            axbonk(axright,xlabel=r'$N$',ylabel='')
            axright.set_ylim( axa.get_ylim())
            axtop.set_xlim( axa.get_xlim())
            axtop.set_xticks([])
            axright.set_yticks([])

        figa.savefig('plots_to_sort/%s_overlaps.pdf'%this_simname)
