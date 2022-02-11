from starter2 import *
import looper2
reload(looper2)
import convex_hull_tools as CHT
reload(CHT)
import hair_dryer
reload(hair_dryer)
import close_tool
reload(close_tool)

loops={}
sim_list=['u501']
loops['u501'] = looper2.load_looper('u501_short_grav.h5')

if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()
if 'ct' not in dir():
    ct = {}
    for this_simname in sim_list:
        ct[this_simname] = close_tool.close_tool( loops[this_simname])
        ct[this_simname].make_distance()

import supersets
reload(supersets)
if 'st' not in dir():
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( loops[this_simname], ht[this_simname])
        st[this_simname].find()

if 1:
    #every pair
    import means_etc
    for this_simname in sim_list:

        htool = ht[this_simname]
        ctool = ct[this_simname]
        stool = st[this_simname]

        hair = hair_dryer.hair_tool( htool.this_looper)
        overlap_matrix = np.zeros( [len(htool.cores_used)]*2) -1
        for nc1,core_id_1 in enumerate(htool.cores_used):
            for nc2,core_id_2 in enumerate(htool.cores_used):
                #if nc2 <= nc1:
                #    continue
                #overlap_matrix[nc1,nc2] = htool.overlaps[core_id_1][nc2]


                if htool.overlaps[core_id_1][nc2] == 0 and htool.overlaps[core_id_2][nc1] > 0:
                    #core_1=core_id_1
                    #core_2=core_id_2
                    #self=htool
#
#                    hull_1 =  self.hulls[core_1]
#                    hull_2 =  self.hulls[core_2]
#                    vert_1 = self.points_3d[core_1][hull_1.vertices,:]
#                    vert_2 = self.points_3d[core_2][hull_2.vertices,:]
#                    points_1 = self.points_3d[core_1]
#                    points_2 = self.points_3d[core_2]
#
#                    in_1_2_3 = CHT.in_hull(points_1, points_2)
#                    fraction =  in_1_2_3.sum()/points_1.shape[0]


                    cdict={core_id_1:'r',core_id_2:'g'}
                    CHT.plot_watershed(htool, core_list=[core_id_1,core_id_2],frames=[0],accumulate=True, 
                                prefix='c%04d_c%04d'%(core_id_1,core_id_2), axis_to_plot=[-1],
                                color_dict=cdict,label_cores=[-1])
                    #hair.run( core_list=[core_id_1,core_id_2], newplots=False, colors=cdict, name = "c%04d_c%04d_"%(core_id_1,core_id_2))
                    #break
            else:
                continue
            break




fig,ax=plt.subplots(1,1)
p1 = overlap_matrix.flatten()
p2 = overlap_matrix.transpose().flatten()
ok = (p1>-1)*(p2>-1)
#ok = slice(None)
p1=p1[ok]
p2=p2[ok]
ax.scatter( p1,p2)
ax.set_xscale('symlog',linthresh=1e-2)
ax.set_yscale('symlog',linthresh=1e-2)
fig.savefig('plots_to_sort/test.png')
