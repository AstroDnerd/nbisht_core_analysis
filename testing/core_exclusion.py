
#run core_overlap on 26 and 27 first.


temp_clump=more_leaf_storage_here['temp_clump']
c1 =   more_leaf_storage_here['stuff']['c1']
c2 =   more_leaf_storage_here['stuff']['c2']
rho1 = more_leaf_storage_here['stuff']['rho1']
rho2 = more_leaf_storage_here['stuff']['rho2']
min1 = more_leaf_storage_here['stuff']['min1']
min2 = more_leaf_storage_here['stuff']['min2']

reload(mountain_top)
which_to_contour= more_leaf_storage_here['stuff']['whicch_to_contour']
this_stack = []
loop_clump = temp_clump
all_excluded=[]
while len(loop_clump.children):
    this_stack.append( mountain_top.leaf_with_center( loop_clump.children))
    loop_clump = this_stack[-1]

point_to_exclude = [c2,c1][which_to_contour]

for nclump,clump in enumerate(this_stack):
    c2_in_new_clump, hull1, point1= mountain_top.is_point_in_hull_kludge( point_to_exclude, clump) 
    print("WAAAA",c2_in_new_clump)

    fig,ax=plt.subplots(1,1)
    ax.scatter( point1[:,1], point1[:,2], c='k', s=1)
    ax.scatter( point_to_exclude[1], point_to_exclude[2],c='g',s=10)
    ax.scatter( c1[1], c1[2],c='r',s=10, marker='*')
    ax.scatter( c2[1], c2[2],c='b',s=10, marker='x')
    ax.scatter( point_to_exclude[1], point_to_exclude[2],c='k',s=15, marker='x')
    fig.savefig('plots_to_sort/exclusion_test%d.png'%nclump)
    plt.close(fig)
