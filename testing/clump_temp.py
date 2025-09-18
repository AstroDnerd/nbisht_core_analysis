temp code. 
Kill later.
Unless its useful.
def find_biggest_contour(this_simname, peak_1, peak_2,temp_clump, directory = None, frame = None, peak_fname=None, work_radius=1e-2, do_projections=False, top_to_bottom=3./4, kludge={},
                    verify=None, leaf_storage=None):
    if directory is None:
        directory = dl.sims[this_simname]
    if frame is None:
        frame = dl.target_frames[this_simname]
    if peak_fname is None:
        peak_fname = dl.peak_list[this_simname]
    ds_name="%s/DD%04d/data%04d"%(directory,frame,frame)
    print("OPENING %s and %s"%(ds_name, peak_fname))
    ds = yt.load(ds_name)
    fptr = h5py.File(peak_fname,'r')
    output={}
    try:
        peaks = fptr['peaks'][()]
    except:
        raise
    finally:
        fptr.close()
    #peaks = nar([ np.array([0.89331055, 0.1159668 , 0.4440918 ])])
    radius = ds.arr(1e-2,'code_length')
    c1 =     ds.arr(peaks[peak_1],'code_length')
    c2 =     ds.arr(peaks[peak_2],'code_length')
    print("C1",c1)
    print("C2",c2)
    sph1 = ds.sphere(c1,radius)
    sph2 = ds.sphere(c2,radius)
    rho1 = get_density(c1, sph1).v
    rho2 = get_density(c2, sph2).v
    min1 = rho1**(top_to_bottom)
    min2 = rho2**(top_to_bottom)
    
    if 1:
        #pw.annotate_clumps([leaf1,leaf2])
        this_stack = []
        loop_clump = temp_clump
        while len(loop_clump.children):
            this_stack.append( leaf_with_center( loop_clump.children))
            loop_clump = this_stack[-1]
        new_clump=None
        #Find the first clump that rejects the other point
        all_excluded=[]
        for clump in this_stack:
            c2_in_new_clump= is_point_in_hull( c2, clump) 
            if not c2_in_new_clump:
                all_excluded.append( clump )
        new_clump = all_excluded[0]
        leaf_storage['this_stack']=this_stack
        #in_the_stack = [ is_point_in_hull(c2, lf, tolerance=1./2048) for lf in this_stack]
        #print("IN", in_the_stack)
        if new_clump is None:
            print("Cannot find a new clump.  Reject one of them.")
            output[peak_1]=-1
            output[peak_2]=-1

        else:
            output[peak_1]=new_clump['density'].min()
            output[peak_2]=output[peak_1]


    if do_projections:
        proj = ds.proj('density',1,center=c1,data_source=sph1)
        pw = proj.to_pw()
        pw.set_cmap('density','Greys')
        pw.set_axes_unit('code_length')
        pw.zoom(0.5/radius.v)

        #pw.annotate_clumps([master_clump]+master_clump.leaves)
        p_size=1
        #pw.annotate_these_particles2(1.0,col='r',positions= leaf1['particle_position'], p_size=p_size)
        #pw.annotate_these_particles2(1.0,col='g',positions= leaf2['particle_position'], p_size=p_size)
        #pw.annotate_clumps([leaf1,leaf2])
        #pw.annotate_text(c1,'*')
        #pw.annotate_text(c2,'x')
        #pw.save('plots_to_sort/%s_peak_split_c%04d_c%04d'%(this_simname,peak_1,peak_2))

        pw.annotate_clear()
        pw.annotate_clumps(this_stack[::-1])
        pw.annotate_text(c1,'*')
        pw.annotate_text(c2,'x')
        pw.save('plots_to_sort/%s_TMP_peak_split_corestack_c%04d_c%04d'%(this_simname,peak_1,peak_2))
        pw.annotate_clear()
        pw.annotate_clumps(all_excluded[::-1])
        pw.annotate_text(c1,'*')
        pw.annotate_text(c2,'x')
        pw.save('plots_to_sort/%s_TMP_peak_split_all_excluded_c%04d_c%04d'%(this_simname,peak_1,peak_2))
        pw.annotate_clear()
        pw.annotate_clumps([new_clump])
        pw.annotate_text(c1,'*')
        pw.annotate_text(c2,'x')
        pw.save('plots_to_sort/%s_TMP_peak_split_newclump_c%04d_c%04d'%(this_simname,peak_1,peak_2))
        #            proj = ds.proj(field,ax,center=center, data_source = sph) 
    output['all_excluded']=all_excluded
    return output
