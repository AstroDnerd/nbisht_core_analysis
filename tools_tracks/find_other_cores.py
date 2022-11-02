from starter2 import *

def centroid_in_box(this_looper,left,right,core_list=None, mini_scrubbers=None, thresh=0.05, last_only=False):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)
    all_cores=np.unique(this_looper.tr.core_ids)

    if mini_scrubbers is None:

        mini_scrubbers={}
        for core_id in all_cores:
            print('in box',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[core_id]=ms
    LLL = left+0
    LLL.shape = LLL.size, 1
    RRR = right+0
    RRR.shape = RRR.size, 1
    if hasattr(LLL,'v'):
        LLL = LLL.v
        RRR = RRR.v
    cores_in_box = []
    shifts = []
    for core in core_list:
        ms=mini_scrubbers[core]
        P1 = np.stack( [ms.mean_x, ms.mean_y, ms.mean_z])
        ok = ((P1>=LLL)*(P1<=RRR)).prod(axis=0)
        if last_only:
            keep = ok[-1]
        else:
            keep=ok.any()
        if keep:
            cores_in_box.append(core)
            shifts.append(0)

    #If the box is outside the domain
    new_extras = []
    if (LLL < 0).any() + (RRR > 1.0).any():
        shift = (LLL<0).astype('int') - (RRR>1).astype('int')
        print("SHIFT",shift)
        individuals = []
        for i in range(3):
            if np.abs( shift[i]):
                sss = np.zeros_like(shift)
                sss[i]=shift[i]
                individuals.append(sss)
        import itertools
        for nshifts in range( len(individuals)):
            for subset in itertools.combinations( individuals, nshifts+1):
                to_shift = nar(subset).sum(axis=0)
                print("TO SHIFT", to_shift)

                Ltest = LLL +to_shift
                Rtest = RRR +to_shift
                for core in core_list:
                    ms=mini_scrubbers[core]
                    P1 = np.stack( [ms.mean_x, ms.mean_y, ms.mean_z])
                    ok = ((P1>=Ltest)*(P1<=Rtest)).prod(axis=0)
                    if last_only:
                        keep = ok[-1]
                    else:
                        keep=ok.any()
                    if keep:
                        new_extras.append(core)
                        shifts.append(-to_shift.flatten())

    print('NEW',new_extras)
    return cores_in_box+new_extras, shifts



def get_other_cores(this_looper,core_list=None, mini_scrubbers=None, thresh=0.05):

    if core_list is None:
        core_list = np.unique(this_looper.tr.core_ids)
    all_cores=np.unique(this_looper.tr.core_ids)

    if mini_scrubbers is None:

        mini_scrubbers={}
        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[core_id]=ms

    other_cores={}
    shift={}
    for main_core in core_list:

        x_ext=extents()
        y_ext=extents()
        z_ext=extents()
        other_cores[main_core]=[]
        ms = mini_scrubbers[main_core]
        mmx = ms.mean_x.mean()
        mmy = ms.mean_y.mean()
        mmz = ms.mean_z.mean()
        left =  nar([mmx-thresh,mmy-thresh,mmz-thresh])
        right = nar([mmx+thresh,mmy+thresh,mmz+thresh])
        for ocore in all_cores:
            if ocore==main_core:
                continue
            oms = mini_scrubbers[ocore]
            P1 = np.stack( [ms.mean_x, ms.mean_y, ms.mean_z])
            P2 = np.stack( [oms.mean_x, oms.mean_y, oms.mean_z])
            delta = P1-P2

            ok = ( ((delta)**2).sum(axis=0) < thresh**2).any()
            if ok:
                other_cores[main_core].append(ocore)


            #look for periodic shifts
            #Anything that's further than 1-delta,  shift it back.
            a_delta = np.abs( P1-P2)
            thresh_inv = 1-thresh
            if ( a_delta >= thresh_inv).any():
                if main_core not in shift:
                    shift[main_core]={}
                shifter = delta+0
                lookit= a_delta >= thresh_inv
                shifter[~lookit] = 0
                shifter = np.sign(shifter)
                ok2 = ((( P1-(P2+shifter))**2).sum(axis=0) < thresh**2).any()
                if ok2:
                    other_cores[main_core].append(ocore)
                shift[main_core][ocore] = shifter[:,-1]
    return other_cores, shift

if 0:
    #Testing
    import three_loopers_u500 as TL
    this_looper=TL.loops['u502']
    thtr=this_looper.tr
    all_cores=np.unique(this_looper.tr.core_ids)
    if 'mini_scrubbers' not in dir():

        mini_scrubbers={}

        for core_id in all_cores:
            print('ms',core_id)
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=False)
            ms.particle_pos(core_id)

            mini_scrubbers[core_id]=ms
    core_list=all_cores
    #core_list=[268]
    #core_list=[373]
    other_cores, shift=get_other_cores( TL.loops['u502'], core_list=core_list, mini_scrubbers=mini_scrubbers)
    print(other_cores)

