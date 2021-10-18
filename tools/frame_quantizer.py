from starter2 import *



def frame_quantizer(sim=None,directory=None, dt=None, thresh=1e-5):
    #Assume that the list of outputs happens every dt
    #AND some extra in between from the cycle-based outputs.
    #Assume that the first non-zero time is dt if not given explicitly
    if directory is None:
        directory = dl.sims[sim]

    OutputLog = "%s/OutputLog"%directory
    fptr = open(OutputLog,'r')
    lines = fptr.readlines()
    fptr.close()
    frames=[]
    times=[]
    for line in lines:
        #DATASET WRITTEN ./DD0001/data0001        0 0.0000000000000000e+00
        sp = line.split()
        frames.append(  int(sp[2][-4:]))
        times.append( float(sp[4]))
    frames=nar(frames)
    times=nar(times)
    if dt is None:
        dt = times[ times>0][0]
        print(dt)
    keep = (times % dt) < thresh
    frames_to_keep = frames[keep]
    times_to_keep = times[keep]
    return frames_to_keep,times_to_keep




