from starter2 import *

#Note: not happy with these peak locations, they seem off.
#Check.

def image_peaks(trackname):
    this_track = track_info.tracks[trackname]
    frame = this_track.target_frame
    simdir = this_track.sim_directory
    outname = 'plots_to_sort/%s_%04d_peak_loc'%(trackname, frame)
    if os.path.exists(outname):
        print("File exists, saving time and exiting. ",outname)
        return 0
    ds = yt.load("%s/DD%04d/data%04d"%(simdir,frame,frame))

    peak_fname = this_track.peak_fname
    fptr = h5py.File(peak_fname,'r')
    peaks = fptr['peaks'][()]
    fptr.close()

    axis = 0
    proj = ds.proj('density',0)
    pw = proj.to_pw()
    for npeak,peak in enumerate(peaks):
        pw.annotate_text(peak,'.')
        pw.annotate_text(peak,r'$%d$'%npeak)
    print("Made image",outname)
    pw.save(outname)
