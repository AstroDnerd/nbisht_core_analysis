
I made a bunch of changes.
% pip install astropy
First and formost, there were two major bugs.
    Bug 1.) Some particles were impossible to get out of some datasets.
    Bug 2.) Some times the track manager would get the particles out of order.

    For Bug 1, it is necessary to pass "bad_particle_list" to "get_target_indices."

Second, I made a way to read and write loopers.
    1.) Make an empty looper like this
        this_looper=looper.core_looper(directory=dl.enzo_directory)
    2.) Import tracks_read_write, and read
        import tracks_read_write as trw
        trw.load_loop(this_looper,"whatever1.h5")
        trw.load_loop(this_looper,"whatever2.h5")
    3.) You should be able to read several files, but watch that you don't
        out-of-memory

Third, one needs to sort the times for the track manager.  
        asort =  np.argsort(thtr.times)
        tsorted = thtr.times[asort]


Look at vel_hist.py for more recent usage.

New tools you'll find useful:
data_puller.py gets particle info and writes hdf5 files.
vel_hist.py shows usage of the reader and making plots with particle data
