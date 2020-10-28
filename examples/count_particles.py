

def count_particles(this_looper,fname='n_particles.txt'):
    fptr = open(fname,'w')
    for core in this_looper.target_indices:
        fptr.write("%d %d\n"%(core, len( this_looper.target_indices[core])))
    fptr.close()
#count_particles(this_looper,'u11_n_particles.txt')
