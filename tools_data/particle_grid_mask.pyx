
import numpy as np
cimport numpy as np
cimport cython

def particle_grid_mask_go_i3(xpos ,
                             ypos ,
                             zpos,
                             grid_left,
                             grid_dx,
                             grid_size,
                             grid_mask,
                             grid_selector,
                             particle_selector, #only for this gri
                             particle_exclude   #for all grids
                          ):
    cdef int nparticles, nx, ny, nz
    nparticles = xpos.shape[0]
    cdef np.int64_t i,j,k,n
    #cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(nparticles, dtype='int32')
    cdef int nselect = 0

    for n in range( nparticles ):
        if particle_exclude[n] > 0:
            continue
        i = <int> ( (xpos[n] - grid_left[0])/grid_dx[0] )
        j = <int> ( (ypos[n] - grid_left[1])/grid_dx[1] )
        k = <int> ( (zpos[n] - grid_left[2])/grid_dx[2] )
        if i  >= 0 and i < grid_size[0] and j  >= 0 and j < grid_size[1] and k  >= 0 and k < grid_size[2] and grid_mask[i][j][k] == 1: 
            particle_selector[n] = 1
            particle_exclude[n] = 1
            grid_selector[0][n]=i
            grid_selector[1][n]=j
            grid_selector[2][n]=k

def particle_grid_mask_go_i2(xpos ,
                             ypos ,
                             zpos,
                             grid_left,
                             grid_dx,
                             grid_size,
                             grid_mask,
                             grid_selector,
                             particle_selector
                          ):
    cdef int nparticles, nx, ny, nz
    nparticles = xpos.shape[0]
    cdef np.int64_t i,j,k,n
    #cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(nparticles, dtype='int32')
    cdef int nselect = 0

    for n in range( nparticles ):
        i = <int> ( (xpos[n] - grid_left[0])/grid_dx[0] )
        j = <int> ( (ypos[n] - grid_left[1])/grid_dx[1] )
        k = <int> ( (zpos[n] - grid_left[2])/grid_dx[2] )
        if i  >= 0 and i < grid_size[0] and j  >= 0 and j < grid_size[1] and k  >= 0 and k < grid_size[2] and grid_mask[i][j][k] == 1: 
            particle_selector[n] = 1
            grid_selector[0][n]=i
            grid_selector[1][n]=j
            grid_selector[2][n]=k

def particle_grid_mask_go_i1(xpos ,
                             ypos ,
                             zpos,
                             grid_left,
                             grid_dx,
                             grid_size,
                             grid_mask,
                             grid_selector,
                             particle_selector
                          ):
    print("THIS CODE DOESN'T WORK STOP USING IT.")










#@particle_mask = particle_clump_mask( xpos, ypos, zpos, grid_left, grid_dx, clump_cut_mask)
