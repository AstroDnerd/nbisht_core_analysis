"""
Particle operations for Lagrangian Volume

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
cimport numpy as np
cimport cython

def mask_particles(np.ndarray[np.int64_t, ndim=1] ids_to_get,
                   np.ndarray[np.int64_t, ndim=1] p_ids,
                   np.ndarray[np.int32_t, ndim=1] mask_to_get):
    cdef int n1, n2, i1, i2
    n1 = ids_to_get.shape[0]
    n2 = p_ids.shape[0]
    cdef np.int64_t id1, id2
    cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(n2, dtype='int32')
    # Assume unsorted
    cdef int found_any = 0
    for i1 in range(n1):
        if mask_to_get[i1] == 1: continue
        id1 = ids_to_get[i1]
        for i2 in range(n2):
            id2 = p_ids[i2]
            if id1 == id2: 
                mask_to_get[i1] = 1
                mask[i2] = 1
                found_any = 1
                break
    return found_any, mask

def mask_particles_sorted_t4(np.ndarray[np.int64_t, ndim=1] ids_to_get,
                   np.ndarray[np.int64_t, ndim=1] p_ids,
                   np.ndarray[np.int32_t, ndim=1] mask_to_get):
    #Assumes p_ids is sorted.
    print('particles take 4')
    cdef int n1, n2, i1, i2
    n1 = ids_to_get.shape[0]
    n2 = p_ids.shape[0]
    cdef np.int64_t id1, id2
    cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(n2, dtype='int32')
    cdef int found_any = 0
    for i1 in range(n1):
        if mask_to_get[i1] == 1: continue
        id1 = ids_to_get[i1]
        for i2 in range(n2):
            id2 = p_ids[i2]
            if id2 > id1:
                break
            if id1 == id2: 
                mask_to_get[i1] = 1
                mask[i2] = 1
                found_any = 1
                break
    return found_any, mask

def mask_particles_sorted_t7(np.ndarray[np.int64_t, ndim=1] ids_to_get,
                   np.ndarray[np.int64_t, ndim=1] p_ids,
                   np.ndarray[np.int32_t, ndim=1] mask_to_get):
    #Assumes p_ids is sorted.
    print('particles take 4')
    cdef int n1, n2, i1, i2
    n1 = ids_to_get.shape[0]
    n2 = p_ids.shape[0]
    cdef np.int64_t id1, id2, start=0
    cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(n2, dtype='int32')
    cdef int found_any = 0
    for i1 in range(n1):
        if mask_to_get[i1] == 1: continue
        id1 = ids_to_get[i1]
        for i2 in range(start,n2):
            id2 = p_ids[i2]
            if id2 > id1:
                break
            if id1 == id2: 
                mask_to_get[i1] = 1
                mask[i2] = 1
                found_any = 1
                start = i2
                break
    return found_any, mask

