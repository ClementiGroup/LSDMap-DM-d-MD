import cython
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)
    double cos(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def cmd(np.ndarray[np.double_t, ndim=2] coord1, np.ndarray[np.double_t, ndim=2] coord2, double r0):

    cdef double value 
    cdef int natoms = coord1.shape[1]
    cdef double r1, r2, c1, c2, sumc1, sumc2, sumc1c2

    sumc1 = sumc2 = sumc1c2 = 0.0

    for i in xrange(natoms):
        for j in xrange(natoms):
            r1 = sqrt((coord1[0,i] - coord1[0,j])**2 + (coord1[1,i] - coord1[1,j])**2 + (coord1[2,i] - coord1[2,j])**2)
            r2 = sqrt((coord2[0,i] - coord2[0,j])**2 + (coord2[1,i] - coord2[1,j])**2 + (coord2[2,i] - coord2[2,j])**2)

            c1 = (1-(r1/r0)**8)/(1-(r1/r0)**12)
            c2 = (1-(r2/r0)**8)/(1-(r2/r0)**12)

            sumc1 += c1
            sumc2 += c2

            if i != j: sumc1c2 += (c1 - c2)**2

    value = sqrt(sumc1c2/sqrt(sumc1*sumc2))

    return value


@cython.boundscheck(False)
@cython.wraparound(False)
def euclidean(np.ndarray[np.double_t, ndim=2] coord1, np.ndarray[np.double_t, ndim=2] coord2):

    cdef double value
    cdef int ndim = coord1.shape[0]
    cdef int natoms = coord1.shape[1]
    cdef double sum

    cdef unsigned int i, j

    sum = 0.0

    for i in xrange(ndim):
        for j in xrange(natoms):
          sum += (coord1[i,j] - coord2[i,j])**2
    value = sqrt(sum)
    return value

@cython.boundscheck(False)
@cython.wraparound(False)
def dihedral(np.ndarray[np.double_t, ndim=1] coord1, np.ndarray[np.double_t, ndim=1] coord2):

    cdef double value
    cdef int ndihedral = coord1.shape[0]
    cdef double sum = 0.0

    for i in xrange(ndihedral):
        sum += 0.5*(1-cos(coord1[i] - coord2[i]))

    value = sqrt(sum/ndihedral)

    return value
