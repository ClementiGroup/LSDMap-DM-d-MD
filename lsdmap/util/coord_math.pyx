import numpy as np
cimport numpy as np
import cython

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def cmd(np.ndarray[np.double_t, ndim=2] coord1, np.ndarray[np.double_t, ndim=2] coord2):

    cdef int natoms = coord1.shape[1]
    cdef double r1, r2, c1, c2, sumc1, sumc2, sumc1c2
    cdef double r0 = 0.01

    sumc1 = sumc2 = sumc1c2 = 0.0

    for i in range(natoms):
        for j in range(natoms):
            r1 = sqrt((coord1[0,i] - coord1[0,j])**2 + (coord1[1,i] - coord1[1,j])**2 + (coord1[2,i] - coord1[2,j])**2)
            r2 = sqrt((coord2[0,i] - coord2[0,j])**2 + (coord2[1,i] - coord2[1,j])**2 + (coord2[2,i] - coord2[2,j])**2)

            c1 = (1-(r1/r0)**8)/(1-(r1/r0)**12)
            c2 = (1-(r2/r0)**8)/(1-(r2/r0)**12)

            sumc1 += c1
            sumc2 += c2

            if i != j: sumc1c2 += (c1 - c2)**2

    value = sqrt(sumc1c2/sqrt(sumc1*sumc2))

    return value
