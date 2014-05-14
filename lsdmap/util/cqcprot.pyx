import numpy as np
cimport numpy as np

cdef extern from "qcprot.h": 
    double calcrmsdrotationalmatrix(const int *plen, double *x1, double *x2, double *x3, double *y1, double *y2, double *y3, double *rot, const double *weight)

def rmsd(np.ndarray[np.double_t,ndim=1] coord1, np.ndarray[np.double_t,ndim=1] coord2):

    cdef int natoms = len(coord1)/3

    cdef np.ndarray[np.double_t,ndim=1] coord1x = coord1[:natoms]
    cdef np.ndarray[np.double_t,ndim=1] coord1y = coord1[natoms:2*natoms]
    cdef np.ndarray[np.double_t,ndim=1] coord1z = coord1[2*natoms:]
    cdef np.ndarray[np.double_t,ndim=1] coord2x = coord2[:natoms]
    cdef np.ndarray[np.double_t,ndim=1] coord2y = coord2[natoms:2*natoms]
    cdef np.ndarray[np.double_t,ndim=1] coord2z = coord2[2*natoms:]

    cdef np.ndarray[np.double_t,ndim=1] rot = np.zeros(natoms)
    cdef np.ndarray[np.double_t,ndim=1] weight = np.ones(natoms)

    cdef double value

    value = calcrmsdrotationalmatrix( &natoms, <double*> coord1x.data, <double*> coord1y.data, <double*> coord1z.data, 
        <double*> coord2x.data, <double*> coord2y.data, <double*> coord2z.data, <double*> rot.data, <double*> weight.data  )

    if np.isnan(value): # sometimes rmsd computation with the same point returns nan
        value = 0.0

    return value
