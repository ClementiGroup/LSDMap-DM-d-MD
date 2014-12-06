cimport numpy as np

cdef public struct BiasedMD:
    int natoms # number of atoms
    int nheavy_atoms # number of heavy atoms
    int step # current step number
    int nsteps # number of steps
    int* heavy_atoms_idxs # indices of heavy atoms
    double vbias # value of biased potential
    double* dcs # values of dcs
    float* coord # coordinates
    float* force # force

cdef public struct DMSConfig:
    int isfirst # 1 if first dmaps iter, 0 if not
    int nstride # number of configs to save
    int ndcs # number of dcs that should be considered
    double kT # kT value

cdef public struct Fit:
    int npoints # number of points used to fit
    char* function # function name
    char* metric # metric name
    double* weights # fitting weights
    double* sigma # sigma values
    double* coords # coordinates of points used to fit 

cdef public struct FEHist:
    int nbins # number of bins along each dimension
    int nnebins # number of non-empty bins (overall)
    int* nebins_idxs # non-empty bins indices
    int* nebins_idxs_s # non-empty bins indices (serial)
    double* steps # steps of histogram along each dimension
    double* bins # bins coordinates
    double* values # values of the free energy (nonempty bins)
    double* gradient # values of the free energy gradient (nonempty bins)

cdef extern from "math.h":
    double exp(double x)

cdef int do_biased_force_low_level(int natoms, np.ndarray[np.float64_t,ndim=2] coord, np.ndarray[np.float64_t,ndim=2] force, double* vbias, double* dcs, DMSConfig *dmsc, Fit *ft, FEHist *feh)

cdef int save_data(np.ndarray[np.float64_t,ndim=2] coord, BiasedMD *bs, DMSConfig *dmsc)
