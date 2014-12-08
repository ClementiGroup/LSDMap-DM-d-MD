from dmaps.kernel.bias cimport FEHist
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double x)
    double exp(double x)

cdef int local_fe_fit_voronoi(double* fit_value, double* fit_gradient, double* data, int ndim, int nneighbors_per_dim, FEHist *feh, double factor_sig)

cdef int get_fe_neighbors(double* data, int ndim, np.ndarray[np.float64_t,ndim=2] neighbors, np.ndarray[np.float64_t,ndim=1] fe_neighbors, int nneighbors_per_dim, FEHist *feh)
