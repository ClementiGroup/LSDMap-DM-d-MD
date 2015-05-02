from math import floor
import itertools as it
import numpy as np
cimport numpy as np

from dmaps.critical.bias cimport FEHist

cdef int local_fe_fit_voronoi(double* fit_value, double* fit_gradient, double* data, int ndim, int nneighbors_per_dim, FEHist *feh, int iter, double factor_sig):

    cdef unsigned int idx, jdx 
    cdef int nneighbors = nneighbors_per_dim**ndim
    cdef double distance, sum_kernel0, sum_kernel1

    cdef np.ndarray[np.float64_t,ndim=1] fe_neighbors = np.zeros(nneighbors)
    cdef np.ndarray[np.float64_t,ndim=1] sigma = np.zeros(ndim)
    cdef np.ndarray[np.float64_t,ndim=1] kernel = np.zeros(nneighbors)

    cdef np.ndarray[np.float64_t,ndim=2] neighbors = np.zeros((nneighbors, ndim))
    cdef np.ndarray[np.float64_t,ndim=2] gradient_kernel = np.zeros((nneighbors, ndim))

    get_fe_neighbors(data, ndim, neighbors, fe_neighbors, nneighbors_per_dim, feh, iter)

    for jdx in xrange(ndim):
        sigma[jdx] = factor_sig*feh.steps[iter*ndim+jdx]

    neighbors_idxs = []
    # compute kernel using only non empty bins
    for idx in xrange(nneighbors):
        if fe_neighbors[idx] != 1.0:
            neighbors_idxs.append(idx)
            # compute distance
            distance = 0.0
            for jdx in xrange(ndim):
                distance += (data[jdx]-neighbors[idx, jdx])**2/sigma[jdx]**2
            distance = sqrt(distance)
            kernel[idx] = exp(-distance**2/2)
            for jdx in xrange(ndim):
                gradient_kernel[idx, jdx] = -(data[jdx]-neighbors[idx, jdx])/sigma[jdx]**2 * kernel[idx]

    # compute the free energy and its gradient
    sum_kernel0 = np.sum(kernel[neighbors_idxs])
    sum_kernel1 = np.sum(fe_neighbors[neighbors_idxs]*kernel[neighbors_idxs])
    fit_value[0] += sum_kernel1/sum_kernel0

    for jdx in xrange(ndim):
        fit_gradient[jdx] += (np.sum(fe_neighbors[neighbors_idxs]*gradient_kernel[neighbors_idxs, jdx])*sum_kernel0 -\
                             sum_kernel1*np.sum(gradient_kernel[neighbors_idxs, jdx]))/sum_kernel0**2

    return 0

cdef int get_fe_neighbors(double* data, int ndim, np.ndarray[np.float64_t,ndim=2] neighbors, np.ndarray[np.float64_t,ndim=1] fe_neighbors, int nneighbors_per_dim, FEHist *feh, int iter):

    cdef int mdx, idx, jdx, kdx, index_s
    cdef int nneighbors, shift_nbins, shift_nnebins

    nneighbors = nneighbors_per_dim**ndim

    shift_nbins = 0
    for mdx in range(iter):
        shift_nbins += feh.nbins[mdx]

    shift_nnebins = 0
    for mdx in range(iter):
        shift_nnebins += feh.nnebins[mdx]

    # recompute the index of the bin in which the point is located
    idxs_center = []
    for jdx in xrange(ndim):
        idxs_center.append(int(floor((data[jdx] - feh.bins[shift_nbins*ndim+feh.nbins[iter]*jdx])/feh.steps[iter*ndim+jdx] + 0.5)))

    indices = []
    if nneighbors_per_dim%2 == 0:
        raise ValueError("number of neighbor bins should not be even")
    for jdx in xrange(ndim):
        indices.append([idxs_center[jdx] + kdx for kdx in range(-(nneighbors_per_dim-1)/2, (nneighbors_per_dim+1)/2)])

    idxs_neighbors = []
    for idx, index_bin in enumerate(it.product(*indices)):
        idxs_neighbors.append(index_bin)
        for jdx in xrange(ndim):
            neighbors[idx, jdx] = feh.bins[shift_nbins*ndim+feh.nbins[iter]*jdx]+index_bin[jdx]*feh.steps[iter*ndim+jdx]

    for idx, index_bin in enumerate(idxs_neighbors):
        # compute index (serial number)
        index_s = 0
        fe_neighbors[idx] = 1.0
        for jdx in xrange(ndim):
            index_s += index_bin[jdx]*feh.nbins[iter]**(ndim-jdx-1)
        # if the point is outside the grid, the serial number does not apply
        if any([idx_bin < 0 or idx_bin >= feh.nbins[iter] for idx_bin in index_bin]):
            fe_neighbors[idx] = 1.0
        else:
            for kdx in xrange(feh.nnebins[iter]):
                if index_s == feh.nebins_idxs_s[shift_nnebins+kdx]:
                    fe_neighbors[idx] = feh.values[shift_nnebins+kdx]
                    break
                elif index_s < feh.nebins_idxs_s[shift_nnebins+kdx]:
                    fe_neighbors[idx] = 1.0
                    break
    return 0
