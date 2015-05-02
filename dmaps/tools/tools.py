import sys
import random
import copy
import math
import numpy as np
import scipy.ndimage as ndimage 
from scipy.sparse import spdiags

def smooth2a(arrayin, nr, nc):

    # Building matrices that will compute running sums.  The left-matrix, eL,
    # smooths along the rows.  The right-matrix, eR, smooths along the
    # columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
    # (2*Nc+1) rectangle centered on element "i".

    row = arrayin.shape[0]
    col = arrayin.shape[1]

    el = spdiags(np.ones((2*nr+1, row)),range(-nr,nr+1), row, row).todense()
    er = spdiags(np.ones((2*nc+1, col)), range(-nc,nc+1), col, col).todense()

    # Setting all "nan" elements of "arrayin" to zero so that these will not
    # affect the summation.  (If this isn't done, any sum that includes a nan
    # will also become nan.)

    a = np.isnan(arrayin)
    arrayin[a] = 0.

    # For each element, we have to count how many non-nan elements went into
    # the sums.  This is so we can divide by that number to get a mean.  We use
    # the same matrices to do this (ie, "el" and "er").

    nrmlize = el.dot((~a).dot(er))
    nrmlize[a] = None

    # Actually taking the mean.

    arrayout = el.dot(arrayin.dot(er))
    arrayout = arrayout/nrmlize

    return arrayout


def do_grid(data, nbins, mins=None, maxs=None, nextrabins=0):
    """
    Performs multi-dimensional histogram (grid containing the indices of each data instead of the local density)
        data should be an array whose columns correspond to the several dimensions
    """

    npoints = data.shape[0] 
    if len(data.shape) == 1:
        ndim = 1
        data = data[:,np.newaxis]
    else:
        ndim = data.shape[1]

    if mins is None or maxs is None:
        if mins is None and maxs is None:
            window = 1e-8 # normally the window should depend on the data
            mins = np.amin(data, axis=0) - window
            maxs = np.amax(data, axis=0) + window
            steps = (maxs - mins)/nbins
            mins = mins - nextrabins*steps
            nbins = nbins + 2*nextrabins
        else:
            raise ValueError('mins and maxs should be both NaN or non NaN!')
    else:
        if nextrabins != 0:
            raise ValueError('nextrabins should set be to 0 when using mins option')
        steps = (maxs - mins)/nbins

    bins = mins + steps/2  + np.array(range(nbins))[:,np.newaxis]*steps[np.newaxis]

    grid = []
    for dim in xrange(ndim):
        grid = [copy.deepcopy(grid) for idx in xrange(nbins)]

    idxs = np.floor((data-mins)/steps[np.newaxis]).astype(int)

    for num, line in enumerate(idxs):
        bin = grid[line[0]]
        if ndim > 1:
            for idx in line[1:]:
                bin = bin[idx]
        bin.append(num)

    return bins, grid

def do_grid_optimized(data, nnebins, nbins_min, nbins_max, mins=None, maxs=None, niters=10, nextrabins=0):
    """
    Performs multi-dimensional histogram (grid containing the indices of each data instead of the local density)
    data should be an array whose columns correspond to the several dimensions, the grid is construct in such
    a way that the number of empty bins is closed to nnebins
    """

    # look for the number of non empty bins using "nbins_min" bins
    bins, grid = do_grid(data, nbins_min, mins, maxs, nextrabins=nextrabins) 
    nnebins_min = len([bin for idxs, bin in nonempty_bins(grid)])
    if nnebins <= nnebins_min:
        return bins, grid

    # look for the number of non empty bins using "nbins_max" bins
    bins, grid = do_grid(data, nbins_max, mins, maxs, nextrabins=nextrabins)
    nnebins_max = len([bin for idxs, bin in nonempty_bins(grid)])
    if nnebins >= nnebins_max:
        return bins, grid

    for idx in xrange(niters):
        nbins_test = (nbins_min + nbins_max)/2
        bins, grid = do_grid(data, nbins_test, mins, maxs, nextrabins=nextrabins)
        nnebins_test = len([bin for idxs, bin in nonempty_bins(grid)])
        if nbins_test in [nbins_min, nbins_max]:
            break
        elif nnebins_test < nnebins:
            nbins_min = nbins_test
        elif nnebins_test > nnebins:
            nbins_max = nbins_test
        else:
            break

    return bins, grid

def compute_free_energy(grid, ndim, weights, cutoff, kT):
    """
    Give the free energy grid from a grid computed using do_grid function
    """
    # get number of bins
    nbins = len(grid)
    free_energy_grid = np.zeros((nbins,)*ndim, dtype='float')
    free_energy_grid.fill(np.nan)

    for idxs, bin in nonempty_bins(grid):
        weight = 0
        for grid_idx in bin:
            weight += weights[grid_idx]
        free_energy_grid[[np.array(idx) for idx in idxs]] = -kT*np.log(weight)

    # smooth the data
    if ndim == 2:
        free_energy_grid = smooth2a(free_energy_grid, 2, 2)

    free_energy_grid = np.copy(free_energy_grid) # without this line it fails

    # rescale so that the maximum value is 0
    free_energy_grid -= np.nanmax(free_energy_grid)
    # rescale if the minimum is < than - cutoff
    min_free_energy_grid = np.nanmin(free_energy_grid)
    if min_free_energy_grid < -cutoff:
        free_energy_grid -= min_free_energy_grid + cutoff
        with np.errstate(invalid='ignore'):
            free_energy_grid[free_energy_grid > 0.0] = np.nan

    return free_energy_grid

def nonempty_bins(grid):
    """
    the first element of the generator is a list containing the coordinates of the non-empty bin in the grid
    whereas the second element is a list with the data located in the bin, they are given by their numbering in "data" of do_grid
    """
    if not grid:
        pass
    elif isinstance(grid[0], list):
        for idx, value in enumerate(grid):
            for subidxs, subvalue in nonempty_bins(value):
                subidxs = [idx] + subidxs 
                yield subidxs, subvalue
    else: 
        yield [], grid

def pick_points_from_grid(grid, npoints):
    """
    Uniformly select "npoints" indices from a grid generated using do_grid function
    """

    idxs_picked_points = []

    # construct list with nonempty bins to save time when drawing if many bins are empty
    nebins = []

    for idxs, bin in nonempty_bins(grid):
        nebins.append(bin)

    nsamples = sum(len(bin) for bin in nebins) 
    if npoints > nsamples:
        raise ValueError("Too few samples in the grid to select %i of them, %i samples detected" %(npoints, nsamples))

    npicked = 0
    while npicked < npoints:
        idx_point = random.choice(random.choice(nebins))
        if idx_point not in idxs_picked_points:
            idxs_picked_points.append(idx_point)
            npicked = npicked + 1

    return idxs_picked_points

def pick_points_optimized(data, npoints, idxs_preselect=None):
    """
    Select "npoints" points from data so as to maximize the euclidean distance between them 
    """

    if idxs_preselect is None:
        nsamples = data.shape[0]
        idxs_preselect = range(nsamples)
    else:
        nsamples = len(idxs_preselect)

    if len(data.shape) == 1:
        ndim = 1
        data = data[:,np.newaxis]
    else:
        ndim = data.shape[1]

    idxs_picked_points = []
    random_index = random.randrange(0, nsamples)
    idxs_picked_points.append(idxs_preselect.pop(random_index))
    for count in xrange(npoints-1):
        max_min_r = 0.
        for i, idx in enumerate(idxs_preselect):
            min_r = 1.e100
            for kdx in idxs_picked_points:
                r = 0
                for jdx in range(ndim):
                    r += (data[idx,jdx] - data[kdx,jdx])**2
                min_r = min(r, min_r)
            if min_r >= max_min_r:
                max_min_r = min_r
                new_i = i
                new_idx = idx
        del idxs_preselect[new_i]
        idxs_picked_points.append(new_idx)

    return idxs_picked_points
