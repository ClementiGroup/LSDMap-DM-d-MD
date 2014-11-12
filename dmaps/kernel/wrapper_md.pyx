# import Python modules
import sys
import ConfigParser
import numpy as np
from math import floor
import itertools as it
import time

# import Cython modules
cimport numpy as np
from libc.stdlib cimport malloc
from libc.string cimport strcpy
from cpython cimport PyObject, Py_INCREF

# import lsdmap modules
from lsdmap.rw import reader, writer
from lsdmap.util import pyqcprot

# import dmaps modules
from dmaps.tools.config import known_prms

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef public struct BiasedMD:
    int natoms # number of atoms
    int nheavy_atoms # number of heavy atoms
    int step # current step number
    int nsteps # number of steps
    int* heavy_atoms_idxs # indices of heavy atoms
    float vbias # value of biased potential
    float* dcs # values of dcs
    float* coord # coordinates
    float* force # force

cdef public struct DMSConfig:
    int isfirst # 1 if first dmaps iter, 0 if not
    int nstride # number of configs to save
    int ndcs # number of dcs that should be considered
    float kT # kT value

cdef public struct Fit:
    int npoints # number of points used to fit
    char* function # function name
    char* metric # metric name
    float* weights # fitting weights
    float* sigma # sigma values
    float* coords # coordinates of points used to fit 

cdef public struct FEHist:
    int nbins # number of bins along each dimension
    int nnebins # number of non-empty bins (overall)
    int* nebins_idxs # non-empty bins indices
    int* nebins_idxs_s # non-empty bins indices (serial)
    float* steps # steps of histogram along each dimension
    float* bins # bins coordinates
    float* values # values of the free energy (nonempty bins)
    float* gradient # values of the free energy gradient (nonempty bins)

cdef DMSConfig dmsc
cdef BiasedMD bs
cdef Fit ft
cdef FEHist feh

cdef extern from "math.h":
    float exp(float x)
    float sqrt(float x)
    float log(float x)

cdef check_parameter(value_found, value, prmname, filename):
    """
    Check parameter from defined value, prmname is the name of the parameter as it is specified in known_prms of file dmaps/tools/config.py
    """
    name, section, tag = known_prms[prmname]
    if value_found != value:
        raise IOError("file " + filename + " " + "should contain %i "%value  + tag  + " according to .ini file (" + \
name + " " + "in section " + section + "), %i detected"%value_found)

cdef public DMSConfig* initDMSConfig(const char* file):

    cdef float kb, temperature

    # initialize parser 
    config = ConfigParser.SafeConfigParser()
    config.read(file)

    # first iteration?
    dmsc.isfirst = config.getint('DMAPS', 'isfirst')

    # number of points saved per replica
    dmsc.nstride = config.getint('DMAPS', 'nstride')

    # number of first dcs used
    dmsc.ndcs = config.getint('DMAPS', 'ndcs')

    kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)

    # temperature
    temperature = config.getint('MD', 'temperature')
    dmsc.kT = kb*temperature

    return &dmsc

cdef public FEHist* initFEHist(DMSConfig *dmsc, const char* file):

    if dmsc.isfirst == 0:

        bins = np.loadtxt('bins.xyz', dtype="f4")
        feh.nbins = bins.shape[0]

        hfile = 'hist.dat'
        nebins_idxs = np.genfromtxt(hfile, usecols=tuple(range(dmsc.ndcs)), dtype="i4")
        nebins_idxs_s = np.genfromtxt(hfile, usecols=dmsc.ndcs, dtype="i4")
        free_energy = np.genfromtxt(hfile, usecols=dmsc.ndcs+1, dtype="f4")
        gradient = np.genfromtxt(hfile, usecols=tuple(range(dmsc.ndcs+2,2*dmsc.ndcs+2)), dtype="f4")
        feh.nnebins = free_energy.shape[0]

        if dmsc.ndcs == 1:
            bins = bins[:,np.newaxis]
            nebins_idxs = nebins_idxs[:,np.newaxis]
            gradient = gradient[:,np.newaxis]

        # allocate memory for the bins
        feh.bins = <float *>malloc(feh.nbins*dmsc.ndcs*sizeof(float))

        # allocate memory for the idxs of nonempty bins
        feh.nebins_idxs = <int *>malloc(feh.nnebins*dmsc.ndcs*sizeof(int))

        # allocate memory for the idxs of nonempty bins
        feh.nebins_idxs_s = <int *>malloc(feh.nnebins*sizeof(int))

        # allocate memory for the values of the free energy (nonempty bins)
        feh.values = <float *>malloc(feh.nnebins*sizeof(float))

        # allocate memory for the gradient (nonempy bins)
        feh.gradient = <float *>malloc(feh.nnebins*dmsc.ndcs*sizeof(float))
        
        for idx in xrange(feh.nbins):
            for jdx in xrange(dmsc.ndcs):
                # update values of the bins
                feh.bins[idx+feh.nbins*jdx] = bins[idx,jdx]

        for idx in xrange(feh.nnebins):
            # update values of the free energy
            feh.values[idx] = free_energy[idx]
            feh.nebins_idxs_s[idx] = nebins_idxs_s[idx]

            for jdx in xrange(dmsc.ndcs):
                # update values of the idxs of nonempty bins
                feh.nebins_idxs[idx+feh.nnebins*jdx] = nebins_idxs[idx,jdx]

                # update values of the gradient
                feh.gradient[idx+feh.nnebins*jdx] = gradient[idx,jdx]

        # allocate memory for the step values of the histogram
        feh.steps = <float *>malloc(dmsc.ndcs*sizeof(float))

        for jdx in xrange(dmsc.ndcs):
            feh.steps[jdx] = feh.bins[feh.nbins*jdx+1] - feh.bins[feh.nbins*jdx]

    return &feh

cdef public Fit* initFit(DMSConfig *dmsc, const char* file):

    cdef unsigned int idx, jdx, kdx
    cdef int ndim, natoms 

    if dmsc.isfirst == 0:

        # initialize parser
        config = ConfigParser.SafeConfigParser()
        config.read(file)

        # number of configs used for the fit 
        npoints = config.getint('FITTING', 'npoints')

        # load configs
        grofile = 'fit.gro'
        rg = reader.open(grofile) # reader .gro file
        coordsfit = rg.readlines().astype('f4')
        check_parameter(coordsfit.shape[0], npoints, 'nfit', grofile) # check number of configs
        ndim = coordsfit.shape[1] # number of spatial dimensions
        natoms = coordsfit.shape[2] # number of atoms

        # allocate memory for the coordinates used for the fit
        ft.coords = <float *>malloc(ft.npoints*3*natoms*sizeof(float))

        # load weights
        wfile = 'fit.w'
        weightsfit = np.loadtxt(wfile, dtype="f4")
        if dmsc.ndcs == 1:
            weightsfit = weightsfit[:,np.newaxis]
        check_parameter(weightsfit.shape[0], npoints, 'nfit', wfile)
        check_parameter(weightsfit.shape[1], dmsc.ndcs, 'ndcs',  wfile)

        # load sigma values
        sigfile = 'fit.sig'
        sigmafit = np.loadtxt(sigfile, dtype="f4")
        if dmsc.ndcs == 1:
            sigmafit = sigmafit[:,np.newaxis]
        check_parameter(sigmafit.shape[0], npoints, 'nfit', sigfile)
        check_parameter(sigmafit.shape[1], dmsc.ndcs, 'ndcs', sigfile)

        ft.npoints = npoints
        # allocate memory for the coordinates used for the fit
        ft.coords = <float *>malloc(ft.npoints*3*natoms*sizeof(float))

        # allocate memory for the weights used for the fit
        ft.weights = <float *>malloc(ft.npoints*dmsc.ndcs*sizeof(float))

        # allocate memory for the values of sigma used for the fit
        ft.sigma = <float *>malloc(ft.npoints*dmsc.ndcs*sizeof(float))

        for idx in xrange(ft.npoints):
            # update coordinates
            for jdx in xrange(3):
                for kdx in xrange(natoms):
                    ft.coords[idx + ft.npoints*(jdx+3*kdx)] = coordsfit[idx,jdx,kdx] # Python-like array

            # update weights
            for jdx in xrange(dmsc.ndcs):
                ft.weights[idx+ft.npoints*jdx] = weightsfit[idx,jdx]

            # update sigma values
            for jdx in xrange(dmsc.ndcs):
                ft.sigma[idx+ft.npoints*jdx] = sigmafit[idx,jdx]

        # get function used for the fit
        function = config.get('FITTING','function')
        ft.function = <char *>malloc(sizeof(char)*len(function))
        strcpy(ft.function, function)

        # get metric's name
        metric = config.get('FITTING', 'metric')
        ft.metric = <char *>malloc(len(metric)*sizeof(char))
        strcpy(ft.metric, metric) # need to copy manually since metric is a local variable

    return &ft

cdef public BiasedMD* initBiasedMD(DMSConfig *dmsc, const char* file):

    cdef unsigned int idx
    cdef int natoms

    # list of heavy atoms
    heavy_atoms = ['C', 'N', 'O', 'P', 'S']
    heavy_atoms_idxs = []

    f = open('tmp.gro', 'r')
    f.next()
    natoms = int(f.next())
    for atom_idx, line in it.izip(xrange(natoms), f):
        atoms = line[8:15].lstrip()
        if not 'W' in atoms: # check if atoms are not from water molecules
            # check if atoms are in the list of heavy atoms
            if any([atom in atoms for atom in heavy_atoms]):
               heavy_atoms_idxs.append(atom_idx)
    f.close()

    # update number of heavy atoms
    bs.nheavy_atoms = len(heavy_atoms_idxs)

    # allocate idxs of heavy atoms
    bs.heavy_atoms_idxs = <int *>malloc(bs.nheavy_atoms*sizeof(int))

    # update idxs of heavy atoms
    for idx in xrange(bs.nheavy_atoms):
        bs.heavy_atoms_idxs[idx] = heavy_atoms_idxs[idx]

    # allocate array of dcs (current iteration)
    bs.dcs = <float *>malloc(dmsc.ndcs*sizeof(float))

    return &bs

cdef public int do_biased_force(BiasedMD *bs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx
    cdef int natoms = bs.natoms, nheavy_atoms = bs.nheavy_atoms
    cdef int nsave, nsavedcs
    cdef float* vbias = <float *>malloc(sizeof(float))
    cdef float* dcs = <float *>malloc(dmsc.ndcs*sizeof(float))

    cdef np.ndarray[np.float32_t,ndim=2] coord = np.zeros((3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] coord_heavy_atoms = np.zeros((3, nheavy_atoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] force_heavy_atoms = np.zeros((3, nheavy_atoms), dtype='f4')

    np.import_array()
    coord = np.copy(c2numpy(bs.coord, 3*natoms).reshape(3, natoms, order='F'))

    # select only heavy atoms to compute the biased force
    heavy_atoms_idxs = [bs.heavy_atoms_idxs[idx] for idx in xrange(bs.nheavy_atoms)]
    coord_heavy_atoms  = coord[:,heavy_atoms_idxs]

    # store configuration and weight after nsave steps
    nsave = max(int(floor(bs.nsteps*1.0/dmsc.nstride)), 1)
    if bs.step%nsave == 0 and bs.step > 0:
        save_data(coord, bs, dmsc)

    # store dcs after nsavedcs steps to compute autocorrelation time 
    nsavedcs = max(int(floor(nsave*0.1)), 1)

    if bs.step%nsavedcs == 0 and bs.step > 0 and dmsc.isfirst == 0:
        with open('autocorr.ev', 'a') as evfile:
            print >> evfile, ' '.join(['%15.7e' % (bs.dcs[idx],) for idx in xrange(dmsc.ndcs)])

    # compute biased force if not first iteration
    if dmsc.isfirst == 0:
        do_biased_force_low_level(nheavy_atoms, coord_heavy_atoms, force_heavy_atoms, vbias, dcs, dmsc, ft, feh)
        for jdx in xrange(dmsc.ndcs):
            bs.dcs[jdx] = dcs[jdx]
        bs.vbias = vbias[0]
        #print bs.force[0], bs.force[1], bs.force[3], force[0][1]
        for idx in xrange(3):
            for jdx, atom_jdx in enumerate(heavy_atoms_idxs):
                bs.force[idx+3*atom_jdx] += force_heavy_atoms[idx][jdx]
        #print bs.force[0], bs.force[1], bs.force[3], force[0][1]
    else:
        bs.vbias = 0.0

    return 0

cdef int do_biased_force_low_level(int natoms, np.ndarray[np.float32_t,ndim=2] coord, np.ndarray[np.float32_t,ndim=2] force, float* vbias, float* dcs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx, ldx
    cdef float* dvbias_ddcs = <float *>malloc(dmsc.ndcs*sizeof(float))
    cdef float* wdf = <float *>malloc(dmsc.ndcs*sizeof(float))

    # quantities used fit the dcs
    cdef fitfunction_type fitfunction, fitderivative
    cdef np.ndarray[np.float32_t,ndim=2] gradient_metric = np.zeros((3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=3] coordsfit = np.zeros((ft.npoints, 3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=3] gradient_dcs = np.zeros((dmsc.ndcs, 3, natoms), dtype='f4')

    # quantities used to fit the free energy and its gradient
    cdef int nneighbors = 5**dmsc.ndcs
    cdef float distance, sum_kernel0, sum_kernel1
    cdef np.ndarray[np.float32_t,ndim=1] sigma = np.zeros(dmsc.ndcs, dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=1] fe_neighbors = np.zeros(nneighbors, dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=1] kernel = np.zeros(nneighbors, dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] dcs_neighbors = np.zeros((nneighbors, dmsc.ndcs), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] gradient_kernel = np.zeros((nneighbors, dmsc.ndcs), dtype='f4')

    #cdef double p1, p2
 
    fitfunction = setfitfunction(ft.function)
    fitderivative = setfitderivative(ft.function)

    #p1 = time.time()
    # load coordinates used for the fit
    coordsfit = np.copy(c2numpy(ft.coords, ft.npoints*3*natoms).reshape(ft.npoints, 3, natoms, order='F'))

    for jdx in xrange(dmsc.ndcs):
        dcs[jdx] = 0.0
        dvbias_ddcs[jdx] = 0.0

    # compute cv and gradient cv
    for idx in xrange(ft.npoints):
        # TODO: make a library with the gradients of each metric
        distance = pyqcprot.CalcRMSDGradient(coord, coordsfit[idx], gradient_metric, None)
        for jdx in xrange(dmsc.ndcs):
            # update dcs
            dcs[jdx] += ft.weights[idx+ft.npoints*jdx] * fitfunction(distance, ft.sigma[idx+ft.npoints*jdx])
            # update gradient dc
            wdf[jdx] = ft.weights[idx+ft.npoints*jdx] * fitderivative(distance, ft.sigma[idx+ft.npoints*jdx])
        for ldx in xrange(natoms):
            for kdx in xrange(3):
                for jdx in xrange(dmsc.ndcs):
                    gradient_dcs[jdx, kdx, ldx] += wdf[jdx] * gradient_metric[kdx, ldx]

    get_fe_neighbors(dcs, dmsc.ndcs, dcs_neighbors, fe_neighbors, feh)

    factor_sig = 0.6
    for jdx in xrange(dmsc.ndcs):
        sigma[jdx] = factor_sig*feh.steps[jdx]

    # compute kernel
    for idx in xrange(nneighbors):
       # compute distance
        distance = 0.0
        for jdx in xrange(dmsc.ndcs):
            distance += (dcs[jdx]-dcs_neighbors[idx, jdx])**2/sigma[jdx]**2
        distance = sqrt(distance)
        kernel[idx] = exp(-distance**2/2)
        for jdx in xrange(dmsc.ndcs):
            gradient_kernel[idx, jdx] = -(dcs[jdx]-dcs_neighbors[idx, jdx])/sigma[jdx]**2 * kernel[idx]

    # compute the free energy and its gradient
    sum_kernel0 = np.sum(kernel)
    sum_kernel1 = np.sum(fe_neighbors*kernel)
    vbias[0] = -sum_kernel1/sum_kernel0
    #print vbias[0]
    for jdx in xrange(dmsc.ndcs):
        dvbias_ddcs[jdx] = -(np.sum(fe_neighbors*gradient_kernel[:,jdx])*sum_kernel0 - sum_kernel1*np.sum(gradient_kernel[:,jdx]))/sum_kernel0**2
       #print -dvbias_ddcs[jdx]

    for ldx in xrange(natoms):
        for kdx in xrange(3):
            for jdx in xrange(dmsc.ndcs):
                force[kdx, ldx] -= dvbias_ddcs[jdx]*gradient_dcs[jdx, kdx, ldx]

    #candidates = range(feh.nnebins)
    #for jdx in xrange(dmsc.ndcs):
    #    new_candidates = candidates
    #    candidates = []
    #    # compute bin numbers of current cunfiguration (the "+ 0.5" is because feh.bins corresponds to the centers of the bins, not the left borders)
    #    bin_idx_sample = int(floor((dcs[jdx] - feh.bins[feh.nbins*jdx])/feh.steps[jdx] + 0.5))
    #    for idx in new_candidates:
    #        bin_idx = feh.nebins_idxs[idx+feh.nnebins*jdx]
    #        if bin_idx == bin_idx_sample:
    #            candidates.append(idx)
    #    if not candidates:
    #        break

    #if candidates:
    #    assert len(candidates) == 1, "more than one candidates have been selected from free energy histogram"
    #    vbias[0] = -feh.values[candidates[0]]
    #    for jdx in xrange(dmsc.ndcs):
    #        dvbias_ddcs[jdx] = -feh.gradient[candidates[0]+feh.nnebins*jdx]
    #    # compute force
    #    for ldx in xrange(natoms):
    #        for kdx in xrange(3):
    #            for jdx in xrange(dmsc.ndcs):
    #                force[kdx, ldx] -= dvbias_ddcs[jdx]*gradient_dcs[jdx, kdx, ldx]
    #else:
    #    vbias[0] = 0.0
    #    for ldx in xrange(natoms):
    #        for kdx in xrange(3):
    #            force[kdx, ldx] = 0.0

    return 0

cdef int get_fe_neighbors(float* dcs, int ndcs, np.ndarray[np.float32_t,ndim=2] dcs_neighbors, np.ndarray[np.float32_t,ndim=1] fe_neighbors, FEHist *feh):

    cdef int idx, jdx, kdx, index_s
    cdef float bin_idx
    cdef int nneighbors = fe_neighbors.shape[0]

    idxs_center = []
    for jdx in xrange(ndcs):
        idxs_center.append(int(floor((dcs[jdx] - feh.bins[feh.nbins*jdx])/feh.steps[jdx] + 0.5)))

    #print idxs_center

    indices = []
    for jdx in xrange(ndcs):
        indices.append([idxs_center[jdx] + kdx for kdx in [-2, -1, 0 , 1, 2]])

    idxs_neighbors = []
    for idx, index_bin in enumerate(it.product(*indices)):
        idxs_neighbors.append(index_bin)
        for jdx in xrange(ndcs): 
            dcs_neighbors[idx, jdx] = feh.bins[feh.nbins*jdx]+index_bin[jdx]*feh.steps[jdx]

    #print idxs_neighbors

    for idx, index_bin in enumerate(idxs_neighbors):
        # compute index (serial)
        index_s = 0
        for jdx in xrange(ndcs):
            index_s += index_bin[jdx]*feh.nbins**(ndcs-jdx-1)
        for kdx in xrange(feh.nnebins):
            if index_s == feh.nebins_idxs_s[kdx]:
                fe_neighbors[idx] = feh.values[kdx]
                break
            elif index_s < feh.nebins_idxs_s[kdx]:
                fe_neighbors[idx] = 0.0
                break
    #print fe_neighbors
    return 0

cdef int save_data(np.ndarray[np.float32_t,ndim=2] coord, BiasedMD *bs, DMSConfig *dmsc):

    w = writer.open('.gro', pattern='tmp.gro')
    w.write(coord, 'confall.gro', mode='a')

    with open('confall.w', 'a') as wfile:
        print >> wfile, '%15.7e' %(exp(bs.vbias/dmsc.kT))

    if dmsc.isfirst == 0:
        with open('confall.ev', 'a') as evfile:
            print >> evfile, ' '.join(['%15.7e' % (bs.dcs[idx],) for idx in xrange(dmsc.ndcs)])

    return 0

cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array
 
        This cannot be done in the constructor as it must recieve C-level
        arguments.
 
        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data            
 
        """
        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]

        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_FLOAT, self.data_ptr)
        return ndarray


cdef np.ndarray c2numpy(void *array, int size):
    """
    Python binding of a C array that does 
    not copy the data allocated in C.
    """
    cdef np.ndarray ndarray

    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(size, <void*> array)
    ndarray = np.array(array_wrapper, copy=False)

    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)
 
    return ndarray

ctypedef float (*fitfunction_type)(float, float)

cdef float _h_multiquadric(float r, float sigma):
    return sqrt((1.0/sigma*r)**2 + 1)

cdef float _h_inverse_multiquadric(float r, float sigma):
    return 1.0/sqrt((1.0/sigma*r)**2 + 1)

cdef float _h_thin_plate(float r, float sigma):
    if r == 0.0:
        return 0.0
    else:
        return r**2 * log(r)

cdef float _d_multiquadric(float r, float sigma):
    return r/(sigma**2*sqrt((r/sigma)**2 + 1))

cdef float _d_inverse_multiquadric(float r, float sigma):
    return r/(sigma**2*((r/sigma)**2 + 1)**(3./2))

cdef float _d_thin_plate(float r, float sigma):
    return r * (1 + 2 *log(r))

cdef fitfunction_type setfitfunction(bytes name):

    name = name.lower()
    if name  == "multiquadric":
        funct = _h_multiquadric
    elif name  == "inverse_multiquadric":
        funct =  _h_inverse_multiquadric
    elif name  == "thin_plate":
        funct =  _h_thin_plate
    else:
        raise ValueError("function " + name + " is not supported for fitting")

    return funct

cdef fitfunction_type setfitderivative(bytes name):

    name = name.lower()
    if name  == "multiquadric":
        deriv = _d_multiquadric
    elif name  == "inverse_multiquadric":
        deriv = _d_inverse_multiquadric
    elif name  == "thin_plate":
        deriv = _d_thin_plate
    else:
        raise ValueError("function " + name + " is not supported for fitting")

    return deriv

