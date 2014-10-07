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
from lsdmap.util import metric as mt
from lsdmap.util import pyqcprot
from lsdmap.rbf import rbf


# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef public struct BiasedMD:
    int natoms
    int nhatoms # number of heavy atoms
    int step
    int nsteps
    int* hatoms_idxs
    float vbias
    float dc1
    float dc2
    float* coord
    float* force

cdef public struct DMSConfig:
    int isfirst
    int is2nddc
    int nstride
    int nsave
    int nsavedcs
    char* metric
    float kT

cdef public struct Fit:
    int npoints
    char* function
    char* metric
    float* weights1
    float* sigma1
    float* weights2
    float* sigma2
    float* coords

cdef public struct FEHist:
    int nbins
    float delta_dc1
    float delta_dc2
    float* bins
    float* values
    float* binsdc1
    float* binsdc2
    float* graddc1
    float* graddc2

cdef DMSConfig dmsc
cdef BiasedMD bs
cdef Fit ft
cdef FEHist feh

cdef extern from "math.h":
    float exp(float x)
    float sqrt(float x)
    float log(float x)

cdef public DMSConfig* initDMSConfig(const char* file):

    cdef bytes metric
    cdef int idx
    cdef float kb, temperature

    # initialize parser 
    config = ConfigParser.SafeConfigParser()
    config.read(file)

    # first iteration?
    dmsc.isfirst = config.getint('DMAPS', 'isfirst')

    # use 2nd dc?
    dmsc.is2nddc = config.getint('DMAPS', 'is2nddc')

    # number of points saved per replica
    dmsc.nstride = config.getint('DMAPS', 'nstride')

    kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
    # temperature
    temperature = config.getint('MD', 'temperature')
    dmsc.kT = kb*temperature

    return &dmsc

cdef public FEHist* initFEHist(DMSConfig *dmsc, const char* file):

    if dmsc.isfirst == 0:

        binsdc = np.loadtxt('bins_fe.xy')
        graddcs = np.loadtxt('gradient_fe.xy')
        free_energy = np.loadtxt('hist_fe.xyz')
        feh.nbins = int(len(binsdc))

        # allocate memory for the bins along DC1
        feh.binsdc1 = <float *>malloc(feh.nbins*sizeof(float))

        if dmsc.is2nddc == 1:
            # allocate memory for the bins along DC2
            feh.binsdc2 = <float *>malloc(feh.nbins*sizeof(float))

        # allocate memory for the gradient along DC1
        feh.graddc1 = <float *>malloc(feh.nbins*feh.nbins*sizeof(float))
        
        # allocate memory for the gradient along DC2
        if dmsc.is2nddc == 1:
            feh.graddc2 = <float *>malloc(feh.nbins*feh.nbins*sizeof(float))

        # allocate memory for the values of the free energy
        feh.values = <float *>malloc(feh.nbins*feh.nbins*sizeof(float))

        for idx in xrange(feh.nbins):
            if dmsc.is2nddc == 0:
                # update values of the bins
                feh.binsdc1[idx] = float(binsdc[idx])
            else:
                # update values of the bins
                feh.binsdc1[idx] = float(binsdc[idx,0])
                feh.binsdc2[idx] = float(binsdc[idx,1])

        for idx in xrange(feh.nbins*feh.nbins):
            if dmsc.is2nddc == 0:
                # update values of the gradient
                feh.graddc1[idx] = float(graddcs[idx])
            else:
                # update values of the gradient
                feh.graddc1[idx] = float(graddcs[idx,0])
                feh.graddc2[idx] = float(graddcs[idx,1])

            # update values of the free energy
            feh.values[idx] = float(free_energy[idx][2])

        feh.delta_dc1 = feh.binsdc1[1] - feh.binsdc1[0]
        feh.delta_dc2 = feh.binsdc2[1] - feh.binsdc2[0]

    return &feh

cdef public Fit* initFit(DMSConfig *dmsc, const char* file):

    cdef unsigned int idx, jdx, kdx
    cdef int natoms

    if dmsc.isfirst == 0:

        rg = reader.open('fit.gro')
        coordsfit = rg.readlines()
        ft.npoints = coordsfit.shape[0]
        natoms = coordsfit.shape[2]

        weightsfit = np.loadtxt('fit.w')        
        sigmafit = np.loadtxt('fit.sig')

        # allocate memory for the coordinates used for the fit
        ft.coords = <float *>malloc(ft.npoints*3*natoms*sizeof(float))

        # allocate memory for the weights used for the fit (1st DC)
        ft.weights1 = <float *>malloc(ft.npoints*sizeof(float))

        # allocate memory for the values of sigma used for the fit (1st DC)
        ft.sigma1 = <float *>malloc(ft.npoints*sizeof(float))

        if dmsc.is2nddc == 1:
            # allocate memory for the weights used for the fit (2nd DC)
            ft.weights2 = <float *>malloc(ft.npoints*sizeof(float))

            # allocate memory for the values of sigma used for the fit (2nd DC)
            ft.sigma2 = <float *>malloc(ft.npoints*sizeof(float))

        for idx in xrange(ft.npoints):
            # update coordinates
            for jdx in xrange(3):
                for kdx in xrange(natoms):
                    ft.coords[idx + ft.npoints*(jdx+3*kdx)] = float(coordsfit[idx][jdx][kdx]) # Python-like array

            if dmsc.is2nddc == 0:
                # update weights
                ft.weights1[idx] = float(weightsfit[idx])

                # update sigma values
                ft.sigma1[idx] = float(sigmafit[idx])
            else:
                # update weights (1st DC)
                ft.weights1[idx] = float(weightsfit[idx,0])

                # update weights (2nd DC)
                ft.weights2[idx] = float(weightsfit[idx,1])

                # update sigma values (1st DC)
                ft.sigma1[idx] = float(sigmafit[idx,0])

                # update sigma values (2nd DC)
                ft.sigma2[idx] = float(sigmafit[idx,1])

        # initialize parser 
        config = ConfigParser.SafeConfigParser()
        config.read(file)

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
    hatoms = ['C', 'O', 'N', 'S']
    hatoms_idxs = []

    f = open('tmp.gro', 'r')
    f.next()
    natoms = int(f.next())
    for atom_idx, line in it.izip(xrange(natoms), f):
        atoms = line[8:15].lstrip()
        if not 'W' in atoms: # check if atoms are not from water molecules
            # check if atoms are in the list of heavy atoms
            if any([hatom in atoms for hatom in hatoms]):
               hatoms_idxs.append(atom_idx)
    f.close()
    # update number of heavy atoms
    bs.nhatoms = len(hatoms_idxs)
    # allocate idxs of heavy atoms
    bs.hatoms_idxs = <int *>malloc(bs.nhatoms*sizeof(int))
    # update idxs of heavy atoms
    for idx in xrange(bs.nhatoms):
        bs.hatoms_idxs[idx] = int(hatoms_idxs[idx])

    return &bs


cdef public int do_biased_force(BiasedMD *bs, DMSConfig *dmsc, Fit *ft, FEHist *feh) except -1:

    cdef unsigned int idx, jdx
    cdef int natoms = bs.natoms, nhatoms = bs.nhatoms
    cdef int nsave, nsavedcs
    cdef float* vbias = <float *>malloc(sizeof(float))
    cdef float* dc1 = <float *>malloc(sizeof(float))
    cdef float* dc2 = <float *>malloc(sizeof(float))

    cdef np.ndarray[np.float32_t,ndim=2] coord = np.zeros((3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] coord_ha = np.zeros((3, nhatoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] force = np.zeros((3, nhatoms), dtype='f4')

    np.import_array()

    coord = np.copy(c2numpy(bs.coord, 3*natoms).reshape(3, natoms, order='F'))

    # select only heavy atoms of the protein to compute the biased force
    hatoms_idxs = [bs.hatoms_idxs[idx] for idx in xrange(bs.nhatoms)]
    coord_ha  = coord[:,hatoms_idxs]

    # store configuration and weight after nsave steps
    nsave = max(floor(bs.nsteps*1.0/dmsc.nstride), 1)
    if bs.step%nsave==0 and bs.step>0:
        save_data(coord, bs, dmsc)

    # store dcs after nsavedcs steps to compute autocorrelation time 
    nsavedcs = max(floor(nsave*0.1), 1)
    if bs.step%nsavedcs==0 and bs.step>0 and dmsc.isfirst == 0:
        with open('autocorr.ev', 'a') as evfile:
            print >> evfile, '%15.7e %15.7e' %(bs.dc1, bs.dc2) 

    if dmsc.isfirst == 0:
        # compute biased force
        do_biased_force_low_level(nhatoms, coord_ha, force, vbias, dc1, dc2, dmsc, ft, feh)
        # update vbias and force
        bs.vbias = vbias[0]
        bs.dc1 = dc1[0]
        bs.dc2 = dc2[0]
        #print bs.force[0], bs.force[1], bs.force[3], force[0][1]
        for idx in xrange(3):
            for jdx, hatom_jdx in enumerate(hatoms_idxs):
                bs.force[idx+3*hatom_jdx] += float(force[idx][jdx])
        #print bs.force[0], bs.force[1], bs.force[3], force[0][1]
    else:
        bs.vbias = 0.0
    return 0

cdef int do_biased_force_low_level(int natoms, np.ndarray[np.float32_t,ndim=2] coord, np.ndarray[np.float32_t,ndim=2] force, float* vbias, float* dc1, float* dc2, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx
    cdef int idx_bindc1, idx_bindc2
    cdef float dc1_value = 0.0, dc2_value = 0.0, dv_ddc1 = 0.0, dv_ddc2 = 0.0
    cdef float wdf_dc1, wdf_dc2

    cdef fitfunction_type fitfunction, fitderivative 
    cdef np.ndarray[np.float32_t,ndim=2] gradient_metric = np.zeros((3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] gradient_dc1 = np.zeros((3, natoms), dtype='f4')
    cdef np.ndarray[np.float32_t,ndim=2] gradient_dc2 = np.zeros((3, natoms), dtype='f4')

    #cdef double p1, p2
 
    fitfunction = setfitfunction(ft.function)
    fitderivative = setfitderivative(ft.function)

    #p1 = time.time()
    # load coordinates used for the fit
    coordsfit = np.copy(c2numpy(ft.coords, ft.npoints*3*natoms).reshape(ft.npoints, 3, natoms, order='F'))

    # compute cv and gradient cv
    for kdx in xrange(ft.npoints):
        # TODO: make a library with the gradients of each metric
        distance = pyqcprot.CalcRMSDGradient(coord, coordsfit[kdx], gradient_metric, None)
        # update dcs
        dc1_value += ft.weights1[kdx] * fitfunction(distance, ft.sigma1[kdx])
        dc2_value += ft.weights2[kdx] * fitfunction(distance, ft.sigma2[kdx])
        # update gradient dc
        wdf_dc1 = ft.weights1[kdx] * fitderivative(distance, ft.sigma1[kdx])
        wdf_dc2 = ft.weights2[kdx] * fitderivative(distance, ft.sigma2[kdx])
        for jdx in xrange(natoms):
            for idx in xrange(3):
                gradient_dc1[idx,jdx] += wdf_dc1 * gradient_metric[idx,jdx]
                gradient_dc2[idx,jdx] += wdf_dc2 * gradient_metric[idx,jdx]

    #p2 = time.time()
    #print p2 - p1

    # store dc values
    dc1[0] = dc1_value  
    dc2[0] = dc2_value

    # compute bin numbers of current cunfiguration
    idx_bindc1 = int(floor((dc1_value - feh.binsdc1[0])/feh.delta_dc1))
    idx_bindc2 = int(floor((dc2_value - feh.binsdc2[0])/feh.delta_dc2))

    if idx_bindc1 < 0 or idx_bindc1 >= feh.nbins or idx_bindc2 < 0 or idx_bindc2 >= feh.nbins:
        vbias[0] = 0.0
        for idx in xrange(3):
            for jdx in xrange(natoms):
                force[idx,jdx] = 0.0
    else:
        idx_hist = idx_bindc1 * feh.nbins + idx_bindc2
        vbias[0] = -feh.values[idx_hist]
        dv_ddc1 = feh.graddc1[idx_hist]
        dv_ddc2 = feh.graddc2[idx_hist]
        # compute force
        for idx in xrange(3):
            for jdx in xrange(natoms):
                force[idx,jdx] = dv_ddc1*gradient_dc1[idx,jdx] + dv_ddc2*gradient_dc2[idx,jdx]
    return 0

cdef int save_data(np.ndarray[np.float32_t,ndim=2] coord, BiasedMD *bs, DMSConfig *dmsc):

    w = writer.open('.gro', pattern='input.gro')
    w.write(coord, 'confall.gro', mode='a')

    with open('confall.w', 'a') as wfile:
        print >> wfile, '%15.7e' %(exp(bs.vbias/dmsc.kT))

    if dmsc.isfirst == 0:
        with open('confall.ev', 'a') as evfile:
            print >> evfile, '%15.7e %15.7e' %(bs.dc1, bs.dc2)

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

