import ConfigParser
from math import floor
import itertools as it
import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc
from libc.string cimport strcpy

from lsdmap.rw import reader, writer
from lsdmap.util import pyqcprot

from dmaps.tools.config import known_prms
from dmaps.tools.rbf cimport *
from dmaps.tools.voronoi cimport *

np.import_array()

cdef DMSConfig dmsc
cdef BiasedMD bs
cdef Fit ft
cdef FEHist feh

cdef check_parameter(value_found, value, prmname, filename):
    """
    Check parameter from defined value, prmname is the name of the parameter as it is specified in known_prms of file dmaps/tools/config.py
    """
    name, section, tag = known_prms[prmname]
    if value_found != value:
        raise IOError("file " + filename + " " + "should contain %i "%value  + tag  + " according to .ini file (" + \
name + " " + "in section " + section + "), %i detected"%value_found)


cdef public DMSConfig* initDMSConfig(const char* file):

    cdef double kb, temperature

    # initialize parser 
    config = ConfigParser.SafeConfigParser()
    config.read(file)

    # first iteration?
    dmsc.isfirst = config.getint('DMAPS', 'isfirst')
    # number of points saved per replica
    dmsc.nstride = config.getint('DMAPS', 'nstride')
    # number of first dcs used
    dmsc.ndcs = config.getint('DMAPS', 'ndcs')
    # temperature
    kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
    temperature = config.getint('MD', 'temperature')
    dmsc.kT = kb*temperature

    return &dmsc

cdef public FEHist* initFEHist(DMSConfig *dmsc, const char* file):

    if dmsc.isfirst == 0:

        bins = np.loadtxt('bins.xyz')
        feh.nbins = bins.shape[0]

        hfile = 'hist.dat'
        nebins_idxs = np.genfromtxt(hfile, usecols=tuple(range(dmsc.ndcs)), dtype="i4")
        nebins_idxs_s = np.genfromtxt(hfile, usecols=dmsc.ndcs, dtype="i4")

        free_energy = np.genfromtxt(hfile, usecols=dmsc.ndcs+1)
        gradient = np.genfromtxt(hfile, usecols=tuple(range(dmsc.ndcs+2,2*dmsc.ndcs+2)))

        feh.nnebins = free_energy.shape[0]

        if dmsc.ndcs == 1:
            bins = bins[:,np.newaxis]
            nebins_idxs = nebins_idxs[:,np.newaxis]
            gradient = gradient[:,np.newaxis]

        if np.any(free_energy > 0.0):
            raise IOError("free energy values in hist.dat should be all negative!")

        # allocate memory for the bins
        feh.bins = <double *>malloc(feh.nbins*dmsc.ndcs*sizeof(double))
        # allocate memory for the idxs of nonempty bins
        feh.nebins_idxs = <int *>malloc(feh.nnebins*dmsc.ndcs*sizeof(int))
        # allocate memory for the idxs of nonempty bins
        feh.nebins_idxs_s = <int *>malloc(feh.nnebins*sizeof(int))

        # allocate memory for the values of the free energy (nonempty bins)
        feh.values = <double *>malloc(feh.nnebins*sizeof(double))
        # allocate memory for the gradient (nonempy bins)
        feh.gradient = <double *>malloc(feh.nnebins*dmsc.ndcs*sizeof(double))
        
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
        feh.steps = <double *>malloc(dmsc.ndcs*sizeof(double))

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
        coordsfit = rg.readlines()

        check_parameter(coordsfit.shape[0], npoints, 'nfit', grofile) # check number of configs
        ndim = coordsfit.shape[1] # number of spatial dimensions
        natoms = coordsfit.shape[2] # number of atoms

        # load weights
        wfile = 'fit.w'
        weightsfit = np.loadtxt(wfile, dtype="f8")
        if dmsc.ndcs == 1:
            weightsfit = weightsfit[:,np.newaxis]

        check_parameter(weightsfit.shape[0], npoints, 'nfit', wfile)
        check_parameter(weightsfit.shape[1], dmsc.ndcs, 'ndcs',  wfile)

        # load sigma values
        sigfile = 'fit.sig'
        sigmafit = np.loadtxt(sigfile, dtype="f8")
        if dmsc.ndcs == 1:
            sigmafit = sigmafit[:,np.newaxis]

        check_parameter(sigmafit.shape[0], npoints, 'nfit', sigfile)
        check_parameter(sigmafit.shape[1], dmsc.ndcs, 'ndcs', sigfile)

        ft.npoints = npoints

        # allocate memory for the coordinates used for the fit
        ft.coords = <double *>malloc(ft.npoints*3*natoms*sizeof(double))
        # allocate memory for the weights used for the fit
        ft.weights = <double *>malloc(ft.npoints*dmsc.ndcs*sizeof(double))
        # allocate memory for the values of sigma used for the fit
        ft.sigma = <double *>malloc(ft.npoints*dmsc.ndcs*sizeof(double))

        for idx in xrange(ft.npoints):
            # update coordinates
            for jdx in xrange(3):
                for kdx in xrange(natoms):
                    ft.coords[idx + ft.npoints*(jdx+3*kdx)] = coordsfit[idx, jdx, kdx] # Python-like array

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
    bs.dcs = <double *>malloc(dmsc.ndcs*sizeof(double))

    return &bs

cdef public int do_biased_force(BiasedMD *bs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx
    cdef int natoms = bs.natoms
    cdef int nheavy_atoms = bs.nheavy_atoms
    cdef int nsave, nsavedcs

    cdef double* vbias = <double *>malloc(sizeof(double))
    cdef double* dcs = <double *>malloc(dmsc.ndcs*sizeof(double))

    cdef np.ndarray[np.float64_t,ndim=2] coord = np.zeros((3, natoms))
    cdef np.ndarray[np.float64_t,ndim=2] coord_heavy_atoms = np.zeros((3, nheavy_atoms))
    cdef np.ndarray[np.float64_t,ndim=2] force_heavy_atoms = np.zeros((3, nheavy_atoms))

    np.import_array()
    for jdx in xrange(3):
        for kdx in xrange(natoms):
            coord[jdx,kdx] = float(bs.coord[jdx+3*kdx])

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
            print >> evfile, ' '.join(['%.18e' % (bs.dcs[idx],) for idx in xrange(dmsc.ndcs)])

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


cdef int do_biased_force_low_level(int natoms, np.ndarray[np.float64_t,ndim=2] coord, np.ndarray[np.float64_t,ndim=2] force, 
                                   double* vbias, double* dcs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx, ldx

    cdef int isempty, bin_idx_s, num_line

    cdef double* fe_value = <double *>malloc(sizeof(double))
    cdef double* fe_gradient = <double *>malloc(dmsc.ndcs*sizeof(double))
    cdef double* wdf = <double *>malloc(dmsc.ndcs*sizeof(double))

    cdef fitfunction_type fitfunction, fitderivative

    cdef np.ndarray[np.float64_t,ndim=1] steps = np.zeros(dmsc.ndcs)
    cdef np.ndarray[np.int32_t,ndim=2] nebins_idxs = np.zeros((feh.nnebins,dmsc.ndcs), dtype='i4')

    cdef np.ndarray[np.float64_t,ndim=2] gradient_metric = np.zeros((3, natoms))
    cdef np.ndarray[np.float64_t,ndim=3] coordsfit = np.zeros((ft.npoints, 3, natoms))
    cdef np.ndarray[np.float64_t,ndim=3] gradient_dcs = np.zeros((dmsc.ndcs, 3, natoms))

    fitfunction = setfitfunction(ft.function)
    fitderivative = setfitderivative(ft.function)

    # load coordinates used for the fit
    for idx in xrange(ft.npoints):
        for jdx in xrange(3):
            for kdx in xrange(natoms):
                coordsfit[idx, jdx, kdx] = ft.coords[idx + ft.npoints*(jdx+3*kdx)]

    for jdx in xrange(dmsc.ndcs):
        dcs[jdx] = 0.0

    # compute DC values and their gradient
    for idx in xrange(ft.npoints):
        distance = pyqcprot.CalcRMSDGradient(coord, coordsfit[idx], gradient_metric, None)
        for jdx in xrange(dmsc.ndcs):
            # update DC's
            dcs[jdx] += ft.weights[idx+ft.npoints*jdx] * fitfunction(distance, ft.sigma[idx+ft.npoints*jdx])
            # update gradient of the DC's
            wdf[jdx] = ft.weights[idx+ft.npoints*jdx] * fitderivative(distance, ft.sigma[idx+ft.npoints*jdx])
        for ldx in xrange(natoms):
            for kdx in xrange(3):
                for jdx in xrange(dmsc.ndcs):
                    gradient_dcs[jdx, kdx, ldx] += wdf[jdx] * gradient_metric[kdx, ldx]

    # estimate the free energy of the current configuration and its derivative
    bin_idx_s = 0 # bin serial number
    bin_idxs = [] # bin indices
    isempty = 1 # by default the bin is empty

    # compute the index of the bin within which the point is located
    for jdx in xrange(dmsc.ndcs):
        index = int(floor((dcs[jdx] - feh.bins[feh.nbins*jdx])/feh.steps[jdx] + 0.5))
        bin_idxs.append(index)
        bin_idx_s += index*feh.nbins**(dmsc.ndcs-jdx-1)

    if any([idx < 0 or idx >= feh.nbins for idx in bin_idxs]):
        # if the point is outside the grid, the serial number does not apply
        isempty = 1
    else:
        # check if the bin is empty or not
        for idx in xrange(feh.nnebins):
            if bin_idx_s == feh.nebins_idxs_s[idx]: # non-empty bin found
                num_line = idx
                isempty = 0
                break
            elif bin_idx_s < feh.nebins_idxs_s[idx]: # the bin is empty
                isempty = 1
                break

    # estimate the value of the free energy and its gradient
    if isempty == 0:
        for jdx in xrange(dmsc.ndcs):
            fe_value[0] = 0.0
            fe_gradient[jdx] = 0.0
        # interpolate the free energy locally
        local_fe_fit_voronoi(fe_value, fe_gradient, dcs, dmsc.ndcs, 5, feh, 0.6)
        # update vbias and the biased force
        vbias[0] = -fe_value[0]
        for ldx in xrange(natoms):
            for kdx in xrange(3):
                for jdx in xrange(dmsc.ndcs):
                    force[kdx, ldx] += fe_gradient[jdx]*gradient_dcs[jdx, kdx, ldx]
    elif isempty == 1:
        vbias[0] = 0.0
        # the force is set as 0.0 
        for ldx in xrange(natoms):
            for kdx in xrange(3):
                force[kdx, ldx] = 0.0

    return 0

cdef int save_data(np.ndarray[np.float64_t,ndim=2] coord, BiasedMD *bs, DMSConfig *dmsc):

    w = writer.open('.gro', pattern='tmp.gro')
    w.write(coord, 'confall.gro', mode='a')

    with open('confall.w', 'a') as wfile:
        print >> wfile, '%.18e' %(exp(bs.vbias/dmsc.kT))

    if dmsc.isfirst == 0:
        with open('confall.ev', 'a') as evfile:
            print >> evfile, ' '.join(['%.18e' % (bs.dcs[idx],) for idx in xrange(dmsc.ndcs)])
    return 0

