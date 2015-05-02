import sys
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

    # iteration number
    dmsc.iter = config.getint('GENERAL', 'iter')
    # number of points saved per replica
    dmsc.nstride = config.getint('GENERAL', 'nstride')
    # number of first dcs used
    dmsc.ndcs = config.getint('GENERAL', 'ndcs')

    # number of MD steps skipped for the computation of the bias potential
    if config.has_option('MD', 'nstepbias'):
        dmsc.nstepbias = config.getint('MD', 'nstepbias')
    else:
        dmsc.nstepbias = 1

    # temperature
    kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
    temperature = config.getfloat('GENERAL', 'temperature')
    dmsc.kT = kb*temperature

    # fraction of the free energy we actually use for the bias potential
    if config.has_option('GENERAL', 'fefrac'):
        dmsc.fefrac = config.getfloat('GENERAL', 'fefrac')
    else:
        dmsc.fefrac = 1.0

    return &dmsc


cdef public FEHist* initFEHist(DMSConfig *dmsc, const char* file):

    cdef int shift, ndx

    if dmsc.iter > 0:

        fedir = '../../fe'
        # allocate memory for the bins
        feh.nbins = <int *>malloc(dmsc.iter*dmsc.ndcs*sizeof(int))
        feh.nnebins = <int *>malloc(dmsc.iter*dmsc.ndcs*sizeof(int))

        # allocate memory for the step values of the histogram
        feh.steps = <double *>malloc(dmsc.iter*dmsc.ndcs*sizeof(double))

        # get the total number of bins and non-empty bins
        for ndx in range(dmsc.iter):
            bins = np.loadtxt(fedir + '/bins.xyz.%i'%ndx)
            feh.nbins[ndx] = bins.shape[0]
            histo = np.loadtxt(fedir + '/hist.dat.%i'%ndx)
            feh.nnebins[ndx] = histo.shape[0]
            
        nbins_tot = 0
        nnebins_tot = 0
        for ndx in range(dmsc.iter):
            nbins_tot += feh.nbins[ndx]
            nnebins_tot += feh.nnebins[ndx]

        feh.nbins_tot = nbins_tot
        feh.nnebins_tot = nnebins_tot

        # allocate memory for the bins
        feh.bins = <double *>malloc(feh.nbins_tot*dmsc.ndcs*sizeof(double))
        # allocate memory for the idxs of nonempty bins
        feh.nebins_idxs_s = <int *>malloc(feh.nnebins_tot*sizeof(int))

        # allocate memory for the values of the free energy (nonempty bins)
        feh.values = <double *>malloc(feh.nnebins_tot*sizeof(double))

        for ndx in range(dmsc.iter):

            bins = np.loadtxt(fedir + '/bins.xyz.%i'%ndx)
            histo = np.loadtxt(fedir + '/hist.dat.%i'%ndx)
            nebins_idxs = histo[:,:dmsc.ndcs].astype(int)
            nebins_idxs_s = histo[:,dmsc.ndcs].astype(int)
            free_energy = histo[:,dmsc.ndcs+1]

            if dmsc.ndcs == 1:
                bins = bins[:,np.newaxis]
                nebins_idxs = nebins_idxs[:,np.newaxis]

            if np.any(free_energy > 0.0):
                raise IOError("free energy values in hist.dat should be all negative!")

            shift = 0
            for mdx in range(ndx):
                shift += feh.nbins[mdx]

            for idx in xrange(feh.nbins[ndx]):
                for jdx in xrange(dmsc.ndcs):
                    # update values of the bins
                    feh.bins[shift*dmsc.ndcs+idx+feh.nbins[ndx]*jdx] = bins[idx,jdx]

            shift = 0
            for mdx in range(ndx):
                shift += feh.nnebins[mdx]

            for idx in xrange(feh.nnebins[ndx]):
                # update values of the free energy
                feh.values[shift+idx] = free_energy[idx]
                feh.nebins_idxs_s[shift+idx] = nebins_idxs_s[idx]

            for jdx in xrange(dmsc.ndcs):
                feh.steps[ndx*dmsc.ndcs+jdx] = bins[1,jdx] - bins[0,jdx]

    return &feh

cdef public Fit* initFit(DMSConfig *dmsc, const char* file):

    cdef unsigned int idx, jdx, kdx
    cdef int ndim, natoms 

    if dmsc.iter > 0:

        # initialize parser
        config = ConfigParser.SafeConfigParser()
        config.read(file)

        # number of configs used for the fit 
        npoints = config.getint('FITTING', 'npoints')

        # load configs
        fitdir = '../../fit'
        grofile = fitdir + '/fit.gro'
        rg = reader.open(grofile) # reader .gro file
        coordsfit = rg.readlines()

        check_parameter(coordsfit.shape[0], npoints, 'nfit', grofile) # check number of configs
        ndim = coordsfit.shape[1] # number of spatial dimensions
        natoms = coordsfit.shape[2] # number of atoms

        # load weights
        wfile = fitdir + '/fit.w'
        weightsfit = np.loadtxt(wfile, dtype="f8")
        if dmsc.ndcs == 1:
            weightsfit = weightsfit[:,np.newaxis]

        check_parameter(weightsfit.shape[0], npoints, 'nfit', wfile)
        check_parameter(weightsfit.shape[1], dmsc.ndcs, 'ndcs',  wfile)

        # load sigma values
        sigfile = fitdir + '/fit.sig'
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
    water_and_ions = ['OW', 'NA', 'CL']
    heavy_atoms_idxs = []
 
    f = open('tmp.gro', 'r')
    f.next()
    natoms = int(f.next())
    for atom_idx, line in it.izip(xrange(natoms), f):
        amino_acid = line[0:9].strip()
        atoms = line[9:15].lstrip()
        if "CXH" not in amino_acid:
        # check if not an hydrogen atom, not water and not an ion
            if atoms[0] != 'H' and atoms not in water_and_ions:
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
    # allocate array for prev. config. and prev. force
    bs.store_coord = <float *>malloc(bs.nheavy_atoms*3*sizeof(float))
    bs.store_biasforce = <float *>malloc(bs.nheavy_atoms*3*sizeof(float))
    bs.HH = 0

    return &bs

cdef public int do_biased_force(BiasedMD *bs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx
    cdef int natoms = bs.natoms
    cdef int nheavy_atoms = bs.nheavy_atoms
    cdef int nsave, nsavedcs
    cdef float DeltaHH

    cdef double* vbias = <double *>malloc(sizeof(double))
    cdef double* dcs = <double *>malloc(dmsc.ndcs*sizeof(double))

    cdef np.ndarray[np.float64_t,ndim=2] coord = np.zeros((3, natoms))
    cdef np.ndarray[np.float64_t,ndim=2] vel = np.zeros((3, natoms))
    cdef np.ndarray[np.float64_t,ndim=2] coord_heavy_atoms = np.zeros((3, nheavy_atoms))
    cdef np.ndarray[np.float64_t,ndim=2] force_heavy_atoms = np.zeros((3, nheavy_atoms))

    np.import_array()
    for jdx in xrange(3):
        for kdx in xrange(natoms):
            coord[jdx,kdx] = float(bs.coord[jdx+3*kdx])

    for jdx in xrange(3):
        for kdx in xrange(natoms):
            vel[jdx,kdx] = float(bs.vel[jdx+3*kdx])

    # select only heavy atoms to compute the biased force
    heavy_atoms_idxs = [bs.heavy_atoms_idxs[idx] for idx in xrange(bs.nheavy_atoms)]
    coord_heavy_atoms  = coord[:,heavy_atoms_idxs]

    # store configuration and weight after nsave steps
    nsave = max(bs.nsteps/dmsc.nstride, 1) # in dms.py, we already check if nsteps is a multiple of nstride 
    if bs.step%nsave == 0 and bs.step > 0:
        save_data(coord, vel, heavy_atoms_idxs, bs, dmsc)

    # compute biased force if not first iteration
    if dmsc.iter > 0:
         # save config. of heavy atoms at time step previous to application of bias force
         if bs.step % dmsc.nstepbias == dmsc.nstepbias-1:
             for idx in xrange(3):
                 for jdx in xrange(len(heavy_atoms_idxs)):
                     bs.store_coord[idx+3*jdx] = coord[idx,heavy_atoms_idxs[jdx]]
         # compute monitor quantity at step subsequent to application of bias force
         if (bs.step % dmsc.nstepbias == 1 and bs.step > 2*dmsc.nstepbias-1):
             for idx in xrange(3):                              
                 for jdx in xrange(len(heavy_atoms_idxs)):
                     bs.HH += (coord_heavy_atoms[idx,jdx]-bs.store_coord[idx+3*jdx])*bs.store_biasforce[idx+3*jdx]/2+(bs.vbias-bs.vbias_prev)
             #print bs.HH

         if bs.step % dmsc.nstepbias == 0:
             do_biased_force_low_level(nheavy_atoms, coord_heavy_atoms, force_heavy_atoms, vbias, dcs, dmsc, ft, feh)
             for jdx in xrange(dmsc.ndcs):
                 bs.dcs[jdx] = dcs[jdx]
             # save vbias before updating it
             bs.vbias_prev =  bs.vbias
             bs.vbias = vbias[0]

             for idx in xrange(3):
                 for jdx, atom_jdx in enumerate(heavy_atoms_idxs):
                     bs.force[idx+3*atom_jdx] += force_heavy_atoms[idx][jdx]*dmsc.nstepbias
                     # store only bias force
                     bs.store_biasforce[idx+3*jdx] = force_heavy_atoms[idx][jdx]*dmsc.nstepbias
    else:
        bs.vbias = 0.0

    return 0

cdef int do_biased_force_low_level(int natoms, np.ndarray[np.float64_t,ndim=2] coord, np.ndarray[np.float64_t,ndim=2] force, 
                                   double* vbias, double* dcs, DMSConfig *dmsc, Fit *ft, FEHist *feh):

    cdef unsigned int idx, jdx, kdx, ldx
    cdef double* wdf = <double *>malloc(dmsc.ndcs*sizeof(double))
    cdef fitfunction_type fitfunction, fitderivative

    cdef np.ndarray[np.float64_t,ndim=2] gradient_metric = np.zeros((3, natoms))
    cdef np.ndarray[np.float64_t,ndim=3] coordsfit = np.zeros((ft.npoints, 3, natoms))
    cdef np.ndarray[np.float64_t,ndim=3] gradient_dcs = np.zeros((dmsc.ndcs, 3, natoms))

    cdef int* isempty = <int *>malloc(sizeof(int))
    cdef double* fe_value = <double *>malloc(sizeof(double))
    cdef double* fe_gradient = <double *>malloc(dmsc.ndcs*sizeof(double))

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

    fe_value[0] = 0.0
    for jdx in range(dmsc.ndcs):
        fe_gradient[jdx] = 0.0

    for ndx in range(dmsc.iter):
        compute_free_energy_dcs(fe_value, fe_gradient, dcs, dmsc.iter-ndx-1, isempty, dmsc, feh)
        if isempty[0] == 1:
            break

    vbias[0] -= dmsc.fefrac*fe_value[0]
    for ldx in xrange(natoms):
        for kdx in xrange(3):
            for jdx in xrange(dmsc.ndcs):
                force[kdx, ldx] += dmsc.fefrac*fe_gradient[jdx]*gradient_dcs[jdx, kdx, ldx]
    return 0

cdef int compute_free_energy_dcs(double* fe_value, double* fe_gradient, double* dcs, int iter, int* isempty, DMSConfig *dmsc, FEHist *feh):

    cdef unsigned int mdx, idx, jdx, kdx, ldx
    cdef int bin_idx_s, num_line
    cdef int shift_nbins, shift_nnebins
    cdef int idxtmp

    cdef np.ndarray[np.int32_t,ndim=2] nebins_idxs = np.zeros((feh.nnebins[iter], dmsc.ndcs), dtype='i4')

    shift_nbins = 0
    for mdx in range(iter):
        shift_nbins += feh.nbins[mdx]

    shift_nnebins = 0
    for mdx in range(iter):
        shift_nnebins += feh.nnebins[mdx]

    # estimate the free energy of the current configuration and its derivative
    bin_idx_s = 0 # bin serial number
    bin_idxs = [] # bin indices
    isempty[0] = 1 # by default the bin is empty

    # compute the index of the bin within which the point is located
    for jdx in xrange(dmsc.ndcs):
        index = int(floor((dcs[jdx] - feh.bins[shift_nbins*dmsc.ndcs+feh.nbins[iter]*jdx])/feh.steps[iter*dmsc.ndcs+jdx] + 0.5))
        bin_idxs.append(index)
        bin_idx_s += index*feh.nbins[iter]**(dmsc.ndcs-jdx-1)

    if any([idxtmp < 0 or idxtmp >= feh.nbins[iter] for idxtmp in bin_idxs]):
        # if the point is outside the grid, the serial number does not apply
        isempty[0] = 1
    else:
        # check if the bin is empty or not
        for kdx in xrange(feh.nnebins[iter]):
            if bin_idx_s == feh.nebins_idxs_s[shift_nnebins+kdx]: # non-empty bin found
                num_line = kdx
                isempty[0] = 0
                break
            elif bin_idx_s < feh.nebins_idxs_s[shift_nnebins+kdx]: # the bin is empty
                isempty[0] = 1
                break

    # estimate the value of the free energy and its gradient
    if isempty[0] == 0:
        # interpolate the free energy locally
        local_fe_fit_voronoi(fe_value, fe_gradient, dcs, dmsc.ndcs, 5, feh, iter, 0.6)
    elif isempty[0] == 1:
        # vbias is set as minus the maximum value of the free energy
        fe_value[0] += 0.0
        # the force is set as 0.0
        for jdx in xrange(dmsc.ndcs):
            fe_gradient[jdx] += 0.0

cdef int save_data(np.ndarray[np.float64_t,ndim=2] coord, np.ndarray[np.float64_t,ndim=2] vel, heavy_atoms_idxs, BiasedMD *bs, DMSConfig *dmsc):

    w = writer.open('.gro', pattern='tmp.gro')
    w.write(np.vstack((coord, vel)), 'confall_aa.gro', mode='a') # write configurations with all atoms
    w.write(np.vstack((coord, vel)), 'confall.gro', idxs_atoms=heavy_atoms_idxs, mode='a') # write configurations with heavy atoms
    w.close()

    with open('confall.w', 'a') as wfile:
        print >> wfile, '%.18e'%(exp(bs.vbias/dmsc.kT))

    if dmsc.iter > 0:
        with open('confall.ev', 'a') as evfile:
            print >> evfile, ' '.join(['%.18e' % (bs.dcs[idx],) for idx in xrange(dmsc.ndcs)])
    return 0
