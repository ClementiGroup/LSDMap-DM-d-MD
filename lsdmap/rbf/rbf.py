import os
import sys
import argparse
import logging
import ConfigParser
import numpy as np
from mpi4py import MPI

from lsdmap.rw import reader
from lsdmap.rw import writer
from lsdmap.mpi import p_index
from lsdmap.util import metric as mt

class RbfLib(object):

    def __init__(self, fit):

        fit_name = fit.lower()
        _mapped = {'inverse': 'inverse_multiquadric',
                       'inverse multiquadric': 'inverse_multiquadric',
                       'thin-plate': 'thin_plate'}

        if fit_name in _mapped:
            self.name = _mapped[fit]
        else:
            self.name = fit_name

        function_name = "_h_" + self.name
        if hasattr(self, function_name):
            self.function = getattr(self, function_name)
            derivative_name = "_d_" + self.name
            if hasattr(self, derivative_name):
                self.derivative = getattr(self, derivative_name)
        else:
            functionlist = [x[3:] for x in dir(self) if x.startswith('_h_')]
            raise ValueError("RbfFit: function must be one of " + ", ".join(functionlist))

    def _h_multiquadric(self, r, sigma):
        return np.sqrt((1.0/sigma*r)**2 + 1)
    
    def _h_inverse_multiquadric(self, r, sigma):
        return 1.0/np.sqrt((1.0/sigma*r)**2 + 1)
    
    def _h_gaussian(self, r, sigma):
        return np.exp(-(1.0/sigma*r)**2)
    
    def _h_linear(self, r, sigma):
        return r
    
    def _h_cubic(self, r, sigma):
        return r**3
    
    def _h_quintic(self, r, sigma):
        return r**5
    
    def _h_thin_plate(self, r, sigma):
        result = r**2 * np.log(r)
        if isinstance(r, (long, float)):
            if r == 0: result = 0
        else:
            result[r == 0] = 0 # the spline is zero at zero
        return result
    
    # derivatives
    def _d_thin_plate(self, r, sigma):
        result = r * (1 + 2 * np.log(r))
        return result

class RbfFit(object):

    def __init__(self, comm, coords, values, distance_matrix=None, metric='rmsd', metric_prms={},\
                  sigma=None, ksigma=None, fit='inverse_multiquadric'):

        self.coords = coords
        self.values = values
 
        self.npoints = self.coords.shape[0]
        self.ndim = self.coords.shape[1]
        self.natoms = self.coords.shape[2]

        self.sigma = sigma
        self.ksigma = ksigma

        if self.values.shape[0] != self.npoints:

            raise ValueError('Number of values provided does not match with the number of coords')

        self.metric = metric
        self.metric_prms = metric_prms
        self.fit = fit

        self.metric_function = mt.Metric(self.metric, ndim=self.ndim, **self.metric_prms).function
        self.fit_function = RbfLib(self.fit).function

        self.distance_matrix = distance_matrix
        self.weights = self.get_weights(comm)

    def initialize_sigma_values(self, sigma, ksigma):

        if sigma is None:
            if ksigma is None:
                return np.empty(self.npoints)
            else:
                sigma = np.zeros(self.npoints, dtype='float')
                for idx, distance in enumerate(self.distance_matrix):
                    idx_neighbors = np.argsort(distance)[1:ksigma+1]  # the first element is the point itself
                    sigma[idx] = self.distance_matrix[idx][idx_neighbors[-1]]
        else:
            if ksigma is not None:
                raise ValueError('RbfFit: conflict sigma values provided but ksigma specified')

        sigma = np.array(sigma)

        if sigma.size == 1:
            sigma = sigma*np.ones(self.npoints)
        if sigma.size !=self.npoints:
            raise ValueError('RbfFit: number of sigma values does not match the number of points') 

        return sigma

    def get_weights(self, comm):

        # compute the distance matrix if needed
        if self.distance_matrix is None:
            logging.info("distance matrix not provided, computing it...")

            idxs_thread = p_index.get_idxs_thread(comm, self.npoints)
            npoints_thread = len(idxs_thread)
            coords_thread = np.array([self.coords[idx] for idx in idxs_thread])
            DistanceMatrix = mt.DistanceMatrix(coords_thread, self.coords, metric=self.metric, metric_prms=self.metric_prms)
            self.distance_matrix = np.vstack(comm.allgather(DistanceMatrix.distance_matrix))
            logging.info("distance matrix computed")
        else: 
            if any(self.distance_matrix.shape[idx] != self.npoints for idx in [0, 1]):
                logging.error("distance matrix provided doesn't match the number of coordinates")

        self.sigma = self.initialize_sigma_values(self.sigma, self.ksigma)

        # invert kernel matrix
        kernel_matrix = self.fit_function(self.distance_matrix, self.sigma)
        # perform singular valu decomposition
        U, s, V = np.linalg.svd(kernel_matrix)
        sinv = 1/s
        # filter noisy singular values values
        sinv[s<0.02*np.max(s)] = 0
        inverse_kernel_matrix = np.dot(np.dot(V.T, np.diag(sinv)), U.T)
        weights = np.dot(inverse_kernel_matrix, self.values)
        logging.info("kernel matrix inverted")

        return weights

    def __call__(self, query):
        distance = np.array([self.metric_function(query, coord) for coord in self.coords])
        return np.sum(self.weights*self.fit_function(distance, self.sigma))


#TODO: create exe class for LSDMap
#TODO: create a function to check if parallel reading is needed
class RbfExe(object):

    def initialize(self, comm, config, args):

        self.config = config
        self.args = args

        # load configurations
        format_struct_file = os.path.splitext(args.struct_file[0])[1]
        if format_struct_file == '.gro': # use lsdmap reader
            struct_file = reader.open(args.struct_file)
            self.npoints = struct_file.nlines

            idxs_thread = p_index.get_idxs_thread(comm, self.npoints)
            coords_thread = struct_file.readlines(idxs_thread)
            self.coords = np.vstack(comm.allgather(coords_thread))
        else: # use numpy
            rank = comm.Get_rank()
            if rank == 0:
                coords = np.loadtxt(args.struct_file[0])
                if len(coords.shape)==0:
                    self.coords = coords[:, np.newaxis, np.newaxis] 
                else:
                    self.coords = coords[:, :, np.newaxis]
            else:
                self.coords = None
            self.coords = comm.bcast(self.coords, root=0)
            self.npoints = self.coords.shape[0]
        logging.info('input coordinates loaded')

        # load file of values
        valfile = reader.open(args.valfile)
        self.values = valfile.readlines()

        format = os.path.splitext(args.valfile)[1]
        if format == '.ev':
            self.fitdcs = True
        else:
            self.fitdcs = False
            if len(self.values.shape) > 2:
                raise ValueError('file of values should contain a single column')

        self.initialize_metric()
        self.function = config.get(self.args.section,'function')

        self.status_sigma = config.get(self.args.section,'status')
        if self.status_sigma == 'constant':
            self.sigma = config.getfloat(self.args.section,'sigma')
            self.ksigma = None
        elif self.status_sigma == 'kneighbor':
            self.sigma = None
            self.ksigma = config.getint(self.args.section,'ksigma')

        if args.embed_file is not None:
            embed_file = reader.open(args.embed_file)
            self.embed_filename = embed_file.filename
            self.npoints_embed = embed_file.nlines

            self.idxs_thread_embed = p_index.get_idxs_thread(comm, self.npoints_embed)

            if hasattr(embed_file, '_skip'): # multi-thread reading
                coords_thread_embed = embed_file.readlines(self.idxs_thread_embed)
                self.coords_embed = np.vstack(comm.allgather(coords_thread_embed))
            else: # serial reading
                if rank == 0:
                    self.coords_embed = embed_file.readlines()
                else:
                    self.coords_embed = None
                self.coords_embed = comm.bcast(self.coords_embed, root=0)


    def initialize_metric(self):

        _known_prms = ['r0']
        config = self.config
        self.metric = config.get(self.args.section, 'metric')

        self.metric_prms = {}
        for prm in _known_prms:
            try:
                self.metric_prms[prm] = config.getfloat(self.args.section, prm)
            except:
                pass

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="Performs radial basis function fit..") # message displayed when typing rbffit -h

        # required options
        # the same config file as the one used for LSDMap can be used
        parser.add_argument("-f",
            type=str,
            dest="config_file",
            required=True,
            help='Configuration file (input): ini')

        parser.add_argument("-c",
            type=str,
            dest="struct_file",
            required=True,
            nargs='*',
            help='Structure file (input): gro, xvg')

        parser.add_argument("-v",
            type=str,
            dest="valfile",
            required=True,
            help='File containing the values of the fitting points (input) : ev, sl, w')

        parser.add_argument("--dc",
            type=int,
            dest="dc_orders",
            nargs='*',
            default=[1],
            help='DC orders (opt., file specified with -v option needs to be of .ev format)')

        parser.add_argument("--embed",
            type=str,
            dest="embed_file",
            help='Structure file with points to be embedded (input, opt.): gro')

        parser.add_argument("--section",
            type=str,
            dest="section",
            default='FITTING',
            help='Section of configuration file where fitting parameters are located')

        return parser

    def embed(self, comm, fit, dc_order=None):

        rank = comm.Get_rank()  # number of the current thread 

        if dc_order is None:
            logging.info('Start embedding procedure')
        else:
            logging.info('Start embedding procedure along DC %i'%dc_order)

        coords_thread_embed = np.array([self.coords_embed[idx] for idx in self.idxs_thread_embed])

        values_thread_embed = []
        for coord in coords_thread_embed:
           values_thread_embed.append(fit(coord))

        values_thread_embed = np.array(values_thread_embed)
        values_embed = np.hstack(comm.allgather(values_thread_embed))

        if dc_order is None:
            logging.info('Embedding done')
        else:
            logging.info('Embedding along DC %i done'%dc_order)

        return values_embed


    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD   # MPI environment
        size = comm.Get_size()  # number of threads
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args() # set argument parser

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file) # set config file parser

        logging.basicConfig(filename='fitting.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        logging.info('intializing Rbf Fitting...')
        self.initialize(comm, config, args)
        logging.info('Fitting initialized')

        if size > self.npoints:
            logging.error("number of threads should be less than the number of frames")
            raise ValueError

        if self.fitdcs: 
            dc_order = args.dc_orders[0]
            logging.info('Start fitting procedure along DC %i'%dc_order) 
            fit = RbfFit(comm, self.coords, self.values[:,dc_order], metric=self.metric, fit=self.function, sigma=self.sigma, ksigma=self.ksigma)        
            logging.info('Fitting along DC %i done'%dc_order)
        else:
            dc_order = None
            logging.info('Start fitting procedure')
            fit = RbfFit(comm, self.coords, self.values, metric=self.metric, fit=self.function, sigma=self.sigma, ksigma=self.ksigma) 
            logging.info('Fitting done')

        weights = fit.weights[np.newaxis]
        sigma = fit.sigma[np.newaxis]

        if args.embed_file is not None:        
            values_embed = self.embed(comm, fit, dc_order=dc_order)
            values_embed = values_embed[np.newaxis]

        if self.fitdcs:
            if len(args.dc_orders) > 1:
                # use distance matrix computed before
                distance_matrix = fit.distance_matrix
                for dc_order in args.dc_orders[1:]:
                    logging.info('Start fitting procedure along DC %i'%dc_order)
                    fit = RbfFit(comm, self.coords, self.values[:,dc_order], distance_matrix=distance_matrix, metric=self.metric,
                                 fit=self.function, sigma=self.sigma, ksigma=self.ksigma)
                    logging.info('Fitting along DC %i done'%dc_order)
                    weights = np.concatenate((weights, fit.weights[np.newaxis]))
                    sigma = np.concatenate((sigma, fit.sigma[np.newaxis]))
                    if args.embed_file is not None:
                        values_embed = np.concatenate((values_embed, self.embed(comm, fit, dc_order)[np.newaxis]))

        if rank == 0:
            wfile = "fit_w.npy"
            np.save(wfile, weights)

            sigfile = "fit_sig.npy"
            np.save(sigfile, sigma)

            if args.embed_file is not None:
                np.save("fit_embed.npy", values_embed)
