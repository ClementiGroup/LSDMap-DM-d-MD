import os
import sys
import numpy as np
import argparse
import ConfigParser
import pickle
import random
import logging
from time import time

from lsdmap.rw import reader
from lsdmap.util import p_arpack
from lsdmap.util import metric as mt

from mpi4py import MPI


class LSDMap(object):

    def initialize(self, config, args):

        self.config = config
        self.args = args

        struct_file = reader.open(args.struct_file)
        self.struct_filename = struct_file.filename
        self.coords = struct_file.readlines()
        logging.info('input coordinates loaded')

        self.npoints = self.coords.shape[0]

        self.initialize_local_scale()
        self.initialize_weights()
        self.initialize_metric()

        self.neigs = 10


    def initialize_local_scale(self):

        config = self.config
        args = self.args

        known_status = ['constant', 'kneighbor', 'user', 'kneighbor_mean']
        _mapped = {'const': 'constant', 'cst': 'constant', 'mean-kneighbor': 'mean_kneighbor'}

        if args.epsfile is None:
            status = config.get('LOCALSCALE','status')
            if status in _mapped:
                status = _mapped[status]
            if not status in known_status:
                raise ValueError("local scale status should be one of "+ ', '.join(known_status))
            if status in ['kneighbor', 'kneighbor_mean']:
                value = None
                self.k = config.getint('LOCALSCALE', 'k')
            if status == 'constant':
                value = config.getfloat('LOCALSCALE', 'epsilon')
                value = value*np.ones(self.npoints)
            if status == 'user':
                raise NameError("status 'user' specified in configuration file but epsfile not provided. Please consider option flag -e")
        else:
            espfile = reader.open(epsfile)
            status ='user'
            value = epsfile.readlines()

        self.status_epsilon = status
        self.epsilon = value


    def initialize_weights(self):

        config = self.config
        args = self.args

        if args.wfile is not None:
            wfile = reader.open(args.wfile)
            weights = wfile.readlines()
        else:
            weights = np.ones(self.npoints, dtype='float')
        self.weights = weights


    def initialize_metric(self):

        _known_prms = ['r0']

        config = self.config
        self.metric = config.get('GENERAL','metric')

        self.metric_prms = {}
        for prm in _known_prms:
            try:
                self.metric_prms[prm] = config.getfloat('GENERAL', prm)
            except:
                pass


    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run LSDMap..") # message displayed when typing lsdmap -h

        # required options
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
            help = 'Structure file (input): gro, xvg')

        # other options
        parser.add_argument("-o",
            type=str,
            dest="output_file",
            help='LSDMap file (output, opt.): lsdmap')

        parser.add_argument("--nn",
            type=str,
            dest="nnfile",
            help='File containing the indices of the nearest neighbors of all configurations of the structure file (output, opt.): nn')

        parser.add_argument("--nneighbors",
            type=int,
            dest="nneighbors",
            help="Number of nearest neighbors stored in .nn file; should be used with option --nn.")

        parser.add_argument("--nneighbors_cutoff",
            type=float,
            dest="nneighbors_cutoff",
            help="Only the indices of neighbors closer than a distance of nneighbors_cutoff are stored in .nn file; should be used with option --nn.")

        parser.add_argument("-w",
            type=str,
            dest="wfile",
            help='File containing the weights of every point in a row (input, opt.): .w')

        parser.add_argument("-e",
            type=str,
            dest="epsfile",
            help='File containing the local scales of every point in a row (input, opt.): .eps')


        return parser


    def get_idxs_thread(self, comm):
        """
        get the indices of the coordinates that will be consider by each processor (use MPI)
        """
        size = comm.Get_size()
        rank = comm.Get_rank()

        npoints_per_thread = np.array([self.npoints/size for idx in range(size)])

        for idx in range(self.npoints%size):
            npoints_per_thread[idx]+=1

        if rank == 0:
            idx_shift = 0
            idxs_threads = []
            for idx in npoints_per_thread:
                idxs_threads.append(range(idx_shift,idx_shift+idx))
                idx_shift += idx
        else:
            idxs_threads = None

        idxs_thread = comm.scatter(idxs_threads, root=0)

        return idxs_thread


    def compute_kernel(self, comm, npoints_thread, distance_matrix_thread, weights_thread, epsilon_thread):
        # for a detailed description of the following operations, see the paper:
            # Determination of reaction coordinates via locally scaled diffusion map
            # Mary A. Rohrdanz, Wenwei Zheng, Mauro Maggioni, and Cecilia Clementi
            # The Journal of Chemical Physics 134, 124116 (2011)

        p_vector_thread = np.zeros(npoints_thread, dtype='float')
        d_vector_thread = np.zeros(npoints_thread, dtype='float')

        # compute LSDMap kernel, Eq. (5) of the above paper
        kernel = np.sqrt((weights_thread[:, np.newaxis]).dot(self.weights[np.newaxis])) * \
                 np.exp(-distance_matrix_thread**2/(2*epsilon_thread[:, np.newaxis].dot(self.epsilon[np.newaxis])))

        p_vector_thread = np.sum(kernel, axis=1)
        p_vector = np.hstack(comm.allgather(p_vector_thread)) # Eq. (6)
        self.p_vector = p_vector

        kernel /= np.sqrt(p_vector_thread[:,np.newaxis].dot(p_vector[np.newaxis])) # Eq. (7)
        d_vector_thread = np.sum(kernel, axis=1)
        d_vector = np.hstack(comm.allgather(d_vector_thread)) # Compute D given between Eqs. (7) and (8)
        self.d_vector = d_vector

        kernel /= np.sqrt(d_vector_thread[:,np.newaxis].dot(d_vector[np.newaxis])) # Eq (8) (slightly modified)

        return kernel


    def save(self, args):
        """
        save LSDMap object in .lsdmap file and eigenvalues/eigenvectors in .eg/.ev files
        """

        lsdmap_filename = args.output_file

        if isinstance(self.struct_filename, list):
            struct_filename = self.struct_filename[0]
        else:
            struct_filename = self.struct_filename

        if lsdmap_filename is None:
            path, ext = os.path.splitext(struct_filename)
            lsdmap_filename = path + '.lsdmap'
            np.savetxt(path + '.eg', np.fliplr(self.eigs[np.newaxis]), fmt='%9.6f')
            np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%15.7e')
        else:
            path, ext = os.path.splitext(lsdmap_filename)
            np.savetxt(path + '.eg', np.fliplr(self.eigs[np.newaxis]), fmt='%9.6f')
            np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%15.7e')
            if ext != '.lsdmap':
                lsdmap_filename = lsdmap_filename + '.lsdmap'

        with open(lsdmap_filename, "w") as file:
            pickle.dump(self, file)


    def save_nneighbors(self, comm, args, DistanceMatrix):

        size = comm.Get_size()  # number of threads
        rank = comm.Get_rank()  # number of the current thread

        if rank == 0:
            try:
                os.remove(args.nnfile)
            except OSError:
                pass
            nnfile = open(args.nnfile, 'a')
            for idx in xrange(size):
                if idx == 0:
                    idx_neighbor_matrix = DistanceMatrix.idx_neighbor_matrix()
                    neighbor_matrix = DistanceMatrix.neighbor_matrix()
                else:
                    idx_neighbor_matrix = comm.recv(source=idx, tag=idx)
                    neighbor_matrix = comm.recv(source=idx, tag=1000+idx) 
                if args.nneighbors is not None: # save first "args.nneighbors" nearest neighbors
                    for idxs in idx_neighbor_matrix:
                        np.savetxt(nnfile, idxs[:args.nneighbors][np.newaxis], fmt='%d')
                elif args.nneighbors_cutoff is not None: # save all nearest neighbors within "args.nneighbors_cutoff"
                    for idxrow, row in enumerate(neighbor_matrix):
                        for idxcol, col in enumerate(row):
                            if col > args.nneighbors_cutoff:
                                break
                        np.savetxt(nnfile, idx_neighbor_matrix[idxrow,:idxcol][np.newaxis], fmt='%d')
                else: # save all nearest neighbors
                    for idxs in idx_neighbor_matrix:
                        np.savetxt(nnfile, idxs[np.newaxis], fmt='%d')
            nnfile.close()
        else:
           idx_neighbor_matrix = DistanceMatrix.idx_neighbor_matrix()
           neighbor_matrix = DistanceMatrix.neighbor_matrix()
           comm.send(idx_neighbor_matrix, dest=0, tag=rank)
           comm.send(neighbor_matrix, dest=0, tag=1000+rank)
         

    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD   # MPI environment
        size = comm.Get_size()  # number of threads
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args() # set argument parser

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file) # set config file parser

        logging.basicConfig(filename='lsdmap.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        # consider calling initialize function right after setting parsers
        logging.info('intializing LSDMap...')
        self.initialize(config, args)
        logging.info('LSDMap initialized')

        if size > self.npoints:
            logging.error("number of threads should be less than the number of frames")
            raise ValueError

        idxs_thread = self.get_idxs_thread(comm)
        npoints_thread = len(idxs_thread)

        coords_thread = np.array([self.coords[idx] for idx in idxs_thread])
        weights_thread = np.array([self.weights[idx] for idx in idxs_thread])

        # compute the distance matrix
        DistanceMatrix = mt.DistanceMatrix(coords_thread, self.coords, metric=self.metric, metric_prms=self.metric_prms)
        distance_matrix_thread = DistanceMatrix.distance_matrix
        logging.info("distance matrix computed")

        # compute kth neighbor local scales if needed
        if self.status_epsilon in ['kneighbor', 'kneighbor_mean']:
            epsilon_thread = []
            idx_neighbor_matrix_thread = DistanceMatrix.idx_neighbor_matrix()
            for idx, line in enumerate(idx_neighbor_matrix_thread):
                cum_weight = 0
                for jdx in line[1:]:
                    cum_weight += self.weights[jdx]
                    if cum_weight >= self.k:
                        break
                epsilon_thread.append(distance_matrix_thread[idx,jdx])

            self.epsilon = np.hstack(comm.allgather(epsilon_thread)) # gather epsilon values
            if self.status_epsilon == 'kneighbor_mean':
                mean_value_epsilon = np.mean(self.epsilon) # compute the mean value of the local scales
                self.epsilon = mean_value_epsilon * np.ones(self.npoints)  # and set it as the new constant local scale

        epsilon_thread = np.array([self.epsilon[idx] for idx in idxs_thread])

        # compute kernel
        kernel = self.compute_kernel(comm, npoints_thread, distance_matrix_thread, weights_thread, epsilon_thread)

        if rank == 0:
            time1 = time()

        # diagonalize kernel
        params= p_arpack._ParallelSymmetricArpackParams(comm, kernel, self.neigs)
        while not params.converged:
            params.iterate()
        eigs, evs = params.extract(return_eigenvectors=True)

        # normalize eigenvectors
        self.evsu = np.copy(evs) # store unormalized eigenvectors
        evs /= np.sqrt(self.d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum(evs**2, axis=0))
        evs /= norm[np.newaxis,:]

        # store eigenvalues/eigenvectors
        self.eigs = eigs
        self.evs = evs

        logging.info("kernel diagonalized")
    
        if rank ==0:    
            self.save(args)
        logging.info("LSDMap object and eigenvalues/eigenvectors saved (.lsdmap and .eg/.ev files)")

        # store nearest neighbors in .nn file if specified via --nn option
        if args.nnfile is not None:
            self.save_nneighbors(comm, args, DistanceMatrix)

        logging.info("LSDMap computation done")
