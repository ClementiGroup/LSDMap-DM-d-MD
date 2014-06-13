import os
import sys
import argparse
import ConfigParser
import pickle
import random
from time import time

import numpy as np
from mpi4py import MPI

from lsdmap.rw import reader
from lsdmap.util import metric as mt


class LSDMapRest(object):

    def initialize(self, args, LSDMap):

        struct_file = reader.open(args.struct_file)
        self.struct_filename = struct_file.filename
        self.coords = struct_file.readlines()
        self.npoints = self.coords.shape[0]

        if args.wfile is not None:
            wfile = reader.open(args.wfile)
            weights = wfile.readlines()
            if npoints != len(weights):
                raise("Number of lines in .w file does not match the number of frames in structure file")
        else:
            weights = np.ones(self.npoints, dtype='float')

        self.weights = weights

        if LSDMap.status_epsilon in ['constant', 'kneighbor_mean']:
            self.epsilon = LSDMap.epsilon[:self.npoints]
        elif LSDMap.status_epsilon == 'user':
            raise NotImplementedError


    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Restriction from physical space to LSDMap variables...") # message displayed when typing rlsdmap -h

        # required options
        parser.add_argument("-c",
           type=str,
           dest="struct_file",
           required=True,
           help = 'Structure file (input): gro')

        parser.add_argument("-s",
           type=str,
           dest="lsdmap_file",
           required=True,
           help ='LSDMap file (input): lsdmap')

        # other options
        parser.add_argument("-e",
            type=str,
            dest="epsfile",
            help='File containing the local scales of each points (input, opt.): txt')
        parser.add_argument("-w",
            type=str,
            dest="wfile",
            help='File containing the weights of each points (input, opt.): txt')
        parser.add_argument("--embed",
            dest="embed",
            action='store_true',
            help='embed LSDMap coordinates to old ones',
            default=False)
        parser.add_argument("-o",
            type=str,
            dest="output_file",
            help='File containing eigenvectors (output, opt.): .ev')

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


    def compute_restriction_kernel(self, comm, LSDMap, npoints_thread, distance_matrix_thread, weights_thread, epsilon_thread):

        p_vector_thread = np.zeros(npoints_thread, dtype='float')
        d_vector_thread = np.zeros(npoints_thread, dtype='float')

        kernel = np.sqrt(weights_thread[:, np.newaxis].dot(LSDMap.weights[np.newaxis])) * \
                 np.exp(-distance_matrix_thread**2/(2*epsilon_thread[:, np.newaxis].dot(LSDMap.epsilon[np.newaxis])))

        p_vector_thread = np.sum(kernel, axis=1)
        kernel /= np.sqrt(p_vector_thread[:,np.newaxis].dot(LSDMap.p_vector[np.newaxis]))
        self.p_vector = p_vector_thread

        d_vector_thread = np.sum(kernel, axis=1)
        kernel /= np.sqrt(d_vector_thread[:,np.newaxis].dot(LSDMap.d_vector[np.newaxis]))
        self.d_vector = np.hstack(comm.allgather(d_vector_thread))

        return kernel


    def save(self, ev_filename):
        """
        save eigenvectors in .ev files
        """

        if isinstance(self.struct_filename, list):
            struct_filename = self.struct_filename[0]
        else:
            struct_filename = self.struct_filename

        if ev_filename is None:
            path, ext = os.path.splitext(struct_filename)
            np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%15.7e')
        else:
            path, ext = os.path.splitext(ev_filename)
            if ext != '.ev':
                ev_filename = ev_filename + '.ev'
            np.savetxt(ev_filename, np.fliplr(self.evs), fmt='%15.7e')


    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD
        size = comm.Get_size()  # number of threads involved
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args()

        with open(args.lsdmap_file, "r") as file:
            LSDMap = pickle.load(file)
        self.initialize(args, LSDMap)

        if size > self.npoints:
            raise ValueError("number of threads should be less than the number of configurations")

        idxs_thread = self.get_idxs_thread(comm)
        npoints_thread = len(idxs_thread)

        coords_thread = np.array([self.coords[idx] for idx in idxs_thread])
        weights_thread = np.array([self.weights[idx] for idx in idxs_thread])

        # compute the distance matrix
        if rank == 0:
            time1 = time()
        DistanceMatrix = mt.DistanceMatrix(coords_thread, LSDMap.coords, metric=LSDMap.metric)
        distance_matrix_thread = DistanceMatrix.distance_matrix
        if rank == 0:
            time2 = time()
            print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

        # compute kth neighbor local scales if needed
        if LSDMap.status_epsilon == 'kneighbor':
            epsilon_thread = DistanceMatrix.neighbor_matrix(k=LSDMap.k+1)[:,-1]
            self.epsilon = np.hstack(comm.allgather(epsilon_thread))

        epsilon_thread = np.array([self.epsilon[idx] for idx in idxs_thread])

        # compute kernel
        kernel = self.compute_restriction_kernel(comm, LSDMap, npoints_thread, distance_matrix_thread, weights_thread, epsilon_thread) 

        # compute eigenvectors from nystrom formula
        evs_thread = kernel.dot(LSDMap.evsu)/LSDMap.eigs
        evs = np.concatenate(comm.allgather(evs_thread), axis=0)

        # normalization
        evs /= np.sqrt(self.d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum((LSDMap.evsu/np.sqrt(LSDMap.d_vector[:,np.newaxis]))**2, axis=0))
        evs /= norm[np.newaxis,:]
        self.evs = evs

        # print eigenvectors and eigenvalues in .ev and .eg files
        if rank == 0:
            self.save(args.output_file)
