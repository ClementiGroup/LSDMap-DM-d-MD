import os
import sys
import argparse
import ConfigParser
import pickle
import random
import logging
from time import time

import numpy as np
from mpi4py import MPI

from lsdmap.rw import reader
from lsdmap.rw import writer
from lsdmap.util import metric as mt


class LSDMapLift(object):

    def initialize(self, args, LSDMap):

        evfile = reader.open(args.evfile)
        targets = evfile.readlines()

        self.ntargets = targets.shape[0]
        self.nevs = min(targets.shape[1], args.nevs)
        self.targets = targets[:,:self.nevs]

        if len(LSDMap.coords.shape) == 2:
            self.ndof = LSDMap.coords.shape[1]
        elif len(LSDMap.coords.shape) == 3:
            self.ndof = LSDMap.coords.shape[1]*LSDMap.coords.shape[2]

        if not all(weight == 1.0 for weight in LSDMap.weights):
            raise ValueError("file %s contains weighted version of LSDMap; cannot lift points from weighted LSDMap")

        self.niters = args.niters
        self.value_b = args.value_b

        temperatures = [args.temperature_max]
        for idx in xrange(args.nsteps-1):
            temperatures.append((1-args.ratio)*temperatures[idx])

        self.temperatures = temperatures

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Lifting from diffusion map to physical variables...") # message displayed when typing llsdmap -h

        # required options
        parser.add_argument("-f",
           type=str,
           dest="evfile",
           required=True,
           help='File containing the eigenvectors components used to lift (input) : ev')

        parser.add_argument("-s",
           type=str,
           dest="lsdmap_file",
           required=True,
           help='LSDMap file (input): lsdmap')

        # other options
        parser.add_argument("-o",
            type=str,
            dest="output_file",
            help='File containing the output configurations (output, opt.): .gro',
            default='out.gro')

        parser.add_argument("-e",
            type=str,
            dest="ev_output_file",
            help='File containing the output eigenvectors (output, opt.): .ev',
            default='out.ev')

        parser.add_argument("-n",
           "--nevs",
           type=int,
           dest="nevs",
           help='Number of eigenvalues used to lift (opt., default=3)',
           default=3)

        parser.add_argument("--temperature_max",
           type=float,
           dest="temperature_max",
           help='Initial (maximum) temperature for simulated annealing (opt., default=1.0)',
           default=1.0)

        parser.add_argument("--nsteps",
           type=int,
           dest="nsteps",
           help='Number of steps for simulated annealing (opt., default=20)',
           default=20)

        parser.add_argument("--niters",
           type=int,
           dest="niters",
           help='Number of iterations for Nelder Mead method (opt., default=100)',
           default=100)

        parser.add_argument("--ratio",
           type=float,
           dest="ratio",
           help='Temperature is decreased by multiplying it by (1-ratio) at each step of simulated annealing (opt., default=0.1)',
           default=0.1)

        parser.add_argument("--value_b",
           type=float,
           dest="value_b",
           help='Best cost function value to start with, (opt., default=1.0)',
           default=1.0)

        return parser


    def initialize_simplex(self, target, LSDMap):

        evs = np.fliplr(LSDMap.evs)[:,:self.nevs]
        diff = np.sqrt(np.sum((evs - target)**2, axis=1))

        simplex_idxs = np.argsort(diff)[:self.ndof+1]
        simplex = np.array([LSDMap.coords[idx] for idx in simplex_idxs])
        evs_simplex = np.array([LSDMap.evs[idx] for idx in simplex_idxs])[:,:self.nevs]
        values = np.sqrt(np.sum((evs_simplex - target)**2, axis=1))

        return simplex, values


    def generate_new_value(self, target, LSDMap, fac):

        fac1 = (1.0 - fac)/self.ndof
        fac2 = fac1 - fac

        coords = self._ssimplex * fac1 - self._simplex[self._idx_hi,:] * fac2
        evs = self.nystrom(coords[np.newaxis], LSDMap, self.nevs)
        value = self.cost_function(evs, target)

        if value[0] <= self._value_b:
            self._coords_b = coords
            self._value_b = value[0]

        value_rand = value[0] - self._temperature*np.log(np.random.rand());

        if value_rand < self._value_hi:
            self._values[self._idx_hi] = value[0]
            self._value_hi = value_rand
            self._ssimplex += coords - self._simplex[self._idx_hi,:]
            self._simplex[self._idx_hi,:] = coords

        return value_rand

    def neldermead(self, target, LSDMap, simplex, values, value_b, temperature):
        """
        Nelder Mead method used to find the minimum of cost function at a given temperature
        """

        niters = self.niters

        self._simplex = simplex
        self._values = values
        self._temperature = temperature
        self._value_b = value_b

        self._ssimplex = np.sum(self._simplex, axis=0) # ssimplex has self.ndof components

        while True:

            values_rand = values - self._temperature*np.log(np.random.rand(self.ndof+1))
            values_rand_idxs = np.argsort(values_rand)

            self._idx_l = values_rand_idxs[0]
            self._idx_hi = values_rand_idxs[-1]
            self._idx_nhi = values_rand_idxs[-2]

            self._value_l = values_rand[self._idx_l]
            self._value_hi = values_rand[self._idx_hi]
            self._value_nhi = values_rand[self._idx_nhi]

            #print rtol
            if niters < 0:
                break

            niters -= 2
            new_value = self.generate_new_value(target, LSDMap, -1.0)

            if new_value <= self._value_l:
                new_value = self.generate_new_value(target, LSDMap, 2.0)
            elif new_value >= self._value_nhi:
                value_save = self._value_hi
                new_value = self.generate_new_value(target, LSDMap, 0.5)
                if new_value >= value_save:
                    for idx in xrange(self.ndof+1):
                        if idx != self._idx_l:
                            self._ssimplex = 0.5*(self._simplex[idx,:]+self._simplex[self._idx_l,:])
                            self._simplex[idx,:] = self._ssimplex
                            evs = self.nystrom(self._ssimplex[np.newaxis], LSDMap, self.nevs)
                            self._values[idx] = self.cost_function(evs, target)
                    niters -= self.ndof
                    self._ssimplex = np.sum(self._simplex, axis=0)
            else:
                niters +=1

        tmp = self._values[0]
        self._values[0] = self._values[self._idx_l]
        self._values[self._idx_l] = tmp

        tmp = self._simplex[0,:]
        self._simplex[0,:] = self._simplex[self._idx_l,:]
        self._simplex[self._idx_l] = tmp

        return self._simplex, self._values, self._value_b

    def nystrom(self, coords, LSDMap, nevs):
        """
        compute eigenvectors of coordinates coords using nystrom formula
        """
        npoints = coords.shape[0]

        # compute the distance matrix
        DistanceMatrix = mt.DistanceMatrix(coords, LSDMap.coords, metric=LSDMap.metric)
        distance_matrix = DistanceMatrix.distance_matrix

        if LSDMap.status_epsilon == 'constant':
            epsilon = LSDMap.epsilon[:npoints]
        elif LSDMap.status_epsilon == 'kneighbor':
            epsilon = DistanceMatrix.neighbor_matrix(k=LSDMap.k+1)[:,-1]

        # compute LSDMap kernel
        kernel = np.exp(-distance_matrix**2/(2*epsilon[:, np.newaxis].dot(LSDMap.epsilon[np.newaxis])))

        p_vector = np.sum(kernel, axis=1)
        kernel /= np.sqrt(p_vector[:,np.newaxis].dot(LSDMap.p_vector[np.newaxis]))

        d_vector = np.sum(kernel, axis=1)
        kernel /= np.sqrt(d_vector[:,np.newaxis].dot(LSDMap.d_vector[np.newaxis]))
        evs = kernel.dot(LSDMap.evsu)/LSDMap.eigs

        # normalize eigenvectors
        evs /= np.sqrt(d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum((LSDMap.evsu/np.sqrt(LSDMap.d_vector[:,np.newaxis]))**2, axis=0))
        evs /= norm[np.newaxis,:]

        evs = np.fliplr(evs)[:,:nevs]

        return evs

    def cost_function(self, evs, target, alpha=1000):
        return alpha*np.sum((evs - target)**2, axis=1)


    def get_idxs_thread(self, comm):
        """
        get the indices of the coordinates that will be consider by each processor (use MPI)
        """
        size = comm.Get_size()
        rank = comm.Get_rank()

        npoints_per_thread = np.array([self.ntargets/size for idx in range(size)])

        for idx in range(self.ntargets%size):
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


    def save(self, args, LSDMap):
        """
        save configurations in output_file and eigenvectors in ev_output_file
        """
        format_struct_file = os.path.splitext(LSDMap.struct_filename)[1]
        struct_file_writer = writer.open(format_struct_file, pattern=LSDMap.struct_filename)
        struct_file_writer.write(self.coords, args.output_file)

        evfile_writer = writer.open('.ev')
        evfile_writer.write(self.evs, args.ev_output_file)


    def run(self):
        #initialize mpi variables
        comm = MPI.COMM_WORLD
        size = comm.Get_size()  # number of threads involved
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args()

        logging.basicConfig(filename='lift.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        # consider calling initialize function right after setting parsers

        logging.info('intializing LSDMap lift...')

        with open(args.lsdmap_file, "r") as file:
            LSDMap = pickle.load(file)

        self.initialize(args, LSDMap)

        if size > self.ntargets:
            logging.error("number of threads should be less than the number of frames")
            raise ValueError

        idxs_thread = self.get_idxs_thread(comm)
        npoints_thread = len(idxs_thread)

        targets_thread = np.array([self.targets[idx] for idx in idxs_thread])
        coords = []

        for target in targets_thread:    
            simplex, values = self.initialize_simplex(target, LSDMap)
            value_b = self.value_b
            for temperature in self.temperatures:
                simplex, values, value_b = self.neldermead(target, LSDMap, simplex, values, value_b, temperature)
            coords.append(simplex[0])
       
        self.coords = np.array(comm.allgather(np.array(coords)))
        self.coords = np.concatenate(self.coords, axis=0)
        self.evs = self.nystrom(self.coords, LSDMap, 10)

        if rank == 0:
            self.save(args, LSDMap)

        logging.info("LSDMap lifting done")
