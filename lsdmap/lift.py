import os
import sys
import argparse
import ConfigParser
import pickle
import random
from time import time

import numpy as np
from mpi4py import MPI

from lsdmap.reader import reader
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
           help='Temperature is decreased by multiplying it by (1-ratio) at each step of simulated annealing (opt., default=0.1) ',
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

        coords = self.ssimplex * fac1 - self.simplex[self.idx_hi,:] * fac2
        evs = self.nystrom(coords[np.newaxis], LSDMap)
        value = self.cost_function(evs, target)

        if value[0] <= self.value_b:
            self.coords_b = coords
            self.value_b = value[0]

        value_rand = value[0] - self.temperature*np.log(np.random.rand());

        if value_rand < self.value_hi:
            self.values[self.idx_hi] = value[0]
            self.value_hi = value_rand
            self.ssimplex += coords - self.simplex[self.idx_hi,:]
            self.simplex[self.idx_hi,:] = coords

        return value_rand

    def neldermead(self, target, LSDMap, simplex, values, value_b, temperature, niters=1000):
        """
        Nelder Mead method used to find the minimum of cost function at a given temperature
        """

        self.simplex = simplex
        self.values = values
        self.temperature = temperature
        self.value_b = value_b

        self.ssimplex = np.sum(self.simplex, axis=0) # ssimplex has self.ndof components

        while True:

            values_rand = values - self.temperature*np.log(np.random.rand(self.ndof+1))
            values_rand_idxs = np.argsort(values_rand)

            self.idx_l = values_rand_idxs[0]
            self.idx_hi = values_rand_idxs[-1]
            self.idx_nhi = values_rand_idxs[-2]

            self.value_l = values_rand[self.idx_l]
            self.value_hi = values_rand[self.idx_hi]
            self.value_nhi = values_rand[self.idx_nhi]

            #print rtol
            if niters < 0:
                break

            niters -= 2
            new_value = self.generate_new_value(target, LSDMap, -1.0)

            if new_value <= self.value_l:
                new_value = self.generate_new_value(target, LSDMap, 2.0)
            elif new_value >= self.value_nhi:
                value_save = self.value_hi
                new_value = self.generate_new_value(target, LSDMap, 0.5)
                if new_value >= value_save:
                    for idx in xrange(self.ndof+1):
                        if idx != self.idx_l:
                            self.ssimplex = 0.5*(self.simplex[idx,:]+self.simplex[self.idx_l,:])
                            self.simplex[idx,:] = self.ssimplex
                            evs = self.nystrom(self.ssimplex[np.newaxis], LSDMap)
                            self.values[idx] = self.cost_function(evs, target)
                    niters -= self.ndof
                    self.ssimplex = np.sum(self.simplex, axis=0)
            else:
                niters +=1

        tmp = self.values[0]
        self.values[0] = self.values[self.idx_l]
        self.values[self.idx_l] = tmp

        tmp = self.simplex[0,:]
        self.simplex[0,:] = self.simplex[self.idx_l,:]
        self.simplex[self.idx_l] = tmp

        return self.simplex, self.values, self.value_b

    def nystrom(self, coords, LSDMap):
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

        evs = np.fliplr(evs)[:,:self.nevs]

        return evs

    def cost_function(self, evs, target, alpha=1000):
        return alpha*np.sum((evs - target)**2, axis=1)

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

        temperature_max = args.temperature_max
        ratio = args.ratio
        nsteps = args.nsteps
        niters = args.niters
        value_b = args.value_b

        temperatures = [temperature_max]
        for idx in xrange(nsteps-1):
            temperatures.append((1-ratio)*temperatures[idx])
        print "temperature max: %.3f" %temperatures[0]
        print "temperature min: %.3f" %temperatures[-1]

        target = self.targets[0]

        simplex, values = self.initialize_simplex(target, LSDMap)
        for temperature in temperatures:
            simplex, values, value_b = self.neldermead(target, LSDMap, simplex, values, value_b, temperature, niters=niters)

        evs = self.nystrom(simplex[0,:][np.newaxis], LSDMap)
        print target, evs[0]
