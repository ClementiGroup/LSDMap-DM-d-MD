import os
import sys
import numpy as np
import argparse
import ConfigParser
import pickle
import random
from time import time

from lsdmap.util import reader, p_arpack
from lsdmap.util import metric as mt

from mpi4py import MPI


known_exts_struct_file = ['.gro', '.txt']

class LSDMapConfig(object):
    """ 
    LSDMapConfig(config, args)
    
    A class used to store all the parameters associated with LSDMap: coordinates of the points, 
    metric, status of the local scale, values of the local scale, weights ....

    Parameters
    ----------

    config: ConfigParser object
        used to read the configuration file

    args: argparse object
        used to get all the arguments that were specified on the command line when running the script.
    """

    def __init__(self, config, args):

        self.config = config
        self.args = args

        self.struct_file = self.struct_file()
        self.npoints = self.struct_file.nframes
        self.coords = self.struct_file.coords
        self.ndof = self.coords.shape[1]

        self.status_epsilon, self.epsilon = self.initialize_local_scale()
        self.weights = self.initialize_weights()

        self.metric = config.get('GENERAL','metric')

        try:
            self.neigs = config.getint('GENERAL','neigs')
        except:
            self.neigs = 10

    def struct_file(self):

        filename, ext = os.path.splitext(self.args.struct_file)
        if ext not in known_exts_struct_file:
            raise TypeError('Structure file extension %s unknown, please use one of ' %ext + ', '.join(known_exts_struct_file))  

        return reader.GroFile(self.args.struct_file)


    def initialize_local_scale(self):

        config = self.config
        args = self.args
        npoints = self.npoints

        known_status = ['constant', 'kneighbor', 'user']
        _mapped = {'const': 'constant', 'cst': 'constant'}

        if args.epsfile is not None:
            espfile = reader.EPSFile(epsfile)
            if npoints!=epsfile.nlines:
                raise("Number of lines in .eps file does not match the number of frames in structure file")
            status ='user'       
            value = epsfile.read()
        else:
            status = config.get('LOCALSCALE','status')
            if status in _mapped:
                status = _mapped[status]
            if not status in known_status:
                raise ValueError("local scale status should be one of "+ ', '.join(known_status))
            if status == 'kneighbor':
                value = None
                self.k = config.getint('LOCALSCALE', 'k')
            if status == 'constant':
                value = config.getfloat('LOCALSCALE', 'epsilon')
                value = value*np.ones(npoints)

            if status == 'user':
                raise NameError("status 'user' specified in configuration file but epsfile not provided. Please consider option flag -e")

        return status, value


    def initialize_weights(self):

        config = self.config
        args = self.args
        npoints = self.npoints

        if args.wfile is not None:
            wfile = reader.WFile(wfile)
            weights = wfile.read()
            if npoints!=wfile.nlines:
                raise("Number of lines in .w file does not match the number of frames in structure file")
        else:
            weights = np.ones(npoints, dtype='float')

        return weights


class LSDMap(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run LSDMap..") # message displayed when typing lsdmap -h

        # required options
        parser.add_argument("-f",
                            type=str,
                            dest="config_file",
                            required=True,
                            help ='Configuration file (input): ini')

        parser.add_argument("-c",
                            type=str,
                            dest="struct_file",
                            required=True,
                            help = 'Structure file (input): gro')

        # other options
        parser.add_argument("-o",
                            type=str,
                            dest="output_file",
                            help='LSDMap file (output)')

        parser.add_argument("-w",
            type=str,
            dest="wfile",
            help='File containing the weights of each points (input, opt.): txt')

        parser.add_argument("-e",
            type=str,
            dest="epsfile",
            help='File containing the local scale (input, opt.): txt')

        return parser

    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD
        size = comm.Get_size()  # number of threads involved
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args()

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file)

        self.LSDMapConfig = LSDMapConfig(config, args) # create an instance of LSDMapConfig class

        if size > self.LSDMapConfig.npoints:
            raise ValueError("number of threads should be less than the number of frames")

        # "split" the coordinates over all the threads.
        npoints_per_thread = np.array([self.LSDMapConfig.npoints/size for idx in range(size)])
        for idx in range(self.LSDMapConfig.npoints%size):
            npoints_per_thread[idx]+=1

        if rank == 0:
            idx_shift = 0; idx_coords_threads = []; coords_threads = []
            for idx in npoints_per_thread:
                idx_coords = range(idx_shift,idx_shift+idx)
                idx_coords_threads.append(idx_coords)
                idx_shift += idx
        else:
            idx_coords_threads = None
 
        idx_coords_thread = comm.scatter(idx_coords_threads, root=0)
        npoints_thread = len(idx_coords_thread)

        coords_thread = np.array([self.LSDMapConfig.coords[idx] for idx in idx_coords_thread])
        weights_thread = np.array([self.LSDMapConfig.weights[idx] for idx in idx_coords_thread])

        # compute the distance matrix
        if rank == 0:
            time1 = time()
        DistanceMatrix = mt.DistanceMatrix(coords_thread, self.LSDMapConfig.coords, metric=self.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix
        if rank == 0:
            time2 = time()
            print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

        # good moment to compute the kth neighbor local scale if needed
        if self.LSDMapConfig.status_epsilon == 'kneighbor':
            epsilon_thread = DistanceMatrix.neighbor_matrix(k=self.LSDMapConfig.k+1)[:,-1]
            self.LSDMapConfig.epsilon = np.hstack(comm.allgather(epsilon_thread)) #used to gather the values of epsilon

        epsilon_thread = np.array([self.LSDMapConfig.epsilon[idx] for idx in idx_coords_thread])

        p_vector_thread = np.zeros(npoints_thread, dtype='float')
        d_vector_thread = np.zeros(npoints_thread, dtype='float')

        if rank == 0: time1 = time()
        # for a detailed description of the following operations, see the paper:
            # Determination of reaction coordinates via locally scaled diffusion map
            # Mary A. Rohrdanz, Wenwei Zheng, Mauro Maggioni, and Cecilia Clementi
            # The Journal of Chemical Physics 134, 124116 (2011)

        # compute LSDMap kernel, Eq. (5) of the above paper
        kernel = np.sqrt(weights_thread[:, np.newaxis].dot(self.LSDMapConfig.weights[np.newaxis])) * \
                 np.exp(-distance_matrix**2/(2*epsilon_thread[:, np.newaxis].dot(self.LSDMapConfig.epsilon[np.newaxis])))

        p_vector_thread = np.sum(kernel, axis=1)   
        p_vector = np.hstack(comm.allgather(p_vector_thread)) # Eq. (6)
        self.p_vector = p_vector

        kernel /= np.sqrt(p_vector_thread[:,np.newaxis].dot(p_vector[np.newaxis])) # Eq. (7)
        d_vector_thread = np.sum(kernel, axis=1)
        d_vector = np.hstack(comm.allgather(d_vector_thread)) # Compute D given between Eqs. (7) and (8)
        self.d_vector = d_vector

        kernel /= np.sqrt(d_vector_thread[:,np.newaxis].dot(d_vector[np.newaxis])) # Eq (8) (slightly modified)

        # diagonalization of LSDMap kernel
        params= p_arpack._ParallelSymmetricArpackParams(comm, kernel, self.LSDMapConfig.neigs) # The last argument refers to the number of eigenvalues to be extracted
        while not params.converged:
            params.iterate()
        eigs, evs = params.extract(return_eigenvectors=True)

        # normalization
        self.evsu = np.copy(evs) # store unormalized evs
        evs /= np.sqrt(d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum(evs**2, axis=0))
        evs /= norm[np.newaxis,:]

        self.eigs = eigs; self.evs = evs

        # print eigenvectors and eigenvalues in .ev and .eg files
        if rank == 0:
            time2 = time()
            print "time estimated to diagonalize the kernel: %.3f " %(time2 - time1)
            path, ext = os.path.splitext(self.LSDMapConfig.struct_file.filename)
            np.savetxt(path + '.eg', np.fliplr(self.eigs[np.newaxis]), fmt='%9.6f')
            np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%15.7e')

            if args.output_file is None:
                args.output_file = path + '.lsdmap'
            else:
                filename, ext = os.path.splitext(args.output_file)
                if ext != '.lsdmap': args.output_file = args.output_file + '.lsdmap'

            with open(args.output_file, "w") as file:
                pickle.dump(self, file)


class LSDMapRest(object):

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

        return parser

    def initialize(self, args, LSDMap):

        filename, ext = os.path.splitext(args.struct_file)
        if ext not in known_exts_struct_file:
            raise TypeError('Structure file extension %s unknown, please use one of ' %ext + ', '.join(_known_exts))

        self.struct_file = reader.GroFile(args.struct_file)
        self.npoints = self.struct_file.nframes
        self.coords = self.struct_file.coords

        if args.wfile is not None:
            wfile = reader.WFile(wfile)
            weights = wfile.read()
            if npoints!=wfile.nlines:
                raise("Number of lines in .w file does not match the number of frames in structure file")
        else:
            weights = np.ones(self.npoints, dtype='float')
    
        self.weights = weights

        if LSDMap.LSDMapConfig.status_epsilon == 'constant':
            self.epsilon = LSDMap.LSDMapConfig.epsilon[:self.npoints]
        elif LSDMap.LSDMapConfig.status_epsilon == 'user':
            raise NotImplementedError 

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

        # "split" the coordinates over all the threads.
        npoints_per_thread = np.array([self.npoints/size for idx in range(size)])
        for idx in range(self.npoints%size):
            npoints_per_thread[idx]+=1

        if rank == 0:
            idx_shift = 0; idx_coords_threads = []; coords_threads = []
            for idx in npoints_per_thread:
                idx_coords = range(idx_shift,idx_shift+idx)
                idx_coords_threads.append(idx_coords)
                idx_shift += idx
        else:
            idx_coords_threads = None

        idx_coords_thread = comm.scatter(idx_coords_threads, root=0)
        npoints_thread = len(idx_coords_thread)

        coords_thread = np.array([self.coords[idx] for idx in idx_coords_thread])
        weights_thread = np.array([self.weights[idx] for idx in idx_coords_thread])    
    
        # compute the distance matrix
        if rank == 0:
            time1 = time()
        DistanceMatrix = mt.DistanceMatrix(coords_thread, LSDMap.LSDMapConfig.coords, metric=LSDMap.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix
        if rank == 0:
            time2 = time(); print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

        # good moment to compute the kth neighbor local scale if needed
        if LSDMap.LSDMapConfig.status_epsilon == 'kneighbor':
            epsilon_thread = DistanceMatrix.neighbor_matrix(k=LSDMap.LSDMapConfig.k+1)[:,-1]
            self.epsilon = np.hstack(comm.allgather(epsilon_thread))

        epsilon_thread = np.array([self.epsilon[idx] for idx in idx_coords_thread])
        kernel = np.sqrt(weights_thread[:, np.newaxis].dot(LSDMap.LSDMapConfig.weights[np.newaxis])) * \
            np.exp(-distance_matrix**2/(2*epsilon_thread[:, np.newaxis].dot(LSDMap.LSDMapConfig.epsilon[np.newaxis])))

        p_vector_thread = np.sum(kernel, axis=1)
        kernel /= np.sqrt(p_vector_thread[:,np.newaxis].dot(LSDMap.p_vector[np.newaxis]))

        d_vector_thread = np.sum(kernel, axis=1)
        kernel /= np.sqrt(d_vector_thread[:,np.newaxis].dot(LSDMap.d_vector[np.newaxis]))
        d_vector = np.hstack(comm.allgather(d_vector_thread))

        evs_thread = kernel.dot(LSDMap.evsu)/LSDMap.eigs
        evs = np.concatenate(comm.allgather(evs_thread), axis=0) 

        # normalization
        evs /= np.sqrt(d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum((LSDMap.evsu/np.sqrt(LSDMap.d_vector[:,np.newaxis]))**2, axis=0))
        evs /= norm[np.newaxis,:]

        # print eigenvectors and eigenvalues in .ev and .eg files
        if rank == 0:
            if args.embed is True:
                LSDMap.evs = np.vstack((LSDMap.evs, evs))
                evs = LSDMap.evs
                with open(args.lsdmap_file, "w") as file:
                    pickle.dump(LSDMap, file)

            path, ext = os.path.splitext(args.struct_file)
            np.savetxt(path + '.ev', np.fliplr(evs), fmt='%15.7e')


class LSDMapLift(object):

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

    def initialize(self, args, LSDMap):

        evfile = reader.EvFile(args.evfile)
        targets = evfile.get_evs()

        self.ntargets = targets.shape[0]
        self.nevs = min(targets.shape[1], args.nevs)
        self.targets = targets[:,:self.nevs]
        self.ndof = LSDMap.LSDMapConfig.ndof

        if not all(weight == 1.0 for weight in LSDMap.LSDMapConfig.weights):
            raise ValueError("file %s contains weighted version of LSDMap; cannot lift points from weighted LSDMap")

    def initialize_simplex(self, target, LSDMap):

        evs = np.fliplr(LSDMap.evs)[:,:self.nevs]
        diff = np.sqrt(np.sum((evs - target)**2, axis=1))

        simplex_idxs = np.argsort(diff)[:self.ndof+1]
        simplex = np.array([LSDMap.LSDMapConfig.coords[idx] for idx in simplex_idxs])
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
        DistanceMatrix = mt.DistanceMatrix(coords, LSDMap.LSDMapConfig.coords, metric=LSDMap.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix

        if LSDMap.LSDMapConfig.status_epsilon == 'constant':
            epsilon = LSDMap.LSDMapConfig.epsilon[:npoints]
        elif LSDMap.LSDMapConfig.status_epsilon == 'kneighbor':
            epsilon = DistanceMatrix.neighbor_matrix(k=LSDMap.LSDMapConfig.k+1)[:,-1]

        # compute LSDMap kernel
        kernel = np.exp(-distance_matrix**2/(2*epsilon[:, np.newaxis].dot(LSDMap.LSDMapConfig.epsilon[np.newaxis])))

        p_vector = np.sum(kernel, axis=1)
        kernel /= np.sqrt(p_vector[:,np.newaxis].dot(LSDMap.p_vector[np.newaxis]))

        d_vector = np.sum(kernel, axis=1)
        kernel /= np.sqrt(d_vector[:,np.newaxis].dot(LSDMap.d_vector[np.newaxis]))
        evs = kernel.dot(LSDMap.evsu)/LSDMap.eigs

        # normalization
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

