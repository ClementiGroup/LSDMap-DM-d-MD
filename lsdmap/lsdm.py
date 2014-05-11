import os, sys
import numpy as np
import argparse
import ConfigParser
import pickle

from time import time

from lsdmap.util import reader, p_arpack
from lsdmap.util import metric as mt

from mpi4py import MPI


_known_exts = ['.gro', '.txt']


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

        self.status_epsilon, self.epsilon = self.initialize_local_scale()
        self.weights = self.initialize_weights()

        self.metric = config.get('GENERAL','metric')

        try: self.neigs = config.getint('GENERAL','neigs')
        except: self.neigs = 10

    def struct_file(self):

        filename, ext = os.path.splitext(self.args.struct_file)
        if ext not in _known_exts: raise TypeError('Structure file extension %s unknown, please use one of ' %ext + ', '.join(_known_exts))  

        return reader.GroFile(self.args.struct_file)


    def initialize_local_scale(self):

        config = self.config
        args = self.args
        npoints = self.npoints

        known_status = ['constant', 'kneighbor', 'user']
        _mapped = {'const': 'constant', 'cst': 'constant'}

        if args.epsfile is not None:
            espfile = reader.EPSFile(epsfile)
            if npoints!=epsfile.nlines: raise("Number of lines in .eps file does not match the number of frames in structure file")
            status ='user'       
            value = epsfile.read()
        else:
            status = config.get('LOCALSCALE','status')
            if status in _mapped: status = _mapped[status]
            if not status in known_status: raise ValueError("local scale status should be one of "+ ', '.join(known_status))
      
            if status == 'kneighbor':
                value = None
                self.k = config.getint('LOCALSCALE', 'k')

            if status == 'constant':
                value = config.getfloat('LOCALSCALE', 'epsilon')
                value = value*np.ones(npoints)

            if status == 'user': raise NameError("status 'user' specified in configuration file but epsfile not provided. Please consider option flag -e")

        return status, value


    def initialize_weights(self):

        config = self.config
        args = self.args
        npoints = self.npoints

        if args.wfile is not None:
            wfile = reader.WFile(wfile)
            weights = wfile.read()
            if npoints!=wfile.nlines: raise("Number of lines in .w file does not match the number of frames in structure file")
        else: weights = np.ones(npoints, dtype='float')

        return weights


class LSDMap(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run LSDMap..") # message displayed when typing lsdmap -h

        # required options
        parser.add_argument("-f", type=str, dest="config_file", required=True, help ='Configuration file (input): ini')
        parser.add_argument("-c", type=str, dest="struct_file", required=True, help = 'Structure file (input): gro')

        # other options
        parser.add_argument("-o", type=str, dest="output_file", help='LSDMap file (output)')
        parser.add_argument("-w", type=str, dest="wfile", help='File containing the weights of each points (input, opt.): txt')
        parser.add_argument("-e", type=str, dest="epsfile", help='File containing the local scale (input, opt.): txt')

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

        if size > self.LSDMapConfig.npoints: raise ValueError("number of threads should be less than the number of frames")

        # "split" the coordinates over all the threads.
        npoints_per_thread = np.array([self.LSDMapConfig.npoints/size for idx in range(size)])
        for idx in range(self.LSDMapConfig.npoints%size): npoints_per_thread[idx]+=1

        if rank == 0:
            idx_shift = 0; idx_coords_threads = []; coords_threads = []
            for idx in npoints_per_thread:
                idx_coords = range(idx_shift,idx_shift+idx)
                idx_coords_threads.append(idx_coords)
                idx_shift += idx
        else: idx_coords_threads = None
 
        idx_coords_thread = comm.scatter(idx_coords_threads, root=0)
        npoints_thread = len(idx_coords_thread)

        coords_thread = np.array([self.LSDMapConfig.coords[idx] for idx in idx_coords_thread])
        weights_thread = np.array([self.LSDMapConfig.weights[idx] for idx in idx_coords_thread])

        # compute the distance matrix
        if rank == 0: time1 = time()
        DistanceMatrix = mt.DistanceMatrix(coords_thread, self.LSDMapConfig.coords, metric=self.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix
        if rank == 0: time2 = time(); print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

        # good moment to compute the kth neighbor local scale if needed
        if self.LSDMapConfig.status_epsilon == 'kneighbor':
            epsilon_thread = DistanceMatrix.neighbor_matrix(k=self.LSDMapConfig.k+1)[:,-1]
            self.LSDMapConfig.epsilon = np.hstack(comm.allgather(epsilon_thread)) #used to gather the values of epsilon

        epsilon_thread = np.array([self.LSDMapConfig.epsilon[idx] for idx in idx_coords_thread])

        p_vector_thread = np.zeros(npoints_thread, dtype='float')
        d_vector_thread = np.zeros(npoints_thread, dtype='float')

        if rank == 0: time1 = time()
        ### for a detailed description of the following operations, see the paper:
            ### Determination of reaction coordinates via locally scaled diffusion map
            ### Mary A. Rohrdanz, Wenwei Zheng, Mauro Maggioni, and Cecilia Clementi
            ### The Journal of Chemical Physics 134, 124116 (2011)

        # compute LSDMap kernel, Eq. (5) of the above paper
        kernel = np.sqrt(weights_thread[:, np.newaxis].dot(self.LSDMapConfig.weights[np.newaxis]))* \
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
        params= p_arpack._ParallelSymmetricArpackParams(comm, kernel, self.LSDMapConfig.neigs) # last argument refers to the number of eigenvalues to be extracted
        while not params.converged:
            params.iterate()
        eigs, evs = params.extract(return_eigenvectors=True)

        # normalization
        self.evsu = np.copy(evs) #store unormalized evs
        evs /= np.sqrt(d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum(evs**2, axis=0))
        evs /= norm[np.newaxis,:]

        self.eigs = eigs; self.evs = evs

        # print eigenvectors and eigenvalues in .ev and .eg files
        if rank == 0:
            time2 = time(); print "time estimated to diagonalize the kernel: %.3f " %(time2 - time1)
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
        parser.add_argument("-c", type=str, dest="struct_file", required=True, help = 'Structure file (input): gro')
        parser.add_argument("-s", type=str, dest="lsdmap_file", required=True, help ='LSDMap file (input): lsdmap')

        # other options
        parser.add_argument("-e", type=str, dest="epsfile", help='File containing the local scales of each points (input, opt.): txt')
        parser.add_argument("-w", type=str, dest="wfile", help='File containing the weights of each points (input, opt.): txt')
        parser.add_argument("--embed", action='store_true', help='embed LSDMap coordinates to old ones')

        return parser

    def initialize(self, args, LSDMap):

        filename, ext = os.path.splitext(args.struct_file)
        if ext not in _known_exts: raise TypeError('Structure file extension %s unknown, please use one of ' %ext + ', '.join(_known_exts))

        self.struct_file = reader.GroFile(args.struct_file)
        self.npoints = self.struct_file.nframes
        self.coords = self.struct_file.coords

        if args.wfile is not None:
            wfile = reader.WFile(wfile)
            weights = wfile.read()
            if npoints!=wfile.nlines: raise("Number of lines in .w file does not match the number of frames in structure file")
        else: weights = np.ones(self.npoints, dtype='float')
    
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

        with open(args.lsdmap_file, "r") as file: LSDMap = pickle.load(file)
        self.initialize(args, LSDMap)

        if size > self.npoints: raise ValueError("number of threads should be less than the number of configurations")

        # "split" the coordinates over all the threads.
        npoints_per_thread = np.array([self.npoints/size for idx in range(size)])
        for idx in range(self.npoints%size): npoints_per_thread[idx]+=1

        if rank == 0:
            idx_shift = 0; idx_coords_threads = []; coords_threads = []
            for idx in npoints_per_thread:
                idx_coords = range(idx_shift,idx_shift+idx)
                idx_coords_threads.append(idx_coords)
                idx_shift += idx
        else: idx_coords_threads = None

        idx_coords_thread = comm.scatter(idx_coords_threads, root=0)
        npoints_thread = len(idx_coords_thread)

        coords_thread = np.array([self.coords[idx] for idx in idx_coords_thread])
        weights_thread = np.array([self.weights[idx] for idx in idx_coords_thread])    
    
        # compute the distance matrix
        if rank == 0: time1 = time()
        DistanceMatrix = mt.DistanceMatrix(self.coords, LSDMap.LSDMapConfig.coords, metric=LSDMap.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix
        if rank == 0: time2 = time(); print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

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

        evs_thread = kernel.dot(LSDMap.evsu)/LSDMap.eigs
        evs = np.hstack(comm.allgather(evs_thread)) 

        # normalization
        evs /= np.sqrt(d_vector_thread[:,np.newaxis])
        norm = np.sqrt(np.sum((LSDMap.evsu/np.sqrt(LSDMap.d_vector[:,np.newaxis]))**2, axis=0))
        evs /= norm[np.newaxis,:]
        self.evs = evs   

        # print eigenvectors and eigenvalues in .ev and .eg files
        if rank == 0:
            path, ext = os.path.splitext(args.struct_file)
            np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%15.7e')


class LSDMapLift(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Lifting from diffusion map to physical variables...") # message displayed when typing llsdmap -h

        # required options
        parser.add_argument("-f", type=str, dest="evfile", required=True, help = 'File containing the eigenvectors components used to lift (input) : ev')
        parser.add_argument("-s", type=str, dest="lsdmap_file", required=True, help ='LSDMap file (input): lsdmap')

        return parser

    def initialize(self, args, LSDMap):

        self.evfile = reader.EvFile(args.evfile)
        self.evs = self.evfile.get_evs()
        self.npoints = self.evs.shape[0]
        self.nevs = self.evs.shape[1]

        if not all(weight == 1.0 for weight in LSDMap.LSDMapConfig.weights):
            raise ValueError("cannot lift points from weighted LSDMap")

        if LSDMap.LSDMapConfig.status_epsilon == 'constant':
            self.epsilon = LSDMap.LSDMapConfig.epsilon[:self.npoints]
        elif LSDMap.LSDMapConfig.status_epsilon == 'user':
            raise ValueError("cannot lift points from LSDMap computed via user-defined local scales")

        if LSDMap.LSDMapConfig.status_epsilon == 'constant':
            self.epsilon = LSDMap.LSDMapConfig.epsilon[:self.npoints]
        else: raise NotImplementedError('LSDMap lift is only implemented when the local scale is constant')


    def rest(self, coords, epsilon, LSDMap):

        # compute the distance matrix
        DistanceMatrix = mt.DistanceMatrix(coords, LSDMap.LSDMapConfig.coords, metric=LSDMap.LSDMapConfig.metric)
        distance_matrix = DistanceMatrix.distance_matrix

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

        return evs

    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD
        size = comm.Get_size()  # number of threads involved
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args()

        with open(args.lsdmap_file, "r") as file: LSDMap = pickle.load(file)
        self.initialize(args, LSDMap)

        coords = LSDMap.LSDMapConfig.coords[0]
        epsilon = np.array([self.epsilon[0]])
        evs = self.rest(coords, epsilon, LSDMap)
