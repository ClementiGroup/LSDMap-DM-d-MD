import os, sys
import numpy as np
import argparse
import ConfigParser

from time import time

from lsdmap.util import reader, p_arpack
from lsdmap.util import metric as mt

from mpi4py import MPI


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
        self.ncoords = self.struct_file.nframes
        self.coords = self.struct_file.coords

        self.status_epsilon, self.epsilon = self.initialize_local_scale()
        self.weights = self.initialize_weights()

        self.metric = config.get('GENERAL','metric')


    def struct_file(self):

        _known_exts = ['.gro', '.txt']
        filename, ext = os.path.splitext(self.args.struct_file)
        if ext not in _known_exts: raise TypeError('Structure file extension %s unknown, please use one of ' %ext + ', '.join(_known_exts))  

        return reader.GroFile(self.args.struct_file)


    def initialize_local_scale(self):

        config = self.config
        args = self.args
        ncoords = self.ncoords

        known_status = ['constant', 'kneighbor', 'user']
        _mapped = {'const': 'constant', 'cst': 'constant'}

        if args.epsfile is not None:
            espfile = reader.EPSFile(epsfile)
            if ncoords!=epsfile.nlines: raise("Number of lines in .eps file does not match the number of frames in structure file")
            status ='user'       
            value = epsfile.read()
        else:
            status = config.get('LOCALSCALE','status')
            if status in _mapped: status = _mapped[status]
            if not status in known_status: raise ValueError("local scale status should be one of "+ ', '.join(known_status))
      
            if status == 'kneighbor':
                value = None

            if status == 'constant':
                value = config.getfloat('LOCALSCALE', 'epsilon')
                value = value*np.ones(ncoords)

            if status == 'user': raise NameError("status 'user' specified in configuration file but epsfile not provided. Please consider option flag -e")

        return status, value


    def initialize_weights(self):

        config = self.config
        args = self.args
        ncoords = self.ncoords

        if args.wfile is not None:
            wfile = reader.WFile(wfile)
            weights = wfile.read()
            if ncoords!=wfile.nlines: raise("Number of lines in .w file does not match the number of frames in structure file")
        else: weights = np.ones(ncoords, dtype='float')

        return weights


def run_lsdmap():

    global LSDMapConfig

    #initialize mpi variables
    comm = MPI.COMM_WORLD
    size = comm.Get_size()  # number of threads involved
    rank = comm.Get_rank()  # number of the current thread 

    parser = argparse.ArgumentParser(description="Run LSDMap..") # message displayed when typing lsdmap -h

    # required options
    parser.add_argument("-f", type=str, dest="config_file", required=True, help ='Configuration file (input): ini')
    parser.add_argument("-c", type=str, dest="struct_file", required=True, help = 'Structure file (input): gro')

    # other options
    parser.add_argument("-w", type=str, dest="wfile", help='File containing the weights of each points (input, opt.): txt')
    parser.add_argument("-e", type=str, dest="epsfile", help='File containing the local scale (input, opt.): txt')

    args = parser.parse_args() # the values of each parameter specified on the command line
                               # can be recovered via args.<dest> where <dest> is the value of dest
                               # specified above when using parser.add_argument(...)

    config = ConfigParser.SafeConfigParser()
    config.read(args.config_file)

    LSDMapConfig = LSDMapConfig(config, args) # create an instance of LSDMapConfig class

    if size > LSDMapConfig.ncoords: raise ValueError("number of threads should be less than the number of frames")

    # "split" the coordinates over all the threads.
    ncoords_per_thread = np.array([LSDMapConfig.ncoords/size for idx in range(size)])
    for idx in range(LSDMapConfig.ncoords%size): ncoords_per_thread[idx]+=1

    if rank == 0:
        idx_shift = 0; idx_coords_threads = []; coords_threads = []
        for idx in ncoords_per_thread:
            idx_coords = range(idx_shift,idx_shift+idx)
            idx_coords_threads.append(idx_coords)
            idx_shift += idx
    else: idx_coords_threads = None
 
    idx_coords_thread = comm.scatter(idx_coords_threads, root=0)
    ncoords_thread = len(idx_coords_thread)

    coords_thread = np.array([LSDMapConfig.coords[idx] for idx in idx_coords_thread])
    weights_thread = np.array([LSDMapConfig.weights[idx] for idx in idx_coords_thread])

    # compute the distance matrix
    if rank == 0: time1 = time()
    DistanceMatrix = mt.DistanceMatrix(coords_thread, LSDMapConfig.coords, metric=LSDMapConfig.metric)
    distance_matrix = DistanceMatrix.distance_matrix
    if rank == 0: time2 = time(); print "time estimated to compute the distance matrix: %.3f " %(time2 - time1)

    # good moment to compute the kth neighbor local scale if needed
    if LSDMapConfig.status_epsilon == 'kneighbor':
        k = config.getint('LOCALSCALE', 'k')
        epsilon_thread = DistanceMatrix.neighbor_matrix(k=k+1)[:,-1]
        LSDMapConfig.epsilon = np.hstack(comm.allgather(epsilon_thread)) #used to gather the values of epsilon

    epsilon_thread = np.array([LSDMapConfig.epsilon[idx] for idx in idx_coords_thread])

    p_vector_thread = np.zeros(ncoords_thread, dtype='float')
    d_vector_thread = np.zeros(ncoords_thread, dtype='float')

    if rank == 0: time1 = time()

    # for a detailed description of the following operations, see the paper:
        # Determination of reaction coordinates via locally scaled diffusion map
        # Mary A. Rohrdanz, Wenwei Zheng, Mauro Maggioni, and Cecilia Clementi
        # The Journal of Chemical Physics 134, 124116 (2011)

    # compute LSDMap kernel, Eq. (5) of the paper mentioned above
    kernel = np.sqrt(weights_thread[:, np.newaxis].dot(LSDMapConfig.weights[np.newaxis]))* \
        np.exp(-distance_matrix**2/(2*epsilon_thread[:, np.newaxis].dot(LSDMapConfig.epsilon[np.newaxis])))

    p_vector_thread = np.sum(kernel, axis=1)   
    p_vector = np.hstack(comm.allgather(p_vector_thread)) # Eq. (6)

    kernel /= np.sqrt(p_vector_thread[:,np.newaxis].dot(p_vector[np.newaxis])) # Eq. (7)
    d_vector_thread = np.sum(kernel, axis=1)
    d_vector = np.hstack(comm.allgather(d_vector_thread)) # Compute D given between Eqs. (7) and (8)

    kernel /= np.sqrt(d_vector_thread[:,np.newaxis].dot(d_vector[np.newaxis])) # Eq (8) (slightly modified)

    # diagonalization of LSDMap kernel
    params= p_arpack._ParallelSymmetricArpackParams(comm, kernel, 10) #10 indicates that 10 eigenvalues will be extracted
    while not params.converged:
        params.iterate()
    eigs, evs = params.extract(return_eigenvectors=True)

    # normalization
    evs /= np.sqrt(d_vector[:,np.newaxis])
    norm = np.sqrt(np.sum(evs**2, axis=0))
    evs /= norm[np.newaxis,:]

    # print eigenvectors and eigenvalues in .ev and .eg files
    if rank == 0:
        time2 = time(); print "time estimated to diagonalize the kernel: %.3f " %(time2 - time1)
        filename, ext = os.path.splitext(LSDMapConfig.struct_file.filename)
        np.savetxt(filename + '.eg', np.fliplr(eigs[np.newaxis]), fmt='%9.6f')
        np.savetxt(filename + '.ev', np.fliplr(evs), fmt='%15.7e')
