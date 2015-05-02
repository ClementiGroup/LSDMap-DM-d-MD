import os
import sys
import numpy as np
import argparse
import ConfigParser
import pickle
import random
import logging

from lsdmap.rw import reader
from lsdmap.util import metric as mt

from scipy.sparse.linalg.eigen.arpack.arpack import eigsh
   
class LSDMap(object):

    def initialize(self, config, args):

        self.config = config
        self.args = args

        struct_file = reader.open(args.struct_file)
        self.struct_filename = struct_file.filename
        self.npoints = struct_file.nlines

        self.coords = struct_file.readlines()
        logging.info('input coordinates loaded')

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

        if args.wfile is None:
            self.weights = np.ones(self.npoints, dtype='float')
        else:
            if os.path.isfile(args.wfile):
                wfile = reader.open(args.wfile)
                weights = wfile.readlines()
                self.weights = self.npoints*weights/np.sum(weights)
            else:
                logging.warning('.w file does not exist, set weights to 1.0.')
                self.weights = np.ones(self.npoints, dtype='float')

    def initialize_metric(self):

        _known_prms = ['r0']

        config = self.config
        self.metric = config.get('LSDMAP', 'metric')

        self.metric_prms = {}
        for prm in _known_prms:
            try:
                self.metric_prms[prm] = config.getfloat('LSDMAP', prm)
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
        parser.add_argument("-w",
            type=str,
            dest="wfile",
            help='File containing the weights of every point in a row (input, opt.): .w')

        parser.add_argument("-e",
            type=str,
            dest="epsfile",
            help='File containing the local scales of every point in a row (input, opt.): .eps')

        return parser

    def compute_kernel(self, distance_matrix):

        # for a detailed description of the following operations, see the paper:
            # Determination of reaction coordinates via locally scaled diffusion map
            # Mary A. Rohrdanz, Wenwei Zheng, Mauro Maggioni, and Cecilia Clementi
            # The Journal of Chemical Physics 134, 124116 (2011)

        # compute LSDMap kernel, Eq. (5) of the above paper
        kernel = np.sqrt(self.weights[:, np.newaxis].dot(self.weights[np.newaxis])) * \
                 np.exp(-distance_matrix**2/(2*self.epsilon[:, np.newaxis].dot(self.epsilon[np.newaxis])))

        p_vector = np.sum(kernel, axis=1) # Eq. (6)
        self.p_vector = p_vector

        kernel /= np.sqrt(p_vector[:,np.newaxis].dot(p_vector[np.newaxis])) # Eq. (7)
        d_vector = np.sum(kernel, axis=1) # Compute D given between Eqs. (7) and (8)
        self.d_vector = d_vector

        kernel /= np.sqrt(d_vector[:,np.newaxis].dot(d_vector[np.newaxis])) # Eq (8) (slightly modified)

        return kernel

    def save(self, config, args):
        """
        save LSDMap object in .lsdmap file and eigenvalues/eigenvectors in .eg/.ev files
        """

        if isinstance(self.struct_filename, list):
            struct_filename = self.struct_filename[0]
        else:
            struct_filename = self.struct_filename

        path, ext = os.path.splitext(struct_filename)
        np.savetxt(path + '.eg', np.fliplr(self.eigs[np.newaxis]), fmt='%9.6f')
        np.savetxt(path + '.ev', np.fliplr(self.evs), fmt='%.18e')

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args() # set argument parser

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file) # set config file parser

        logging.basicConfig(filename='lsdmap.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        logging.info('intializing LSDMap...')
        self.initialize(config, args)
        logging.info('LSDMap initialized')

        # compute the distance matrix
        DistanceMatrix = mt.DistanceMatrix(self.coords, self.coords, metric=self.metric, metric_prms=self.metric_prms)

        distance_matrix = DistanceMatrix.distance_matrix
        logging.info("distance matrix computed")

        # compute kth neighbor local scales if needed
        if self.status_epsilon in ['kneighbor', 'kneighbor_mean']:
            #epsilon = []
            self.epsilon = np.zeros(self.npoints,float)
            idx_neighbor_matrix = DistanceMatrix.idx_neighbor_matrix()
            for idx, line in enumerate(idx_neighbor_matrix):
                cum_weight = 0
                for jdx in line[1:]:
                    cum_weight += self.weights[jdx]
                    if cum_weight >= self.k:
                        break
                #epsilon.append(distance_matrix[idx,jdx])
                self.epsilon[idx] = distance_matrix[idx,jdx]

            if self.status_epsilon == 'kneighbor_mean':
                mean_value_epsilon = np.mean(self.epsilon) # compute the mean value of the local scales
                self.epsilon = mean_value_epsilon * np.ones(self.npoints)  # and set it as the new constant local scale

            logging.info("kneighbor local scales computed")

        # compute kernel
        kernel = self.compute_kernel(distance_matrix)

        # diagonalize kernel
        eigs, evs = eigsh(kernel, k=self.neigs)

        # normalize eigenvectors
        self.evsu = np.copy(evs) # store unormalized eigenvectors
        evs /= np.sqrt(self.d_vector[:,np.newaxis])
        norm = np.sqrt(np.sum(evs**2, axis=0))
        evs /= norm[np.newaxis,:]
        
        # store eigenvalues/eigenvectors
        self.eigs = eigs
        self.evs = evs

        logging.info("kernel diagonalized")
        self.save(config, args)

        logging.info("Eigenvalues/eigenvectors saved (.eg/.ev files)")
        logging.info("LSDMap computation done")

if __name__ == '__main__':
    LSDMap().run()
