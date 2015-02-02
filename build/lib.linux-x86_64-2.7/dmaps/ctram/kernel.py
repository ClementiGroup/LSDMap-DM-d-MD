import os
import sys
import imp
import ConfigParser
import argparse
from math import floor
import logging
import subprocess
import shutil
import numpy as np
from mpi4py import MPI
from lsdmap.rw import reader
from lsdmap.rw import writer
from dmaps.tools import tools
from dmaps.ctram import wrapper
from dmaps.ctram import ctramfe

class CTRAM4DMapSConfig(object):

    def __init__(self, settings, config):

        self.iter = settings.iter # current iteration
        self.nreplicas = settings.nreplicas
        self.nstride = config.nstride
        self.nstates = min(self.iter+1, config.nstates_ctram)
        self.ntrajs = self.nreplicas*self.nstates
        self.ntau = config.ntau_ctram

        self.iters = range(self.iter-self.nstates+1, self.iter+1)

        idxs_coords = []
        for idx in xrange(self.nreplicas):
            for jdx in xrange(0, self.nstride, self.ntau):
                idxs_coords.append(jdx+idx*self.nstride)
        self.nsamples_per_state = len(idxs_coords)
        self.nsamples = self.nsamples_per_state*self.nstates

        logging.info("cTRAM: preparing input files for embedding")
        logging.info("cTRAM: merging .gro files and copying files fit.gro and fit.ev")
        subprocess.check_call('rm -rf ctram; mkdir ctram', shell=True)
        for kdx, iter in enumerate(self.iters):
            # get configurations to embed
            if kdx == self.nstates-1:
                dir = "."
            else:
                dir = "iter%i"%iter
                if not os.path.isdir(dir):
                    raise IOError('folder ' + dir + ' does not exist')
            gr = reader.open(dir + '/fit/embed.gro')
            coords = gr.readlines()
            gr.close()
            gw = writer.open('.gro', pattern=dir + '/fit/embed.gro')
            gw.write(coords[idxs_coords], 'ctram/embed.gro', mode='a')
            # copy fit.gro and fit.ev files of old maps if iter > 0
            if iter > 0:
                olddir = "iter%i"%(iter-1)
                if os.path.exists(olddir):
                    shutil.copyfile(olddir + '/fit/fit.gro', 'ctram/fit%i.gro'%iter)
                    shutil.copyfile(olddir + '/fit/fit.ev', 'ctram/fit%i.ev'%iter)
                else:
                    raise IOError('folder ' + olddir + ' does not exist')

        # copy fit.gro and fit.ev files of the current map
        shutil.copyfile("fit/fit.gro", "ctram/fit.gro")
        shutil.copyfile("fit/fit.ev", "ctram/fit.ev")

        logging.info("cTRAM: writing script for the fitting procedure")
        self.write_ctram_script("ctram/ctram.sh", settings, config)
      
        logging.info("cTRAM: all the files were successfully written")

    def write_ctram_script(self, filename, settings, config):

        # directories 
        iters_embed = " ".join([str(idx) for idx in self.iters if idx > 0])
        dcs_options = '--dc ' + ' '.join([str(num+1) for num in xrange(config.ndcs)])
        cores = settings.cores
        # write PBS script
        with open(filename, 'w') as file:
            script = """#!/bin/bash

# this script is used to compute the DCs needed to apply cTRAM

# embed points on the new map
mpiexec -n %(cores)s rbffit -f config.ini -c fit.gro -v fit.ev --embed embed.gro %(dcs_options)s
mv fit.embed confall.ev.embed

# embed points on previous maps
for idx in %(iters_embed)s; do
    mpiexec -n %(cores)s rbffit -f config.ini -c fit${idx}.gro -v fit${idx}.ev --embed embed.gro %(dcs_options)s
    mv fit.embed confall.ev.embed${idx}.old
done
rm -rf fit.w fit.sig
rm -rf fitting.log
rm -rf embed.gro
        """ %locals()
            file.write(script)

    def embed_configurations(self, umgr, settings):

        import radical.pilot

        logging.info('cTRAM: starting embedding...')
        print "Setting up..."
        cu = radical.pilot.ComputeUnitDescription()

        iters_embed = [str(idx) for idx in self.iters if idx > 0]
        
        cu.input_staging = [settings.inifile, 'ctram/ctram.sh', 'ctram/fit.gro'] + ['ctram/fit'+idx+'.gro' for idx in iters_embed] + \
            ['ctram/fit.ev'] + ['ctram/fit'+idx+'.ev' for idx in iters_embed] + ['ctram/embed.gro']
        cu.executable = 'bash ctram.sh' 
        cu.output_staging = ['confall.ev.embed > ctram/confall.ev.embed'] + ['confall.ev.embed'+idx+'.old > ctram/confall.ev.embed'+idx+'.old' for idx in iters_embed]
        cu.cleanup = True
        cu.cores = settings.cores

        unit = umgr.submit_units(cu)
        unit.wait()

        logging.info('cTRAM: embeddding done')
        print 'Setup time:', (unit.stop_time - unit.start_time).total_seconds()

    def prepare_inputs(self, config):

        nbins = config.nbins_ctram
        nsamples = self.nsamples
        nstates = self.nstates
        ntrajs = self.ntrajs
        nstride = self.nstride
        nreplicas = self.nreplicas
        nsamples_per_state = self.nsamples_per_state
        iters = self.iters
        kT = config.kT

        dcs = np.loadtxt('ctram/confall.ev.embed')
        ndcs = config.ndcs
        if ndcs == 1:
            dcs = dcs[:,np.newaxis]
        self.dcs = dcs

        fefrac = config.fefrac

        logging.info("cTRAM: construct histogram from points embedded on the latest map")
        if nbins is None:
            nnebins = int(5*nsamples**(1./2))
            nbins_min = int(nsamples**(1./(2*ndcs)))
            nbins_max = int(nsamples**(1./2))
            # construct the histogram with an optimal number of bins
            bins, grid = tools.do_grid_optimized(dcs, nnebins, nbins_min, nbins_max)
        else:
            bins, grid = tools.do_grid(dcs, nbins, nextrabins=0)

        steps = bins[1,:] - bins[0,:]
        nbins = bins.shape[0]

        traj_total_mem = np.zeros((nsamples, nstates+1))
        nclusters = 0
        for index_bin, bin in tools.nonempty_bins(grid):
            nclusters += 1

        self.nclusters = nclusters
        self.grid = grid
        self.bins = bins
        self.steps = steps
        logging.info("cTRAM: the histrogram has been built using %i bins per dim (%i non-empty bins)"%(nbins, nclusters))

        N_mem = np.zeros((nclusters, nstates))
        nebins = np.zeros(nclusters)
        logging.info("cTRAM: constructing the first column of traj_total_mem and N_mem")
        for idx, (index_bin, bin) in enumerate(tools.nonempty_bins(grid)):
            index_bin_s = 0
            for jdx in range(ndcs):
                index_bin_s += index_bin[jdx]*nbins**(ndcs-jdx-1)
            # store the indices corresponding to nonempty bins
            nebins[idx] = index_bin_s
            for index_sample in bin:
                traj_total_mem[index_sample, 0] = idx
                # find the thermodynamic state associated with the sample
                N_mem[idx, index_sample/nsamples_per_state] += 1

        self.N_mem = N_mem
        logging.info("cTRAM: constructing the rest of traj_total_mem")
        if nstates > 1:
            for kdx in xrange(nstates):
                # check if the thermodynamic state corresponds to iter0 or not
                if iters[kdx] > 0:
                    # load dc values
                    olddcs = np.loadtxt('ctram/confall.ev.embed%i.old'%iters[kdx])
                    free_energy_k = ctramfe.get_free_energy_ctram(olddcs, "iter%i"%(iters[kdx]-1) + '/fe')
                    traj_total_mem[:, kdx+1] = -fefrac*free_energy_k/kT
                else: # the state corresponds to normal MD simulations
                    for idx in range(nsamples):
                        traj_total_mem[idx, kdx+1] = 0.0
        else:
            # for configurations associated with the kth thermodynamic state, simply read the reduced potential from the weights
            dir = "iter%i"%iters[0]
            weights = np.loadtxt(dir + "/confall.w") 
            for num_line, idx in enumerate(idxs_state):
                traj_total_mem[idx, 1] = log(weights[num_line])

        self.traj_total_mem = traj_total_mem

        logging.info("cTRAM: constructing log_gamma_0")
        log_gamma_0 = np.zeros((nclusters, nstates))

        nnebins_0 = int(5*nsamples_per_state**(1./2))
        nbins_0_min = int(nsamples_per_state**(1./(2*ndcs)))
        nbins_0_max = int(nsamples_per_state**(1./2))
        # construct the histogram with an optimal number of bins
        for kdx in xrange(nstates):
            dcs_k = dcs[kdx*nsamples_per_state:(kdx+1)*nsamples_per_state]
            bins_0, grid_0 = tools.do_grid_optimized(dcs_k, nnebins_0, nbins_0_min, nbins_0_max, mins=bins[0,:]-steps/2, maxs=bins[-1,:]+steps/2)
            steps_0 = bins_0[1,:] - bins_0[0,:]
            free_energy_k = tools.compute_free_energy(grid_0, ndcs, np.ones(nsamples_per_state), np.inf, 1)
            # log_gamma_0
            log_gamma_0_k = np.zeros(nclusters)
            log_gamma_0_k.fill(np.nan)
            # stretch array
            for index_bin, bin in tools.nonempty_bins(grid):
                # compute coordinates of the non-empty bins of the new grid in the old grid 
                old_idxs = [np.array(int(floor((bins[idx,num]-bins_0[0,num])/steps_0[num] + 0.5))) for num, idx in enumerate(index_bin)]
                index_bin_s = 0
                for jdx in range(ndcs):
                    index_bin_s += index_bin[jdx]*nbins**(ndcs-jdx-1)
                # find non-empty bin corresponding to index_bin_s
                for ldx in xrange(nclusters):
                    nebin = nebins[ldx]
                    if index_bin_s == nebin:
                       log_gamma_0_k[ldx] = free_energy_k[old_idxs]
                       break
                    elif index_bin_s < nebin:
                       break
            log_gamma_0_k[np.isnan(log_gamma_0_k)] = np.nanmax(log_gamma_0_k) # set all the nan values as the global maximum
            log_gamma_0_k = log_gamma_0_k - np.amin(log_gamma_0_k)
            log_gamma_0[:, kdx] = -log_gamma_0_k

        self.log_gamma_0 = log_gamma_0

        logging.info("cTRAM: constructing count_matrices")
        count_matrices = np.zeros((nclusters, nclusters, nstates))
        for idx in xrange(ntrajs): 
            for kdx in range(0,nstride-1):
                count_matrices[traj_total_mem[idx*nstride+kdx, 0], traj_total_mem[idx*nstride+kdx+1, 0], idx/nreplicas] += 1

        self.count_matrices = count_matrices

        logging.info("cTRAM: count_matrices has %i x %i elements for each thermodynamic state (%i)"%(nclusters, nclusters, nstates))
        nzcounts = np.nonzero(count_matrices)
        np.savetxt('ctram/count_matrices.dat', np.vstack(nzcounts + (count_matrices[nzcounts],)).T,
                   fmt=" ".join(["%10.1i" for idx in range(3)]) + " " + "%10.1i")

        # save other data
        np.savetxt('ctram/traj_total_mem.dat', traj_total_mem, fmt="%6.1i" + " ".join(["%30.18e" for jdx in range(nstates)]))
        np.savetxt('ctram/N_mem.dat', np.hstack((N_mem, nebins[:,np.newaxis])), fmt=" ".join(["%6.1i" for jdx in range(nstates)]) + " %6.1i")
        np.savetxt('ctram/bins.xyz', bins, fmt='%.18e')
        np.savetxt('ctram/log_gamma_0.dat', log_gamma_0, fmt='%.18e')

        # remove the file embed.gro as it can be quite heavy
        subprocess.check_call('rm -rf ctram/embed.gro', shell=True)

    def do_ctram(self, umgr, settings):

        logging.info('cTRAM: computing cTRAM...')
        print "Performing cTRAM..."
        cu = radical.pilot.ComputeUnitDescription()

        cu.input_staging = [settings.inifile, 'ctram/traj_total_mem.dat', 'ctram/N_mem.dat', 'ctram/count_matrices.dat', 'ctram/log_gamma_0.dat']
        cu.executable = 'ctram4dmaps -f ' + settings.inifile + ' -t traj_total_mem.dat -m N_mem.dat -c count_matrices.dat -l log_gamma_0.dat'
        cu.output_staging = ['log_gamma.dat > ctram/log_gamma.dat', 'log_weight.dat > ctram/log_weight.dat']
        cu.mpi = True
        cu.cores = settings.cores
        cu.cleanup = True

        unit = umgr.submit_units(cu)
        unit.wait()

        self.log_gamma = np.loadtxt('ctram/log_gamma.dat')
        self.log_weight = np.loadtxt('ctram/log_weight.dat')

        logging.info('cTRAM: cTRAM done')
        print 'cTRAM Execution time:', (unit.stop_time - unit.start_time).total_seconds()


# The following class should be used on the remote machine
class CTRAM4DMapSExe(object):

    def create_arg_parser(self):
        parser = argparse.ArgumentParser(description="Run cTRAM for Diffusion Map Sampling...") # message displayed when typing ctram4dmaps -h

        # required options
        parser.add_argument("-f",
            type=str,
            dest="config_file",
            required=True,
            help='Configuration file for Dmap Sampling (input): ini')

        parser.add_argument("-t",
            type=str,
            dest="trajfile",
            required=True,
            help='file containing traj_total_mem (input), traj_total_mem is a N*(K+1) matrix. \
            N is the total number of samples and K is the number of thermodynamic states. \
            Each row of traj_total_mem represents a sample point. Column 1 is the index of \
            the discrete cluster, and Columns 2 ~ K+1 are reduced potential values \
            (in unit of kT) at K thermodynamic states')

        parser.add_argument("-m",
            type=str,
            dest="memfile",
            required=True,
            help='file containing N_mem (input), N_mem is m*K matrix. m is the number of \
            discrete clusters. N_mem(i,k) is the number of occurrence of the i-th cluster at the thermodynamic state k.')

        parser.add_argument("-c",
            type=str,
            dest="countfile",
            required=True,
            help='file containing non zeros elements of count_matrices (input), count_matrices is a 3D matrix with size m*m*K.\
            For each k, count_matrices(:,:,k) is the count matrix of the thermodynamic state k. The file should contain 4 columns \
            the first three columns indicate the matrix indices of the non-zero elements and the fourth column indicates the values')

        parser.add_argument("-l",
            type=str,
            dest="log_gamma_0_file",
            required=True,
            help='file containing log_gamma_0, the initial value of log_gamma (input), log_gamma is a m*K matrix. \
            The value of exp(log_gamma(:,k)) is equal to the stationary distribution of the Markov model at \
            thermodynamic state k after normalization, and sum(exp(log_gamma(:,k))) is proportional to the partition function of thermodynamic state k.')

        return parser

    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD   # MPI environment
        size = comm.Get_size()  # number of threads
        rank = comm.Get_rank()  # number of the current thread 

        if rank == 0:
            parser = self.create_arg_parser() 
            args = parser.parse_args()

            config = ConfigParser.SafeConfigParser()
            config.read(args.config_file) # set config file parser
            niters = config.getint('CTRAM', 'niters')

            traj_total_mem = np.loadtxt(args.trajfile)
            traj_total_mem[:,0] += 1

            N_mem = np.loadtxt(args.memfile)
            N_mem = N_mem[:,:-1]

            state_num, temperature_num = N_mem.shape
            data_num = traj_total_mem.shape[0]

            count_matrices = np.zeros((state_num, state_num, temperature_num))
            idxs = np.genfromtxt(args.countfile, usecols=(0,1,2), dtype="i4")
            value = np.genfromtxt(args.countfile, usecols=3, dtype="i4")
            count_matrices[(idxs[:,0], idxs[:,1], idxs[:,2])] = value

            log_gamma_0 = np.loadtxt(args.log_gamma_0_file)
            log_gamma, log_weight = wrapper.ctram_wrapper(traj_total_mem, N_mem, count_matrices, niters, log_gamma_0)

            np.savetxt('log_gamma.dat',log_gamma, fmt='%.18e')
            np.savetxt('log_weight.dat',log_weight, fmt='%.18e')
