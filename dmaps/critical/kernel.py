import os
import sys
import time

import logging
import random
import numpy as np
import subprocess
from math import sqrt, floor

from dmaps.tools import tools
from lsdmap.rw import reader
from lsdmap.rw import writer
import ConfigParser

class DMapSamplingWorker(object):

    def do_preprocessing_md(self, config, args):

        print "Preparing .gro files..."
        subprocess.check_call('rm -rf md; mkdir md', shell=True)

        iter = config.iter
        nreplicas = config.nreplicas
        ncpus = args.ncpus

        if iter == 0:
            with open('input.gro', 'w') as ifile:
                for idx in xrange(nreplicas):
                    with open(config.startgro, 'r') as sfile:
                       for line in sfile:
                          print >> ifile, line.replace("\n", "")

        grofile = open('input.gro', 'r')
        grofile.next()
        natoms = int(grofile.next())
        for idx, line in enumerate(grofile):
            pass
        nlines = idx + 3
        ncoords = nlines/(natoms+3)
        grofile.close()

        if ncoords < ncpus:
            raise ValueError("the number of runs should be greater or equal to the number of CPUs.")

        ncoords_per_thread = [ncoords/ncpus for _ in xrange(ncpus)]
        nextra_coords = ncoords%ncpus

        for idx in xrange(nextra_coords):
            ncoords_per_thread[idx] += 1

        with open('input.gro', 'r') as grofile:
            for idx in xrange(ncpus):
                os.mkdir('md/core%i'%idx)
                grofile_thread = 'md/core%i/input.gro'%idx
                with open(grofile_thread, 'w') as grofile_t:
                    nlines_per_thread = ncoords_per_thread[idx]*(natoms+3)
                    for jdx in xrange(nlines_per_thread):
                        line = grofile.readline()
                        if line:
                            line = line.replace("\n", "")
                            print >> grofile_t, line
                        else:
                            break

    def do_postprocessing_md(self, config, args):

        iter = config.iter
        ncpus = args.ncpus

        files_md_iter0 = ['confall.gro', 'confall_aa.gro', 'confall.w']
        if iter == 0:
            files_md = files_md_iter0
        else:
            files_md = files_md_iter0 + ['confall.ev']
          
        for filename in files_md:
            with open(filename, 'w') as newfile:
                for idx in range(ncpus):
                    oldfilename = 'md/core%i/'%idx+filename
                    with open(oldfilename, 'r') as oldfile:
                        for line in oldfile:
                            print >> newfile, line.replace('\n', '')
        if iter > 0:
            subprocess.check_call(['mv', 'confall.ev', 'confall.ev.embed.old'])

    def run_md(self, config, args):

        print "(1) Run biased MD simulations"

        print "Preprocessing..."
        logging.info('Preprocessing MD...')
        self.do_preprocessing_md(config, args)

        print 'Running...'
        logging.info('Running MD...')

        tcpu1 = time.time()
        ncpus = args.ncpus

        subprocess.check_call("mpiexec -n %i p_mdrun_d"%ncpus, shell=True)

        logging.info('MD done')
        print "Postprocessing..."
        logging.info('Postprocessing MD')
        self.do_postprocessing_md(config, args)

        tcpu2 = time.time()
        print 'Simulation Execution Time: ', tcpu2 - tcpu1
        print

    def do_preprocessing_lsdmap(self, config):

        npoints = config.npoints
        nlsdmap = config.nlsdmap
        nbins= int(sqrt(nlsdmap))

        logging.info('Load configurations in confall.gro')
        # read .gro file containing the configurations
        gr = reader.open('confall.gro')
        self.coords_all = gr.readlines()
        gr.close()

        logging.info('Load weights in confall.w')
        # read .w file containing the weights
        self.weights_all = np.loadtxt('confall.w')
        if self.coords_all.shape[0] != npoints:
            logging.error('Number of coordinates in confall.gro (%i) and number of coordinates expected from config file (%i) do no match' \
                          %(self.coords_all.shape[0], npoints))
       
        if config.iter == 0:
            idxs_lsdmap = random.sample(range(npoints), nlsdmap)
        else:
            logging.info('Select configurations for lsdmap according to the previous map...')
            logging.info('Build histogram')
            dcs = np.loadtxt('confall.ev.embed.old')
            bins, grid = tools.do_grid(dcs, nbins)
            logging.info('Select points for lsdmap uniformly')
            idxs_lsdmap = tools.pick_points_from_grid(grid, nlsdmap)

        self.coords_lsdmap = self.coords_all[idxs_lsdmap]
        self.weights_lsdmap = self.weights_all[idxs_lsdmap]

        subprocess.check_call('rm -rf lsdmap; mkdir lsdmap', shell=True)

        gw = writer.open('.gro', pattern='confall.gro')
        logging.info('Write configurations in lsdmap/lsdmap.gro')
        gw.write(self.coords_lsdmap, 'lsdmap/lsdmap.gro')
        logging.info('Write weights in lsdmap/lsdmap.w')
        np.savetxt('lsdmap/lsdmap.w', self.weights_lsdmap, fmt='%.18e')
        if config.metric=='pca':
          subprocess.check_call("cp grompp_pca.mdp lsdmap", shell=True)
          subprocess.check_call("cp topol.top lsdmap", shell=True)

    def do_postprocessing_lsdmap(self, config):

        ndcs = config.ndcs
        logging.info('Store DCs computed')
        dcs = np.loadtxt('lsdmap/lsdmap.ev')
        if ndcs == 1:
            self.dcs_lsdmap = dcs[:,1][:,np.newaxis]
        else:
            self.dcs_lsdmap = dcs[:,1:ndcs+1]

    def run_lsdmap(self, config, args):

        print "(2) Run LSDMap"
        tcpu1 = time.time()

        print 'Preprocessing...'
        logging.info('Preprocessing LSDMap...')
        self.do_preprocessing_lsdmap(config)
        logging.info('LSDMap preprocessing done')

        print 'Running...'
        logging.info('Starting LSDMap')

        curdir = os.getcwd()
        os.chdir('lsdmap')

        ncpus = args.ncpus
        inifile = args.config_file

        p = subprocess.check_call('mpiexec -n ' + str(ncpus) + ' lsdmap -f ../' + inifile \
             + ' -c lsdmap.gro -w lsdmap.w', shell=True)
        os.chdir(curdir)

        logging.info("LSDMap done")
        tcpu2 = time.time()

        logging.info('Postprocessing LSDMap...')
        self.do_postprocessing_lsdmap(config)
        logging.info('LSDMap postprocessing done')

        print 'Total Analysis time : ', tcpu2 - tcpu1
        print

    def do_preprocessing_fit(self, config):

        dcs_lsdmap = self.dcs_lsdmap
        nfit = config.nfit
        nbins= int(sqrt(config.nlsdmap))

        logging.info('Select configurations used for the fitting...')
        logging.info('Build histogram')
        bins, grid = tools.do_grid(dcs_lsdmap, nbins)
        logging.info('Select fitting points uniformly along the DCs')

        nlsdmap = dcs_lsdmap.shape[0]
        npreselect = min(nlsdmap, 1000)
        if nfit < npreselect:
            # preselection
            idxs_preselect = tools.pick_points_from_grid(grid, npreselect)
            # selection
            logging.info('Optimize selection of fitting points')
            idxs_fit = tools.pick_points_optimized(dcs_lsdmap, nfit, idxs_preselect=idxs_preselect)
        else:
            idxs_fit = tools.pick_points_from_grid(grid, nfit)

        coords_fit = self.coords_lsdmap[idxs_fit]
        dcs_fit = dcs_lsdmap[idxs_fit]

        subprocess.check_call('rm -rf fit; mkdir fit', shell=True)

        gw = writer.open('.gro', pattern='confall.gro')
        logging.info('Write configurations in fit/fit.gro')
        # write gro file used for the fit
        gw.write(coords_fit, 'fit/fit.gro')

        # write ev file used for the fit
        logging.info('Write DCs in fit/fit.ev')
        np.savetxt('fit/fit.ev', np.hstack((np.ones((nfit,1)), dcs_fit)), fmt='%.18e')

    def run_fit(self, config, args):

        print '(3) Fit LSDMap coordinates...'
        tcpu1 = time.time()

        print 'Preprocessing...'
        logging.info("Preprocessing Fitting...")
        self.do_preprocessing_fit(config)
        logging.info("Fit preprocessing done")

        dcs_options = '--dc ' + ' '.join([str(num+1) for num in xrange(config.ndcs)])

        print 'Fitting...'
        logging.info('Starting Fitting')
        curdir = os.getcwd()
        os.chdir('fit')

        ncpus = args.ncpus
        inifile = args.config_file

        for idx in range(config.iter+1):
            if idx < config.iter:
                embed_dir = '../iter%i'%idx
            else:
                embed_dir = '..'
            subprocess.check_call('mpiexec -n ' + str(ncpus) + ' rbffit -f ../' + inifile + \
                ' -c fit.gro -v fit.ev --embed ' + embed_dir + '/confall.gro ' + dcs_options, shell=True)
            subprocess.check_call('mv fit.embed confall.ev.embed.%i'%idx, shell=True)

        os.chdir(curdir)

        # copy the embedded configurations of the current step 
        subprocess.check_call('cp fit/confall.ev.embed.%i confall.ev.embed'%config.iter, shell=True)

        logging.info("Fitting done")
        tcpu2 = time.time()

        print 'Total Analysis time : ', tcpu2 - tcpu1
        print

    def do_free_energy(self, config):

        ndcs = config.ndcs
        cutoff = config.cutoff
        kT = config.kT

        print '(4) Compute the free energy (histogram)'
        subprocess.check_call('rm -rf fe; mkdir fe', shell=True)

        for idx in range(config.iter+1):
            dcs = np.loadtxt('fit/confall.ev.embed.%i'%idx)
            if ndcs == 1:
                dcs = dcs[:,np.newaxis]
            if config.nbinsfe is None:
                nvalues = dcs.shape[0]
                nnebins = int(5*nvalues**(1./2))
                nbins_min = int(nvalues**(1./(2*ndcs)))
                nbins_max = int(nvalues**(1./2))
                bins, grid = tools.do_grid_optimized(dcs, nnebins, nbins_min, nbins_max, nextrabins=0)
            else:
                bins, grid = tools.do_grid(dcs, config.nbinsfe, nextrabins=0)
    
            nbins = bins.shape[0]
            logging.info("Build free energy histogram with %i nbins along each dimension..."%bins.shape[0])
            free_energy_grid = tools.compute_free_energy(grid, ndcs, np.ones(config.npoints), config.cutoff, kT)
    
            nebins = np.where(~np.isnan(free_energy_grid))
            nebins_s = np.sum([nebins[jdx]*nbins**(ndcs-jdx-1) for jdx in xrange(ndcs)], axis=0)
    
            steps = bins[1,:] - bins[0,:]
            grads = np.gradient(free_energy_grid,*steps.tolist())
            if ndcs == 1:
                grads = [grads]
            grads = [np.nan_to_num(grad) for grad in grads]
    
            logging.info("Save free energy histogram")
            np.savetxt('fe/hist.dat.%i'%idx, np.vstack(nebins + (np.array(nebins_s), free_energy_grid[nebins],) + tuple([grads[jdx][nebins] for jdx in range(ndcs)])).T,
                       fmt=" ".join(["%6.1i" for jdx in range(ndcs)]) + " " + "%10.1i %.18e " + " ".join(["%.18e" for jdx in range(ndcs)]))
            np.savetxt('fe/bins.xyz.%i'%idx, bins, fmt='%.18e')
        print

    def select_new_points(self, config):

        print '(5) Select new points'
        # take only the dcs of the last iteration
        dcs = np.loadtxt('confall.ev.embed')
        ndcs = config.ndcs
        if ndcs == 1:
            dcs = dcs[:,np.newaxis]
        logging.info("Pick new configurations for the next iteration.")

        if config.uniform_sampling == 1: # select uniformly the configurations in the DC space
            nvalues = dcs.shape[0]
            nnebins = int(5*nvalues**(1./2))
            nbins_min = int(nvalues**(1./(2*ndcs)))
            nbins_max = int(nvalues**(1./2))
            bins, grid = tools.do_grid_optimized(dcs, nnebins, nbins_min, nbins_max, nextrabins=0)

            npreselect = min(dcs.shape[0], 500)
            if config.nreplicas < npreselect:
                # preselection
                idxs_preselect = tools.pick_points_from_grid(grid, npreselect)
                # selection
                idxs_new_coords = tools.pick_points_optimized(dcs, config.nreplicas, idxs_preselect=idxs_preselect)
            else:
                idxs_new_coords = tools.pick_points_from_grid(grid, config.nreplicas)
        elif config.uniform_sampling == 0: # take the last configurations of each traj as the new starting points
            idxs_new_coords = [(idx + config.nstride -1) for idx in range(0,config.nstride*config.nreplicas,config.nstride)]

        # sort elements of the list
        idxs_new_coords.sort()

        # get the number of atoms from startgro
        f = open(config.startgro, 'r')
        f.next()
        natoms = int(f.next())
        f.close()

        shift = 0
        logging.info('Save new configurations in output.gro')
        with open('output.gro', 'w') as outfile:
            with open('confall_aa.gro', 'r') as infile:
                for idx in idxs_new_coords:
                    # skip lines before the next configuration selected
                    for jdx in range((natoms+3)*(idx-shift)):
                        a=infile.next()
                    # print the configuration selected
                    for jdx in range((natoms+3)):
                        print >> outfile, infile.next().replace("\n", "")
                    shift = idx+1

        # save dcs of new points (check)
        np.savetxt('output.ev', dcs[idxs_new_coords], fmt='%.18e')
