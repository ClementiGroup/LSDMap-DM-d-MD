import os
import sys
import time

import logging
import random
import numpy as np

import radical.pilot
from dmaps.tools import tools
from dmaps.tools.config import known_pre_exec
from lsdmap.rw import reader
from lsdmap.rw import writer

class DMapSamplingWorker(object):

    def do_preprocessing_md(self, settings):

        size = settings.cores

        print "Preparing .gro files..."
        os.system('rm -rf md; mkdir md')

        grofile = open('input.gro', 'r')
        grofile.next()
        natoms = int(grofile.next())
        for idx, line in enumerate(grofile):
            pass
        nlines = idx + 3
        ncoords = nlines/(natoms+3)
        grofile.close()

        if ncoords < size:
            raise ValueError("the number of runs should be greater or equal to the number of threads.")

        ncoords_per_thread = [ncoords/size for _ in xrange(size)]
        nextra_coords = ncoords%size

        for idx in xrange(nextra_coords):
            ncoords_per_thread[idx] += 1

        with open('input.gro', 'r') as grofile:
            for idx in xrange(size):
                grofile_thread =  'md/input%i.gro'%idx
                with open(grofile_thread, 'w') as grofile_t:
                    nlines_per_thread = ncoords_per_thread[idx]*(natoms+3)
                    for jdx in xrange(nlines_per_thread):
                        line = grofile.readline()
                        if line:
                            line = line.replace("\n", "")
                            print >> grofile_t, line
                        else:
                            break

    def do_postprocessing_md(self, settings):

        files_md_iter0 = ['confall.gro', 'confall.w']
        if settings.iter == 0:
            files_md = files_md_iter0
        else:
            files_md = files_md_iter0 + ['confall.ev', 'autocorr.ev']
          
        for filename in files_md:
            name, ext = os.path.splitext(filename)
            with open(filename, 'w') as newfile:
                for idx in range(settings.cores):
                    with open('md/%s%i%s'%(name, idx, ext), 'r') as oldfile:
                        for line in oldfile:
                            print >> newfile, line.replace('\n', '')

        if settings.iter > 0:
            os.system('mv confall.ev confall.ev.embed.old')

    def run_md(self, umgr, settings):

        print "Preprocessing..."
        logging.info('Preprocessing MD...')
        self.do_preprocessing_md(settings)

        print 'Starting simulation...'
        logging.info('Running MD...')
        cud_list = []

        tcpu1 = time.time()

        for idx in xrange(settings.cores):

            cu = radical.pilot.ComputeUnitDescription()
            cu.executable = "/bin/bash"
            cu.arguments = "run.sh"
            cu.pre_exec = known_pre_exec[settings.remote_host]
            
            cu.input_staging = ['run_md.sh > run.sh', 'md/input%i.gro > input.gro'%idx, settings.mdpfile, settings.topfile, settings.inifile]
            if settings.iter > 0:
               cu.input_staging.extend(['fe/bins.xyz', 'fe/hist.dat', 'fit/fit.gro', 'fit/fit.w', 'fit/fit.sig'])
            cu.output_staging = ['confall.gro > md/confall%i.gro'%idx, 'confall.w > md/confall%i.w'%idx]
            if settings.iter > 0:
               cu.output_staging.extend(['confall.ev > md/confall%i.ev'%idx, 'autocorr.ev > md/autocorr%i.ev'%idx])
            cu.cores = 1
            cu.cleanup = True

            cud_list.append(cu)
        units = umgr.submit_units(cud_list)

        start_times = []
        end_times = []

        # Wait for all compute units to finish.
        umgr.wait_units()

        logging.info('MD done')
        print "Postprocessing..."
        logging.info('Postprocessing MD')
        self.do_postprocessing_md(settings)

        tcpu2 = time.time()

        print 'Total Simulation Time : ', tcpu2 - tcpu1

        for unit in units:
            start_times.append(unit.start_time)
            end_times.append(unit.stop_time)
        print 'Simulation Execution Time : ', (max(end_times) - min(start_times)).total_seconds()

    def do_preprocessing_lsdmap(self, settings, npoints, nlsdmap, nbins):

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
       
        if settings.iter == 0:
            idxs_lsdmap = random.sample(range(npoints), nlsdmap)
        else:
            logging.info('Select configurations for lsdmap according to the previous map...')
            logging.info('Build histogram')
            dcs = np.loadtxt('confall.ev.embed.old')
            bins, grid = tools.do_grid(dcs, nbins)
            logging.info('Select fitting points uniformly')
            idxs_lsdmap = tools.pick_points_from_grid(grid, nlsdmap)

        self.coords_lsdmap = self.coords_all[idxs_lsdmap]
        self.weights_lsdmap = self.weights_all[idxs_lsdmap]

        os.system('rm -rf lsdmap; mkdir lsdmap')

        logging.info('Write configurations in lsdmap/lsdmap_aa.gro')
        gw = writer.open('.gro', pattern=settings.startgro)
        gw.write(self.coords_lsdmap, 'lsdmap/lsdmap_aa.gro')
        logging.info('Write weights in lsdmap/lsdmap.w')
        np.savetxt('lsdmap/lsdmap.w', self.weights_lsdmap, fmt='%15.7e')

        logging.info('Create file containing only heavy atoms')
        os.system('cd lsdmap; echo 2 | trjconv -f lsdmap_aa.gro -s lsdmap_aa.gro -o lsdmap.gro &>/dev/null; cd ..')

    def do_postprocessing_lsdmap(self, settings, ndcs):

        logging.info('Store DCs computed')
        dcs = np.loadtxt('lsdmap/lsdmap.ev')
        if ndcs == 1:
            self.dcs_lsdmap = dcs[:,1][:,np.newaxis]
        else:
            self.dcs_lsdmap = dcs[:,1:ndcs+1]

    def run_lsdmap(self, umgr, settings, npoints, nlsdmap, nbins, ndcs):

        print 'Starting LSDMap'
        tcpu1 = time.time()

        logging.info('Preprocessing LSDMap...')
        self.do_preprocessing_lsdmap(settings, npoints, nlsdmap, nbins)
        logging.info('LSDMap preprocessing done')

        logging.info('Starting LSDMap')

        cu = radical.pilot.ComputeUnitDescription()
        cu.pre_exec = known_pre_exec[settings.remote_host]

        cu.input_staging = [settings.inifile, 'lsdmap/lsdmap.gro', 'lsdmap/lsdmap.w']

        cu.executable = 'lsdmap -f ' + settings.inifile + ' -c lsdmap.gro -w lsdmap.w'
        cu.output_staging = ['lsdmap.ev > lsdmap/lsdmap.ev', 'lsdmap.eg > lsdmap/lsdmap.eg', 'lsdmap.log > lsdmap/lsdmap.log']
        cu.mpi = True
        cu.cores = settings.cores
        cu.cleanup = True

        unit = umgr.submit_units(cu)
        unit.wait()

        logging.info("LSDMap done")
        tcpu2 = time.time()

        logging.info('Postprocessing LSDMap...')
        self.do_postprocessing_lsdmap(settings, ndcs)
        logging.info('LSDMap postprocessing done')

        print 'LSDMap Execution time : ',(unit.stop_time - unit.start_time).total_seconds()
        print 'Total Analysis time : ', tcpu2 - tcpu1

    def do_preprocessing_fit_dcs(self, settings, nfit, nbins):

        dcs_lsdmap = self.dcs_lsdmap

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

        os.system('rm -rf fit; mkdir fit')

        logging.info('Write configurations in fit/fit_aa.gro')
        # write gro file used for the fit
        gw = writer.open('.gro', pattern=settings.startgro)
        gw.write(coords_fit, 'fit/fit_aa.gro')

        logging.info('Create file containing only heavy atoms')
        os.system('cd fit; echo 2 | trjconv -f fit_aa.gro -s fit_aa.gro -o fit.gro &>/dev/null; cd ..')

        logging.info('Write DCs in fit/fit.ev')

        # write ev file used for the fit
        np.savetxt('fit/fit.ev', np.hstack((np.ones((nfit,1)), dcs_fit)), fmt='%15.7e')

        os.system('echo 2 | trjconv -f confall.gro -s confall.gro -o fit/embed.gro &>/dev/null')


    def run_fit_dcs(self, umgr, settings, nfit, nbins, ndcs):

        print 'Starting Fitting'
        tcpu1 = time.time()

        logging.info("Preprocessing Fitting...")
        self.do_preprocessing_fit_dcs(settings, nfit, nbins)
        logging.info("Fit preprocessing done")

        dcs_options = '--dc ' + ' '.join([str(num+1) for num in xrange(ndcs)])

        logging.info('Starting Fitting')
        cu = radical.pilot.ComputeUnitDescription()
        cu.pre_exec = known_pre_exec[settings.remote_host]
        cu.input_staging = [settings.inifile, 'fit/fit.gro', 'fit/fit.ev', 'fit/embed.gro']
        cu.executable = 'rbffit -f ' + settings.inifile + ' -c fit.gro -v fit.ev --embed embed.gro ' + dcs_options
        cu.output_staging = ['fit.w > fit/fit.w', 'fit.sig > fit/fit.sig', 'fit.embed > confall.ev.embed']
        cu.mpi = True
        cu.cleanup = True
        cu.cores = settings.cores

        unit = umgr.submit_units(cu)
        unit.wait()

        logging.info("Fitting done")
        tcpu2 = time.time()

        print 'Fitting Execution time : ',(unit.stop_time - unit.start_time).total_seconds()
        print 'Total Analysis time : ', tcpu2 - tcpu1


    def do_free_energy(self, umgr, settings, nbins, cutoff, kT, ndcs):

        print 'Starting Free Energy Estimate'
        os.system('rm -rf fe; mkdir fe')

        dcs = np.loadtxt('confall.ev.embed')
        if ndcs == 1:
            dcs = dcs[:,np.newaxis]
        weights = np.loadtxt('confall.w')

        nvalues = dcs.shape[0]
        nnebins = int(nvalues**(3./4))
        nbins_min = int(nvalues**(1./(2*ndcs)))
        nbins_max = int(nvalues**(1./2))
        bins, grid = tools.do_grid_optimized(dcs, nnebins, nbins_min, nbins_max, nextrabins=1)
 
        logging.info("Build free energy histogram with %i nbins along each dimension..."%len(bins))
        free_energy_grid = tools.compute_free_energy(grid, ndcs, weights, cutoff, kT)

        steps = bins[1,:] - bins[0,:]
        grads = np.gradient(free_energy_grid,*steps.tolist())
        if ndcs == 1:
            grads = [grads]

        # considering only bins where the free energy is non zero
        #nebins = np.nonzero(free_energy_grid)

        # one may want to include bins where the free energy is zero but not the force
        nzgrads = (grads[0] != 0.0)
        if ndcs > 1:
            for dim in xrange(1,ndcs):
                nzgrads = np.logical_or(nzgrads, grads[dim] != 0.0)
        nebins = np.where(nzgrads)

        nbins = bins.shape[0]
        nebins_s = np.sum([nebins[idx]*nbins**(ndcs-idx-1) for idx in xrange(ndcs)], axis=0)

        logging.info("Save free energy histogram")
        np.savetxt('fe/hist.dat', np.vstack(nebins + (np.array(nebins_s), free_energy_grid[nebins],) + tuple([grads[idx][nebins] for idx in range(ndcs)])).T,
                   fmt=" ".join(["%6.1i" for idx in range(ndcs)]) + " " + "%10.1i %15.7e " + " ".join(["%15.7e" for idx in range(ndcs)]))

        np.savetxt('fe/bins.xyz', bins, fmt='%15.7e')

        logging.info("Pick new configurations for the next iteration.")
        npreselect = min(dcs.shape[0], 500)
        if settings.nreplicas < npreselect:
            # preselection
            idxs_preselect = tools.pick_points_from_grid(grid, npreselect)
            # selection
            idxs_new_coords = tools.pick_points_optimized(dcs, settings.nreplicas, idxs_preselect=idxs_preselect)
        else:
            idxs_new_coords = tools.pick_points_from_grid(grid, settings.nreplicas)

        new_coords = self.coords_all[idxs_new_coords]

        logging.info('Save new configurations in output.gro')
        # save new coordinates
        gw = writer.open('.gro', pattern=settings.startgro)
        gw.write(new_coords, 'output.gro')

        # save dcs of new points (check)
        np.savetxt('output.ev', dcs[idxs_new_coords], fmt='%15.7e')
