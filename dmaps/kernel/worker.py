import os
import sys
import time

import logging
import random
import numpy as np

import radical.pilot
from dmaps.tools import tools
from dmaps.tools.config import known_pre_exec, tmpfiles
from lsdmap.rw import reader
from lsdmap.rw import writer

class DMapSamplingWorker(object):

    def do_preprocessing_md(self, settings, tmpdir):

        size = settings.cores

        print "Preparing .gro files..."

        os.system('rm -rf ' + tmpdir)
        os.makedirs(tmpdir)

        grofile = open(tmpfiles['ingro'][0], 'r')
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

        with open(tmpfiles['ingro'][0], 'r') as grofile:
            for idx in xrange(size):
                grofile_thread = tmpdir + '/' + 'input%s.gro' %idx
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

        with open('tmp.gro', 'w') as output_grofile:
            for idx in range(settings.cores):
                with open('tmp/out%i.gro' %idx, 'r') as output_file:
                    for line in output_file:
                        print >> output_grofile, line.replace("\n", "")
                #os.remove('tmp/out%i.gro' %idx)

        with open('confall.gro', 'w') as output_grofile:
            for idx in range(settings.cores):
                with open('tmp/confall%i.gro' %idx, 'r') as output_file:
                    for line in output_file:
                        print >> output_grofile, line.replace("\n", "")
                #os.remove('tmp/confall%i.gro' %idx)

        with open('confall.w', 'w') as output_wfile:
            for idx in range(settings.cores):
                with open('tmp/confall%i.w' %idx, 'r') as output_file:
                    for line in output_file:
                        print >> output_wfile, line.replace("\n", "")
                #os.remove('tmp/confall%i.w' %idx)

        if settings.iter > 0:
            with open('confall.ev.embed.old', 'w') as output_evfile:
                for idx in range(settings.cores):
                    with open('tmp/confall%i.ev' %idx, 'r') as output_file:
                        for line in output_file:
                            print >> output_evfile, line.replace("\n", "")
                    #os.remove('tmp/confall%i.ev' %idx)
            with open('autocorr.ev', 'w') as output_evfile:
                for idx in range(settings.cores):
                    with open('tmp/autocorr%i.ev' %idx, 'r') as output_file:
                        for line in output_file:
                            print >> output_evfile, line.replace("\n", "")
                    #os.remove('tmp/autocorr%i.ev' %idx)


    def run_md(self, umgr, settings):

        curdir = os.getcwd()
        tmpdir = curdir + '/' + 'tmp'
        fitdir = curdir + '/' + 'fit'
        fedir = curdir + '/' + 'fe'
     
        print "Preprocessing..."
        logging.info('Preprocessing MD...')
        self.do_preprocessing_md(settings, tmpdir)

        print 'Starting simulation...'
        logging.info('Running MD...')
        cud_list = []

        p1 = time.time()

        for idx in xrange(settings.cores):
            cu = radical.pilot.ComputeUnitDescription()
            cu.executable = "/bin/bash"
            cu.arguments = "run.sh"
            cu.pre_exec = known_pre_exec[settings.remote_host]
            cu.input_staging = ['run_md.sh > run.sh', tmpdir + '/' + 'input%i.gro > input.gro' %idx,\
                settings.mdpfile, settings.topfile, settings.inifile]
            if settings.iter > 0:
               cu.input_staging.extend([fedir + '/' + 'bins.xy > bins_fe.xy', fedir + '/' + 'gradient.xy > gradient_fe.xy', fedir + '/' + 'hist.xyz > hist_fe.xyz',\
                                  fitdir + '/' + 'fit.gro', fitdir + '/' + 'fit.w', fitdir + '/' + 'fit.sig'])
            cu.output_staging = ['confall.gro > tmp/confall%i.gro' %idx, 'out.gro > tmp/out%i.gro' %idx, 'confall.w > tmp/confall%i.w' %idx]
            if settings.iter > 0:
               cu.output_staging.extend(['confall.ev > tmp/confall%i.ev' %idx, 'autocorr.ev > tmp/autocorr%i.ev' %idx])

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

        p2 = time.time()

        print 'Total Simulation Time : ', (p2-p1)

        for unit in units:
            start_times.append(unit.start_time)
            end_times.append(unit.stop_time)

        print 'Simulation Execution Time : ', (max(end_times)-min(start_times)).total_seconds()

    # TODO: do the preprocessing in parallel to optimize the loading of all the configurations
    def do_preprocessing_lsdmap(self, settings, npoints, nlsdmap, nbins, lsdmapdir):

        logging.info('Load coordinates in confall.gro')

        # read .gro file containing all the configurations
        grofile = reader.open("confall.gro")
        self.coords_all = grofile.readlines()

        logging.info('Load weights in confall.w')

        # read .w file containing all the weights
        wfile = reader.open("confall.w")
        self.weights_all = wfile.readlines()
        if self.coords_all.shape[0] != npoints:
            logging.error("Number of coordinates in confall.gro (%i) and number of coordinates implied in config file (%i) do no match" \
                          %(self.coords_all.shape[0], npoints))
       
        if settings.iter == 0:
            self.idxs_lsdmap = random.sample(range(npoints), nlsdmap)
        else:
            # select lsdmap points according to the previous map
            old_dcs = np.loadtxt('confall.ev.embed.old')
            logging.info('Select configurations for lsdmap according to the previous map...')
            logging.info('Build 2D histogram')
            old_dc1s = old_dcs[:,0]
            old_dc2s = old_dcs[:,1]
            # build histogram
            bins1, bins2, hist_idxs = tools.do_hist2D(old_dc1s, old_dc2s, nbins)
            # select the fitting points uniquely and uniformly along the first two DCs
            logging.info('Select fitting points uniformly')
            self.idxs_lsdmap = tools.pick_points_from_hist2D(hist_idxs, nbins, nlsdmap)

        self.coords_lsdmap = self.coords_all[self.idxs_lsdmap]
        self.weights_lsdmap = self.weights_all[self.idxs_lsdmap]

        os.system('rm -rf ' + lsdmapdir)
        os.makedirs(lsdmapdir)

        logging.info('Write configurations in ' + lsdmapdir + '/' + 'lsdmap_aa.gro')
        # write gro file used for lsdmap
        grofile_w = writer.open('.gro', pattern='confall.gro')
        grofile_w.write(self.coords_lsdmap, lsdmapdir + '/' + 'lsdmap_aa.gro')

        logging.info('Write weights in ' + lsdmapdir + '/' + 'lsdmap.w')
        # write gro file used for lsdmap
        np.savetxt(lsdmapdir + '/' + 'lsdmap.w', self.weights_lsdmap, fmt='%15.7e')

        logging.info('Create file lsdmap.gro containing only heavy atoms')
        os.system("echo 2 | trjconv -f " + lsdmapdir + '/' + "lsdmap_aa.gro -s " + lsdmapdir + '/' + "lsdmap_aa.gro \
                   -o " + lsdmapdir + '/' + "lsdmap.gro &>/dev/null")

    def do_postprocessing_lsdmap(self, settings, lsdmapdir):

        logging.info('Store DCs computed')

        evfile = reader.open(lsdmapdir + '/' + 'lsdmap.ev')
         
        dcs = evfile.readlines()
        self.dc1s = dcs[:,1]
        self.dc2s = dcs[:,2]


    def run_lsdmap(self, umgr, settings, npoints, nlsdmap, nbins):

        curdir = os.getcwd()
        lsdmapdir = curdir + '/' + 'lsdmap'

        print 'Starting LSDMap'
        p1=time.time()

        logging.info('Preprocessing LSDMap...')
        self.do_preprocessing_lsdmap(settings, npoints, nlsdmap, nbins, lsdmapdir)
        logging.info('LSDMap preprocessing done')

        logging.info('Starting LSDMap')

        cu = radical.pilot.ComputeUnitDescription()
        cu.pre_exec = known_pre_exec[settings.remote_host]
        cu.input_staging = [curdir + '/' + settings.inifile, lsdmapdir + '/' + 'lsdmap.gro', lsdmapdir + '/' +  'lsdmap.w']
        cu.executable = 'lsdmap' +  ' -f ' + settings.inifile + ' -c ' + 'lsdmap.gro ' + ' -w ' + 'lsdmap.w'
        cu.output_staging = ['lsdmap.ev > ' + lsdmapdir + '/' + 'lsdmap.ev', 'lsdmap.eg > ' + lsdmapdir + '/' + 'lsdmap.eg', "lsdmap.log > " + lsdmapdir + '/' + "lsdmap.log"]
        cu.mpi = True
        cu.cores = settings.cores
        cu.cleanup = True

        unit = umgr.submit_units(cu)
        unit.wait()

        logging.info("LSDMap done")
        p2=time.time()

        logging.info('Postprocessing LSDMap...')
        self.do_postprocessing_lsdmap(settings, lsdmapdir)
        logging.info('LSDMap postprocessing done')

        print 'LSDMap Execution time : ',(unit.stop_time - unit.start_time).total_seconds()
        print 'Total Analysis time : ',p2 - p1


    def do_preprocessing_fit_dcs(self, settings, nfit, nbins, fitdir):

        logging.info('Select configurations used for the fitting...')
        logging.info('Build 2D histogram')
        # build histogram
        bins1, bins2, hist_idxs = tools.do_hist2D(self.dc1s, self.dc2s, nbins)
        # select the fitting points uniquely and uniformly along the first two DCs
        logging.info('Select fitting points uniformly along the DCs')
        ndcs = self.dc1s.shape[0]
        ndcs_preselect_max = min(ndcs, 3000)
        if nfit < ndcs_preselect_max:
            # preselection
            idxs_preselect_dcs = tools.pick_points_from_hist2D(hist_idxs, nbins, ndcs_preselect_max)
            # selection
            idxs_fit = tools.pick_points_2D_optimized(self.dc1s, self.dc2s, nfit, idxs_preselect=idxs_preselect_dcs)
        else:
            idxs_fit = tools.pick_points_from_hist2D(hist_idxs, nbins, nfit)

        self.idxs_fit = idxs_fit

        self.coords_fit = self.coords_lsdmap[idxs_fit]
        self.dc1s_fit = self.dc1s[idxs_fit]
        self.dc2s_fit = self.dc2s[idxs_fit]

        os.system('rm -rf ' + fitdir)
        os.makedirs(fitdir)

        logging.info('Write configurations in ' + fitdir + '/' +'fit_aa.gro')
        # write gro file used for the fit
        grofile_w = writer.open('.gro', pattern='confall.gro')
        grofile_w.write(self.coords_fit, fitdir + '/' + 'fit_aa.gro')

        logging.info('Create file fit.gro containing only heavy atoms')
        os.system("echo 2 | trjconv -f " + fitdir + '/' + "fit_aa.gro -s " + fitdir + '/' + "fit_aa.gro \
                   -o " + fitdir + '/' + "fit.gro &>/dev/null")

        logging.info('Write DCs in ' + fitdir + '/' + 'fit.ev')
        # write ev file used for the fit
        np.savetxt(fitdir + '/' + 'fit.ev', np.array([np.ones(nfit), self.dc1s_fit, self.dc2s_fit]).T, fmt='%15.7e')

    def run_fit_dcs(self, umgr, settings, nfit, nbins):

        curdir = os.getcwd()
        fitdir = curdir + '/' + 'fit'

        print 'Starting Fitting'
        p1=time.time()

        logging.info("Preprocessing Fitting...")
        self.do_preprocessing_fit_dcs(settings, nfit, nbins, fitdir)
        os.system("echo 2 | trjconv -f " + curdir + '/' + "confall.gro -s " + curdir + '/' + "confall.gro \
                   -o " + fitdir + '/' + "embed.gro &>/dev/null")
        logging.info("Fit preprocessing done")

        logging.info('Starting Fitting')

        cu = radical.pilot.ComputeUnitDescription()
        cu.pre_exec = known_pre_exec[settings.remote_host]
        cu.input_staging = [curdir + '/' + settings.inifile, fitdir + '/' + 'fit.gro', fitdir + '/' + 'fit.ev', fitdir + '/' + 'embed.gro']
        cu.executable = 'rbffit' + ' -f ' + settings.inifile + ' -c fit.gro -v fit.ev --embed embed.gro  --dc 1 2'
        cu.output_staging = ["fit.w > fit/fit.w", "fit.sig > fit/fit.sig", "fit.embed > confall.ev.embed"]
        cu.mpi = True
        cu.cleanup = True
        cu.cores = settings.cores

        unit = umgr.submit_units(cu)
        unit.wait()

        logging.info("Fitting done")
        p2 = time.time()

        print 'Fitting Execution time : ',(unit.stop_time - unit.start_time).total_seconds()
        print 'Total Analysis time : ', p2 - p1


    def do_free_energy(self, nbins, cutoff, kT):

        curdir = os.getcwd()
        fedir = curdir + '/' + 'fe'

        os.system('rm -rf ' + fedir)
        os.makedirs(fedir)

        bins1, bins2, free_energy_grid = self.compute_free_energy_hist(nbins, cutoff, kT, fedir)
        grad1, grad2 = np.gradient(free_energy_grid, bins1[1]-bins1[0], bins2[1]-bins2[0])

        nbinstot = bins1.shape[0]

        # construct dc1 and dc2 grids
        bins1_grid = bins1[:, np.newaxis].dot(np.ones((1, nbinstot)))
        bins2_grid = np.ones((nbinstot,1)).dot(bins2[np.newaxis])

        logging.info("Save free energy histogram")

        xyzfile = open(fedir + '/' + 'hist.xyz', 'w')

        for line in np.dstack((bins1_grid, bins2_grid,free_energy_grid)):
            for bin1, bin2, fe in line:
                print >> xyzfile, '%15.7e %15.7e %15.7e' %(bin1, bin2, fe)

        xyzfile.close()

        np.savetxt(fedir + '/' + 'bins.xy', np.array([bins1, bins2]).T, fmt='%15.7e')

        logging.info("Save free energy gradient")
        slfile = open(fedir + '/' + 'gradient.xy', 'w')

        for line in np.dstack((grad1, grad2)):
            for g1, g2 in line:
                print >> slfile, '%15.7e %15.7e' %(g1, g2)

        slfile.close()
        return

    def compute_free_energy_hist(self, nbins, cutoff, kT, fedir):

        dcs = np.loadtxt("confall.ev.embed")
        wfile = reader.open("confall.w")
        weights = wfile.readlines()

        nextrabins = 1
        nsmooth = 2

        nbinstot = nbins+2*nextrabins

        self.dc1s_embed = dcs[:,0]
        self.dc2s_embed = dcs[:,1]

        logging.info("Build free energy histogram...")
        # build histogram
        bins1, bins2, hist_idxs = tools.do_hist2D(self.dc1s_embed, self.dc2s_embed, nbins, nextrabins=nextrabins)

        # build free energy grid
        free_energy_grid = np.zeros((nbinstot, nbinstot))
        for idx, row in enumerate(hist_idxs):
            for jdx, col in enumerate(row):
                npoints = len(col)
                weight = 0
                for grid_idx in col:
                    weight += weights[grid_idx]
                if npoints == 0:
                    free_energy_grid[idx, jdx] = None
                else:
                    free_energy_grid[idx, jdx] = -kT*np.log(weight)

        # give values to all NaN
        free_energy_grid[np.isnan(free_energy_grid)] = np.nanmax(free_energy_grid) + 0.1

        # smooth the grid
        free_energy_grid = tools.smooth2a(free_energy_grid, nsmooth, nsmooth)

        # rescale so that the maximum value is 0
        free_energy_grid -= np.max(free_energy_grid)

        # rescale if the minimum is < than - cutoff
        min_free_energy_grid = np.min(free_energy_grid)
        if min_free_energy_grid < -cutoff :
            free_energy_grid += -min_free_energy_grid - cutoff
            free_energy_grid[free_energy_grid>0.0] = 0.0

        free_energy_grid = np.copy(free_energy_grid) # without this line it fails

        # rescale the bins so that it takes the middle of each old bin
        bins1 = (bins1[1]-bins1[0])/2 + bins1[:-1]
        bins2 = (bins2[1]-bins2[0])/2 + bins2[:-1]

        self.nbinsfe = nbins
        self.nextrabinsfe = nextrabins
        self.hist_idxs_fe = hist_idxs

        return bins1, bins2, free_energy_grid

    def pick_new_points(self, settings):

        curdir = os.getcwd()

        logging.info("Pick new configurations for the next iteration.")
        ndcs = self.dc1s_embed.shape[0]
        ndcs_preselect_max = min(ndcs, 1000)
        if settings.nreplicas < ndcs_preselect_max:
            # preselection
            idxs_preselect_dcs = tools.pick_points_from_hist2D(self.hist_idxs_fe, self.nbinsfe+2*self.nextrabinsfe, ndcs_preselect_max)
            # selection
            idxs_new_dcs = tools.pick_points_2D_optimized(self.dc1s_embed, self.dc2s_embed, settings.nreplicas, idxs_preselect=idxs_preselect_dcs)
        else:
            idxs_new_dcs = tools.pick_points_from_hist2D(self.hist_idxs_fe, self.nbinsfe+2*self.nextrabinsfe, settings.nreplicas)

        new_coords = self.coords_all[idxs_new_dcs]

        logging.info("Save new configurations in " + curdir + '/' + 'output.gro')
        # save new coordinates
        grofile_w = writer.open('.gro', pattern='confall.gro')
        grofile_w.write(new_coords, 'output.gro')

        # save dcs of new points (check)
        np.savetxt("output.ev", np.array([self.dc1s_embed[idxs_new_dcs], self.dc2s_embed[idxs_new_dcs]]).T, fmt='%15.7e')
