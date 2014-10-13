import os
import sys
import time

import shutil
import logging
import random
import numpy as np
import radical.pilot

from dmaps.tools import tools
from lsdmap.rw import reader
from lsdmap.rw import writer


stampede_pre_exec = ["module load -intel intel/14.0.1.106", "module load python",
"PYTHONPATH=$PYTHONPATH:/opt/apps/intel14/mvapich2_2_0/python/2.7.6/lib/python2.7",
"PYTHONPATH=$PYTHONPATH:/opt/apps/intel14/mvapich2_2_0/python/2.7.6/lib/python2.7/site-packages",
"PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages","export PYTHONPATH"]

davinci_pre_exec = ["PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.7",
"PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.7/site-packages/numpy",
"export PYTHONPATH"]

biou_pre_exec = ["PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.7",
"PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.7/site-packages/numpy",
"export PYTHONPATH"]


known_pre_exec = {"stampede.tacc.utexas.edu": stampede_pre_exec, "xsede.stampede": stampede_pre_exec,
"davinci.rice.edu": davinci_pre_exec, "rice.davinci": davinci_pre_exec, "rice.biou": biou_pre_exec}

class DMapSamplingWorker(object):

    def do_preprocessing_md(self, settings, tmpdir):

        size = settings.cores

        print "Preparing .gro files..."

        shutil.rmtree(tmpdir, ignore_errors=True)
        os.makedirs(tmpdir)

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
    def do_preprocessing_lsdmap(self, settings, npoints, nlsdmap, nbins, border_frac, is2nddc, lsdmapdir):

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
            if is2nddc == 1:
                logging.info('Select configurations for lsdmap according to the previous map...')
                logging.info('Build 2D histogram')
                old_dc1s = old_dcs[:,0]
                old_dc2s = old_dcs[:,1]
                # build histogram
                bins1, bins2, hist_idxs = tools.do_hist2D(old_dc1s, old_dc2s, nbins)
                # select the fitting points uniquely and uniformly along the first two DCs
                logging.info('Select fitting points uniformly')
                self.idxs_lsdmap = tools.draw_points_hist2D(hist_idxs, nbins, nlsdmap, border_frac=border_frac)
            else:
                raise NotImplementedError("is2nddc should be equal to 1!")

        self.coords_lsdmap = self.coords_all[self.idxs_lsdmap]
        self.weights_lsdmap = self.weights_all[self.idxs_lsdmap]

        shutil.rmtree(lsdmapdir, ignore_errors=True)
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


    def run_lsdmap(self, umgr, settings, npoints, nlsdmap, nbins, border_frac, is2nddc):

        curdir = os.getcwd()
        lsdmapdir = curdir + '/' + 'lsdmap'

        print 'Starting LSDMap'
        p1=time.time()

        logging.info('Preprocessing LSDMap...')
        self.do_preprocessing_lsdmap(settings, npoints, nlsdmap, nbins, border_frac, is2nddc, lsdmapdir)
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


    def do_preprocessing_fit_dcs(self, settings, nfit, nbins, border_frac, is2nddc, fitdir):

        if is2nddc == 1:
            logging.info('Select configurations used for the fitting...')
            logging.info('Build 2D histogram')
            # build histogram
            bins1, bins2, hist_idxs = tools.do_hist2D(self.dc1s, self.dc2s, nbins)
            # select the fitting points uniquely and uniformly along the first two DCs
            logging.info('Select fitting points uniformly')
            idxs_fit = tools.draw_points_hist2D(hist_idxs, nbins, nfit, border_frac=border_frac)
        else:
            raise NotImplementedError("is2nddc should be equal to 1!")

        self.idxs_fit = idxs_fit

        self.coords_fit = self.coords_lsdmap[idxs_fit]
        self.dc1s_fit = self.dc1s[idxs_fit]
        self.dc2s_fit = self.dc2s[idxs_fit]

        shutil.rmtree(fitdir, ignore_errors=True)
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

    def run_fit_dcs(self, umgr, settings, nfit, nbins, border_frac, is2nddc):

        curdir = os.getcwd()
        fitdir = curdir + '/' + 'fit'

        print 'Starting Fitting'
        p1=time.time()

        logging.info("Preprocessing Fitting...")
        self.do_preprocessing_fit_dcs(settings, nfit, nbins, border_frac, is2nddc, fitdir)
        os.system("echo 2 | trjconv -f " + curdir + '/' + "confall.gro -s " + curdir + '/' + "confall.gro \
                   -o " + fitdir + '/' + "embed.gro &>/dev/null")
        logging.info("Fit preprocessing done")

        logging.info('Starting Fitting')

        if is2nddc == 1:
            dc_options = "--dc 1 2"
        else:
            dc_options = "--dc 1"

        cu = radical.pilot.ComputeUnitDescription()
        cu.pre_exec = known_pre_exec[settings.remote_host]
        cu.input_staging = [curdir + '/' + settings.inifile, fitdir + '/' + 'fit.gro', fitdir + '/' + 'fit.ev', fitdir + '/' + 'embed.gro']
        cu.executable = 'rbffit' + ' -f ' + settings.inifile + ' -c fit.gro -v fit.ev --embed embed.gro' +  ' ' + dc_options
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


#    def do_preprocessing_fit_free_energy(self, nbins, nfit, cutoff, kT, is2nddc, fedir):
#    
#        bins1, bins2, free_energy_grid = self.compute_free_energy_hist(nbins, cutoff, kT, is2nddc)
#
#        ninter = nbinstot**2/nfit
#        shutil.rmtree(fedir, ignore_errors=True)
#        os.makedirs(fedir)
#
#        logging.info("Save results free energy histogram")
#
#        xyfile = open(fedir + '/' + 'free_energy.xy', 'w')
#        slfile = open(fedir + '/' + 'free_energy.sl', 'w')
#
#        # construct dc1 and dc2 grids
#        bins1_grid = bins1[:, np.newaxis].dot(np.ones((1, nbinstot)))
#        bins2_grid = np.ones((nbinstot,1)).dot(bins2[np.newaxis])
#
#        free_energy_grid = np.copy(free_energy_grid) # without this line dstack fails
#
#        for idx, line in enumerate(np.dstack((bins1_grid, bins2_grid, free_energy_grid))):
#            for jdx, (bin1, bin2, fe) in enumerate(line):
#                if (idx*nbinstot+jdx)%ninter == 0:
#                    print >> xyfile, '%15.7e %15.7e' %(bin1, bin2)
#                    print >> slfile, '%15.7e' %fe
#
#        xyfile.close()
#        slfile.close()
#
#        # save boundaries
#        np.savetxt(fedir + '/' + 'bound.txt', np.array([bins1[0], bins1[-1], bins2[0], bins2[-1]])[np.newaxis], fmt='%15.7e')
#
#        self.nbins_fe = nbins
#        self.nextrabins_fe = nextrabins
#        self.hist_fe = hist
#
#    def run_fit_free_energy(self, umgr, settings, nbins, nfit, cutoff, kT, is2nddc):
#       
#        curdir = os.getcwd()
#        fedir = curdir + '/' + 'fe'
#
#        print 'Starting Free Energy Estimate'
#        p1 = time.time()
#
#        logging.info("Preprocessing free energy...")
#        self.do_preprocessing_free_energy(nbins, nfit, cutoff, kT, is2nddc, fedir)
#        logging.info("Preprocessing done")
#
#        logging.info("Fit free energy...")
#        cu = radical.pilot.ComputeUnitDescription()
#        cu.executable = 'rbffit' + ' -f ' +  curdir + '/' + settings.inifile + ' -c ' + fedir + '/' + 'free_energy.xy' + \
#                        ' -v ' + fedir + '/' + 'free_energy.sl' + ' --section ' + 'FE_FITTING'
#        cu.output_data = ["fit.w > " + fedir + "/fit.w", "fit.sig > " + fedir + "/fit.sig"]
#        cu.mpi = True
#        cu.cores = settings.cores
#
#        unit = umgr.submit_units(cu)
#        unit.wait()
#
#        logging.info("Fit free energy done...")
#        p2 = time.time()
#
#        print 'Fitting Execution time : ',(unit.stop_time - unit.start_time).total_seconds()
#        print 'Total Analysis time : ', p2 - p1
#
#
    def do_free_energy(self, nbins, cutoff, kT, is2nddc):

        curdir = os.getcwd()
        fedir = curdir + '/' + 'fe'

        shutil.rmtree(fedir, ignore_errors=True)
        os.makedirs(fedir)

        bins1, bins2, free_energy_grid = self.compute_free_energy_hist(nbins, cutoff, kT, is2nddc, fedir)
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

    def compute_free_energy_hist(self, nbins, cutoff, kT, is2nddc, fedir):

        dcs = np.loadtxt("confall.ev.embed")
        wfile = reader.open("confall.w")
        weights = wfile.readlines()

        nextrabins = 1
        nsmooth = 2

        nbinstot = nbins+2*nextrabins

        if is2nddc == 1:
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

        else:
            raise NotImplementedError("is2nddc should be equal to 1!")

        # rescale the bins so that it takes the middle of each old bin
        bins1 = (bins1[1]-bins1[0])/2 + bins1[:-1]
        bins2 = (bins2[1]-bins2[0])/2 + bins2[:-1]

        self.nbinsfe = nbins
        self.nextrabinsfe = nextrabins
        self.hist_idxs_fe = hist_idxs

        return bins1, bins2, free_energy_grid

    def pick_new_points(self, settings, is2nddc):

        curdir = os.getcwd()

        logging.info("Pick new configurations for the next iteration.")
        if is2nddc == 1:
            ndcs_pre_max = 500
            ndcs_pre = max(settings.nreplicas, ndcs_pre_max) # ndcs used for the preselection
            # preselection
            idxs_dcs_pre = tools.draw_points_hist2D(self.hist_idxs_fe, self.nbinsfe+2*self.nextrabinsfe, ndcs_pre)
            # selection
            if settings.nreplicas <= ndcs_pre_max:
                idxs_new_dcs = []
                random_index = random.randrange(0, ndcs_pre)
                idxs_new_dcs.append(idxs_dcs_pre.pop(random_index))
                for count in xrange(settings.nreplicas-1):
                    max_min_r = 0.
                    for i, idx in enumerate(idxs_dcs_pre):
                        min_r = 1.e100
                        for jdx in idxs_new_dcs:
                            r = (self.dc1s_embed[idx] - self.dc1s_embed[jdx])**2 + (self.dc2s_embed[idx] - self.dc2s_embed[jdx])**2
                            min_r = min(r, min_r)
                        if min_r >= max_min_r:
                            max_min_r = min_r
                            new_i = i
                            new_idx = idx
                    del idxs_dcs_pre[new_i]
                    idxs_new_dcs.append(new_idx)
            else:
                idxs_new_dcs = idxs_dcs_pre
        else:
            raise NotImplementedError("is2nddc should be equal to 1!")

        new_coords = self.coords_all[idxs_new_dcs]

        logging.info("Save new configurations in " + curdir + '/' + 'output.gro')
        # save new coordinates
        grofile_w = writer.open('.gro', pattern='confall.gro')
        grofile_w.write(new_coords, 'output.gro')

        # save dcs of new points (check)
	np.savetxt("output.ev", np.array([self.dc1s_embed[idxs_new_dcs], self.dc2s_embed[idxs_new_dcs]]).T, fmt='%15.7e')
