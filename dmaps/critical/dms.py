import os
import sys
import time
import imp

import stat
import glob
import shutil

import argparse
import ConfigParser
import logging
import subprocess
import numpy as np

from dmaps.tools.config import platforms
from dmaps.critical import kernel as dmsk

class DMapSamplingConfig(object):

    def __init__(self, iter, args):

        if sys.platform in platforms['mac']:
            self.sedarg = " '' "
        else:
            self.sedarg = " "

        config = ConfigParser.SafeConfigParser()
        config.read(args.config_file)

        self.iter = iter

        self.load_parameters(config, args)
        self.initialize_md(config, args)
        self.metric = config.get('LSDMAP', 'metric')
        
    def check_md_file(self, config, option, default=None):

        if config.has_option('MD', option):
            filename = config.get('MD', option)
        elif default is not None:
            filename = default
        else:
            raise KeyError('Option ' + option + ' in the config file is mandatory!')

        if not os.path.isfile(filename):
            logging.error("file does not exist:" + filename)
            raise IOError("file does not exist:" + filename)

        return filename

    def initialize_md(self, config, args):

        # get the name of the starting MD file
        self.startgro = self.check_md_file(config, 'startgro')

        if self.startgro == 'input.gro':
            raise ValueError("name input.gro is already allocated for md input file and thus can not be assigned to variable startgro! \
                              Please use a different name")

        topfile = self.check_md_file(config, 'topfile', default='topol.top')
        mdpfile = self.check_md_file(config, 'mdpfile', default='grompp.mdp')

        # check grompp options
        grompp_options = ''
        if config.has_option('MD', 'ndxfile'):
            grompp_options += ' -n ../../' + config.get('MD', 'ndxfile')
        else:
            grompp_options += ''

        if config.has_option('MD', 'grompp_options'):
            grompp_options += ' ' + config.get('MD', 'grompp_options')
        else:
            grompp_options += ''

        mdrun_options = ''
        if config.has_option('MD', 'mdrun_options'):
            mdrun_options += ' ' + config.get('MD', 'mdrun_options')
        else:
            mdrun_options += ''

        self.write_md_script('run_md.sh', args.config_file, topfile, mdpfile, grompp_options=grompp_options, mdrun_options=mdrun_options)

        # set number of steps in MD simulations
        self.nsteps = config.getint('MD', 'nsteps')
        if self.nsteps%self.nstride != 0:
            raise ValueError("nstride (number of configurations saved per replica) should be a multiple of nsteps, please update file " + args.config_file)

        subprocess.check_call("sed -i" + self.sedarg + "'s/nsteps.*/nsteps                   = %i/g' "%self.nsteps + mdpfile, shell=True)
        logging.info("nsteps for MD simulations was set to %i "%self.nsteps)

        # number of MD steps skipped for the computation of the bias potential
        if config.has_option('MD', 'nstepbias'):
            self.nstepbias = config.getint('MD', 'nstepbias')
            nsave = self.nsteps/self.nstride
            if nsave%self.nstepbias != 0:
                raise ValueError("the number of steps skipped when applying the biased force must divide exactly\
                    (with no remainder) the number of steps in between two consecutive config savings")
        else:
            self.nstepbias = 1

        try:
            self.solvation = config.get('MD','solvation').lower()
        except:
            self.solvation = 'none'

    def write_md_script(self, filename, inifile, topfile, mdpfile, grompp_options=None, mdrun_options=None):

        # write script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash

# this script was automatically created when using DMap Sampling

startgro=input.gro
tmpstartgro=tmp.gro

natoms=$(sed -n '2p' $startgro)
nlines_per_frame=$((natoms+3))

nlines=`wc -l $startgro| cut -d' ' -f1`
nframes=$((nlines/nlines_per_frame))

for idx in `seq 1 $nframes`; do

  start=$(($nlines_per_frame*(idx-1)+1))
  end=$(($nlines_per_frame*idx))
  sed "$start"','"$end"'!d' $startgro > $tmpstartgro

  # gromacs preprocessing & MD
  grompp -f ../../%(mdpfile)s -c $tmpstartgro -p ../../%(topfile)s %(grompp_options)s &> grompp.log
  mdrun -nt 1 -dms ../../%(inifile)s %(mdrun_options)s &> mdrun.log

done

# remove temporary files
rm -f $tmpstartgro
        """ % locals()
            file.write(script)

        os.chmod(filename, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

    def load_parameters(self, config, args):

        inifile = args.config_file

        # number of configurations saved per replica
        self.nstride = config.getint('GENERAL', 'nstride')

        # total number of configurations
        self.nreplicas = args.nreplicas
        self.npoints = self.nreplicas * self.nstride

        # number of first dcs used
        self.ndcs = config.getint('GENERAL', 'ndcs')

        # cutoff when computing the free energy
        ncutoff =  config.getint('GENERAL', 'ncutoff')

        # temperature in Kelvins
        temperature = config.getint('GENERAL', 'temperature')
        kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
        self.kT = kb*temperature

        self.cutoff = ncutoff*self.kT

        # fraction of the free energy we actually use for the bias potential
        if config.has_option('GENERAL', 'fefrac'):
            self.fefrac = config.getfloat('GENERAL', 'fefrac')
        else:
            self.fefrac = 1.0

        # number of bins used to compute the free energy histogram
        if config.has_option('GENERAL', 'nbinsfe'):
            self.nbinsfe = config.getint('GENERAL', 'nbinsfe')
        else:
            self.nbinsfe = None

        # check if uniform sampling is disabled when restarting
        if config.has_option('DMAPS', 'uniform_sampling'):
            self.uniform_sampling = config.getint('DMAPS', 'uniform_sampling')
            if self.uniform_sampling not in [0,1]:
                raise ValueError("option uniform_sampling should be 0 or 1!")
        else:
            self.uniform_sampling = 1

        # number of configurations used to compute lsdmap
        nlsdmap = config.getint('LSDMAP', 'nlsdmap')
        if nlsdmap > self.npoints:
            logging.warning("the number of configs required for LSDMap (%i) is larger than the total number of configs expected every iteration (%i)"%(nlsdmap, self.npoints))
            logging.warning("set the number of configs used for LSDMap to %i" %self.npoints)
            self.nlsdmap = self.npoints
        else:
            self.nlsdmap = nlsdmap

        # number of points used for the fitting of the dc values
        self.nfit = config.getint('FITTING', 'npoints')

        # get ctram parameters
        if config.has_section('CTRAM'):
            self.isctram = config.getint('CTRAM', 'isctram')
            if self.isctram == 1:
                self.ntau_ctram = config.getint('CTRAM', 'ntau')
                self.nstates_ctram = config.getint('CTRAM', 'nstates')
                self.niters_ctram = config.getint('CTRAM', 'niters')
                if config.has_option('CTRAM', 'nbins'):
                    self.nbins_ctram = config.getint('CTRAM', 'nbins')
                else:
                    self.nbins_ctram = None
                # check if ntau is larger than nstride
                if self.ntau_ctram > self.nstride-1:
                    raise ErrorValue("ntau should be less than nstride, please update config file!")
                else:
                    nsamples = (self.nstride-1)/self.ntau_ctram*self.nreplicas*self.nstates_ctram
                    logging.info("cTRAM: max number of sample points: %i"%nsamples)
        else:
            self.isctram = 0

class DMapSamplingExe(object):

    def update(self, iter, args, config):

        # before updating save data
        subprocess.check_call('rm -rf iter%i'%iter, shell=True)
        subprocess.check_call('mkdir iter%i'%iter, shell=True)
        subprocess.check_call('cp ' + ' '.join(['confall.gro', 'confall.w', 'confall.ev.embed', 'output.gro', 'output.ev']) + ' iter%i'%iter, shell=True)
        subprocess.check_call('cp -r ' + ' '.join(['lsdmap', 'fit', 'fe']) + ' iter%i'%iter, shell=True)

        if iter > 0:
            subprocess.check_call('cp confall.ev.embed.old iter%i'%iter, shell=True)
            if config.isctram == 1:
                subprocess.check_call('cp -r ctram iter%i'%iter, shell=True)

        # at the end of the iteration, the output .gro file becomes the input .gro file of the new iteration
        subprocess.check_call('mv output.gro input.gro', shell=True)
        # change the number of the current iteration in the config file
        subprocess.check_call('sed -i' + config.sedarg + "'s/iter=.*/iter=%i/g' "%(iter+1) + args.config_file, shell=True)

        return iter+1

    def restart_from_iter(self, iter, args):
    
        if sys.platform in platforms['mac']:
            sedarg = " '' "
        else:
            sedarg = " "

        logging.info("restarting from iteration %i" %iter)
        # remove iter folders with iter number >= iter
        for dirname in glob.glob("iter*"):
            num = int(dirname[4:])
            if num >= iter:
                shutil.rmtree(dirname)
        # remove the folders of the current iteration
        for dirname in ['md', 'lsdmap', 'fit', 'fe', 'ctram']:
            if os.path.exists(dirname):
                shutil.rmtree(dirname)

        # copy the files needed to continue from the iteration num_iter
        shutil.copyfile("iter%i/output.gro"%(iter-1), "input.gro")
        for folder in ['lsdmap', 'fit', 'fe', 'ctram']:
            dirname = "iter%i/"%(iter-1)+folder
            if os.path.exists(dirname):
                shutil.copytree(dirname, folder)
    
        # update iter in the config file
        subprocess.check_call('sed -i' + sedarg + "'s/iter=.*/iter=%i/g' "%iter + args.config_file, shell=True)
    
    def restart(self, args):

        if sys.platform in platforms['mac']:
            sedarg = " '' "
        else:
            sedarg = " "

        # update iter in the config file
        subprocess.check_call('sed -i' + sedarg + "'s/iter=.*/iter=0/g' " + args.config_file, shell=True)

        # remove iter folders with iter number >= num_iter
        for dirname in glob.glob("iter*"):
            shutil.rmtree(dirname)

        for dirname in ['md', 'lsdmap', 'fit', 'fe', 'ctram']:
            if os.path.exists(dirname):
                shutil.rmtree(dirname)

    def run(self):

        parser = argparse.ArgumentParser(description="Run Diffusion Map Sampling...")

        parser.add_argument('ncpus',
            metavar='ncpus',
            type=int,
            help='total number of cpus used')

        parser.add_argument('nreplicas',
            metavar='nreplicas',
            type=int,
            help='total number of replicas')

        parser.add_argument("-f",
            type=str,
            dest="config_file",
            default='config.ini',
            help='File containing parameters (input): -')

        parser.add_argument("-n",
            type=int,
            dest='niters',
            default=1e20,
            help='number of iterations')

        parser.add_argument("--restart",
            action="store_true",
            dest="restart",
            default=False,
            help='restart from scratch')

        parser.add_argument("--checkpoint",
            type=int,
            dest="iter_cp",
            help='restart from a given iteration')

        parser.add_argument("--skipmd",
            action="store_true",
            dest="skipmd",
            default=False,
            help='skip first MD simulations')

        args = parser.parse_args()

        logging.basicConfig(filename='dmaps.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        if args.restart:
            self.restart(args)
            if args.iter_cp is not None:
                logging.error("checkpoint option can not be set together with restart option")
                raise ValueError("checkpoint option can not be set together with restart option")
            if args.skipmd:
                logging.error("skipmd option can not be set together with restart option")
                raise ValueError("skipmd option can not be set together with restart option")
            iter = 0

        elif args.iter_cp is not None:
            if args.iter_cp == 0:
                logging.error("checkpoint option can not be set to 0, use restart option instead")
            elif args.iter_cp > 0:
                self.restart_from_iter(args.iter_cp, args)
                if args.skipmd:
                    logging.error("skipmd option can not be set together with checkpoint option")
                    raise ValueError("skipmd option can not be set together with checkpoint option")
            else:
                logging.error("argument of checkpoint option should be a positive integer (iteration number to restart from)")
                raise ValueError("argument of checkpoint option should be a positive integer (iteration number to restart from)")
            iter = args.iter_cp

        else:
            config = ConfigParser.SafeConfigParser()
            config.read(args.config_file)
            # if restart or checkpoint options are disabled, restart from the iteration specified in the config file
            iter = config.getint('GENERAL', 'iter')
            if iter == 0:
                self.restart(args)
            else:
                self.restart_from_iter(iter, args)
        config = DMapSamplingConfig(iter, args)
        t_start=time.time()
        # main loop
        for idx in xrange(args.niters):
            logging.info("START ITERATION %i"%config.iter)
            print 'Iteration %i\n'%config.iter
            # run biased MD
            dmapsworker = dmsk.DMapSamplingWorker()
            if idx > 0 or not args.skipmd:
                dmapsworker.run_md(config, args)
            # run LSDMap and fit
            dmapsworker.run_lsdmap(config, args)
            dmapsworker.run_fit(config, args) # fit configurations of the current iteration
            # compute the free energy
            dmapsworker.do_free_energy(config)
            # select the new configurations for the next iteration
            dmapsworker.select_new_points(config)
            # update for next iteration
            config.iter = self.update(config.iter, args, config)
            t_current=time.time()
            subprocess.call("echo "+str(config.iter)+" " +str((t_current-t_start)*args.ncpus/3600) +" >> dmaps_cputime.log", shell=True) 
            print '--------------------------------------------------'
