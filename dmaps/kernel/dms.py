import os
import sys
import time
import imp
import argparse
import subprocess
import ConfigParser

from math import sqrt, floor
import logging
import radical.pilot
import numpy as np

from dmaps.tools import pilot
from dmaps.tools.config import platforms, tmpfiles
from dmaps.kernel import worker

class DMapSampling(object):

    def initialize(self, settings):

        for key, value in tmpfiles.iteritems():
            if settings.startgro == value[0]:
                raise ValueError("name " + settings.startgro + " is already allocated for " + value[1] + " and thus can not be assigned to variable startgro! \
	                      Please use a different name")

        if hasattr(settings, "inifile"):
            inifile = settings.inifile
        else:
            inifile = "config.ini"
            settings.inifile = inifile

        if not os.path.isfile(inifile):
            logging.error(".ini file does not exist:" + inifile)
            raise IOError(".ini file does not exist:" + inifile)

        if settings.iter == 0:
            with open(tmpfiles['ingro'][0], 'w') as ifile:
                for idx in xrange(settings.nreplicas):
                    with open(settings.startgro, 'r') as sfile:
                       for line in sfile:
                          print >> ifile, line.replace("\n", "")
            os.system("sed -i" + self.sedarg + "'s/isfirst=.*/isfirst=1/g' " + settings.inifile)
        else:
            os.system("sed -i" + self.sedarg + "'s/isfirst=.*/isfirst=0/g' " + settings.inifile)

        # load parameters inside inifile
        self.load_parameters(settings)

        # create executables
        self.write_md_script("run_md.sh", settings)

    def load_parameters(self, settings):

        config = ConfigParser.SafeConfigParser()
        config.read(settings.inifile)

        if hasattr(settings, "mdpfile"):
            mdpfile = settings.mdpfile
        else:
            mdpfile = "grompp.mdp"
            settings.mdpfile = mdpfile

        if not os.path.isfile(mdpfile):
            logging.error(".mdp file does not exist:" + mdpfile)
            raise IOError(".mdp file does not exist:" + mdpfile)

        # number of configurations saved per replica
        self.nstride = config.getint('DMAPS', 'nstride')

        # set number of steps in MD simulations
        self.nsteps_min = 10000

        if config.has_option('MD', 'nsteps'):
            self.nsteps = config.getint('MD', 'nsteps')
            os.system("sed -i" + self.sedarg + "'s/nsteps.*/nsteps                   = %i/g' "%self.nsteps + settings.mdpfile)
            logging.info("nsteps for MD simulations was set to %i "%self.nsteps)
            self.nsteps_guess = False
        else:
            if settings.iter == 0:
                # if first iteration, the number of steps is set to 10000 or so
                self.nsteps = self.nstride*(self.nsteps_min/self.nstride)
                # update nsteps line in .mdp file
                os.system("sed -i" + self.sedarg + "'s/nsteps.*/nsteps                   = %i/g' "%self.nsteps + settings.mdpfile)
                logging.info("nsteps for MD simulations was set to %i (default value)"%self.nsteps)
            else:
                self.nsteps = int(subprocess.check_output("cat " + settings.mdpfile + " | sed -n -e 's/^.*nsteps.*=//p' | tr -d ' '", shell=True))
            self.nsteps_guess = True

        # first iteration?
        self.isfirst = config.getint('DMAPS', 'isfirst')

        # total number of configurations
        self.npoints = settings.nreplicas * self.nstride

        # number of configurations used to compute lsdmap
        nlsdmap = config.getint('LSDMAP', 'nlsdmap')
        if nlsdmap > self.npoints:
            logging.warning("number of configs required for LSDMap (%i) is larger than the total number of configs expected every iteration (%i)"%(nlsdmap, self.npoints))
            logging.warning("set the number of configs used for LSDMap to %i" %self.npoints)
            self.nlsdmap = self.npoints
        else:
            self.nlsdmap = nlsdmap

        # number of bins used to build the histogram for the free energy
        self.nbinsfe = int(floor(self.npoints**(1./3)))

        # temperature in Kelvins
        temperature = config.getint('MD', 'temperature')
        kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
        self.kT = kb*temperature

        try:
            # solvation?
            self.solvation = config.get('MD','solvation').lower()
        except:
            self.solvation = 'none'

        if self.solvation in ['none', 'implicit']: 
            os.system("sed -i" + self.sedarg + "'s/ref_t.*/ref_t               = %i/g' "%temperature + settings.mdpfile)
        elif self.solvation == 'explicit':
            os.system("sed -i" + self.sedarg + "'s/ref_t.*/ref_t               = %i %i/g' "%(temperature, temperature) + settings.mdpfile)
        else:
	    raise ValueError("solvation option in configuration file not understood, should one among none, implicit and explicit")

        os.system("sed -i" + self.sedarg + "'s/gen_temp.*/gen_temp                = %i/g' "%temperature + settings.mdpfile)
 
        # cutoff when computing the free energy
        ncutoff =  config.getint('DMAPS', 'ncutoff')
        self.cutoff = ncutoff*self.kT

        # number of points used for the fitting of the dc values
        self.nfitdcs = config.getint('FITTING', 'npoints')

        # number of bins used to build the histogram to select the fitting points
        self.nbinsdcs = int(floor(self.nlsdmap**(1./3)))

        # number of bins used to build the histogram to select the lsdmap points
        self.nbins_lsdmap = int(floor(self.npoints**(1./3)))


    def write_md_script(self, filename, settings):

        inifile = settings.inifile
        mdpfile = settings.mdpfile

        # check topfile
        if hasattr(settings, "topfile"):
            topfile = settings.topfile
        else:
            topfile = "topol.top"
            settings.topfile = topfile

        if not os.path.isfile(topfile):
            logging.error(".top file does not exist:" + topfile)
            raise IOError(".top file does not exist:" + topfile)

        # check grompp and mdrun options
        if hasattr(settings, "grompp_options"):
            grompp_options = settings.grompp_options
        else:
            grompp_options = ""

        if hasattr(settings, "mdrun_options"):
            mdrun_options = settings.mdrun_options
        else:
            mdrun_options = ""

        # write script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash

# this script was automatically created when using DMap Sampling

startgro=input.gro
tmpstartgro=tmp.gro
outgro=out.gro

natoms=$(sed -n '2p' $startgro)
nlines_per_frame=$((natoms+3))

nlines=`wc -l $startgro| cut -d' ' -f1`
nframes=$((nlines/nlines_per_frame))

rm -rf $outgro

for idx in `seq 1 $nframes`; do

  start=$(($nlines_per_frame*(idx-1)+1))
  end=$(($nlines_per_frame*idx))
  sed "$start"','"$end"'!d' $startgro > $tmpstartgro

  # gromacs preprocessing & MD
  grompp -f %(mdpfile)s -c $tmpstartgro -p %(topfile)s %(grompp_options)s &> /dev/null
  mdrun -nt 1 -dms %(inifile)s -s topol.tpr %(mdrun_options)s &> mdrun.log

  # store data
  cat confout.gro >> $outgro

done

# remove temporary files
rm -f $tmpstartgro
        """ % locals()
            file.write(script)


    def update(self, args, settings):

        # before updating save data
        os.system("rm -rf iter%i"%settings.iter)
        os.system("mkdir iter%i"%settings.iter)
        os.system("cp confall.gro confall.w confall.ev.embed output.ev iter%i"%settings.iter)
        os.system("cp output.ev " + tmpfiles['outgro'][0] + " iter%i" %settings.iter)
        os.system("cp -r fit fe lsdmap iter%i"%settings.iter)       
 
        if settings.iter > 0:
	    os.system("cp confall.ev.embed.old iter%i"%settings.iter)

        autocorrelation_time_dc1, autocorrelation_time_dc2 = self.get_dcs_autocorrelation_time(settings)

        if self.nsteps_guess:
            min_autocorrelation_time_dc = min(autocorrelation_time_dc1, autocorrelation_time_dc2)
            self.nsteps = max(self.nstride*min_autocorrelation_time_dc/settings.nreplicas, self.nstride*(self.nsteps_min/self.nstride))
            os.system("sed -i" + self.sedarg + "'s/nsteps.*/nsteps                   = %i/g' "%self.nsteps + settings.mdpfile)
            logging.info("nsteps for next MD simulations has been set to %i"%self.nsteps)

        os.system("mv " + tmpfiles['outgro'][0] + ' ' + tmpfiles['ingro'][0])
        os.system("sed -i" + self.sedarg + "'s/iter=.*/iter=%i/g' "%(settings.iter+1) + args.setfile)
        os.system("sed -i" + self.sedarg + "'s/isfirst=.*/isfirst=0/g' " + settings.inifile)

        return settings.iter+1


    def get_dcs_autocorrelation_time(self, settings):

        # if first iteration, compute autocorrelation time from computed dcs 
        if settings.iter == 0:
            os.system('cp confall.ev.embed autocorr.ev')

        dcs = np.loadtxt('autocorr.ev')
        ndcs = dcs.shape[0]
        ndcs_per_replica = ndcs/settings.nreplicas
       
        step = np.linspace(0, self.nsteps, ndcs_per_replica+1).astype(int)
        step = step[:-1]
        autocorrelation_dc1 = np.zeros(ndcs_per_replica)
        autocorrelation_dc2 = np.zeros(ndcs_per_replica)

        logging.info("Compute DCs autocorrelation time")

        for idx in xrange(settings.nreplicas):
       
            # load dc1s values
            dc1s = dcs[idx*ndcs_per_replica:(idx+1)*ndcs_per_replica, 0]
        
            vardc1s = dc1s.var()
            meandc1s = dc1s.mean()
            dc1s -= meandc1s
        
            autocorrelation_tmp = np.correlate(dc1s, dc1s, mode='full')[-ndcs_per_replica:]
            autocorrelation_tmp /= vardc1s*(np.arange(ndcs_per_replica, 0, -1))
            autocorrelation_dc1 += autocorrelation_tmp
        
            # load dc2s values
            dc2s = dcs[idx*ndcs_per_replica:(idx+1)*ndcs_per_replica, 1]
       
            vardc2s = dc2s.var()
            meandc2s = dc2s.mean()
            dc2s -= meandc2s
        
            autocorrelation_tmp = np.correlate(dc2s, dc2s, mode='full')[-ndcs_per_replica:]
            autocorrelation_tmp /= vardc2s*(np.arange(ndcs_per_replica, 0, -1))
            autocorrelation_dc2 += autocorrelation_tmp
       
        for idx in xrange(ndcs_per_replica):
            if autocorrelation_dc1[idx] <= autocorrelation_dc1[0]/2:
                break
        autocorrelation_time_dc1 = step[idx]

        for idx in xrange(ndcs_per_replica):
            if autocorrelation_dc2[idx] <= autocorrelation_dc2[0]/2:
                break
        autocorrelation_time_dc2 = step[idx]

        logging.info("Autocorrelation time along DC1: %i steps"%autocorrelation_time_dc1)
        logging.info("Autocorrelation time along DC2: %i steps"%autocorrelation_time_dc2)

        autocorrdir = 'autocorr'
        os.system("rm -rf " + autocorrdir)
        os.makedirs(autocorrdir)

        np.savetxt(autocorrdir + '/' + 'dcs.xyz', np.array([step, autocorrelation_dc1, autocorrelation_dc2]).T, fmt='%15.7e')

        return autocorrelation_time_dc1, autocorrelation_time_dc2


    def restart_from_iter(self, num_iter, args):

        logging.info("restarting from iteration %i" %num_iter)
        # remove iter folders with iter number >= num_iter
        os.system('for dir in iter*; do num=`cut -d "r" -f 2 <<< $dir`; if [ "$num" -ge %i ]; then rm -rf $dir; fi ; done'%num_iter)
        os.system("cp iter%i/"%(num_iter-1) + tmpfiles['outgro'][0] + " " + tmpfiles['ingro'][0])
        os.system("rm -rf fit lsdmap autocorr tmp fe")
        os.system("cp -r iter%i/{fit,fe,lsdmap} ."%(num_iter-1))

        # update iter in settings file
        os.system("sed -i" + self.sedarg + "'s/iter=.*/iter=%i/g' "%num_iter + args.setfile)

        return

    def restart(self, args):
       
        os.system("sed -i" + self.sedarg + "'s/iter=.*/iter=0/g' settings")
        os.system("rm -rf iter* fit lsdmap autocorr tmp fe")

        return

    def run(self):

        parser = argparse.ArgumentParser(description="Run Diffusion Map Driven Adaptive Sampling...")
        parser.add_argument("-f", type=str, dest="setfile", required=True, help='File containing settings (input): -')
        parser.add_argument("--restart", action="store_true", dest="restart", default=False, help='restart from scratch')
        parser.add_argument("--checkpoint", type=int, dest="num_iter", help='restart from a given iteration')

        args = parser.parse_args()

        logging.basicConfig(filename='dmaps.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)

        if sys.platform in platforms['mac']:
            self.sedarg = " '' "
        else:
            self.sedarg = " "

        if args.restart:
            self.restart(args)
            if args.num_iter is not None:
	        logging.error("checkpoint option can not be set together with restart option")

        if args.num_iter is not None:
            if args.num_iter == 0:
                logging.error("checkpoint option can not be set to 0, use restart option instead")
            elif args.num_iter > 0:
                self.restart_from_iter(args.num_iter, args)
            else:
                logging.error("argument of checkpoint option should be a positive integer (iteration number to restart from)")

        settings = imp.load_source('setfile', args.setfile)
	umgr, session = pilot.startPilot(settings)

        # initialize dmap sampling
        self.initialize(settings)

        # main loop
        for idx in xrange(settings.niters):
            logging.info("START ITERATION %i"%settings.iter)
            # run biased MD
            dmapsworker = worker.DMapSamplingWorker()
            dmapsworker.run_md(umgr, settings)
            # run LSDMap
            dmapsworker.run_lsdmap(umgr,settings, self.npoints, self.nlsdmap, self.nbins_lsdmap)
            # run fit
            dmapsworker.run_fit_dcs(umgr, settings, self.nfitdcs, self.nbinsdcs)
            # estimate free energy
            dmapsworker.do_free_energy(self.nbinsfe, self.cutoff, self.kT)
            # pick new points
            dmapsworker.pick_new_points(settings)
            # update for next iteration
            settings.iter = self.update(args, settings)
           
        session.close()
