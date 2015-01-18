import os
import sys
import time
import imp
import argparse
import ConfigParser
import logging
import radical.pilot
import numpy as np

from dmaps.tools import pilot
from dmaps.tools.config import platforms
from dmaps.critical import kernel as dmsk

class DMapSamplingConfig(object):

    def __init__(self, settings):

        if settings.startgro == 'input.gro':
            raise ValueError("name input.gro is already allocated for md input file and thus can not be assigned to variable startgro! \
	                      Please use a different name")

        if hasattr(settings, "inifile"):
            inifile = settings.inifile
        else:
            inifile = "config.ini"
            settings.inifile = inifile

        if not os.path.isfile(inifile):
            logging.error(".ini file does not exist:" + inifile)
            raise IOError(".ini file does not exist:" + inifile)

        if sys.platform in platforms['mac']:
            self.sedarg = " '' "
        else:
            self.sedarg = " "

        if settings.iter == 0:
            with open('input.gro', 'w') as ifile:
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
        self.nsteps = config.getint('MD', 'nsteps')
        if self.nsteps%self.nstride != 0:
            raise ValueError("nstride (number of configurations saved per replica) should be a multiple of nsteps, please update file " + settings.inifile)
        os.system("sed -i" + self.sedarg + "'s/nsteps.*/nsteps                   = %i/g' "%self.nsteps + settings.mdpfile)
        logging.info("nsteps for MD simulations was set to %i "%self.nsteps)

        # first iteration?
        self.isfirst = config.getint('DMAPS', 'isfirst')

        # number of first dcs used
        self.ndcs = config.getint('DMAPS', 'ndcs')

        # total number of configurations
        self.npoints = settings.nreplicas * self.nstride

        # number of configurations used to compute lsdmap
        nlsdmap = config.getint('LSDMAP', 'nlsdmap')
        if nlsdmap > self.npoints:
            logging.warning("the number of configs required for LSDMap (%i) is larger than the total number of configs expected every iteration (%i)"%(nlsdmap, self.npoints))
            logging.warning("set the number of configs used for LSDMap to %i" %self.npoints)
            self.nlsdmap = self.npoints
        else:
            self.nlsdmap = nlsdmap

        # temperature in Kelvins
        temperature = config.getint('MD', 'temperature')
        kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
        self.kT = kb*temperature

        try:
            # solvation?
            self.solvation = config.get('MD','solvation').lower()
        except:
            self.solvation = 'none'

        # update temperature in .mdp file
        if self.solvation in ['none', 'implicit']: 
            os.system("sed -i" + self.sedarg + "'s/ref_t.*/ref_t               = %i/g' "%temperature + settings.mdpfile)
        elif self.solvation == 'explicit':
            os.system("sed -i" + self.sedarg + "'s/ref_t.*/ref_t               = %i %i/g' "%(temperature, temperature) + settings.mdpfile)
        else:
	    raise ValueError("solvation option in configuration file not understood, should one among none, implicit and explicit")
        os.system("sed -i" + self.sedarg + "'s/gen_temp.*/gen_temp                = %i/g' "%temperature + settings.mdpfile)
 
        # number of points used for the fitting of the dc values
        self.nfit = config.getint('FITTING', 'npoints')

        # cutoff when computing the free energy
        ncutoff =  config.getint('DMAPS', 'ncutoff')
        self.cutoff = ncutoff*self.kT

        # get ctram parameters
        if config.has_section("CTRAM"):
            self.isctram = config.getint("CTRAM", "isctram")
            if self.isctram == 1:
                self.ntau_ctram = config.getint("CTRAM", "ntau")
                self.nstates_ctram = config.getint("CTRAM", "nstates")
                self.niters_ctram = config.getint("CTRAM", "niters")
                if config.has_option("nbins", "CTRAM"):
                    self.nbins_ctram = config.getint("nbins", "CTRAM")
                else:
                    self.nbins_ctram = None
                # check if ntau is larger than nstride
                if self.ntau_ctram > self.nstride-1:
                    raise ErrorValue("ntau should be less than nstride, please update config file!")
                else:
                    nsamples = (self.nstride-1)/self.ntau_ctram*settings.nreplicas*self.nstates_ctram
                    logging.info("cTRAM: max number of sample points: %i"%nsamples)
        else:
            self.isctram = 0

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

natoms=$(sed -n '2p' $startgro)
nlines_per_frame=$((natoms+3))

nlines=`wc -l $startgro| cut -d' ' -f1`
nframes=$((nlines/nlines_per_frame))

for idx in `seq 1 $nframes`; do

  start=$(($nlines_per_frame*(idx-1)+1))
  end=$(($nlines_per_frame*idx))
  sed "$start"','"$end"'!d' $startgro > $tmpstartgro

  # gromacs preprocessing & MD
  grompp -f %(mdpfile)s -c $tmpstartgro -p %(topfile)s %(grompp_options)s &> /dev/null
  mdrun -nt 1 -dms %(inifile)s -s topol.tpr %(mdrun_options)s &> mdrun.log

done

# remove temporary files
rm -f $tmpstartgro
        """ % locals()
            file.write(script)


class DMapSamplingExe(object):

    def update(self, args, settings, config):
    
        # before updating save data
        os.system('rm -rf iter%i'%settings.iter)
        os.system('mkdir iter%i'%settings.iter)
        os.system('cp ' + ' '.join(['confall.gro', 'confall.w', 'confall.ev.embed', 'output.gro', 'output.ev']) + ' iter%i'%settings.iter)
        os.system('cp -r ' + ' '.join(['lsdmap', 'fit', 'fe', 'ctram']) + ' iter%i'%settings.iter)

        if settings.iter > 0:
            os.system('cp confall.ev.embed.old iter%i'%settings.iter)
    
        os.system('mv output.gro input.gro')
        os.system('sed -i' + config.sedarg + "'s/iter=.*/iter=%i/g' "%(settings.iter+1) + args.setfile)
        os.system('sed -i' + config.sedarg + "'s/isfirst=.*/isfirst=0/g' " + settings.inifile)
    
        return settings.iter+1


    def restart_from_iter(self, num_iter, args):
    
        if sys.platform in platforms['mac']:
            sedarg = " '' "
        else:
            sedarg = " "

        logging.info("restarting from iteration %i" %num_iter)
        # remove iter folders with iter number >= num_iter
        os.system('for dir in iter*; do num=`cut -d "r" -f 2 <<< $dir`; if [ "$num" -ge %i ]; then rm -rf $dir; fi ; done'%num_iter)
        os.system('rm -rf ' + ' '.join(['md','lsdmap', 'fit', 'fe', 'ctram']))
        os.system('cp iter%i/output.gro input.gro'%(num_iter-1))
        os.system('cp -r iter%i/{'%(num_iter-1) + ','.join(['lsdmap', 'fit', 'fe', 'ctram']) + '} .')
    
        # update iter in settings file
        os.system('sed -i' + sedarg + "'s/iter=.*/iter=%i/g' "%num_iter + args.setfile)
    
    def restart(self, args):

        if sys.platform in platforms['mac']:
            sedarg = " '' "
        else:
            sedarg = " "

        # update iter in settings file
        os.system('sed -i' + sedarg + "'s/iter=.*/iter=0/g' " + args.setfile)
        os.system('rm -rf iter* ' + ' '.join(['lsdmap', 'fit', 'fe', 'md', 'ctram']))
    
    def run(self):

        parser = argparse.ArgumentParser(description="Run Diffusion Map Sampling...")
        parser.add_argument("-f", type=str, dest="setfile", required=True, help='File containing settings (input): -')
        parser.add_argument("--restart", action="store_true", dest="restart", default=False, help='restart from scratch')
        parser.add_argument("--checkpoint", type=int, dest="num_iter", help='restart from a given iteration')
        
        args = parser.parse_args()
        
        logging.basicConfig(filename='dmaps.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)
        
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
        config = DMapSamplingConfig(settings) 
        umgr, session = pilot.startPilot(settings)
        
        # main loop
        for idx in xrange(settings.niters):
            logging.info("START ITERATION %i"%settings.iter)
            # run biased MD
            dmapsworker = dmsk.DMapSamplingWorker()
            dmapsworker.run_md(umgr, settings)
            # run LSDMap and fit
            dmapsworker.run_lsdmap(umgr, settings, config)
            dmapsworker.run_fit(umgr, settings, config) # fit configurations of the current iteration
            # compute the free energy
            dmapsworker.do_free_energy(umgr, settings, config)
            # select the new configurations for the next iteration
            dmapsworker.select_new_points(settings, config)
            # update for next iteration
            settings.iter = self.update(args, settings, config)
           
        session.close()
