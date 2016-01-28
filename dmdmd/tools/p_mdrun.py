import os
import sys
import time
import glob
import shutil
import argparse
import subprocess
import ConfigParser
import logging

class ParallelMDruns(object):

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run parallel mdrun..")

        # required options

        parser.add_argument("-f",
            type=str,
            dest="inifile",
            required=True)
        
        parser.add_argument("-c",
            type=str, 
            dest="grofile",
            required=True)
        
        parser.add_argument("-o",
            type=str,
            dest="output_file",
            default="out.gro")

        return parser

    def run(self):
        logging.basicConfig(filename='md.log',
                            filemode='w',
                            format="%(levelname)s:%(name)s:%(asctime)s: %(message)s",
                            datefmt="%H:%M:%S",
                            level=logging.DEBUG)
        
      
        parser = self.create_arg_parser()
        args = parser.parse_args() # set argument parser

        config = ConfigParser.SafeConfigParser()
        config.read('lsdmap.ini')


        tcpu0 = time.time()

        # preprocessing
        print "Preprocessing MD..."
        nthreads=config.getint('MD','nthreads')
        nthreads_gromacs=config.getint('MD','nthreads_gromacs')
        nnodes=nthreads/nthreads_gromacs
        ##gmx_suffix
        ##mpi_command
        if config.has_option('MD','additional'):
          additional=config.get('MD','additional')[1:-1]
        else:
          additional=''
        logging.info('Starting with #: '+additional)
        set_pid=set()
        p2=subprocess.call("rm -r tmp",shell=True)
        p2=subprocess.call("mkdir tmp",shell=True)

        print "Creating subdirectories"
        for index in range(nnodes):    
          p2=subprocess.call("mkdir tmp/core"+str(index) ,shell=True)


        grofile = open(args.grofile, 'r')            
        grofile.next()
        natoms = int(grofile.next())
        for idx, line in enumerate(grofile):
                pass
        nlines = idx + 3
        ncoords = nlines/(natoms+3) 
        grofile.close()

        tcpu1 = time.time()
        print "Time used for preprocessing (parallelization): %.2fs" %(tcpu1 - tcpu0)
        print "Run GROMACS preprocessing and MD..."  
        current_n_coord = 0    
        max_index=nnodes
        p2=subprocess.call('rm '+str(args.output_file),shell=True)
        while(current_n_coord < ncoords):
          set_pid=set()
          logging.info('Starting with #: '+str(current_n_coord))
          for index in range(nnodes): 
            if current_n_coord < ncoords:
               start=(natoms+3)*(current_n_coord)+1
               end=(natoms+3)*(current_n_coord+1)
               p=subprocess.Popen("sed '"+str(start)+','+str(end)+"!d' "+str(args.grofile)+' > tmp/core'+str(index)+"/tmp.gro; cd tmp/core"+str(index) +" ; aprun -n  1 -N  1 -d "+str(nthreads_gromacs)+" grompp_mpi "+additional+" -f ../../grompp.mdp -c tmp.gro -p ../../topol.top -maxwarn 2; aprun -n 1 -N 1 -d "+str(nthreads_gromacs)+" mdrun_mpi 1>/dev/null 2>/dev/null; rm \#*",shell=True)
               #logging.info("sed '"+str(start)+','+str(end)+"!d' "+str(args.grofile)+' > tmp/core'+str(index)+"/tmp.gro; cd tmp/core"+str(index) +"; aprun -n 1 -N 1 -d "+str(nthreads_gromacs)+" grompp_mpi "+additional+" -f ../../grompp.mdp -c tmp.gro -p ../../topol.top -maxwarn 2 1>/dev/null 2>/dev/null; aprun -n 1 -N 1 -d "+str(nthreads_gromacs)+" mdrun_mpi 1>/dev/null 2>/dev/null; rm \#*")
               set_pid.add(p)
               current_n_coord=current_n_coord+1
            else:
              max_index=index
              break

          for pid in set_pid:
            pid.communicate()

          for index in range(nnodes): 
            if index<max_index:
              p2=subprocess.call('cat tmp/core'+str(index)+'/confout.gro >> '+str(args.output_file),shell=True)

        # post processing
        tcpu2 = time.time() 
        print "Time used for MD: %.2fs" %(tcpu2 - tcpu1)
        
       
        return
if __name__ == '__main__':
    ParallelMDruns().run()
