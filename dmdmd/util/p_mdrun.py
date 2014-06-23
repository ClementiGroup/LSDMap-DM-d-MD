import os
import sys
import time
import glob
import shutil
import argparse

from mpi4py import MPI

class ParallelMDruns(object):

    def write_script(self, script, directory, rank, mdpfile, grofile, topfile, output_grofile, tprfile='topol.tpr', trrfile='traj.trr', edrfile='ener.edr', shebang='/bin/bash', ndxfile='', grompp_options='', mdrun_options=''):

        file_path = directory + '/' + script
    
        if not ndxfile:
            ndxfile_option=''
        else:
            ndxfile_option='-n '+ndxfile

        with open(file_path, 'w') as file:
            script ="""#!%(shebang)s

# this script was generated automatically by thread %(rank)i

startgro=%(grofile)s
tmpstartgro=tmpstart.gro
outgro=%(output_grofile)s

natoms=$(sed -n '2p' $startgro)
nlines_per_frame=$((natoms+3))

nlines=`wc -l %(directory)s/$startgro| cut -d' ' -f1`
nframes=$((nlines/nlines_per_frame))

rm -rf $outgro

for idx in `seq 1 $nframes`; do

  start=$(($nlines_per_frame*(idx-1)+1))
  end=$(($nlines_per_frame*idx))
  sed "$start"','"$end"'!d' $startgro > $tmpstartgro

  # gromacs preprocessing & MD
  grompp %(grompp_options)s -f %(mdpfile)s -c $tmpstartgro -p %(topfile)s %(ndxfile_option)s -o %(tprfile)s 1>/dev/null 2>/dev/null
  mdrun -nt 1 %(mdrun_options)s -s %(tprfile)s -o %(trrfile)s -e %(edrfile)s 1>/dev/null 2>/dev/null

  # store data
  cat confout.gro >> $outgro

done

# remove temporary files
rm -f $tmpstartgro
        """ % locals()
            file.write(script)


    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description="Run parallel mdrun..")

        # required options

        parser.add_argument("-f",
            type=str,
            dest="mdpfile",
            required=True)
        
        parser.add_argument("-c",
            type=str, 
            dest="grofile",
            required=True)
        
        parser.add_argument("-p",
            type=str,
            dest="topfile",
            required=True)
        
        # other options
        
        parser.add_argument("-n",
            type=str,
            dest="ndxfile")
        
        parser.add_argument("-o",
            type=str,
            dest="output_file",
            default="out.gro")
        
        parser.add_argument("--grompp_options",
            type=str, 
            dest="grompp_options",
            default="")
        
        parser.add_argument("--mdrun_options",
            type=str,
            dest="mdrun_options",
            default="")
        
        parser.add_argument("-t",
            action="store",
            type=str,
            dest="tmpdir")
        
        return parser

    def run(self):

        #initialize mpi variables
        comm = MPI.COMM_WORLD   # MPI environment
        size = comm.Get_size()  # number of threads
        rank = comm.Get_rank()  # number of the current thread 

        parser = self.create_arg_parser()
        args = parser.parse_args() # set argument parser

        tcpu0 = time.time()

        # preprocessing
        curdir = os.getcwd()
        if args.tmpdir is None:
            args.tmpdir = curdir + '/tmp'
        rundir = args.tmpdir + '/' + 'thread%i'%rank

        if rank == 0:
            print "Creating subdirectories %s/threadX (1<=X<=%i)..." %(args.tmpdir, size)

        shutil.rmtree(rundir, ignore_errors=True)
        os.makedirs(rundir)
        comm.Barrier()

        if rank==0:
            print "Preparing .gro files..." 
            rundirs=[rundir]

            for idx in range(1,size):
                rundirs.append(comm.recv(source=idx, tag=idx))

            grofile = open(args.grofile, 'r')            
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

            with open(args.grofile, 'r') as grofile:
                for idx in xrange(size):
                    grofile_thread = rundirs[idx] + '/' + 'start.gro'
                    with open(grofile_thread, 'w') as grofile_t:
                        nlines_per_thread = ncoords_per_thread[idx]*(natoms+3)
                        for jdx in xrange(nlines_per_thread):
                            line = grofile.readline()
                            if line:
                                line = line.replace("\n", "")
                                print >> grofile_t, line
                            else:
                                break
        else:
            comm.send(rundir, dest=0, tag=rank)


        if rank==0:
            print "copying .mdp and .top files..."

        shutil.copy(args.mdpfile, rundir+ '/' + 'grompp.mdp')
        shutil.copy(args.topfile, rundir + '/' + 'topol.top')
       
        if args.ndxfile is not None:
            shutil.copy(args.ndxfile, rundir + '/' + args.ndxfile) # copy ndxfile if given
       
        # copying .itp files supposing that those are located in the same directory as the .top file
        topdir = os.path.split(args.topfile)[0]
       
        for itpfile in glob.glob(topdir + '*.itp'):
            shutil.copy(topdir + itpfile, rundir + '/' + itpfile) # copy .itp files
       
        comm.Barrier()
        if rank == 0:
            tcpu1 = time.time()
            print "Time used for preprocessing (parallelization): %.2fs" %(tcpu1 - tcpu0)
       
        script = 'run.sh'
        self.write_script(script, rundir, rank, 'grompp.mdp', 'start.gro', 'topol.top', 'out.gro',\
                     ndxfile=args.ndxfile, grompp_options=args.grompp_options, mdrun_options=args.mdrun_options)
       
        comm.Barrier()
       
        if rank == 0:
            print "Run GROMACS preprocessing and MD..."
       
        os.chdir(rundir)
        os.system('chmod +x' + ' ' + script)
        os.system('./' + script)
        os.chdir(curdir)
       
        comm.Barrier()
       
        # post processing
        if rank == 0:
            tcpu2 = time.time() 
            print "Time used for MD: %.2fs" %(tcpu2 - tcpu1)
       
            shutil.rmtree(args.output_file, ignore_errors=True) 
       
            with open(args.output_file, 'w') as output_file:
                for rundir in rundirs:
                    with open(rundir + '/' + 'out.gro', 'r') as output_file_thread:
                        for line in output_file_thread:
                            print >> output_file, line.replace("\n", "")
            
            print "Output data have been saved in %s" %args.output_file
