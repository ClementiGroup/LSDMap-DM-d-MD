import os
import sys

#platforms
platforms = {'mac': ['darwin', 'os2', 'os2emx', 'riscos', 'atheos'], 'windows': ['win32', 'cygwin'], 'linux': ['linux2']}

# tmpfiles dictionary, each elements is a tuple with the filename and a brief description
tmpfiles = {'ingro': ('input.gro', 'md input .gro file'), 'outgro' : ('output.gro', 'md output .gro file')}

# known pre exec for radical pilot
stampede_pre_exec = ["module load -intel intel/14.0.1.106", "module load python",
"PYTHONPATH=$PYTHONPATH:/opt/apps/intel14/mvapich2_2_0/python/2.7.6/lib/python2.7",
"PYTHONPATH=$PYTHONPATH:/opt/apps/intel14/mvapich2_2_0/python/2.7.6/lib/python2.7/site-packages",
"PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages","export PYTHONPATH"]

davinci_pre_exec = [""]
biou_pre_exec = [""]

known_pre_exec = {"stampede.tacc.utexas.edu": stampede_pre_exec, "rice.davinci": davinci_pre_exec, "rice.biou": biou_pre_exec}
