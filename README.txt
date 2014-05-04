======
LSDMap
======

LSDMap package is used to compute Locally Scaled Diffusion Map. 
Typical usage is to call "lsdmap" script:

    lsdmap -f <configuration_file> -c <structure_file> <other_options>

or using MPI:

    mpiexec -n <number_of_processors> lsdmap -f <configuration_file> -c <structure_file> <other_options>


A typical example of configuration file can be found in ./examples. The
structure file should contain all the configurations needed to compute
LSDMap. 


Prerequisites
=============

Before installing LSDMap package, make sure that you have the following
packages installed:

* NumPy; version 1.4.1 or larger

* Scipy; version 0.10.0 or larger

* Cython; version 0.15 or larger

* mpi4py; version 1.0 or larger

Version 2.6.x or 2.7.x of python should be used. 


Installation
============

The Python Distutils are used to build and install LSDMap, so it is
fairly simple to get things ready to go. Following are very simple
instructions on how to proceed:

1. First, make sure that you have NumPy, SciPy, Cython and mpi4py 
   installed. If not, get them from http://numpy.scipy.org/,
   http://cython.org/, http://mpi4py.scipy.org/. Compile/install them.

2. From the main lsdmap distribution directory run this command,
   (plus any extra flags like --prefix to specify the installation
   directory)::

        python setup.py install

After installation, make sure that the folder bin inside your 
installation directory is included in your PATH. It contains the 
executable "lsdmap" that is used to compute LSDMap. Tests can be run
in the folder ./examples. This folder contains the structure file 
aladip_1000.gro which contains 1000 configurations of alanine dipeptide 
in vacuum and an example of configuration file that should be used to
compute LSDMap. To test the program, simply type in this folder:

    lsdmap -f config.ini -c aladip_1000.gro

After execution, a '.ev' and a '.eg' file must have been generated. They contain
the eigenvectors and eigenvalues of the Fokker-Planck matrix, respectively. In 
the '.ev' file, the first column corresponds to the values of the first eigenvector
(largest eigenvalue), the second column corresponds to the values of the second one,
and so on. The first line corresponds to the values of the eigenvectors for the
first configuration given in the structure file, the second line corresponds to
the second configuration, and so on.

LSDMap can be computed using MPI using a command similar to:

    mpiexec -n <number_of_processors> lsdmap -f <configuration_file> -c <structure_file>

For more information on lsdmap command, simply type:

    lsdmap -h

