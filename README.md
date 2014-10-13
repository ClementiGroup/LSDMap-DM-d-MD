LSDMap
======

LSDMap package is used to compute Locally Scaled Diffusion Map. 
Typical usage is to call "lsdmap" script:

    lsdmap -f <configuration_file> -c <structure_file> <other_options>

or using MPI:

    mpiexec -n <number_of_processors> lsdmap -f <configuration_file> -c <structure_file> <other_options>

A typical example of configuration file is ./examples/lsdmap/config.ini
The structure file should contain all the configurations needed to compute
LSDMap. After execution of the script, a .ev and a .eg file should have
been generated containing the eigenvectors and eigenvalues, respectively. 


Prerequisites
-------------

Before installing LSDMap package, make sure that you have the following
packages installed:

* NumPy; version 1.4.1 or larger

* Scipy; version 0.10.0 or larger

* mpi4py; version 1.0 or larger

* cython; version 0.21 or later

Version 2.6.x or 2.7.x of python should be used. 


Installation
------------

The Python Distutils are used to build and install LSDMap, so it is
fairly simple to get things ready to go. Following are very simple
instructions on how to proceed:

1. First, make sure that you have NumPy, SciPy, mpi4py and cython
   installed. If not, get them from http://numpy.scipy.org/,
   http://mpi4py.scipy.org/. Compile/install them.

2. From the main lsdmap distribution directory run this command,
   (plus any extra flags like --prefix or --user to specify the
   installation directory)::

        python setup.py install

After installation, make sure that the folder bin inside your 
installation directory is included in your PATH (normally it should). 
It contains the executable "lsdmap" that is used to compute LSDMap. 
Tests can be run in the folder ./examples/lsdmap. This folder contains
the structure file aladip_1000.gro which contains 1000 configurations of
alanine dipeptide in vacuum and an example of configuration file (.ini)
that should be used to compute LSDMap. To test the program, simply type
in this folder:

    lsdmap -f config.ini -c aladip_1000.gro

After execution, a '.ev' and a '.eg' file must have been generated. They contain
the eigenvectors and eigenvalues of the Fokker-Planck operator, respectively. In 
'.ev' file, the first column corresponds to the values of the first eigenvector
(largest eigenvalue), the second column corresponds to the values of the second
one, and so on. The first line corresponds to the values of the eigenvectors for
the first configuration given in the structure file, the second line corresponds
to the second configuration, and so on.

LSDMap can be computed using MPI using a command similar to:

    mpiexec -n <number_of_processors> lsdmap -f <configuration_file> -c <structure_file>

For more information on lsdmap command, simply type:

    lsdmap -h


DM-d-MD
=======

Executing

        python setup.py install

will also install the DM-d-MD (Diffusion-Map-directed molecular dynamics) package.
A typical usage of DM-d-MD is to call:

    dmdmd -f <configuration_file>


Prerequisites
-------------

In order to use DM-d-MD, it is required that GROMACS has been correctly
installed and that "grompp" and "mdrun" commands are working properly.


Testing
-------

Folder ./examples/dmdmd contains an example of DM-d-MD configuration
file (dmdmd.ini) as well as files required to run GROMACS MD simulations
for the photoactive yellow protein (PYP). DM-d-MD can be launched by 
executing the command:

    dmdmd -f dmdmd.ini

within the specified folder.


DMap Sampling
=============

Installation steps to use Diffusion Map Sampling:

1. Install LSDMap package on the remote machine (see section LSDMap above)

2. Install radical pilot

3. Install the modified (serial, simple precision) version of GROMACS needed 
   to run diffusion map sampling on the remote machine. The package can be
   git-cloned using the following command:

	git clone git://git.code.sf.net/p/dmaps-gromacs/code dmaps-gromacs-code

   Once cloned, follow the following steps to install GROMACS: 

   a. build GROMACS Makefile

	cd dmaps-gromacs-code
	mkdir cmake-build
	cd cmake-build
	cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/path/to/local/directory -DGMX_GPU=OFF -DGMX_MPI=OFF

   b. bind libdms.so to GROMACS
    First, find the location of the library libdms.so that was generated when
    installing lsdmap, it should be /path/to/site-packages/dmaps/kernel/libdms.so, 
    where path/to/site-packages is the path to the python folder site-packages
    where lsdmap has been installed. If you specified the local directory with
    the prefix option when installing lsdmap, i.e., --prefix=/path/to/local/directory,
    then the location of libdms.so is likely to be something like
    /path/to/local/directory/lib/pythonX.Y/site-packages/dmaps/kernel/libdms.so
    where X and Y in pythonX.Y corresponds to the version of python that is used.

  Then do the following still from the folder "cmake-build" that you created
in 2) a):

        vi ../src/kernel/CMakeLists.txt

  At line 77-78 of this file, you can find something like:

  # DMap sampling
  set(LIBDMS path/to/libdms.so)

  change path/to/libdms.so to the path to libdms.so you found above, save and
quit.

  Then compile GROMACS by typing (from cmake-build folder):

        make install

  You are ready to use LSDMap driven adaptive sampling!
