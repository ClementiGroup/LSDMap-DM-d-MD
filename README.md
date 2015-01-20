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

See the paper W. Zheng, M. A. Rohrdanz, M. Maggioni and C. Clementi, J. Chem. Phys., 2011, 134, 144109 for more information on how LSDMap works.


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
   http://mpi4py.scipy.org/, http://cython.org/. Compile/install them.

2. From the main lsdmap distribution directory run this command
   (plus any extra flags, e.g., --prefix or --user to specify the
   installation directory):

        python setup.py install

After installation, make sure that the folder bin inside your installation
directory is included in your PATH. It contains the executable "lsdmap"
that is used to compute LSDMap. Tests can be run in the folder
examples/lsdmap. This folder contains the structure file aladip_1000.gro
which contains 1000 configurations of alanine dipeptide in vacuum and
an example of configuration file (.ini) that should be used to compute
LSDMap. To test the program, simply type in this folder:

	lsdmap -f config.ini -c aladip_1000.gro

After execution, a file ".ev" and a file ".eg" must have been generated.
They contain the eigenvectors and eigenvalues of the Fokker-Planck operator,
respectively. In the file ".ev", the first column corresponds to the values
of the first eigenvector (largest eigenvalue), the second column corresponds
to the values of the second one, and so on. The first line corresponds to
the values of the eigenvectors for the first configuration given in the
structure file, the second line corresponds to the second configuration, 
and so on.

LSDMap can be computed using MPI using a command similar to:

	mpiexec -n <number_of_processors> lsdmap -f <configuration_file> -c <structure_file>

For more information on lsdmap command, simply type:

	lsdmap -h


DM-d-MD
=======

DM-d-MD (Diffusion-Map directed Molecular Dynamics) is an adaptive sampling
algorithm based on LSDMap. For an introduction to DM-d-MD, see the paper
J.Preto and C. Clementi, Phys. Chem. Chem. Phys., 2014, 16, 19181-19191.

Besides LSDMap, DM-d-MD requires GROMACS to be correctly installed.

DM-d-MD is automatically installed when installing LSDMap via the command:

        python setup.py install

A typical usage of DM-d-MD is to call:

	dmdmd -f <configuration_file>

Prerequisites
-------------

In order to use DM-d-MD, it is required that GROMACS has been correctly
installed and that "grompp" and "mdrun" commands are working properly
for a serial utilization. If not, please visit http://www.gromacs.org/.

Testing
-------

The Folder examples/dmdmd contains an example of DM-d-MD configuration
file (dmdmd.ini) as well as files required to run GROMACS MD simulations
for the photoactive yellow protein (PYP). DM-d-MD can be launched by 
executing the command:

	dmdmd -f dmdmd.ini

within the specified folder.


Diffusion-Map Sampling
======================

DMap Sampling (Diffusion Map Sampling) is our most recent adaptive sampling
algorithm. It combines LSDMap and techniques from Umbrella Sampling, metadynamics
and more recently, TRAM (Transition-based Reweighting Analysis Method) to offer
the best way of exploring configuration spaces of macromolecular systems and
estimating their free energy landscapes.

Besides LSDMap, DMap Sampling requires RADICAL-Pilot, a pilot-job system, to be
installed properly on the local workstation. Using RADICAL-Pilot implies that
the user prepares the necessary input files and DMap Sampling configuration
file on their local workstation, and launches the job from there, but the
calculations are then performed on the execution host, which is typically an
HPC resource. To install RADICAL-Pilot, see the paragraph Installation below.

DMap Sampling consists in conducting biased Molecular Dynamics simulations (MD)
with the biased potential is a local estimate of the free energy of the 
system along Diffusion Map coordinates. Since DMap coordinates are associated
with slow time scales of MD simulations, it becomes easier to explore a wider 
region of the configuration space without remaining in local minima. To run
the biased MD, an appropriate hacked version of GROMACS is needed. To install
it, see the paragraph Installation below.

Installation
------------

As mentioned above, DMap Sampling uses RADICAL-Pilot -- see the website
http://radicalpilot.readthedocs.org/en/latest/ for more information -- implying
that a typical DMap Sampling job is launched from a local workstation whereas
heavy computations which require the utilization of many CPUs (e.g., to conduct
biased MD simulation, to compute LSDMap, ...) will be performed on a remote host.
The above procedure implies:

1. Password-less ssh login to the remote host.

2. Installing the LSDMap package on both the local workstation and the remote 
machine -- see above regarding the installation of LSDMap.

3. Installing RADICAL-Pilot and make sure it is well configured to work on the
remote host. For example, RADICAL-Pilot can be installed using pip, see the website
http://radicalpilot.readthedocs.org/en/latest/installation.html for more information.

RADICAL-Pilot has been designed to work with many different HPC clusters, located
in many countries (e.g. “archer”, located at EPSRC, UK, "davinci", located at Rice
University, Houston, USA). To consult the list of pre-configured ressources, see:
http://radicalpilot.readthedocs.org/en/latest/resources.html#chapter-resources.
If your remote machine is not in the list, you will need to write a custom Resource 
Configuration File. See the instructions on the following webpage:
http://radicalpilot.readthedocs.org/en/latest/machconf.html#writing-a-custom-resource-configuration-file
To make sure RADICAL-Pilot is working well, see: http://radicalpilot.readthedocs.org/en/latest/testing.html

4. Installing the standard version of GROMACS on the local workstation. See:
   http://www.gromacs.org/ for more information.

5. Installing the modified (serial, simple precision) version of GROMACS on the
   remote machine. The package can be git-cloned using the following command:

	git clone git://git.code.sf.net/p/dmaps-gromacs/code dmaps-gromacs-code

   Note that the package includes a modified version of GROMACS 4.6.1. Once cloned,
   follow the following steps to install the package: 

   a. Create the GROMACS Makefile. The modified version of GROMACS should
      be configured similarly to the standard version of GROMACS 4.6, with
      single precision and GPU and MPI options disabled. Please refer to
      http://www.gromacs.org/Documentation/Installation_Instructions_4.6
      for more information on how to build GROMACS. A typical way of
      configuring our modified version of GROMACS could be the following
      once the installation directory has been git cloned:

	cd dmaps-gromacs-code
	mkdir cmake-build
	cd cmake-build
	cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/path/to/local/directory -DGMX_GPU=OFF -DGMX_MPI=OFF

      Beware that the above commands is a simple example. Some machines may
      require additional options to have GROMACS work properly. Again, refer
      to the GROMACS webpage above for more information.
      
   b. Before installing GROMACS, link the DMap Sampling library to it. Before
      using make or make install, you should tell GROMACS that it should
      consider a library called libdms.so at the compilation. The library
      was created when the LSDMap package was installed. First, you
      should find the location of this particular library. If the installation
      of LSDMap went properly, the corresponding path should be 
      /path/to/site-packages/lsdmap/dmaps/critical/libdms.so, where 
      path/to/site-packages is the path to the python "site-packages" folder where
      LSDMap has been installed. Then, starting from the folder 
      dmaps-gromacs-code that was git cloned earlier, you should edit the file
      src/kernel/CMakeLists.txt. At line 78 in this file, "path/to/libdms.so" in 
      "set(LIBDMS path/to/libdms.so)" should be modified to correspond to the
      absolute path to the library "libdms.so" that you found.

   c. Compile GROMACS by typing (from the folder cmake-build that was created in
      step a):

	make install

You have installed all the prerequisites needed to use Diffusion Map Sampling!

Testing
-------

The folder examples/dmaps contains a set of files that can be used to run
Diffusion Map Sampling on the Rice University cluster BlueBiou. If you want
to run DMap Sampling on a different remote host, please change the parameter
"remote_host" in the file settings of the example folder. Note that the parameters
uname and queue in this file correspond to the username and the name of the
queue on the remote host, respectively, and should be updated by the user.

To run Diffusion Map Sampling, use the following command from the examples/dmaps
directory:

	dmaps -f settings

As the job progresses, logging messages will be written to the console.
If Diffusion Map Sampling finished correclty it should have created 3 folders,
corresponding to the first three iterations of the DMap Sampling procedure (in
the file settings, niters=3 by default), called iterX, where X is the number
of the iteration.

To analyze the results of a particular iteration, go inside the corresponding
"iter" folder. From there, you can use the script phipsi.sh which is located in
the parent directory, by typing:

	../phipsi.sh

This will compute collective variables that are commonly used to analyze results
of the system tested here, alanine dipeptide. A bunch of python scripts are available
in the folder examples/dmaps/scripts to plot specific figures. For example, the
file examples/dmaps/scripts/plotfephipsi.py can be used to plot the free energy
landscape as a a function of the coordinates phi and psi computed after running
the script phipsi.sh. From the same "iter" folder, using:

	python ../scripts/plotfephipsi.py

will plot the desired figure. Note that it requires the modules matplotlib and
densplot to be installed on the local host. The former can be downloaded from
http://matplotlib.org/, whereas the latter can be git cloned 
from https://github.com/jp43/densplot and installed the usual way using 

	python setup.py install
