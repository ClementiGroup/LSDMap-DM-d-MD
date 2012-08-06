#Settings for use with Cray XE6 HECToR
#using the Cray Fortran compiler.
#Optimisation is turned on by default

# mpi fortran 90 compiler
f90=ftn
mpif90=ftn
fflags=-em 

# C compiler
cc=cc
cflags=

# library path
libpath=-L/path/to/ARPACK 
# lapack and blas
lib_lapack_blas=
# parpack and arpack
lib_parpack_arpack=-lparpack -larpack 
