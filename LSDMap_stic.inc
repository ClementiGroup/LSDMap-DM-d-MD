# mpi fortran 90 compiler
f90=pgf90
mpif90=mpif90
fflags=-O3 -byteswapio

# C compiler
cc=gcc
cflags=-O3

# library path
libpath=-L${HOME}/lib

# lapack and blas
lib_lapack_blas=-llapack_stic -lblas_stic

# parpack and arpack
lib_parpack_arpack=-lparpack_stic -larpack_stic
