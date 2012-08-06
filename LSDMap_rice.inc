# mpi fortran 90 compiler
f90=xlf90
mpif90=mpif90
fflags=-O3 -q64 -qextname -qmaxmem=30000 

# C compiler
cc=xlc
cflags=-O3

# library path
libpath=-L${HOME}/lib -L/opt/apps/ibm/lib64
# lapack and blas
lib_lapack_blas=-llapack_biou -lessl
# parpack and arpack
lib_parpack_arpack=-lparpack_biou -larpack_biou
