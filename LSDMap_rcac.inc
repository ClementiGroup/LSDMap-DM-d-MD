# MPI Fortran 90 compiler
f90=ifort
mpif90=mpif90
fflags=-O -i8

# C compiler
cc=gcc
cflags=-O3

# library path (-L<path1> -L<path2> ...)
libpath=-L/opt/intel/mkl/10.2.5.035/lib/em64t -L${HOME}/lib

# lapack and blas (-l<lib_name1> -l<lib_name2> ...)
lib_lapack_blas=-lmkl_blas95_ilp64 -lmkl_lapack95_ilp64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -limf

# parpack and arpack (-l<lib_name1> -l<lib_name2> ...)
lib_parpack_arpack=-lparpack_steele -larpack_steele
