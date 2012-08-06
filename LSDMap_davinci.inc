# MPI Fortran 90 compiler
f90=ifort
mpif90=mpif90
fflags=-O3 
##-assume realloc_lhs -check all -traceback -warn all -fstack-protector -assume protect_parens

# C compiler
cc=gcc
cflags=-O3

# library path (-L<path1> -L<path2> ...)
libpath=-L/opt/apps/intel/2011.0.013/mkl/lib/intel64/ -L${HOME}/lib

# lapack and blas (-l<lib_name1> -l<lib_name2> ...)
lib_lapack_blas=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

# parpack and arpack (-l<lib_name1> -l<lib_name2> ...)
lib_parpack_arpack=-lparpack_davinci -larpack_davinci
