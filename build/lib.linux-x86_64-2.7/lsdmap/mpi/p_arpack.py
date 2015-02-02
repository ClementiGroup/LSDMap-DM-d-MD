import os, sys
import numpy as np

from scipy.sparse.linalg.eigen.arpack import _arpack
from scipy.sparse.linalg.eigen.arpack.arpack import _ArpackParams, _type_conv, eigsh, ArpackError, _SAUPD_ERRORS

from mpi4py import MPI

class _ParallelSymmetricArpackParams(_ArpackParams):

    def __init__(self, comm, matrix_thread, k, **kargs):

        self.comm = comm
        self.matrix_thread = matrix_thread
        self.nrows_thread = np.shape(self.matrix_thread)[0]
        n = comm.allreduce(self.nrows_thread, op=MPI.SUM)
        tp = matrix_thread.dtype.char

        _ArpackParams.__init__(self, n, k, tp, **kargs) # arpack.py l.311

        self.bmat = 'I'
        self.B = lambda x: x

        self.workd = np.zeros(3 * n, self.tp)
        self.workl = np.zeros(self.ncv * (self.ncv + 8), self.tp)

        ltr = _type_conv[self.tp]
        if ltr not in ["s", "d"]:
            raise ValueError("Input matrix is not real-valued.")

        self._arpack_solver = _arpack.__dict__[ltr + 'saupd']
        self._arpack_extract = _arpack.__dict__[ltr + 'seupd']

        self.ipntr = np.zeros(11, "int")
        self.iterate_infodict = _SAUPD_ERRORS[ltr]

    def iterate(self):

        self.ido, self.tol, self.resid, self.v, self.iparam, self.ipntr, self.info = \
            self._arpack_solver(self.ido, self.bmat, self.which, self.k,
                                self.tol, self.resid, self.v, self.iparam,
                                self.ipntr, self.workd, self.workl, self.info)

        xslice = slice(self.ipntr[0] - 1, self.ipntr[0] - 1 + self.n)
        yslice = slice(self.ipntr[1] - 1, self.ipntr[1] - 1 + self.n)

        if self.ido in [-1, 1]:
            _workdy_thread = self.matrix_thread.dot(self.workd[xslice])
            self.workd[yslice] = np.hstack(self.comm.allgather(_workdy_thread))
        elif self.ido == 2:
            self.workd[yslice] = self.B(self.workd[xslice])
        elif self.ido == 3:
            raise ValueError("ARPACK requested user shifts.  Assure ISHIFT==0")
        else:
            self.converged = True

            if self.info == 0:
                pass
            elif self.info == 1:
                self._raise_no_convergence()
            else:
                raise ArpackError(self.info, infodict=self.iterate_infodict)


    def extract(self, return_eigenvectors=True):
        rvec = return_eigenvectors
        ierr = 0
        howmny = 'A'  # return all eigenvectors
        sselect = np.zeros(self.ncv, 'int')  # unused
        d, z, ierr = self._arpack_extract(rvec, howmny, sselect, self.sigma,
                                          self.bmat, self.which, self.k,
                                          self.tol, self.resid, self.v,
                                          self.iparam[0:7], self.ipntr,
                                          self.workd[0:2 * self.n],
                                          self.workl, ierr)
        if ierr != 0:
            raise ArpackError(ierr, infodict=self.extract_infodict)
        k_ok = self.iparam[4]
        d = d[:k_ok]
        z = z[:, :k_ok]

        if return_eigenvectors:
            return d, z
        else:
            return d


def run_test():

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if size!=3:
        if rank==0: print "###Error: the number of processors should be 3. Please use: mpiexec -n 3 python p_arpack.py"
        sys.exit(1)

    if rank == 0: a = np.array([[1,2,0]], dtype='float');
    if rank == 1: a = np.array([[2,0,3]], dtype='float');
    if rank == 2: a = np.array([[0,3,0]], dtype='float');

    slice_thread=slice(rank,rank+1)
    params= _ParallelSymmetricArpackParams(comm, a, 1)

    while not params.converged:
        params.iterate()

    vals, vecs = params.extract(return_eigenvectors=True)
    if rank == 0:
        print 'parallel arpack'
        print vecs, '\n'

    A = np.array([[1,2,0],[2,0,3],[0,3,0]], dtype='float')
    vals, vecs = eigsh(A, k=1)

    if rank==0:
        print 'serial arpack'
        print vecs

if __name__=="__main__":
    run_test()

