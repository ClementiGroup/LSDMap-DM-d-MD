import numpy as np
from mpi4py import MPI


def get_idxs_thread(comm, npoints):
    """ Get indices for processor using Scatterv

    Note: 
    -----
        Uppercase mpi4py functions require everything to be in C-compatible
    types or they will return garbage!
    """

    size = comm.Get_size()
    rank = comm.Get_rank()

    npoints_thread = np.zeros(size,dtype=np.intc)
    offsets_thread = np.zeros(size,dtype=np.intc)

    for idx in range(size):
        npoints_thread[idx] = npoints/size
        offsets_thread[idx] = sum(npoints_thread[:idx])

    for idx in range(npoints % size):
        npoints_thread[idx] += 1
        offsets_thread[idx + 1:] += 1

    npoints_thread = tuple(npoints_thread)
    offsets_thread = tuple(offsets_thread)

    idxs_thread = np.zeros(npoints_thread[rank],dtype=np.intc)
    idxs = np.arange(npoints,dtype=np.intc)

    comm.Scatterv((idxs, npoints_thread, offsets_thread, MPI.INT), idxs_thread, root=0)
    return idxs_thread, npoints_thread, offsets_thread

def get_ravel_offsets(npoints_thread,natoms):
    """ Get lengths and offsets for gathering trajectory fragments """
    size = len(npoints_thread)
    ravel_lengths = np.zeros(size,dtype=np.intc)
    ravel_offsets = np.zeros(size,dtype=np.intc)

    for i in range(size):
        ravel_lengths[i] = npoints_thread[i]*3*natoms
        ravel_offsets[i] = sum(ravel_lengths[:i])

    ravel_lengths = tuple(ravel_lengths)
    ravel_offsets = tuple(ravel_offsets)

    return ravel_lengths, ravel_offsets
