import numpy as np

def get_idxs_thread(comm, npoints):
    """
    get the indices of the coordinates that will be consider by each CPU
    """
    size = comm.Get_size()
    rank = comm.Get_rank()

    npoints_per_thread = np.array([npoints/size for idx in xrange(size)])

    for idx in range(npoints%size):
        npoints_per_thread[idx] += 1

    if rank == 0:
        idx_shift = 0
        idxs_threads = []
        for idx in npoints_per_thread:
            idxs_threads.append(range(idx_shift,idx_shift+idx))
            idx_shift += idx
    else:
        idxs_threads = None

    idxs_thread = comm.scatter(idxs_threads, root=0)

    return idxs_thread

def get_idxs_thread_v(comm, npoints):
    """ Get indices for processor using Scatterv

    Note: 
    -----
        Uppercase mpi4py functions require everything to be in C-compatible
    types or they will return garbage!
    """

    size = comm.Get_size()
    rank = comm.Get_rank()

    npoints_per_thread = np.zeros(size,dtype=np.int32)
    offsets_per_thread = np.zeros(size,dtype=np.int32)

    for idx in range(size):
        npoints_per_thread[idx] = npoints/size
        offsets_per_thread[idx] = sum(npoints_per_thread[:idx])

    for idx in range(npoints % size):
        npoints_per_thread[idx] += 1
        offsets_per_thread[idx + 1:] += 1

    npoints_per_thread = tuple(npoints_per_thread)
    offsets_per_thread = tuple(offsets_per_thread)

    idxs_thread = np.zeros(npoints_per_thread[rank],dtype=np.int32)
    idxs = np.arange(npoints,dtype=np.int32)

    comm.Scatterv([idxs, npoints_per_thread, offsets_per_thread, MPI.INT], idxs_thread, root=0)
    return idxs_thread
