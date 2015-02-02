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
