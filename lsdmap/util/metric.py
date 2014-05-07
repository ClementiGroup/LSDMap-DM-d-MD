import os
import cqcprot

import numpy as np

class Metric(object):
    """
    Metric(metric)

    A class providing the function of two variables used to compute the distance between two points.
    
    Parameters
    ----------

    metric: string or callable
        if metric is a string, metric is the name of the metric and should be associated with one
        of the _h_ functions below. For example, if metric='rmsd', the function used to compute 
        the distance will be given by the function _h_rmsd. If metric is a callable, it should be
        used directly to compute the distance.

    """
    def __init__(self, metric):

        if type(metric)==str:
            self.metric_name = metric
            self.function = self.get_function_from_name()

        elif callable(metric):
            self.metric_name = str(metric)
            self.function = metric

    def _h_rmsd(self, coord1, coord2):
       return cqcprot.rmsd(coord1, coord2)


    def _h_dihedral(self, coord1, coord2):
       raise NotImplementedError


    def get_function_from_name(self):

        metric_name = self.metric_name.lower()
        _mapped = {'dihedral':'phi_psi'}

        if metric_name in _mapped:
            metric_name = _mapped[metric_name]

        metric_name = "_h_" + metric_name

        if hasattr(self, metric_name):
                self._metric_func = getattr(self, metric_name)
        else:
                metriclist = [x[3:] for x in dir(self) if x.startswith('_h_')]
                raise ValueError("function must be one of " + ", ".join(metriclist))

        return self._metric_func


class DistanceMatrix(object):
    """
    DistanceMatrix(coords1, coords2, metric)

    A class to compute the distance matrix and neighbor distance matrix


    Parameters
    ----------

    coords1, coords2: lists or arrays
        coords1, coords2 are the two sets of coordinates used to compute the distance matrix
        coords1, coords2 should be lists or arrays given in the format: 
            coordsX = [coords_config1_setX, coords_config2_setX, ..., coords_congfignX_setX]
        where coords_configY_setX is a 1D numpy array containing the coordinates of 
        configuration Y of set X. 
              
    metric: string or callable, optional
        The metric used to compute the distance matrix. If metric is a callable, it is used
        directly to compute the distance. If metric is a string, it should be one of the 
        metric names provided by the Metric class above (rmsd, euclidean,...). The default is 
        'rmsd'.
    
        Note: If metric='rmsd', the coords_configY_setX mentioned above are supposed 
        to be given in format [x1, x2,...,xm, y1, y2, ..., ym, z1, z2, ..., zm].


    Examples
    --------

    >>> import numpy as np # import numpy library
    >>> coords1 = np.array([[1.0, 3.5, 1.7], [4.8, 5.2, 3.7]]) # length: n1
    >>> coords2 = np.array([[5.5, 2.4, 1.1], [1.8, 3.2, 0.9]]) # length: n2
    >>> DistanceMatrix = DistanceMatrix(coords1, coords2, metric='rmsd') # create instance of DistanceMatrix
    >>> distance_matrix = DistanceMatrix.distance_matrix()  # returns the distance matrix (dimensions: n1 * n2)
    >>> neighbor_distance_matrix = DistanceMatrix.neighbor_matrix(k=10) # returns the matrix of distances 
                         from each point of coords1 and its k nearest neighbors among coords2 (dimensions: n1 * k)
    >>> idx_neighbor_matrix = DistanceMatrix.idx_neighbor_matrix(k=10) # same as the neighbor distance matrix 
                         but contains the index of the neighbors inside coords2 array (dimensions: n1 * k)
    """

    def __init__(self, coords1, coords2, metric='rmsd'):
       
       self.coords1 = self.check_coords(coords1)
       self.ncoords1 = len(self.coords1)

       self.nxyz_coords= len(self.coords1[0]) 
       self.metric = Metric(metric).function

       self.coords2 = self.check_coords(coords2)
       self.ncoords2 = len(self.coords2)

       self.maxsize = 3E8;


    def check_coords(self, coords):

        if isinstance(coords, np.ndarray) or isinstance(coords, list): 
            for coord in coords:
                if not isinstance(coord, np.ndarray): raise TypeError("Not all coordinates are numpy arrays")
        else:
            raise TypeError("Coordinates should be given as a list or numpy array")

        return coords


    def __getattr__(self, name):
        if name == "distance_matrix":
            if hasattr(self, '_distance_matrix'):
                return self._distance_matrix
            else:
                return self.get_distance_matrix()

        return DistanceMatrix.__getattribute__(self, name)


    def get_distance_matrix(self):

        if (self.ncoords1*self.ncoords2) > self.maxsize: raise ValueError("Large distance matrix expected! use more threads to avoid too much memory")

        matrix = np.zeros((self.ncoords1, self.ncoords2), dtype='float')
        paircoords = ((coord1, coord2) for coord1 in self.coords1 for coord2 in self.coords2)

        for kdx, (coord1, coord2) in enumerate(paircoords):
            idx = kdx/self.ncoords2; jdx = kdx%self.ncoords2;
            matrix[idx][jdx] = self.metric(coord1, coord2) # compute distance matrix

        self._distance_matrix = matrix

        return matrix


    def neighbor_matrix(self, **kargs):

        neighbor_matrix, idx_neighbor_matrix = self.get_neighbor_matrix(**kargs)
        return neighbor_matrix


    def idx_neighbor_matrix(self, **kargs):

        neighbor_matrix, idx_neighbor_matrix = self.get_neighbor_matrix(**kargs)
        return idx_neighbor_matrix


    def get_neighbor_matrix(self, k=None):

        if k is not None:
            if k >= self.ncoords2:
                print "Warning: k > = number of data points "
        else: k = self.ncoords2

        if (self.ncoords1*k) > self.maxsize: raise ValueError("Large distance matrix expected! use more threads to avoid too much memory")

        neighbor_matrix = np.zeros((self.ncoords1, k), dtype='float')
        idx_neighbor_matrix = np.zeros((self.ncoords1, k), dtype='int')

        if hasattr(self, '_distance_matrix'):
            for idx, distance in enumerate(self.matrix):
                idx_neighbors = np.argsort(distance)[:k]  # the first element is the point itself
                idx_neighbor_matrix[idx] = idx_neighbors
                neighbor_matrix[idx] = [distance[idx_neighbor] for idx_neighbor in idx_neighbors]
        else:
            for idx, coord1 in enumerate(self.coords1):
                distance = np.array([self.metric(coord1, coord2) for coord2 in self.coords2])    
                idx_neighbors = np.argsort(distance)[:k]  # the first element is the point itself
                idx_neighbor_matrix[idx] = idx_neighbors
                neighbor_matrix[idx] = [distance[idx_neighbor] for idx_neighbor in idx_neighbors]

        return neighbor_matrix, idx_neighbor_matrix
