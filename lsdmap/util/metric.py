import sys
import numpy as np
import cython
import copy
import pyqcprot
import util

global MAXSIZE
MAXSIZE = 5E8

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
    def __init__(self, metric, ndim=3, **kwargs):

        self.ndim = ndim
        self.prms = kwargs

        if type(metric)==str:
            self.metric_name = copy.copy(metric)
            self.function = self.get_function()

        elif callable(metric):
            self.metric_name = str(metric)
            self.function = metric

    def _h_rmsd(self, coord1, coord2):
        return pyqcprot.CalcRMSDRotationalMatrix(coord1, coord2, None, None)

    def _h_cmd(self, coord1, coord2):
        return util.cmd(coord1, coord2, self.r0)

    def _h_dihedral(self, coord1, coord2):
        return util.dihedral(coord1, coord2)

    def _h_euclidean(self, coord1, coord2):
        return util.euclidean(coord1, coord2)

    def get_function(self):

        metric_name = self.metric_name.lower()
        _mapped = {}

        if metric_name in _mapped:
            metric_name = _mapped[metric_name]

        metric1D_names = ['dihedral']
        metric3D_names = ['rmsd', 'cmd']

        if metric_name in metric1D_names:
            if self.ndim != 1:
                raise ValueError('%s is a 1D metric. Please check the spatial dimensions of your coordinates'%metric_name)

        if metric_name in metric3D_names:
            if self.ndim != 3:
                raise ValueError('%s is a 3D metric. Please check the spatial dimensions of your coordinates'%metric_name)

        if metric_name == 'cmd':
            if 'r0' in self.prms:
                self.r0 = self.prms['r0']
            else:
                self.r0 = 0.75

        metric_func_name = "_h_" + metric_name

        if hasattr(self, metric_func_name):
            self._metric_func = getattr(self, metric_func_name)
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
        where coords_configY_setX is a 1D or 2D numpy array containing the coordinates of 
        configuration Y of set X. 
              
    metric: string or callable, optional
        The metric used to compute the distance matrix. If metric is a callable, it is used
        directly to compute the distance. If metric is a string, it should be one of the 
        metric names provided by the Metric class above (rmsd, euclidean,...). The default is 
        'rmsd'.
    
        Note: If metric is a 3D-metric name, the coords_configY_setX mentioned above are supposed 
        to be given in format [[x1, x2,...,xm], [y1, y2, ..., ym], [z1, z2, ..., zm]].


    Examples
    --------

    >>> import numpy as np # import numpy library
    >>> coords1 = np.array([[1.0, 6.0], [3.5, 1.5], [1.7, 1.2]], [[4.8, 1.5], [5.2, 1.7], [3.7, 2.1]]]) # length: n1
    >>> coords2 = np.array([[1.0, 6.0], [3.5, 1.5], [1.7, 1.2]], [[4.8, 1.5], [5.2, 1.7], [3.7, 2.1]], [[1.0, 6.0], [3.5, 1.5], [1.7, 1.2]]])# length: n2
    >>> DistanceMatrix = DistanceMatrix(coords1, coords2, metric='rmsd') # create instance of DistanceMatrix
    >>> distance_matrix = DistanceMatrix.distance_matrix  # returns the distance matrix (dimensions: n1 * n2)
    >>> neighbor_distance_matrix = DistanceMatrix.neighbor_matrix(k=10) # (dimensions: n1 * k)
    >>> idx_neighbor_matrix = DistanceMatrix.idx_neighbor_matrix(k=10) # (dimensions: n1 * k)
    """

    def __init__(self, coords1, coords2, metric='rmsd', metric_prms={}):
       
       self.coords1 = np.array(coords1)
       self.coords2 = np.array(coords2)

       shape_coords1 = self.coords1.shape
       shape_coords2 = self.coords2.shape

       if len(shape_coords1) != len(shape_coords2):
           raise TypeError('coords1 and coords2 should have the same number of dimensions!')


       if len(shape_coords1) == 3:
           if shape_coords1[1] == shape_coords2[1]:
               self.ndim = shape_coords1[1]
               self.natoms = shape_coords1[2]
               if self.ndim > 3:
                   raise TypeError('the number of rows of each coordinate should be less than 4 (number of spatial dimensions)')
               #if self.ndim == 1:
               #    self.coords1 = np.squeeze(self.coords1, axis=(1,))
               #    self.coords2 = np.squeeze(self.coords2, axis=(1,))
           else:
               raise TypeError('coords1 and coords2 have not the same number of spatial dimensions')
       elif len(shape_coords1) == 2:
           self.ndim = 1
           self.natoms = shape_coords1[1]
       else:
           raise TypeError('coords1 and coords2 should have a number of dimensions 1 < ndim < 4;                                                                                                      if only one coordinate is used, consider using coords1[np.newaxis] or coords2[np.newaxis]')
           
       self.metric = Metric(metric, ndim=self.ndim, **metric_prms).function
       self.ncoords1 = self.coords1.shape[0]
       self.ncoords2 = self.coords2.shape[0]
       self.maxsize = MAXSIZE

    def __getattr__(self, name):
        if name == "distance_matrix":
            if hasattr(self, '_distance_matrix'):
                return self._distance_matrix
            else:
                return self.get_distance_matrix()
        return DistanceMatrix.__getattribute__(self, name)

    def get_distance_matrix(self):
        if (self.ncoords1*self.ncoords2) > self.maxsize:
            raise ValueError("Large distance matrix expected! use more threads to avoid too much memory")

        matrix = np.zeros((self.ncoords1, self.ncoords2))
        for idx, coord1 in enumerate(self.coords1):
            for jdx, coord2 in enumerate(self.coords2):
                matrix[idx, jdx] = self.metric(coord1, coord2)

        #matrix = np.array([[self.metric(coord1, coord2) for coord2 in self.coords2]
        #    for coord1 in self.coords1])

        matrix[np.isnan(matrix)] = 0.0
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
        else:
            k = self.ncoords2

        if (self.ncoords1*k) > self.maxsize:
            raise ValueError("Large distance matrix expected! use more threads to avoid too much memory")

        neighbor_matrix = np.zeros((self.ncoords1, k), dtype='float')
        idx_neighbor_matrix = np.zeros((self.ncoords1, k), dtype='int')

        if hasattr(self, '_distance_matrix'):
            for idx, distance in enumerate(self._distance_matrix):
                idx_neighbors = np.argsort(distance)[:k]
                idx_neighbor_matrix[idx] = idx_neighbors
                neighbor_matrix[idx] = [distance[idx_neighbor] for idx_neighbor in idx_neighbors]
        else:
            for idx, coord1 in enumerate(self.coords1):
                distance = np.array([self.metric(coord1, coord2) for coord2 in self.coords2])    
                idx_neighbors = np.argsort(distance)[:k]
                idx_neighbor_matrix[idx] = idx_neighbors
                neighbor_matrix[idx] = [distance[idx_neighbor] for idx_neighbor in idx_neighbors]

        return neighbor_matrix, idx_neighbor_matrix


def get_neighbor_matrix(distance_matrix,k=None):

    ncoords1 = distance_matrix.shape[0]
    ncoords2 = distance_matrix.shape[1]

    if k is not None:
        if k >= ncoords2:
            print "Warning: k > = number of data points "
    else:
        k = ncoords2

    if (ncoords1*k) > MAXSIZE:
        raise ValueError("Large distance matrix expected! use more threads to avoid too much memory")

    neighbor_matrix = np.zeros((ncoords1, k), dtype='float')
    idx_neighbor_matrix = np.zeros((ncoords1, k), dtype='int')

    for idx, distance in enumerate(distance_matrix):
        idx_neighbors = np.argsort(distance)[:k]
        idx_neighbor_matrix[idx] = idx_neighbors
        neighbor_matrix[idx] = [distance[idx_neighbor] for idx_neighbor in idx_neighbors]

    return neighbor_matrix, idx_neighbor_matrix
