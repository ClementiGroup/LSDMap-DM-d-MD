import numpy as np
import itertools as it


class XyzError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):

        self.filename = filename
        self.file = open(filename, 'r')

    @property
    def natoms(self):
        f = open(self.filename, 'r')
        natoms = int(f.next())
        f.close()
        return natoms

    @property
    def nlines(self):
        f = open(self.filename, 'r')
        for idx, line in enumerate(f):
            pass
        nlines = (idx+1)/(self.natoms+2)
        f.close()
        return nlines

    def read(self):

        config = self.next()
        while config is not None:
            yield config
            config = self.next()

    def next(self):

        try:
            natoms = int(self.file.next())
        except StopIteration:
            return None

        try:
            self.file.next()
        except StopIteration:
            raise XyzError("File ended unexpectedly when reading 2nd line.")

        firstline = self.file.next().split()
        ndim = len(firstline)
        config = np.zeros((ndim, natoms), dtype='float')

        for idx in xrange(ndim):
            config[idx, 0] = firstline[idx]

        for idx_atom, line in it.izip(xrange(natoms-1), self.file):
            firstline = self.file.next().split()
            for idx in xrange(ndim):
                config[idx, 0] = firstline[idx]

        return config

    def readlines(self, *args):
        coords = []
        coord = self.next()
        while coord is not None:
            coords.append(coord)
            coord = self.next()
        return np.array(coords)

    def close(self):
        self.file.close()

    def __iter__(self):
        return self.read()

    readline = next

