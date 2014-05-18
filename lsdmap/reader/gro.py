import numpy as np
import itertools as it


class GroError(Exception):
    pass

class Reader(object):

    def __init__(self, filename):

        self.filename = filename
        self.file = open(filename, 'r')

    def read(self):

        config = self.next()
        while config is not None:
            yield config
            config = self.next()

    def next(self):

        try:
            self.file.next()
        except StopIteration:
            return None

        try:
            natoms = int(self.file.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")

        config = np.zeros((3, natoms), dtype='float')

        for idx_atom, line in it.izip(xrange(natoms), self.file):
            x, y, z = map(float, line.split()[3:6])
            config[0, idx_atom] = x
            config[1, idx_atom] = y
            config[2, idx_atom] = z

        try:
            box_line = self.file.next()
        except StopIteration:
            raise GroError("File ended unexpectedly when reading box line.")

        return config

    def readlines(self):
        configs = []
        config = self.next()
        while config is not None:
            configs.append(config)
            config = self.next()
        return np.array(configs)


    def close(self):
        self.file.close()

    def __iter__(self):
        return self.read()

    readline = next
