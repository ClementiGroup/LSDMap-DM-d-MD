"""Reader for .xvg files."""

import numpy as np
import itertools as it


class XvgError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):

        if isinstance(filename, basestring):
            filename = [filename]
        else:
            filename = filename
        self.filename = filename
        self.nfiles = len(self.filename)
        self.file = [open(name, 'r') for name in self.filename]

        if 'col' in kwargs:
            self._col = kwargs['col']
        else:
            self._col = None

    def read(self):

        coords = self.next()
        while coords is not None:
            yield coords
            coords = self.next()

    def next(self):
        try:
            coords = []
            for file in self.file:
                line = file.next().split()
                if self._col is not None:
                    coords.extend(map(float, line)[self._col])
                else: #if self._col is not specified read column 1
                    coords.append(map(float, line)[1])
            return np.array(coords)
        except StopIteration:
            return None


    def readlines(self):
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
