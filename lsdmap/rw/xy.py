"""Reader/Writer for xy files (two columns)."""

import numpy as np
import itertools as it

class XyError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):

        _known_types = {'float': float, 'int': int}

        self.filename = filename       
        self.file = open(self.filename, 'r')
        try:
            type = kwargs['type']
            if type not in _known_types:
                raise XyError('type specified %s unknown'%type)
            self.str2num = _known_types[type]
        except KeyError:
            self.str2num = float

    @property
    def nlines(self):
        f = open(self.filename, 'r')
        for idx, line in enumerate(f):
            pass
        nlines = idx+1
        f.close()
        return nlines


    def read(self):
        value = self.next()
        while value is not None:
            yield value
            value = self.next()

    def next(self):
        try:
            config = np.zeros((2, 1), dtype='float')
            x, y = map(float, self.file.next().split())
            config[0, 0] = x
            config[1, 0] = y
            return config
        except StopIteration:
            return None

    def close(self):
        self.file.close()

    def readlines(self):
        values = []
        value = self.next()
        while value is not None:
            values.append(value)
            value = self.next()
        return np.array(values)

    def __iter__(self):
        return self.read()

    readline = next


class Writer(object):

    def __init__(self, **kwargs):
        pass

    def write(self, values, filename):
        np.savetxt(filename, values, fmt='%15.7e')
