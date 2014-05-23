"""Reader/Writer for .ev files"""

import numpy as np
import itertools as it
import sl

class EvError(Exception):
    pass

class Reader(sl.Reader):

    def read(self):
        evs = self.next()
        while evs is not None:
            yield evs
            evs = self.next()

    def next(self):
        try:
            values = map(float, self.file.next().split()) 
            return np.array(values)
        except StopIteration:
            return None

    def readlines(self):
        evss = []
        evs = self.next()
        while evs is not None:
            evss.append(evs)
            evs = self.next()
        return np.array(evss)

    readline = next

class Writer(object):

    def write(self, evs, filename, mode='w'):
        np.savetxt(filename, evs, fmt='%15.7e')
