"""Reader for .ev files"""

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
            return map(float, self.file.next().split())
        except StopIteration:
            return None

    def readlines(self):
        evss = []
        evs = self.next()
        while evs is not None:
            evss.append(evs)
            evs = self.next()
        return np.array(evss)

