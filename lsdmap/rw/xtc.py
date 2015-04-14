import numpy as np
import itertools as it

#from xdrfile.libxdrfile2 import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc_numframes, read_xtc, DIM, exdrOK

class XtcError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):
        try:
            import xdrfile.libxdrfile2 as libxdr
        except ImportError:
            raise ImportError("libxdrfile2 library needed for xtc trajectories!")

        self.filename = filename
        self.velocities = False
        self.file = libxdr.xdrfile_open(filename, 'r')
        self.natoms = libxdr.read_xtc_natoms(self.filename)
        self.nframes, self.offsets = libxdr.read_xtc_numframes(self.filename)
        self.box = np.zeros((3,3),dtype=np.float32)
        self.config = np.zeros((self.natoms,3),dtype=np.float32)

    @property
    def box(self):
        raise NotImplementedError

    def read(self):
        config,status = self.next()
        while status == libxdr.exdrOK:
            yield config
            config,status = self.next()

    def next(self):
        status = libxdr.read_xtc(self.file, self.box, self.config)[0]
        return np.array(self.config,copy=True,dtype=np.float64), status

    def _skip(self):
        ## NOT DONE.
        try:
            self.file.next()
        except StopIteration:
            return None

        try:
            self.file.next()
        except StopIteration:
            raise XtcError("File ended unexpectedly when reading box line.")
        return None


    def readlines(self, *args):
        if len(args) == 0:
            configs = np.zeros((self.nframes,self.natoms,3),dtype=np.float64)
            config, status = self.next()
            i = 0
            while status == libxdr.exdrOK:
                configs[i,:,:] = config
                config, status = self.next()
                i += 1
        elif len(args) == 1:
            lines = args[0]
            if isinstance(lines, int):
                lines = [lines]
            else:
                lines = list(set(lines))
                lines.sort()
            lines = np.array(lines)
            lines = np.hstack((-1, lines))
            sklines = np.diff(lines) - 1
            configs = []
            for skline in sklines:
                for idx in xrange(skline):
                    self._skip()
                config = self.next()
                configs.append(config)
        else:
            raise XtcError("invalid number of arguments to readlines")

        return configs

    def close(self):
        libxdr.xdrfile_close(self.file)

    def __iter__(self):
        return self.read()

    readline = next
