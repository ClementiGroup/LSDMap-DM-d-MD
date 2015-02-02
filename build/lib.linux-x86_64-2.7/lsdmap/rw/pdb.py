import numpy as np
import itertools as it

_known_entries = [ "ATOM  ", "HETATM", "ANISOU", "CRYST1",
    "COMPND", "MODEL", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK",
    "CONECT"]


class PDBError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):

        self.filename = filename
        self.file = open(filename, 'r')
        try:
            self.velocities = kwargs['velocities']
        except KeyError:
            self.velocities = False

    @property
    def firstline(self):
        f = open(self.filename, 'r')
        firstline = f.next()
        f.close()
        return firstline

    @property
    def natoms(self):
        f = open(self.filename, 'r')
        f.next()
        natoms = int(f.next())
        f.close()
        return natoms

    @property
    def nlines(self):
        f = open(self.filename, 'r')
        for idx, line in enumerate(f):
            pass
        nlines = (idx+1)/(self.natoms+3)
        f.close()
        return nlines

    @property
    def residues(self):
        return self._read_column(0, 8)

    @property
    def atoms(self):
        return self._read_column(8, 15)

    @property
    def atoms_nums(self):
        return self._read_column(15, 20)

    def _read_column(self, start, end):
        f = open(self.filename, 'r')
        try:
            f.next()
            natoms = int(f.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")
        column = []

        for idx_atom, line in it.izip(xrange(natoms), f):
            column.append(line[start:end].lstrip())

        f.close()
        return column


    def read(self):

        config = self.next()
        while config is not None:
            yield config
            config = self.next()

    def next(self):

        for line in self.file:
            entry = line[:6].strip()
            if record == 'END':
                break


        try:
            self.file.next()
        except StopIteration:
            return None

        try:
            natoms = int(self.file.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")

        if self.velocities:
            config = np.zeros((6, natoms), dtype='float')

            for idx_atom, line in it.izip(xrange(natoms), self.file):
                x, y, z, vx, vy, vz = map(float, line[20:].split())
                config[0, idx_atom] = x
                config[1, idx_atom] = y
                config[2, idx_atom] = z
                config[3, idx_atom] = vx
                config[4, idx_atom] = vy
                config[5, idx_atom] = vz
 
        else:
            config = np.zeros((3, natoms), dtype='float')

            for idx_atom, line in it.izip(xrange(natoms), self.file):
                x, y, z = map(float, line[20:].split()[:3])
                config[0, idx_atom] = x
                config[1, idx_atom] = y
                config[2, idx_atom] = z

        try:
            box_line = self.file.next()
        except StopIteration:
            raise GroError("File ended unexpectedly when reading box line.")

        return config


    def _skip(self):
        try:
            self.file.next()
        except StopIteration:
            return None

        try:
            natoms = int(self.file.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")

        for atom in it.izip(xrange(natoms), self.file):
            pass

        try:
            self.file.next()
        except StopIteration:
            raise GroError("File ended unexpectedly when reading box line.")
        return None


    def readlines(self, *args):
        if len(args) == 0:
            configs = []
            config = self.next()
            while config is not None:
                configs.append(config)
                config = self.next()
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
            raise GroError("invalid number of arguments to readlines")

        return np.array(configs)

    def close(self):
        self.file.close()

    def __iter__(self):
        return self.read()

    readline = next


class Writer(object):
    pass
