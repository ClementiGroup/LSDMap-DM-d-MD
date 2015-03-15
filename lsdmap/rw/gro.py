import numpy as np
import itertools as it


class GroError(Exception):
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
        return self._read_column(0, 9)

    @property
    def atoms(self):
        return self._read_column(9, 15)

    @property
    def atoms_nums(self):
        return self._read_column(15, 20)

    @property
    def box(self):
        f = open(self.filename, 'r')
        try:
            f.next()
            natoms = int(f.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")

        for idx_atom, line in it.izip(xrange(natoms), f):
            pass
        boxline = f.next()
        f.close()
        return map(float, boxline.split())


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

    fmt = {'numatoms': "%5d\n",   # number of atoms
               # coordinates output format, see http://chembytes.wikidot.com/g-grofile
               'xyz': "%8s%7s%5s%8.3f%8.3f%8.3f\n",                  # coordinates only
               'xyz_v': "%8s%7s%5s%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", # coordinates and velocities
               # unitcell
               }

    def __init__(self, **kwargs):
        try:
            self.pattern = Reader(kwargs['pattern'])
        except KeyError:
            raise NotImplementedError('a file pattern is needed to write .gro files, please use writer.open(format, pattern=...)')

    def write(self, coords, filename, idxs_atoms=None, mode='w'):

        coords = np.array(coords)

        if len(coords.shape) not in (2, 3):
            raise ValueError('coordinates should be an array 2 or 3 dimensions') 
        elif len(coords.shape) == 2:
            coords = np.array([coords])

        ncoords = coords.shape[0]
        ndim = coords.shape[1]
        if idxs_atoms is None:
            natoms = coords.shape[2]
        else:
            natoms = len(idxs_atoms)

        if ndim not in (3, 6) :
            raise ValueError('coordinates should be an array with 3 or 6 columns (number of spatial dimensions) ')

        residues = self.pattern.residues
        atoms = self.pattern.atoms
        atoms_nums = self.pattern.atoms_nums
        if idxs_atoms is not None:
            residues = [residues[idx] for idx in idxs_atoms]
            atoms = [atoms[idx] for idx in idxs_atoms]
            atoms_nums = ['%i'%(idx+1) for idx in range(len(idxs_atoms))]
        box = self.pattern.box

        with open(filename, mode) as file:
            # Atom descriptions and coords
            if ndim == 3:
                for coord in coords:
                    file.write(self.pattern.firstline)
                    file.write(self.fmt['numatoms'] %natoms)
                    for idx, [residue, atom, atom_num] in enumerate(it.izip(residues, atoms, atoms_nums)):
                        if idxs_atoms is None:
                           output_line = self.fmt['xyz'] % \
                                (residue, atom, atom_num, coord[0,idx], coord[1,idx], coord[2,idx])
                        else:
                           output_line = self.fmt['xyz'] % \
                                (residue, atom, atom_num, coord[0,idxs_atoms[idx]], coord[1,idxs_atoms[idx]], coord[2,idxs_atoms[idx]])
                        file.write(output_line)
                    for edge in box:
                        file.write('%10.5f'% edge)
                    file.write('\n')
            elif ndim == 6:
                for coord in coords:
                    file.write(self.pattern.firstline)
                    file.write(self.fmt['numatoms'] %natoms)
                    for idx, [residue, atom, atom_num] in enumerate(it.izip(residues, atoms, atoms_nums)):
                        if idxs_atoms is None:
                            output_line = self.fmt['xyz_v'] % \
                                (residue, atom, atom_num, coord[0,idx], coord[1,idx], coord[2,idx], coord[3,idx], coord[4,idx], coord[5,idx])
                        else:
                            output_line = self.fmt['xyz_v'] % \
                                (residue, atom, atom_num, coord[0,idxs_atoms[idx]], coord[1,idxs_atoms[idx]], coord[2,idxs_atoms[idx]], \
                                    coord[3,idxs_atoms[idx]], coord[4,idxs_atoms[idx]], coord[5,idxs_atoms[idx]])
                        file.write(output_line)
                    for edge in box:
                        file.write('%10.5f'% edge)
                    file.write('\n')

    def close(self):
        self.pattern.close()
