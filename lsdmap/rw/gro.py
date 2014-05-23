import numpy as np
import itertools as it


class GroError(Exception):
    pass

class Reader(object):

    def __init__(self, filename):

        self.filename = filename
        self.file = open(filename, 'r')

    @property
    def residues(self):
        return self._read_column(0)

    @property
    def atoms(self):
        return self._read_column(1)

    @property
    def atoms_nums(self):
        return self._read_column(2)

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

        return map(float, f.next().split())

    def _read_column(self, idx):
        f = open(self.filename, 'r')
        try:
            f.next()
            natoms = int(f.next())
        except StopIteration:
            raise GroError("File ended unexpectedly when reading number of atoms.")
        column = []

        for idx_atom, line in it.izip(xrange(natoms), f):
            column.append(line.split()[idx])

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


class Writer(object):

    fmt = {'numatoms': "%5d\n",   # number of atoms
               # coordinates output format, see http://chembytes.wikidot.com/g-grofile
               'xyz': "%8s%7s%5s%8.3f%8.3f%8.3f\n",                  # coordinates only
               'xyz_v': "%8s%7s%5s%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", # coordinates and velocities
               # unitcell
               }

    def __init__(self, pattern):
        self.pattern = Reader(pattern)

    def write(self, coords, filename, format='xyz', mode='w'):

        coords = np.array(coords)

        if len(coords.shape) not in (2, 3):
            raise ValueError('coordinates should be a array 2 or 3 dimensions') 
        elif len(coords.shape) == 2:
            coords = np.array([coords])

        ncoords = coords.shape[0]
        ndim = coords.shape[1]
        natoms = coords.shape[2]

        if ndim not in (3, 6) :
            raise ValueError('coordinates should be a array with 3 or 6 rows (number of spatial dimensions) ')

        residues = self.pattern.residues
        atoms = self.pattern.atoms
        atoms_nums = self.pattern.atoms_nums
        box = self.pattern.box

        if natoms != len(atoms):
            raise GroError('pattern used to create writer object has not the same number of atoms that the coordinates provided')
        
        with open(filename , 'w') as file:
            # Atom descriptions and coords
            for coord in coords:
                file.write('Protein\n')
                file.write(self.fmt['numatoms'] %natoms)
                for idx, [residue, atom, atom_num] in enumerate(it.izip(residues, atoms, atoms_nums)):
                    output_line = self.fmt['xyz'] % \
                        (residue, atom, atom_num, coord[0,idx], coord[1,idx], coord[2,idx])
                    file.write(output_line)
                for edge in box:
                    file.write('%10.5f'% edge)
                file.write('\n')
