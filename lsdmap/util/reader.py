import os, sys
import numpy as np
import linecache


class GroFile(object):
    """
    Grofile(filename)

    Class used to read .gro files (GROMACS).

    Parameters
    ----------
        filename: string
            name of .gro file

    Attributes
    ----------
        natoms: number of atoms per configuration
        nlines_per_frame: number of lines in the file per configuration (usually natoms+3)
        nlines: total number of lines in the file
        nframes: total number of configurations
        coords: list containing the coordinates of all the configurations in the file
            coords = [coords_config1, coords_config2, ..., coords_congfign]
        where coords_configX is a 1D numpy array containing the coordinates of 
        configuration X given in the format:
            [x1, x2, ..., xn, y1, y2, ...,yn , z1, z2, ...,zn] 

    """
    def __init__(self, filename):

        self.filename = filename
        self.natoms = self.get_natoms()
        self.nlines_per_frame = self.natoms+3
        self.nlines = self.get_nlines()
        self.nframes = self.get_nframes()
        self.coords = self.get_coords()


    def get_natoms(self):
        natoms=int(linecache.getline(self.filename, 2))
        return natoms


    def get_nlines(self):
        with open(self.filename, 'r') as file:
            nlines = sum(1 for line in file)
        return nlines


    def get_nframes(self):
        assert self.nlines%(self.nlines_per_frame)==0, "number of lines in %s is not a multiple of natoms+3" %self.filename
        nframes=self.nlines/self.nlines_per_frame
        return nframes


    def get_nframes(self):
        assert self.nlines%(self.nlines_per_frame)==0, "number of lines in %s is not a multiple of natoms+3" %self.filename
        nframes=self.nlines/self.nlines_per_frame
        return nframes


    def get_coords(self):

        _natoms=self.natoms

        coords=[]
        shift=0

        with open(self.filename, 'r') as file:
            for frame in xrange(self.nframes):
                start=frame*self.nlines_per_frame+2-shift
                end=(frame+1)*self.nlines_per_frame-1-shift
                idx_lines=range(start,end)

                coord=np.zeros(3*self.natoms, dtype='float')
                idx_atom=0

                for idx, line in enumerate(file):
                    if idx in idx_lines:
                        x,y,z = map(float,line.split()[3:6])
                        coord[idx_atom]=x
                        coord[_natoms+idx_atom]=y
                        coord[2*_natoms+idx_atom]=z
                        idx_atom+=1
                        if idx_atom==self.natoms:
                            coords.append(coord)
                            shift+=idx+1
                            break 

        
        return np.array(coords)


# The following classes and functions are not directly used to compute LSDMap

    def get_coords_from_cvs(self, cvs, filenames, idx_column=2):
        
        if not isinstance(cvs, np.ndarray): cvs = np.array(cvs)
        ndim = len(cvs.shape)

        if not isinstance(filenames, list): filenames = [filenames]
        nfiles = len(filenames)

        if ndim == 1: ncvs = len(cvs);
        if ndim == 2: ncvs = cvs.shape[1];

        if nfiles != ncvs: raise ValueError("the number of cvs and filenames are not the same!")

        cvs_from_files = np.empty((self.nframes, ncvs))
        for idx_file, filename in enumerate(filenames):
            with open(filename, 'r') as file:
                for idx_line, line in enumerate(file):
                    line = [float(cv) for cv in line.split() if cv]
                    cvs_from_files[idx_line, idx_file] = line[idx_column-1]

        if ndim == 1:
            dist = np.sqrt((cvs-cvs_from_files)**2)
            idx = dist.argmin()
            return [idx, self.coords[idx]]


        if ndim == 2:
            idx_configs = []
            for cv in cvs:
                dist = np.sqrt(np.sum((cv-cvs_from_files)**2, axis=1))
                idx = dist.argmin()
                idx_configs.append([idx, self.coords[idx]])

            return idx_configs
                

    def save_coords(coords_idxs, filename):
        raise NotImplementedError
        


class OneColumnFile(object):
    """
    OneColumnFile(filename)

    Class used to read text files containing a single column of values.

    """
    def __init__(self, filename):

        self.filename = filename
        self.nlines = self.get_nlines()

    def get_nlines(self):
        with open(self.filename, 'r') as file:
            nlines = sum(1 for line in file)
        return nlines

    def read(self):
 
        with open(self.filename, "r") as file:
            value = np.array(map(float, file.read().splitlines()), dtype='float')

        return value

class EpsFile(OneColumnFile):
    pass

class WFile(OneColumnFile):
    pass


class EvFile(object):
    """
    EvFile(filename)

    Class used to read .ev files (LSDMap output file containing the eigenvectors).

    Parameters
    ----------
        filename: string, name of .ev file.

    """


    def __init__(self, filename):

        self.filename = filename
        self.nlines = self.get_nlines()


    def get_nlines(self):
        with open(self.filename, 'r') as file:
            nlines = sum(1 for line in file)
        return nlines


    def get_evs(self, nums=None):

        with open(self.filename,'r') as evfile:
            first_line = map(float, evfile.readline().split())
            self.nevs = len(first_line)
            evs = np.empty((self.nlines, self.nevs))
            evs[0][:]= np.array(first_line)
            for idx, line in enumerate(evfile):
                evs[idx+1][:] = np.array(map(float, line.split()))
        return evs

