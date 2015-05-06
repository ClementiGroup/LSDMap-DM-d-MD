""" Get input coordinates


Goal:
- Support commonly used formats: xtc, trr, gro, pdb, xyz
- Code that is easy to maintain and read, with minimal
  outside dependencies.

"""

import numpy as np
import linecache

global supported_formats
supported_formats = ["gro","xtc"]

#############################################################################
# General utility functions
#############################################################################
def get_coordinates(filename,idxs=None):
    ext = filename.split(".")[-1]
    reader = {"gro":gro_coords,"xtc":xtc_coords}

    if ext not in supported_formats:
        raise IOError("Extension %s is not supported. Use %s" % (ext,supported_formats.__str__()))
    else:
        coords = reader[ext](filename,idxs=idxs)
    return coords

def get_nframes_natoms(filename):
    ext = filename.split(".")[-1]
    reader = {"gro":gro_nframes,"xtc":xtc_nframes}

    if ext not in supported_formats:
        raise IOError("Extension %s is not supported. Use %s" % (ext,supported_formats.__str__()))
    else:
        nframes, natoms = reader[ext](filename)[:2]
    return nframes,natoms

def supports_parallel_reading(filename):
    """ Query if filetype can be read in parallel """
    parallel_reading = ["gro","xtc"]
    ext = filename.split(".")[-1]
    if ext in parallel_reading:
        return True
    else:
        return False

#############################################################################
# Format specific helper functions: xtc format
#############################################################################
def xtc_coords(filename,idxs=None):
    """ Get the coordinates """

    try:
        import xdrfile.libxdrfile2 as libxdr
    except ImportError:
        raise ImportError("libxdrfile2 library needed for xtc trajectories!")

    file = libxdr.xdrfile_open(filename, "r")
    natoms = libxdr.read_xtc_natoms(filename)
    nframes, offsets = libxdr.read_xtc_numframes(filename)
    box = np.zeros((3,3),dtype=np.float32)
    config = np.zeros((natoms,3),dtype=np.float32)

    if idxs is not None:
        # Seek to read frames specified by idxs.
        nframes = len(idxs)
        coords = np.zeros((nframes,3,natoms),dtype="float")
        for i in xrange(nframes):
            libxdr.xdr_seek(file,offsets[idxs[i]],0)
            status,step,timestmp,prec = libxdr.read_xtc(file, box, config)
            coords[i,:,:] = config.astype("float").T
    else:
        coords = np.zeros((nframes,3,natoms),dtype=np.float32)
        status = libxdr.exdrOK
        i = 0
        while status == libxdr.exdrOK:
            status,step,timestmp,prec = libxdr.read_xtc(file, box, config)
            coords[i,:,:] = config.astype("float").T
            i += 1

    libxdr.xdrfile_close(file)

    return coords

def xtc_nframes(filename):
    """ Get the coordinates """

    try:
        import xdrfile.libxdrfile2 as libxdr
    except ImportError:
        raise ImportError("libxdrfile2 library needed for xtc trajectories!")

    file = libxdr.xdrfile_open(filename, "r")
    natoms = libxdr.read_xtc_natoms(filename)
    nframes, offsets = libxdr.read_xtc_numframes(filename)
    return nframes, natoms

#############################################################################
# Format specific helper functions: gro format
#############################################################################
def gro_coords(filename,idxs=None):
    """ Get the coordinates. gro files have fixed number of lines per frame """
    
    nframes, natoms, framesize, lines_per_frame = gro_nframes(filename)
    
    grabxyz = lambda x: [float(x.split()[3]),float(x.split()[4]),float(x.split()[5])]

    if idxs is not None:
        nframes = len(idxs)
        coords = np.zeros((nframes,3,natoms),dtype="float")
        with open(filename,"rb") as f:
            for i in xrange(nframes):
                for j in xrange(1,natoms + 1):
                    # Parse coordinates from frame
                    line_num = (lines_per_frame*idxs[i]) + 2 + j
                    framexyz = grabxyz(linecache.getline(filename,line_num))
                    coords[i,:,j - 1] = np.array(framexyz,dtype="float")
    else:
        coords = np.zeros((nframes,3,natoms),dtype="float")
        with open(filename,"rb") as f:
            for i in xrange(nframes):
                for j in xrange(1,natoms + 1):
                    # Parse coordinates from frame
                    line_num = (lines_per_frame*i) + 2 + j
                    framexyz = grabxyz(linecache.getline(filename,line_num))
                    coords[i,:,j - 1] = np.array(framexyz,dtype="float")

    return coords

def gro_nframes(filename):
    """ Get the number atoms and frames """

    with open(filename) as f:
        f.readline()
        natoms = int(f.readline())
    lines_per_frame = natoms + 3

    # Get line count without loading whole file into memory
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    nlines = i + 1

    if (nlines % lines_per_frame) != 0:
        raise IOError("Number of lines does not evenly divide ")
    else:
        nframes = nlines/lines_per_frame

    # Get framesize in bytes
    framesize = 0
    with open(filename) as f:
        for i in xrange(lines_per_frame):
            framesize += len(f.readline())
    return nframes, natoms, framesize, lines_per_frame

#############################################################################
# Format specific helper functions: xyz format
#############################################################################
def xyz_coords(filename,idxs=None):
    """ Get the coordinates """
    
    raise NotImplementedError

    nframes, natoms, framesize, lines_per_frame = gro_nframes(filename)
    
    grabxyz = lambda x: [float(x.split()[3]),float(x.split()[4]),float(x.split()[5])]

    if idxs is not None:
        nframes = len(idxs)
        coords = np.zeros((nframes,3,natoms),dtype="float")
        with open(filename,"rb") as f:
            for i in xrange(nframes):
                # Seek to frame of interest
                f.seek(idxs[i]*framesize)
                frame = f.read(framesize).split("\n")[2:2+natoms]
                print frame
                # Parse coordinates from frame
                coords[i,:,:] = np.array(map(grabxyz,frame),dtype="float").T
    else:
        coords = np.zeros((nframes,3,natoms),dtype="float")
        with open(filename,"rb") as f:
            for i in xrange(nframes):
                frame = f.read(framesize).split("\n")[2:2+natoms]
                # Parse coordinates from frame
                coords[i,:,:] = np.array(map(grabxyz,frame),dtype="float").T

    return coords

def xyz_nframes(filename):
    """ Get the number atoms and frames """

    with open(filename) as f:
        f.readline()
        firstline = f.readline()
        firstbytes = len(firstline)
        natoms = int(firstline)
    lines_per_frame = natoms

    # Get line count without loading whole file into memory
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    nlines = i + 1

    if ((nlines - 1) % lines_per_frame) != 0:
        raise IOError("Number of lines does not evenly divide ")
    else:
        nframes = (nlines - 1)/lines_per_frame

    # Get framesize in bytes
    framesize = 0
    with open(filename) as f:
        firstline = f.readline()
        for i in xrange(lines_per_frame):
            framesize += len(f.readline())

    return nframes, natoms, framesize, lines_per_frame

