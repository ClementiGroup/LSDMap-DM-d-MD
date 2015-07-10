import sys
import os
from math import exp
import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
from dmaps.tools import tools
from densplot import hist2d

config = ConfigParser.SafeConfigParser()
config.read('../config.ini')
ndcs = config.getint('DMAPS', 'ndcs')

iterfolder = os.path.split(os.getcwd())[1]
iteration = int(iterfolder[4:])

bfile = "../iter" + str(iteration-1) + "/fe/bins.xyz"
bins = np.loadtxt(bfile, dtype="f4")
nbins = bins.shape[0]

hfile = "../iter" + str(iteration-1) + "/fe/hist.dat"
nebins_idxs = np.genfromtxt(hfile, usecols=tuple(range(ndcs)), dtype="i4")
free_energy = np.genfromtxt(hfile, usecols=ndcs+1, dtype="f4")
gradient = np.genfromtxt(hfile, usecols=tuple(range(ndcs+2,2*ndcs+2)), dtype="f4")

if ndcs > 1:
    x = bins[:,0][:, np.newaxis].dot(np.ones((1, nbins)))
    y = np.ones((nbins,1)).dot(bins[:,1][np.newaxis])

    fe = np.zeros((nbins,)*ndcs)
    fe[tuple([nebins_idxs[:,jdx] for jdx in xrange(ndcs)])] = free_energy
    fe = fe - np.amin(fe) # make the free energy positive to avoid dashed plots
    cp = plt.contour(x, y, fe, 13, colors='k', hold='on')

    plt.xlabel('DC1', size=16)
    plt.ylabel('DC2', size=16)

    dcs = np.loadtxt("confall.ev.embed.old")
    nstride = config.getint('DMAPS', 'nstride')

    ntraj = [1, 2, 3]
    pltstyle = ['-or', '-ob', '-og'] 

    for idx in xrange(len(ntraj)):
        dcs = dcs[:,0:2]
        dc1s = dcs[ntraj[idx]*nstride:(ntraj[idx]+1)*nstride,0]
        dc2s = dcs[ntraj[idx]*nstride:(ntraj[idx]+1)*nstride,1]
        plt.plot(dc1s, dc2s, pltstyle[idx])
    plt.show()

else:

    x = bins
    fe = np.zeros((nbins,)*ndcs)
    fe[nebins_idxs] = free_energy
    fe = (fe - np.amin(fe))/(8.314e-3*300)
    plt.plot(x, fe, '-ro')

    x = bins
    fe = np.zeros((nbins,)*ndcs)
    fe[nebins_idxs] = gradient/3000
    plt.plot(x, fe, '-go')

    dc = np.loadtxt("confall.ev.embed.old")
    bins, grid = tools.do_grid(dc, 100)
    weights = np.ones(dc.shape[0])
    fe = tools.compute_free_energy(grid, 1, weights, 20, 1)
    plt.plot(bins, np.exp(-fe), 'o', color=(0,0,0), markersize= 5)

    plt.show()