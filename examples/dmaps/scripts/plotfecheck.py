import sys
import os
import numpy as np
import ConfigParser
from scipy import ndimage
from dmaps.tools import tools

import matplotlib.pyplot as plt

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

    if ndcs == 3:
        fe = np.sum(fe, axis=2)
    cp = plt.contour(x, y, fe, 13, colors='k', hold='on')

    dcs = np.loadtxt("confall.ev.embed.old")
    dc1s = dcs[:,0]
    dc2s = dcs[:,1]

    weights = np.loadtxt("confall.w")
    fe = -np.log(weights)
    fe = fe - np.amin(fe)

    cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)

    plt.xlabel('DC1', size=16)
    plt.ylabel('DC2', size=16)

    cb = plt.colorbar(cp)
    pad = 10
    cb.set_label(r'Free Energy (kT units)',labelpad=pad)

    plt.show()

else:
    x = bins
    fe = np.zeros((nbins,)*ndcs)
    fe[nebins_idxs] = free_energy
    fe = (fe - np.amin(fe))/(8.314e-3*300)

    plt.plot(x, fe, '-ro')

    dc = np.loadtxt("confall.ev.embed.old")

    weights = np.loadtxt("confall.w")
    fe = -np.log(weights)
    fe = (fe - np.amin(fe))

    plt.plot(dc, fe, 'o', color=(0,0,0), markersize= 5)

    plt.show()
