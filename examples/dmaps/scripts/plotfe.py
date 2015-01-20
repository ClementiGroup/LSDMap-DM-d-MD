import sys
import os
from math import exp, floor
import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
from dmaps.tools import tools
from densplot import hist2d

config = ConfigParser.SafeConfigParser()
config.read('../config.ini')
ndcs = config.getint('DMAPS', 'ndcs')

bfile = "fe/bins.xyz"
bins = np.loadtxt(bfile, dtype="f4")
nbins = bins.shape[0]
steps = bins[1,:] - bins[0,:]

hfile = "fe/hist.dat"
nebins_idxs = np.genfromtxt(hfile, usecols=tuple(range(ndcs)), dtype="i4")
nebins_idxs_s = np.genfromtxt(hfile, usecols=ndcs, dtype="i4")
free_energy = np.genfromtxt(hfile, usecols=ndcs+1, dtype="f4")
gradient = np.genfromtxt(hfile, usecols=tuple(range(ndcs+2,2*ndcs+2)), dtype="f4")

evfile = "ctram/confall.ev.embed"
dcs = np.loadtxt(evfile)
if ndcs == 1:
    dcs = dcs[:,np.newaxis]

if ndcs >= 2:
    idxs_bins = tuple([np.zeros(dcs.shape[0], dtype="i4") for jdx in xrange(ndcs)])
    for jdx in range(ndcs):
        for idx in range(dcs.shape[0]):
            idxs_bins[jdx][idx] = int(floor((dcs[idx,jdx] - bins[0][jdx])/steps[jdx] - 0.5))
    fe = np.zeros((nbins,)*ndcs)
    fe.fill(np.amax(free_energy))
    fe[tuple([nebins_idxs[:,jdx] for jdx in xrange(ndcs)])] = free_energy/(8.314e-3*300)
    fe = fe[idxs_bins] - np.amin(fe)

    if ndcs == 3:
        fe = np.sum(fe,axis=2)
    cp = plt.scatter(dcs[:,0], dcs[:,1], s=10, c=fe, marker='o', linewidth=0.)

    plt.xlabel('DC1', size=16)
    plt.ylabel('DC2', size=16)

    cb = plt.colorbar(cp)
    pad = 10
    cb.set_label(r'Free Energy (kT)',labelpad=pad,size=16)

    plt.show()

elif ndcs == 1:
    x = bins
    fe = np.zeros((nbins,)*ndcs)
    fe[nebins_idxs] = free_energy
    fe = (fe - np.amin(fe))/(8.314e-3*300)
    plt.plot(x, fe, '-ro')

    plt.show()
