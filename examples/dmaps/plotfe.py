import sys
import argparse
import numpy as np

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-u",
    default=False,
    action='store_true',
    dest='unweighted')

args = parser.parse_args()

if args.unweighted:
    xyzfile = np.loadtxt("fe/hist_unweighted.xyz")
else:
    xyzfile = np.loadtxt("fe/hist.xyz")

nbins = np.sqrt(xyzfile.shape[0])

dc1s = xyzfile[:,0].reshape((nbins,nbins))
dc2s = xyzfile[:,1].reshape((nbins,nbins))
fe = xyzfile[:,2].reshape((nbins,nbins))

#cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)
cp = plt.contourf(dc1s, dc2s, fe)
plt.contour(dc1s, dc2s, fe, cp.levels, colors='k', hold='on')

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10

plt.show()
