import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

# required options
parser.add_argument("--lsdmap", action='store_true', default=False)
parser.add_argument("--fit", action='store_true', default=False)
parser.add_argument("--embed", action='store_true', default=True)

args = parser.parse_args()


if args.lsdmap:
    evfile = np.loadtxt("lsdmap/lsdmap.ev")
    dc1s = evfile[:,1]
    dc2s = evfile[:,2]
    phi = np.loadtxt('lsdmap/phiALA2.xvg')
elif args.fit:
    evfile = np.loadtxt("fit/fit.ev")
    dc1s = evfile[:,1]
    dc2s = evfile[:,2]
    phi = np.loadtxt('fit/phiALA2.xvg')
elif args.embed:
    evfile = np.loadtxt("confall.ev.embed")
    dc1s = evfile[:,0]
    dc2s = evfile[:,1]
    phi = np.loadtxt('psiALA2.xvg')

phi = phi[:,1]
phi = np.mod(phi,360)

cp = plt.scatter(dc1s, dc2s, s=10, c=phi, marker='o', linewidth=0.)

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10

cb.set_label(r'Psi',labelpad=pad, size=16)

plt.show()





