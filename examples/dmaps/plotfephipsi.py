import os
import numpy as np
import matplotlib.pyplot as plt
from densplot import hist2d

# load MD data
tphi = np.loadtxt('/Users/jordane/work/MD/aladip/phiALA2.xvg')
tpsi = np.loadtxt('/Users/jordane/work/MD/aladip/psiALA2.xvg')

x = tphi[:,1]
y = tpsi[:,1]

x = np.mod(x, 360)
y = np.mod(y+120, 360)

x, y, z = hist2d.make(x, y, 100, 100, plot_style='contour', free_energy_plot=True, idx_smoothing=3)
plt.contour(x, y, z, 15, colors='k', hold='on')

# load data procedures
tphi = np.loadtxt('phiALA2.xvg')
tpsi = np.loadtxt('psiALA2.xvg')
weight = np.loadtxt('confall.w')
#weight = None

x = tphi[:,1]
y = tpsi[:,1]

x = np.mod(x, 360)
y = np.mod(y+120, 360)

nbins = 100
x, y, z = hist2d.make(x, y, nbins, nbins, plot_style='scatter', weight=weight, free_energy_plot=True, idx_smoothing=3)
cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0.)

plt.xlabel('Phi', size=16)
plt.ylabel('Psi', size=16)

cb = plt.colorbar(cp)
pad = 10
cb.set_label(r'Free Energy (kT units)',labelpad=pad)

plt.show()


