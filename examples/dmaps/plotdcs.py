import sys
import numpy as np

from densplot import hist2d
import matplotlib.pyplot as plt


# load MD data
tphi = np.loadtxt('/Users/jordane/work/MD/aladip/phiALA2.xvg')
tpsi = np.loadtxt('/Users/jordane/work/MD/aladip/psiALA2.xvg')

x = tphi[:,1]
y = tpsi[:,1]

x = np.mod(x, 360)
y = np.mod(y+120, 360)

x, y, z = hist2d.make(x, y, 100, 100, plot_style='contour', free_energy_plot=True, idx_smoothing=3)
plt.contour(x, y, z, 15, colors='k', hold='on')

phifile = np.loadtxt("phiALA2.xvg")
phi = phifile[:,1]

psifile = np.loadtxt("psiALA2.xvg")
psi = psifile[:,1]

#evfile = np.loadtxt("fit.ev")
#evfile = np.loadtxt("lsdmap.ev")
evfile = np.loadtxt("confall.ev.embed")
dc1s = evfile[:,0]

phi = np.mod(phi,360)
psi = np.mod(psi+120, 360)

cp = plt.scatter(phi, psi, s=10, c=dc1s, marker='o', linewidth=0.)

plt.xlabel('Phi', size=16)
plt.ylabel('Psi', size=16)

cb = plt.colorbar(cp)
pad = 10

plt.show()





