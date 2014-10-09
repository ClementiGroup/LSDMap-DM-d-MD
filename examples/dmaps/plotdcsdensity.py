import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from densplot import hist2d

# load MD data

iterfolder = os.path.split(os.getcwd())[1]
iteration = int(iterfolder[4:])

xyzfile = np.loadtxt("../iter" + str(iteration-1) + "/fe/hist.xyz")
nbins = np.sqrt(xyzfile.shape[0])

dc1s = xyzfile[:,0].reshape((nbins,nbins))
dc2s = xyzfile[:,1].reshape((nbins,nbins))
fe = xyzfile[:,2].reshape((nbins,nbins))
fe = fe - np.amin(fe)

cp = plt.contour(dc1s, dc2s, fe, 13, colors='k', hold='on')

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

dcs = np.loadtxt("confall.ev.embed.old")
dc1s = dcs[:,0]
dc2s = dcs[:,1]

z = np.loadtxt("phiALA2.xvg")
z = z[:,1]

z = np.mod(z, 360)

#nbins = 100
#dc1s, dc2s, z = hist2d.make(dc1s, dc2s, 100, 100, plot_style='scatter', free_energy_plot=True, idx_smoothing=3)
cp = plt.scatter(dc1s, dc2s, s=10, c=z, marker='o', linewidth=0.)

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10
#cb.set_label(r'-log(population)',labelpad=pad)
cb.set_label(r'phi',labelpad=pad,size=16)

plt.show()
