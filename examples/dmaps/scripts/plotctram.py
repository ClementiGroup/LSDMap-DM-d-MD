import os
import numpy as np
import matplotlib.pyplot as plt
from dmaps.tools import tools
from densplot import hist2d

dcs = np.loadtxt("confall.ev.embed")
dc1s = dcs[:,0]
dc2s = dcs[:,1]

nbins = 100
weight = np.loadtxt('confall.w')
dc1s, dc2s, z = hist2d.make(dc1s, dc2s, nbins, nbins, plot_style='scatter', weight=weight, free_energy_plot=True, idx_smoothing=0)

plt.figure(1)
cp = plt.scatter(dc1s, dc2s, s=10, c=z, marker='o', linewidth=0., vmax=8)

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10
cb.set_label(r'Free Energy (kT)',labelpad=pad,size=16)

dc1s, dc2s, z = hist2d.make(dc1s, dc2s, nbins, nbins, plot_style='contour', weight=weight, free_energy_plot=True, idx_smoothing=8)

plt.figure(2)
levels = range(0,int(np.amax(z)),2)
cs = plt.contour(dc1s, dc2s, z, levels, colors='k', hold='on', origin='lower')
plt.clabel(cs, fmt="%2.1f", colors = 'w', fontsize=14)

dcs = np.loadtxt("ctram/confall.ev.embed")
dc1s = dcs[:,0]
dc2s = dcs[:,1]

nbins = 100
weight = np.exp(np.loadtxt('ctram/log_weight.dat'))
dc1s, dc2s, z = hist2d.make(dc1s, dc2s, nbins, nbins, plot_style='scatter', weight=weight, free_energy_plot=True, idx_smoothing=0)

cp = plt.scatter(dc1s, dc2s, s=10, c=z, marker='o', linewidth=0., vmax=8)

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10
cb.set_label(r'Free Energy (kT)',labelpad=pad,size=16)

plt.show()
