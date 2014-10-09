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


iterfolder = os.path.split(os.getcwd())[1]
iteration = int(iterfolder[4:])

dcs = np.loadtxt("confall.ev.embed.old")
dc1s = dcs[:,0]
dc2s = dcs[:,1]

weights = np.loadtxt("confall.w")
kb = 8.31451070e-3 #kJ.mol^(-1).K^(-1)
kT = kb*300

fe = -kT*np.log(weights)
fe = fe - np.amin(fe)

nbins = 100
cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)

plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10
cb.set_label(r'Free Energy (kT units)',labelpad=pad)

plt.show()
