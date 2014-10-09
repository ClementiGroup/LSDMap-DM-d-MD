import os
import numpy as np
import matplotlib.pyplot as plt
from densplot import hist2d

xyzfile = np.loadtxt("../confall.ev.embed")
nbins = 100

weight = -np.log(np.loadtxt("../confall.w"))

dc1s = xyzfile[:,0]#.reshape((nbins,nbins))
dc2s = xyzfile[:,1]#.reshape((nbins,nbins))
#fe = xyzfile[:,2].reshape((nbins,nbins))

#cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)
#plt.contour(dc1s, dc2s, fe, 10, colors='k', hold='on')
cp=plt.scatter(dc1s, dc2s, s=13, c=weight, marker='o', linewidth=0.)


plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

cb = plt.colorbar(cp)
pad = 10
cb.set_label(r'Free Energy last iteration (kT units)',labelpad=pad)

plt.show()


