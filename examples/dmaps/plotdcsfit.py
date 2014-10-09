import os
import numpy as np
import matplotlib.pyplot as plt
from densplot import hist2d

xyzfile = np.loadtxt("../confall.ev.embed")
nbins = 100

dc1s = xyzfile[:,0]#.reshape((nbins,nbins))
dc2s = xyzfile[:,1]#.reshape((nbins,nbins))
#fe = xyzfile[:,2].reshape((nbins,nbins))

#cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)
#plt.contour(dc1s, dc2s, fe, 10, colors='k', hold='on')
plt.plot(dc1s, dc2s, '.', color=(1,0,0))


xyzfile = np.loadtxt("../lsdmap/lsdmap.ev")
nbins = 100

dc1s = xyzfile[:,1]#.reshape((nbins,nbins))
dc2s = xyzfile[:,2]#.reshape((nbins,nbins))
#fe = xyzfile[:,2].reshape((nbins,nbins))

#cp = plt.scatter(dc1s, dc2s, s=10, c=fe, marker='o', linewidth=0.)
#plt.contour(dc1s, dc2s, fe, 10, colors='k', hold='on')
plt.plot(dc1s, dc2s, '.', color=(1,0.5,0.3))

# load data procedures
evs = np.loadtxt('fit.ev')
x = evs[:,1]
y = evs[:,2]

nbins = 10
#x, y, z = hist2d.make(x, y, nbins, nbins, plot_style='contour')
#cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0.)
#cp = plt.contour(x, y, z, 15, colors='k', hold='on')

plt.plot(x,y,'.',color=(0,0.5,0), markersize=10)

plt.title('Distribution of configurations used for the fit',size=16)
plt.xlabel('DC1', size=16)
plt.ylabel('DC2', size=16)

#cb = plt.colorbar(cp)
#pad = 10
#cb.set_label(r'density',labelpad=pad,size=16)

plt.show()


