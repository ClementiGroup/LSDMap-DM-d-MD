import sys
import os
from math import exp
import ConfigParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from densplot import hist2d


a=np.arange(1,20,1)
b=np.arange(20,50,10)
b2=np.arange(50,700,50)
c=np.append(a,b)
c=np.append(c,b2)
for i in c:
        print i
        plt.clf()
        plt.xlabel('DC1', size=16)
        plt.ylabel('DC2', size=16)

        dcs = np.loadtxt("results/iter"+str(i)+".ev")
        dcs = dcs[:,0:3]

        #prevent flipping
        if max(dcs[:,1])<min(dcs[:,1])*-1:
          dc1s = -dcs[:,1]
        else:
          dc1s = dcs[:,1]
        if max(dcs[:,2])<min(dcs[:,2])*-1:
          dc2s = -dcs[:,2]
        else:
          dc2s = dcs[:,2]

        dc1s, dc2s, z = hist2d.make(dc1s, dc2s, 100, 100, plot_style='scatter', free_energy_plot=True, idx_smoothing=3)
        cp = plt.scatter(dc1s, dc2s, s=10, c=z, marker='o', linewidth=0.)

        plt.xlabel('DC1', size=16)
        plt.ylabel('DC2', size=16)

        cb = plt.colorbar(cp)
        pad = 10
        cb.set_label(r'population',labelpad=pad,size=16)
        #cb.set_label(r'phi',labelpad=pad,size=16)
        plt.savefig("plot_lsdmap_iter"+str(i)+'.png', bbox_inches='tight')
