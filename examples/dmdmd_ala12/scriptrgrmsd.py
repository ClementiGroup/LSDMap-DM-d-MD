import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from densplot import hist2d
ref_file='start'
#subprocess.call("echo 2 | trjconv -s " + ref_file +".gro -f " +ref_file+".gro -o "+ref_file+"_ha.gro", shell=True)
a=np.arange(1,20,1)
b=np.arange(20,50,10)
b2=np.arange(50,700,50)
c=np.append(a,b)
c=np.append(c,b2)
for i in c:
        plt.clf()
        w_file='results/iter'+str(i)+'.w'
        gro_file='results/iter'+str(i)+'.gro'


        subprocess.call("echo 0 0 | g_rms -xvg none -f " +gro_file+ " -s "+ref_file+".gro 1>/dev/null 2>&1", shell=True)
        subprocess.call("echo 0 | g_gyrate -xvg none -f " +gro_file+ " -s " +gro_file+ " 1>/dev/null 2>&1", shell=True)

        tphi = np.loadtxt('rmsd.xvg')
        tpsi = np.loadtxt('gyrate.xvg')
        weight = np.loadtxt(w_file)

        x = tphi[:,1]
        y = tpsi[:,1]

        nbins = 100
        x, y, z = hist2d.make(x, y, nbins, nbins, plot_style='scatter', weight=weight, free_energy_plot=True, idx_smoothing=3)
        cp = plt.scatter(x, y, s=10, c=z, marker='o', linewidth=0.)#, vmax=3.5)
        print "Free energy in range:", min(z), max(z)
        plt.xlabel('RMSD (nm)', size=16)
        plt.ylabel('Rg (nm)', size=16)

        cb = plt.colorbar(cp)
        pad = 10
        cb.set_label(r'Free Energy (kcal/mol units)',labelpad=pad)
        plt.axis([0.0,1.0, 0.4,1.4])
        plt.savefig('plot_rgrmsd_iter'+str(i)+'.png', bbox_inches='tight')
