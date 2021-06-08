import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os
import pylab

signals = ["s11", "s15", "s15S", "s15G", "s20", "s20S", "s25", "s40"]
detectors = ["","aLIGO", "AdV", "CE1", "CE2", "ET_B", "ET_C", "ET_D",""]

filename='distance.txt'
a= np.loadtxt(filename, dtype='f', delimiter=' ',usecols=range(8))

plt.figure()
index=range(1,8,1)
dd=np.zeros((8,7))

for i in range(0,8,1):
    dd[i,:]=index

lineObjects=plt.semilogy(dd.T, a, marker='d',markersize=10,linestyle='')
plt.grid()
plt.xlim((0,8))
plt.ylim((5,3000))

ax=plt.gca()
labels = [item.get_text() for item in ax.get_xticklabels()]
for i in range(0,9):
    labels[i] = detectors[i]
ax.set_xticklabels(labels)

plt.ylabel('$d_r$ [kpc]')
plt.legend(iter(lineObjects), signals,loc='upper left',ncol=5,numpoints=1)

fig_name = 'dist_allwvfs_2G3G.png';
plt.savefig(fig_name)
fig_name = 'dist_allwvfs_2G3G.eps';
plt.savefig(fig_name)
fig_name = 'dist_allwvfs_2G3G.pdf';
plt.savefig(fig_name)

plt.show()
