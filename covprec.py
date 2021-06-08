import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os
import pylab

def minmax(a, ind_mean, ind_var):
    N,Y=a.shape
    
    out=np.zeros((N,2))    
    for i in range(N):
        out[i,0]=max(0, a[i,ind_mean]-a[i,ind_var])
        out[i,1]=min(1, a[i,ind_mean]+a[i,ind_var])
    return out

def buildbox(a,index):
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    #print('Number of distances: ', dist_nb)
    
    data = []
    box_name = []
    q1 = []
    q2 = []
    q3 = []
    for i in range(dist_nb):
        ext=np.where((a[:,0] == dist[i]))
        data_boxi = a[ext[0][0]:ext[0][-1],index]
        data.append(data_boxi)
        if (i%5 == 0):
            box_name.append(str(dist[i]))
        else:
            box_name.append('')
        q1.append(np.percentile(data_boxi, 25))
        q2.append(np.percentile(data_boxi, 50))
        q3.append(np.percentile(data_boxi, 75))
    return dist,q1,q2,q3


# 0 dist
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

detectors = ["aLIGO", "CE1", "CE2", "ET_B", "ET_C", "ET_D"]
signals = ["s11.2--LS220", "s15.0--GShen", "s15.0--SFHo", "s15.0--LS220", "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220"]

#signals = ["s11.2--LS220", "s15.0--GShen", "s15.0--SFHo", "s15.0--LS220", "s20.0--LS220", "s20.0--SFHo"]




det="aLIGO"

name="s20.0--SFHo"
#name="s11.2--LS220"
#name="s15.0--GShen"
#name="s15.0--LS220"
#name="s15.0--SFHo"
#name="s20.0--LS220"
#name="s25.0--LS220"
#name="s40.0--LS220"

try:
    os.stat(name)
except:
    os.mkdir(name)       

dist_nb=38
sig_nb=8
qq1=np.zeros((dist_nb,sig_nb))
qq2=np.zeros((dist_nb,sig_nb))
qq3=np.zeros((dist_nb,sig_nb))
dd=np.zeros((dist_nb,sig_nb))

qq_p1=np.zeros((dist_nb,sig_nb))
qq_p2=np.zeros((dist_nb,sig_nb))
qq_p3=np.zeros((dist_nb,sig_nb))

ind=0

for sig in signals:
    title=name + ' ' + det
    filename='results/results_AA_prewhiten_f2_' + sig + '_' + det + '.txt'

    a= np.loadtxt(filename, dtype='f', delimiter=' ')
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    print("Distance nb:",dist_nb, dist[10])

    # coverage
    dist,q1,q2,q3=buildbox(a,1)
    qq1[:,ind]=q1
    qq2[:,ind]=q2
    qq3[:,ind]=q3
    dd[:,ind]=dist

    dist,q1,q2,q3=buildbox(a,5)
    qq_p1[:,ind]=q1
    qq_p2[:,ind]=q2
    qq_p3[:,ind]=q3
    ind=ind+1
    
#plt.figure(figsize=(4, 3))
plt.figure()
colormap = plt.cm.gist_rainbow
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, sig_nb))))
lineObjects = plt.plot(qq_p2, qq2)
plt.legend(iter(lineObjects), signals, loc='best')

#plt.xlim((1,19))
#plt.ylim((0,1))
plt.xlabel('Precision')
plt.ylabel('Coverage')
plt.grid(True)

#fig_name = det + '_allwvfs.png';
#plt.savefig(fig_name)
#fig_name = det + '_allwvfs.eps';
#plt.savefig(fig_name)
#fig_name = det + '_allwvfs.pdf';
#plt.savefig(fig_name)

plt.show()


    

