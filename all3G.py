import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os
import pylab
from scipy.ndimage.filters import gaussian_filter1d

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

def find_chardist(dist,q2, threshold, type):
    x=-1
    hs=(dist[1]-dist[0])/2
    dist_nb=len(q2)
    if type == 1:
        test=0
        for i in range(dist_nb):
            if q2[i]>threshold:
                test=1
            if (q2[i]<threshold and test==1) :
                x=dist[i]-hs
                test=0
    else:
        test=0
        for i in range(dist_nb):
            if q2[i]<threshold:
                test=1
            if (q2[i]>threshold and test==1) :
                x=dist[i]-hs
                test=0
        
    return x

def smooth(y, sigma):
    n_rows=len(y)
    n_cols=len(y[0])
    y_smooth=np.zeros((n_rows,n_cols))    
    for i in range(n_cols):
        y_smooth[:,i] = gaussian_filter1d(y[:,i], sigma=sigma)
    return y_smooth

# 0 dist
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

detectors = ["CE1", "CE2", "ET_B", "ET_C", "ET_D"]
signals = ["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220"]

sig=signals[5]

try:
    os.stat(sig)
except:
    os.mkdir(sig)       

dist_nb=60
det_nb=5
qq1=np.zeros((dist_nb,det_nb+1))
qq2=np.zeros((dist_nb,det_nb+1))
qq3=np.zeros((dist_nb,det_nb+1))
dd=np.zeros((dist_nb,det_nb+1))

#quantity="coverage"
quantity="delta"

if quantity == "coverage":
    index=1
else:
    index=5

qq1[1,det_nb]=100
qq3[1,det_nb]=0
    
ind=0
for det in detectors:

    # Noise results    
    filename='results_intercept_false/results_AA_prewhiten_f2_noise_' + det + '.txt'
    a= np.loadtxt(filename, dtype='f', delimiter=' ')
    
    # percentile
    q50=np.percentile(a[:,index], 50)
    q95=np.percentile(a[:,index], 95)
    q5=np.percentile(a[:,index], 5)
    
    qq1[:,det_nb]=min(q5, qq1[1,det_nb])
    qq2[:,det_nb]=q50
    qq3[:,det_nb]=max(q95, qq3[1,det_nb])
    dd[:,det_nb]=0
    print(qq1[1,det_nb])
    filename='results_intercept_false/results_AA_prewhiten_f2_' + sig + '_' + det + '.txt'

    a= np.loadtxt(filename, dtype='f', delimiter=' ')
    dist=np.unique(a[:,0]) 
    dist_nb=len(dist)
    dist_max=dist[dist_nb-1]    
    print("Distance nb:",dist_nb, dist_max)

    
    dist,q1,q2,q3=buildbox(a,index)
    qq1[:,ind]=q1
    qq2[:,ind]=q2
    qq3[:,ind]=q3
    dd[:,ind]=dist
    ind=ind+1
    
#plt.figure(figsize=(4, 3))
plt.figure()
colormap = plt.cm.gist_rainbow
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, det_nb))))
lineObjects = plt.plot(dd[:,0:det_nb], smooth(qq2[:,0:det_nb],3))
lineObjects = plt.plot(dd[:,0:det_nb], qq2[:,0:det_nb], marker="+",linestyle="")
plt.fill_between(dist, qq1[:,det_nb], qq3[:,det_nb], alpha=0.05, facecolor='b')

plt.legend(iter(lineObjects), detectors, loc='best')

plt.xlim((1,dist_max))

if quantity == "coverage":
    plt.ylim((0,1.02))
    plt.ylabel('Coverage')
    plt.legend(iter(lineObjects), detectors, loc='lower left')
else:
    plt.ylim((0.,0.6))
    plt.ylabel('$\Delta$')
    plt.legend(iter(lineObjects), detectors, loc='upper left')

plt.xlabel('Distance [kpc]')

plt.grid(True)
plt.title('s20S')

fig_name = sig + '/' + sig + '_' + quantity + '_all3G.png';
plt.savefig(fig_name)
fig_name = sig + '/' + sig + '_' + quantity + '_all3G.eps';
plt.savefig(fig_name)
fig_name = sig + '/' + sig + '_' + quantity + '_all3G.pdf';
plt.savefig(fig_name)

plt.show()


    

