import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os


# 0 fcut
# 1 dist

# 2 covpbb              Fixed error
# 3 medbandwidth
# 4 absolute residual (mean)
# 5 MSE (mean)
# 6 precision (mean)

# 7 covpbb              Variable variance
# 8 medbandwidth
# 9 absolute residual (mean)
# 10 MSE (mean)
# 11 precision (mean)

name="s20-gw_10kpc16384"
#name="s11.2--LS220--GravA"
#name="s15.0--GShen--GravA"
#name="s15.0--LS220--GravA"
#name="s15.0--SFHo--GravA"
#name="s20.0--LS220--GravA"
#name="s25.0--LS220--GravA"
name="s40.0--LS220"

method='median'
type='prewhiten'
filename='results_' + type + '_fcut1000Hz_AA_' + name + '.txt'
a= np.loadtxt(filename, dtype='f', delimiter=' ')
dist=np.unique(a[:,1]) 
N=len(dist)
print('Number of distances:' , N)

try:
    os.stat(name)
except:
    os.mkdir(name)       

covpbb=np.zeros((N,4),dtype=float)
medbw=np.zeros((N,4),dtype=float)
MSE=np.zeros((N,4),dtype=float)
ratio=np.zeros((N,4),dtype=float)

col=['b', 'r']
lab=['lm', 'lmvar']
colmap=['viridis', 'jet']
mar=['o', 'x']

offset=5
for j in range(2):
    offset=j*5
    for i in range(N):
#        ext=np.where((a[:,1] == dist[i]) & (a[:,offset+6]<1))
        ext=np.where((a[:,1] == dist[i]))
        if method == 'median':
            covpbb[i,0]=np.median(a[ext,offset+2])
            covpbb[i,1]=np.median(abs(a[ext,offset+2]-covpbb[i,0]))
        else:
            covpbb[i,0]=np.mean(a[ext,offset+2])
            covpbb[i,1]=np.std(a[ext,offset+2])
            
        covpbb[i,2]=max(0,covpbb[i,0]-covpbb[i,1])
        covpbb[i,3]=min(1,covpbb[i,0]+covpbb[i,1])    
        
        if method == 'median':
            medbw[i,0]=np.median(a[ext,offset+3])
        else:
            medbw[i,0]=np.mean(a[ext,offset+3])

        if method == 'median':
            MSE[i,0]=np.median(a[ext,offset+5])
            MSE[i,1]=np.median(abs(a[ext,offset+5]-MSE[i,0]))
        else:
            MSE[i,0]=np.mean(a[ext,offset+5])
            MSE[i,1]=np.std(a[ext,offset+5])
            
        MSE[i,2]=max(0,MSE[i,0]-MSE[i,1])
        MSE[i,3]=MSE[i,0]+MSE[i,1]
        
        if method == 'median':
            ratio[i,0]=np.median(a[ext,offset+6])
            ratio[i,1]=np.median(abs(a[ext,offset+6]-ratio[i,0]))
        else:
            ratio[i,0]=np.mean(a[ext,offset+6])
            ratio[i,1]=np.std(a[ext,offset+6])

        ratio[i,2]=max(0,ratio[i,0]-ratio[i,1])
        ratio[i,3]=ratio[i,0]+ratio[i,1]
        
        sample=np.where((a[:,1] == dist[i]))
        sample_size=len(sample[0])
        ind=2+j*5
#        rec_ext=np.where((a[:,1] == dist[i]) & (a[:,ind]<0))
#        rec_frac[i]=(sample_size-len(rec_ext[0]))/float(sample_size)     

        #if i==10:
        #    print(len(ext), a[ext,offset+6])
        
        #print(len(ext[0]), dist[i], ratio[i,0],ratio[i,1])        
        # Plotting

    
    plt.figure(1)
    plt.plot(dist, covpbb[:,0], col[j], label=lab[j])
    plt.fill_between(dist, covpbb[:,2], covpbb[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('Coverage probability')
    plt.ylim([0, 1.05])
    plt.title(name)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name=name + '/covpbb_' + type + '.png';
    plt.savefig(fig_name)

#    plt.figure(4)
#    plt.scatter(a[:,offset+3], a[:,offset+2], c=a[:,1],marker=mar[j], label=lab[j],s=50)
#    if j==1:
#        plt.colorbar()
#    plt.xlim([0.000, 0.0015])
#    plt.xlabel('Estimate bandwidth median')
#    plt.ylabel('Coverage probability')
#    plt.title(name)
#    plt.grid(True)
#    plt.legend(loc='lower right')
#    fig_name=name + '/covpbb_bandwidth' + type + '.png';
#    plt.savefig(fig_name)

    plt.figure(5)
    plt.scatter(medbw[:,0], covpbb[:,0],c=dist,marker=mar[j], label=lab[j],s=50)
    if j==1:
        cbar=plt.colorbar()
        cbar.set_label("Distance (kpc)")
    plt.xlim([0.000, 0.0015])
    plt.xlabel('Estimate bandwidth median')
    plt.ylabel('Coverage probability')
    plt.title(name)
    plt.grid(True)
    plt.legend(loc='lower right')
    fig_name=name + '/covpbb_bandwidth_median_' + type + '.png';
    plt.savefig(fig_name)

    
    
    plt.figure(2)
    plt.plot(dist, ratio[:,0], col[j], label=lab[j])
    plt.fill_between(dist, ratio[:,2], ratio[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('|true ratio - ratio|/true ratio')
    plt.ylim([0, 1.05])
    plt.title(name)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name=name + '/precision_' + type + '.png';
    plt.savefig(fig_name)
   
    plt.figure(3)
    plt.plot(dist, MSE[:,0], col[j], label=lab[j])
    plt.fill_between(dist, MSE[:,2], MSE[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('MSE')
    plt.ylim([0, 1e-7])
    plt.title(name)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name=name + '/MSE_' + type + '.png';
    plt.savefig(fig_name)

#    plt.figure(4)
#    plt.plot(dist, rec_frac, col[j], label=lab[j])
#    plt.xlabel('Distance [kpc]')
#    plt.ylabel('Fraction of reconstructed g-mode')
#    plt.grid(True)
#    plt.ylim([0, 1.05])    
#    plt.legend()



plt.show()




