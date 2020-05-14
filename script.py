import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import math

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

type='prewhiten'
filename='results_' + type + '_fcut1000Hz_AA.txt'
a= np.loadtxt(filename, dtype='f', delimiter=' ')
dist=np.unique(a[:,1]) 
N=len(dist)
print('Number of distances:' , N)
print(len(a))

covpbb=np.zeros((N,4),dtype=float)
medbw=np.zeros((N,4),dtype=float)
MSE=np.zeros((N,4),dtype=float)
ratio=np.zeros((N,4),dtype=float)

col=['b', 'r']
lab=['lm', 'lmvar']

offset=5
for j in range(2):
    offset=j*5
    for i in range(N):
#        ext=np.where((a[:,1] == dist[i]) & (a[:,offset+6]<1))
        ext=np.where((a[:,1] == dist[i]))
        covpbb[i,0]=np.mean(a[ext,offset+2])
        covpbb[i,1]=np.std(a[ext,offset+2])
        covpbb[i,2]=max(0,covpbb[i,0]-covpbb[i,1])
        covpbb[i,3]=min(1,covpbb[i,0]+covpbb[i,1])    
        
        medbw[i,0]=np.mean(a[ext,offset+3])
        
        MSE[i,0]=np.mean(a[ext,offset+5])
        MSE[i,1]=np.std(a[ext,offset+5])
        MSE[i,2]=max(0,MSE[i,0]-MSE[i,1])
        MSE[i,3]=MSE[i,0]+MSE[i,1]
        
        ratio[i,0]=np.median(a[ext,offset+6])
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
        
        print(len(ext[0]), dist[i], ratio[i,0],ratio[i,1])        
        # Plotting


    
    plt.figure(1)
    plt.plot(dist, covpbb[:,0], col[j], label=lab[j])
    plt.fill_between(dist, covpbb[:,2], covpbb[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('Coverage probability')
    plt.ylim([0, 1.05])
    plt.title(type)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name='covpbb_' + type + '.png';
    plt.savefig(fig_name)

    
    plt.figure(2)
    plt.plot(dist, ratio[:,0], col[j], label=lab[j])
    plt.fill_between(dist, ratio[:,2], ratio[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('|true ratio - ratio|/true ratio')
    plt.ylim([0, 1.05])
    plt.title(type)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name='precision_' + type + '.png';
    plt.savefig(fig_name)
    
    plt.figure(3)
    plt.plot(dist, MSE[:,0], col[j], label=lab[j])
    plt.fill_between(dist, MSE[:,2], MSE[:,3], alpha=0.5, facecolor=col[j])
    plt.xlabel('Distance [kpc]')
    plt.ylabel('MSE')
    plt.ylim([0, 3e-7])
    plt.title(type)
    plt.grid(True)
    plt.legend(loc='best')
    fig_name='MSE_' + type + '.png';
    plt.savefig(fig_name)

#    plt.figure(4)
#    plt.plot(dist, rec_frac, col[j], label=lab[j])
#    plt.xlabel('Distance [kpc]')
#    plt.ylabel('Fraction of reconstructed g-mode')
#    plt.grid(True)
#    plt.ylim([0, 1.05])    
#    plt.legend()



plt.show()




