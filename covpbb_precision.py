import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os


# 0 dist
# 1 covpbb           
# 2 medbandwidth
# 3 absolute residual
# 4 MSE
# 5 precision


name="s20-gw"
#name="s11.2--LS220--GravA"
#name="s15.0--GShen--GravA"
#name="s15.0--LS220--GravA"
#name="s15.0--SFHo--GravA"
#name="s20.0--LS220--GravA"
#name="s25.0--LS220--GravA"
#name="s40.0--LS220"

detectors = ["aLIGO", "CE1", "CE2", "ET_B", "ET_C", "ET_D"]
signals = ["s20-gw", "s11.2--LS220", "s15.0--GShen", "s15.0--LS220", "s20.0--LS220", "s25.0--LS220", "s40.0--LS220"]

detectors=["aLIGO"]
signals=["s20-gw"]
method='median'
type='prewhiten_f2'

# Plotting options
col='b'
lab='median'
colmap='viridis'
mar='+'

sig=signals[0]
for sig in signals:
    for det in detectors:
        filename='results/results_AA_' + type + '_' + sig + '_' + det + '.txt'
        a= np.loadtxt(filename, dtype='f', delimiter=' ')
        dist=np.unique(a[:,1]) 
        N=len(dist)
        title=sig + ' ' + det
        print('Number of distances:' , N)

        try:
            os.stat(sig)
        except:
            os.mkdir(sig)       

        covpbb=np.zeros((N,4),dtype=float)
        medbw=np.zeros((N,4),dtype=float)
        prec=np.zeros((N,4),dtype=float)


# We compute the median and the median absolute deviation (MAD).
# MAD = median ( X_i - median(X_i))
# Another deviation measure is the standard deviation of the median (SDM)
# SDM = sqrt(pi/2) * sigma / sqrt(n)
# where sigma is the standard deviation of X_i
# The median absolute deviation(MAD) is a robust measure of how spread
# out a set of data is. The variance and standard deviation are also measures
# of spread, but they are more affected by extremely high or extremely low values
# and non normality. If your data is normal, the standard deviation is usually
# the best choice for assessing spread. However, if your data is not normal,
# the MAD is one statistic you can use instead.


        for i in range(N):
            ext=np.where((a[:,1] == dist[i]))
            if method == 'median':
                covpbb[i,0]=np.median(a[ext,7])
                MAD=np.median(abs(a[ext,7]-covpbb[i,0]))
                covpbb[i,1]=MAD
            else:
                covpbb[i,0]=np.mean(a[ext,7])
                covpbb[i,1]=np.std(a[ext,7])
        
            covpbb[i,2]=max(0,covpbb[i,0]-covpbb[i,1])
            covpbb[i,3]=min(1,covpbb[i,0]+covpbb[i,1])
        
            if method == 'median':
                medbw[i,0]=np.median(a[ext,8])
            else:
                medbw[i,0]=np.mean(a[ext,8])
        
            if method == 'median':
                prec[i,0]=np.median(a[ext,11])
                MAD=np.median(abs(a[ext,11] - prec[i,0]))
                prec[i,1]=MAD
            else:
                prec[i,0]=np.mean(a[ext,11])
                prec[i,1]=np.std(a[ext,11])
                
            prec[i,2]=max(0,prec[i,0]-prec[i,1])
            prec[i,3]=prec[i,0]+prec[i,1]
            
            sample=np.where((a[:,1] == dist[i]))
            sample_size=len(sample[0])
            ind=2

        # plt.figure()
        # plt.plot(dist, covpbb[:,0], col, label=lab)
        # plt.fill_between(dist, covpbb[:,2], covpbb[:,3], alpha=0.5, facecolor=col)
        # plt.xlabel('Distance [kpc]')
        # plt.ylabel('Coverage probability')
        # plt.ylim([0, 1.05])
        # plt.title(title)
        # plt.grid(True)
        # plt.legend(loc='best')
        # fig_name=sig + '/' + sig + '_covpbb_' + det + '.png';
        # plt.savefig(fig_name)
        
        # plt.figure()
        # plt.scatter(medbw[:,0], covpbb[:,0],c=dist,marker=mar,s=50)
        # cbar=plt.colorbar()
        # cbar.set_label("Distance (kpc)")
        
        # plt.xlim([0.00015, 0.0003])
        # plt.xlabel('Estimate bandwidth median')
        # plt.ylabel('Coverage probability')
        # plt.title(title)
        # plt.grid(True)
        # plt.legend(loc='lower right')
        # fig_name=sig + '/' + sig + '_covpbb_bandwidth_median_' + det + '.png';
        # plt.savefig(fig_name)
        
        
        # plt.figure()
        # plt.plot(dist, prec[:,0], col, label=lab)
        # plt.fill_between(dist, prec[:,2], prec[:,3], alpha=0.5, facecolor=col)
        # plt.xlabel('Distance [kpc]')
        # plt.ylabel('precision')
        # plt.ylim([0, 1.05])
        # plt.title(title)
        # plt.grid(True)
        # plt.legend(loc='best')
        # fig_name=sig + '/' + sig + '_precision_' + det + '.png';
        # plt.savefig(fig_name)
        
        plt.figure(figsize=(4, 3))
        #plt.figure()
        plt.plot(dist, covpbb[:,0], 'b', label='Coverage probability')
        plt.fill_between(dist, covpbb[:,2], covpbb[:,3], alpha=0.5, facecolor='b')
        
        plt.plot(dist, prec[:,0], 'r', label='Precision')
        plt.fill_between(dist, prec[:,2], prec[:,3], alpha=0.5, facecolor='r')
        
        plt.xlabel('Distance [kpc]')
        plt.ylabel('Fraction')
        plt.ylim([0, 1.05])
        plt.title(title)
        plt.grid(True)
        plt.legend(loc='best')
        fig_name=sig + '/' + sig + '_covpbb_prec_' + det + '_b200.png';
        plt.savefig(fig_name)

    plt.show()




