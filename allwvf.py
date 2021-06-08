import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy import stats
import math
import os,sys,getopt
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


def main(argv):
   folder = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:",["ifolder="])
   except getopt.GetoptError:
      print('test.py -i <inputfolder>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfolder>')
         sys.exit()
      elif opt in ("-i", "--ifolder"):
         inputfolder = arg

   outputfile = inputfolder + '/distance.txt'
   print('Input folder is ', inputfolder)
   print('Output file is ', outputfile)

   return inputfolder,outputfile

if __name__ == "__main__":
   inputfolder,outputfile = main(sys.argv[1:])

# 0 dist
# 1 covpbb                
# 2 medbandwidth
# 3 absolute residual (mean)
# 4 MSE (mean)
# 5 precision (mean)

detectors = ["aLIGO", "ADV", "KAGRA", "CE1", "CE2", "ET_B", "ET_C", "ET_D"]
signals = ["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220"]

#signal_names=["s11.2--LS220", "s15.0--LS220", "s15.0--SFHo", "s15.0--GShen", "s20.0--LS220", "s20.0--SFHo", "s25.0--LS220", "s40.0--LS220", "no signal"]

signal_names=["s11", "s15", "s15S", "s15G", "s20", "s20S", "s25", "s40", "no signal"]

fid=open("distance_KAGRA.txt","w+")

for det in detectors:

    dist_nb=60
    sig_nb=8+1
    qq1=np.zeros((dist_nb,sig_nb))
    qq2=np.zeros((dist_nb,sig_nb))
    qq3=np.zeros((dist_nb,sig_nb))
    dd=np.zeros((dist_nb,sig_nb))
    
#    quantity="coverage"
    quantity="delta"
    
    if quantity == "coverage":
        index=1
    else:
        index=5

    # Noise results    
    filename=inputfolder + '/results_AA_prewhiten_f2_noise_' + det + '.txt'
    a= np.loadtxt(filename, dtype='f', delimiter=' ')

    # percentile
    q50=np.percentile(a[:,index], 50)
    q95=np.percentile(a[:,index], 95)
    q5=np.percentile(a[:,index], 5)
    
    qq1[:,sig_nb-1]=q5
    qq2[:,sig_nb-1]=q50
    qq3[:,sig_nb-1]=q95
    dd[:,sig_nb-1]=0
    
    if quantity == "coverage":
        threshold=q95
    else:
        threshold=q5

    ind=0
    for sig in signals:
        filename=inputfolder + '/results_AA_prewhiten_f2_' + sig + '_' + det + '.txt'
        
        a= np.loadtxt(filename, dtype='f', delimiter=' ')
        dist=np.unique(a[:,0]) 
        dist_nb=len(dist)
        
        dist,q1,q2,q3=buildbox(a,index)
        qq1[:,ind]=q1
        qq2[:,ind]=q2
        qq3[:,ind]=q3
        dd[:,ind]=dist

        #    x=find_chardist(dist,q2,threshold,index)
        #    print(signals[ind],x)
        
        sqq2=smooth(qq2,3)
        x=find_chardist(dist,sqq2[:,ind],threshold,index)
        print(signals[ind],x)
        
        ind=ind+1
        fid.write("%d " % x)

    fid.write("\n")    
    # set up the distance for the "no signal" results    
    dd[:,sig_nb-1]=dist
    dist_max=dist[dist_nb-1]
    
    #plt.figure(figsize=(4, 3))
    plt.figure()
    colormap = plt.cm.gist_rainbow
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, sig_nb))))
    lineObjects = plt.plot(dd, smooth(qq2,3))
    plt.plot(dd[:,0:sig_nb-1], qq2[:,0:sig_nb-1], marker="+",linestyle="")
    plt.fill_between(dist, qq1[:,sig_nb-1], qq3[:,sig_nb-1], alpha=0.05, facecolor='b')
    

    #plt.xlim((1,30))

    if quantity == "coverage":
        plt.ylim((0,1.02))
        plt.xlim((1,dist_max))
        plt.xlim((1,50))
        plt.ylabel('Coverage')
        plt.legend(iter(lineObjects), signal_names, loc='lower right')
    else:
        plt.ylim((0.05,1.05*qq2[0,sig_nb-1]))
        plt.xlim((1,dist_max))
        plt.xlim((1,50))
        plt.ylabel('$\Delta$')
        plt.legend(iter(lineObjects), signal_names, loc='upper right')

    plt.xlabel('Distance [kpc]')

    plt.grid(True)

    plt.title(det)
    fig_name = det + '/' + det + '_' + quantity + '_allwvfs.png';
    plt.savefig(fig_name)
    fig_name = det + '/' + det + '_' + quantity + '_allwvfs.eps';
    plt.savefig(fig_name)
    fig_name = det + '/' + det + '_' + quantity + '_allwvfs.pdf';
    plt.savefig(fig_name)

    plt.show()

fid.close()
    

