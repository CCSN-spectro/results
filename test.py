import matplotlib.pyplot as plt
import numpy as np

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

    return data,box_name,q1,q2,q3    


# Fixing random state for reproducibility
np.random.seed(19680801)

# fake up some data
spread = np.random.rand(50) * 100
center = np.ones(25) * 50
flier_high = np.random.rand(10) * 100 + 100
flier_low = np.random.rand(10) * -100
data = np.concatenate((spread, center, flier_high, flier_low))
#print(data)
q1=(np.percentile(data, 25))
q2=(np.percentile(data, 50))
q3=(np.percentile(data, 75))

IQR=q3-q1
print(q1-1.5*IQR)
print(q3+1.5*IQR)
print(q1)
print(q3)
print(IQR)

fig1, ax1 = plt.subplots()
ax1.set_title('Basic Plot')
dict=ax1.boxplot(data)
plt.plot(np.ones(95),data,marker="o")
print('boxes')
for item in dict['boxes']:
    print(item.get_ydata())
print('whiskers')
for item in dict['whiskers']:
    print(item.get_ydata())
plt.grid()
plt.show()

