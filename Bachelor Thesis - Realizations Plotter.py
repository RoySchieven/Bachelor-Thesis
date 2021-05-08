"""
This code samples a realization K_n(M,p) of the two-dimensional fractal percolation process 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
import time
import matplotlib.cm as cm

## Plotting single realization

M=2
n=7
p=0.88

start_time = time.time()

K=np.array([[1]])

for k in range (1, n+1):
    X=np.zeros((M**(k),M**(k)))
    for i in range (0,M**(k-1)):
        for j in range (0,M**(k-1)):
            if K[i,j]==1:
                for m in range (i*M,(i+1)*M):
                    for l in range (j*M,(j+1)*M):
                        X[m,l]=ber.rvs(p)
    K=X

print("--- %s seconds ---" % (time.time() - start_time))

plot=plt.imshow(K, cmap=cm.binary)
plot.axes.get_xaxis().set_visible(False)
plot.axes.get_yaxis().set_visible(False)
plt.show()   

## Plotting multiple steps

M=2
n=8
p=0.80

rows=2
cols=2

start_time = time.time()

K=np.array([[1]])
plot=[2,4,6,8]
f,axarr=plt.subplots(rows, cols)

c=0

for k in range (1, n+1):
    X=np.zeros((M**(k),M**(k)))
    for i in range (0,M**(k-1)):
        for j in range (0,M**(k-1)):
            for m in range (i*M,(i+1)*M):
                for l in range (j*M,(j+1)*M):
                    X[m,l]=K[i,j]*ber.rvs(p)
    K=X
    if k in plot:
        axarr[int(c/cols),c%cols].imshow(K, cmap=cm.binary) 
        axarr[int(c/cols),c%cols].set_title("n = %i" %k)
        axarr[int(c/cols),c%cols].axes.get_xaxis().set_visible(False)
        axarr[int(c/cols),c%cols].axes.get_yaxis().set_visible(False)
        c+=1

print("--- %s seconds ---" % (time.time() - start_time))

plt.show()   