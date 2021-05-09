"""
This code takes the results of sampling classifications without confidence intervals and determines its bootstrap p-value and plots a histogram.
"""

import numpy as np
import matplotlib.pyplot as plt

##
l=[]
u=np.array([0.46673,0.46259,0.5651,0.49448,0.47741,0.46403,0.40618,0.4783,0.55168,0.50543,0.50115,0.44819,0.54658,0.51449,0.45118,0.54519])
for i in range (0,10000000):
    x=np.random.choice(u,size=len(u))
    l.append(np.mean(x)<=1-0.556)
print(np.mean(l))

##

u=np.array([0.46673,0.46259,0.5651,0.49448,0.47741,0.46403,0.40618,0.4783,0.55168,0.50543,0.50115,0.44819,0.54658,0.51449,0.45118,0.54519])
b=np.linspace(0,1,30)

plt.hist(u, bins=b)
plt.axvline(x=0.444, ymin=0, ymax=6,lw=1.5, color='r', linestyle='--')
plt.xlabel(r'$Q^{54}_{min}(0.885)$')
plt.ylabel('Amount')

plt.show()
