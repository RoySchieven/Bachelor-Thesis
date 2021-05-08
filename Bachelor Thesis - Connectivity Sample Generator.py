"""
This code estimates the connectivity of the two-dimensional fractal percolation process after n steps for certain M and p
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
from scipy.stats import binom_test
import time

## Sampling
# First run the functions defined below!
# Be warned! Calculation time can be large. To track process, enable the print comment in the function Steekproef.

M=2
p=0.775
n=9
N=3000  # Sample size

start_time = time.time()
e,sketch=steekproef(M,p,n,N)
print("--- %s seconds ---" % (time.time() - start_time))
print(e/N, M, p, n)

d=binom_test(e,N,0.556,'less')
print(d)

steps=np.arange(N)
plt.plot(steps, sketch)
plt.xlabel('k')
plt.ylabel('p')
plt.show()

## Range
# This part samples for a range of probabilities p (to determine suitable p for deeper analysis)

M=2
pstart=0.5
pend=0.55
K=10  # Amount of different values of p
n=3
N=500
e=np.zeros(K)
d=np.zeros(K)

steps=np.arange(N)
start_time = time.time()
mesh=(pend-pstart)/K
for k in range (0,K):
    p=pstart+mesh*k
    e[k],sketch=steekproef(M,p,n,N)
    print(e[k]/N, M, p, n)
    d[k]=binom_test(e[k],N,0.556,'less')
    print(d[k])
    plt.plot(steps, sketch, label=p)
print("--- %s seconds ---" % (time.time() - start_time))

plt.xlabel('k')
plt.ylabel('p')
plt.legend()
plt.show()


## Functions

# Determine which subsquares on the egde of K are in the same component as (beginx, beginy):
def component(K,beginx,beginy,enterdirection,C,X):
    Sx=beginx
    Sy=beginy
    dir=enterdirection
    d=len(K[0])-1
    while True:
        ban=[]
        if Sx==0:
            X[0,Sy]=C
            ban.append(2)
        if Sy==0:
            X[1,Sx]=C
            ban.append(3)
        if Sx==d:
            X[2,Sy]=C
            ban.append(0)
        if Sy==d:
            X[3,Sx]=C           
            ban.append(1)
        moves=np.array([[Sx+1,Sy],[Sx,Sy+1],[Sx-1,Sy],[Sx,Sy-1]])     
        for i in range (-1,3):
            if (dir+i)%4 not in ban:
                if K[moves[(dir+i)%4,1],moves[(dir+i)%4,0]]==1:
                    Sx=moves[(dir+i)%4,0]
                    Sy=moves[(dir+i)%4,1]
                    dir=(dir+i)%4
                    break
        if Sx==beginx and Sy==beginy:
            if dir==enterdirection or dir==(enterdirection+2)%4:
                break
            moves=np.array([[Sx+1,Sy],[Sx,Sy+1],[Sx-1,Sy],[Sx,Sy-1]])  
            if K[moves[enterdirection,1],moves[enterdirection,0]]==0:
                break
    return(X)
        
# Given K compile a list which keeps track which subsquares on the edge are connected with eachother
def equivalence(K):
    z=len(K[0])-1
    X=np.array([np.zeros(z+1,int),np.zeros(z+1,int),np.zeros(z+1,int),np.zeros(z+1,int)])
    for j in range (0,z+1):
        if K[j,0]==1:
            X[0,j]=-1
        if K[0,j]==1:
            X[1,j]=-1
        if K[j,z]==1:
            X[2,j]=-1
        if K[z,j]==1:
            X[3,j]=-1
    C=1
    for j in range (z,-1,-1):
        if X[0,j]==-1:
            X=component(K,0,j,0,C,X)  
            C+=1          
    for j in range (0,z+1):
        if X[1,j]==-1:
            X=component(K,j,0,1,C,X) 
            C+=1 
    for j in range (0,z+1):
        if X[2,j]==-1:
            X=component(K,z,j,2,C,X) 
            C+=1 
    for j in range (z,-1,-1):
        if X[3,j]==-1:
            X=component(K,j,z,3,C,X) 
            C+=1 
    return(X)

# Given the edge connectivity information X1 and X2, make sure connected components of the two together have the same number:
def equaling(X1,X2,glueside,begin,m):
    for i in range (0,len(X2[glueside-2])):
        if X1[glueside][i+begin]!=0:
            x=X1[glueside][i+begin]
            y=X2[glueside-2][i]
            if y>m:
                for j in range (0,4):
                    X2[j]=[x if s==y else s for s in X2[j]]
            elif y!=0 and x!=y:
                for j in range (0,4):
                    X1[j]=[y if s==x else s for s in X1[j]]  
    return(X1,X2)

# Given the edge connectivity information X1 and X2, glue the two together to obtain edge connectivity information of a bigger matrix:
def glueing(X1,X2,beginx,beginy):
    m=max(map(max,X1))
    for i in range (0,4):
        X2[i]=[s+m if s!=0 else s for s in X2[i]]
    if beginy==0:
        X1,X2=equaling(X1,X2,2,beginy,m)
        X=np.array([X1[0],np.append(X1[1],X2[1]),X2[2],np.append(X1[3],X2[3])])        
    elif beginx==0:
        X1,X2=equaling(X1,X2,3,beginx,m)
        for i in range (0,len(X2[3])):
            X1[3][i]=X2[3][i]
        X=np.array([np.append(X1[0],X2[0]),X1[1],np.append(X1[2],X2[2]),X1[3]])  
    else:
        X1,X2=equaling(X1,X2,2,beginy,m)
        X1,X2=equaling(X1,X2,3,beginx,m)
        for i in range (0,len(X2[3])):
            X1[2][i+beginy]=X2[2,i]
        for i in range (0,len(X2[3])):
            X1[3][i+beginx]=X2[3,i]
        X=X1
    return(X)

# Generate a starting realization K1:
def realisation(M,p):   
    X=np.zeros((M,M),int)
    for m in range (0,M):
        for l in range (0,M):
            X[m,l]=ber.rvs(p)
    return(X)

# Sample the edge connectivity information corresponding to K_n(M,p):
def simulation(M,p,n):
    if n==2:
        if ber.rvs(p) == 1:
            X=equivalence(realisation(M,p))
        else:
            X=np.reshape(np.zeros(4*M,int),(4,M))
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
            else:
                Y=np.reshape(np.zeros(4*M,int),(4,M))        
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
            else:
                Y=np.reshape(np.zeros(4*M,int),(4,M))    
            X=glueing(X,Y,0,-1)             
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=equivalence(realisation(M,p))
                else:
                    Y=np.reshape(np.zeros(4*M,int),(4,M))
                X=glueing(X,Y,i*M,j*M)
    else:
        if ber.rvs(p) == 1:
            X=simulation(M,p,n-1)
        else:
            X=np.reshape(np.zeros(4*M**(n-1),int),(4,M**(n-1)))
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1)
                X=glueing(X,Y,-1,0) 
            else:
                X=np.array([X[0],np.append(X[1],np.zeros(M**(n-1),int)),np.zeros(M**(n-1),int),np.append(X[3],np.zeros(M**(n-1),int))])                 
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1)
                X=glueing(X,Y,0,-1) 
            else:
                for k in range (0,M**(n-1)):
                    X[3][k]=0
                X=np.array([np.append(X[0],np.zeros(M**(n-1),int)),X[1],np.append(X[2],np.zeros(M**(n-1),int)),X[3]])                  
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=simulation(M,p,n-1)
                    X=glueing(X,Y,i*(M**(n-1)),j*(M**(n-1)))
                else:
                    for k in range (0,M**(n-1)):
                        X[2][k+j*(M**(n-1))]=0
                    for k in range (0,M**(n-1)):
                        X[3][k+i*(M**(n-1))]=0     
    return(X)

# Given edge connectivity information X, check whether there are connected sides:
def checker(X):
    path=False
    for i in range (0,3):
        for j in range (0,len(X[0])):
            if X[i][j]!=0:
                for k in range (i+1,4):
                    for m in range (0,len(X[0])):
                        if X[i][j]==X[k][m]:
                            path=True
                            break
    return(path)

# Take a sample of size N of edge connectivity information and count how many times there are connected sides:
def steekproef(M,p,n,N):
    c=0
    sketch=np.zeros(N)
    for i in range (0,N):
        if checker(simulation(M,p,n)) is True:
            c+=1
        sketch[i]=c/(i+1)
        #print(i,c)
    return(c,sketch)