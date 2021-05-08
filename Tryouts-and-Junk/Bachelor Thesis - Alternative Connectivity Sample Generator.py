"""
This code estimates the connectivity of the two-dimensional fractal percolation process after n steps for certain M and p. The computation time is somewhat higher than the regular code.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
from scipy.stats import binom_test
import time

## Sampling

M=2
p=0.81
n=11
N=300

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

## Functions

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

def realisation(M,p):   
    X=np.zeros((M,M),int)
    for m in range (0,M):
        for l in range (0,M):
            X[m,l]=ber.rvs(p)
    return(X)

def simulation(M,p,n):
    if n==1:
        if ber.rvs(p) == 1:
            X=np.array([[1],[1],[1],[1]])
        else:
            X=np.array([[0],[0],[0],[0]])
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=np.array([[1],[1],[1],[1]])
                X=glueing(X,Y,-1,0)  
            else:
                X=np.array([X[0],np.append(X[1],[0]),[0],np.append(X[3],[0])])       
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=np.array([[1],[1],[1],[1]])
                X=glueing(X,Y,0,-1)   
            else:
                X[3][0]=0
                X=np.array([np.append(X[0],[0]),X[1],np.append(X[2],[0]),X[3]])                     
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=np.array([[1],[1],[1],[1]])
                    X=glueing(X,Y,i,j)
                else:
                    Y=np.array([[0],[0],[0],[0]])
                    X[2][j]=0
                    X[3][i]=0                 
    else:
        if ber.rvs(p) == 1:
            X=simulation(M,p,n-1)
        else:
            X=np.array([np.zeros(M**(n-1),int),np.zeros(M**(n-1),int),np.zeros(M**(n-1),int),np.zeros(M**(n-1),int)])
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
                    
def steekproef(M,p,n,N):
    c=0
    sketch=np.zeros(N)
    for i in range (0,N):
        if checker(simulation(M,p,n)) is True:
            c+=1
        sketch[i]=c/(i+1)
        print(i,c)
    return(c,sketch)