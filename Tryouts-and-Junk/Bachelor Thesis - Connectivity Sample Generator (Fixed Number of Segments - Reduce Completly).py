
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
from scipy.stats import binom_test
import time

## Sampling

M=2
p=0.8
n=11
N=10
fixed=3  # M^fixed = Amount of side segments in final result

b=n-fixed+1

start_time = time.time()
e=steekproef(M,p,n,N,b)
print("--- %s seconds ---" % (time.time() - start_time))
print(e/N)

d=binom_test(e,N,0.556,'less')
print(d)


## Functions

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

def equaling(X1,X2,glueside,begin,m):
    for i in range (0,len(X2[0])):
        if X1[glueside][i+begin]!=0:
            if X2[glueside-2,i]>m:
                x=X1[glueside][i+begin]
                y=X2[glueside-2][i]
                for j in range (0,4):
                    X2[j]=[x if s==y else s for s in X2[j]]
            elif X2[glueside-2][i]!=0 and X1[glueside][i+begin]!=X2[glueside-2][i]:
                x=X1[glueside][i+begin]
                y=X2[glueside-2][i]
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

def steekproef(M,p,n,N,buildup):
    c=0
    for i in range (0,N):
        if checker(simulation(M,p,n,buildup)) is True:
            c+=1
        #print(i,c)
    return(c)

def reduce(X):
    Y=np.array([[1],[2],[3],[4]])
    for i in range (0,4):
        for j in range (0,len(X[0])):
            if X[i,j]!=0:
                x=X[i,j]
                for k in range (i+1,4):
                    if x in X[k]:
                        if Y[k]==k+1:
                             Y[k]=Y[i]
                        else:
                            Y[i]=Y[k]
                        X[k]=[0 if s==x else s for s in X[k]]
                X[i]=[0 if s==x else s for s in X[i]]
    return(Y)
            
def simulation(M,p,n,buildup):
    if n==2:
        empty=np.array([np.zeros(M,int),np.zeros(M,int),np.zeros(M,int),np.zeros(M,int)])
        if ber.rvs(p) == 1:
            X=equivalence(realisation(M,p))
        else:
            X=empty
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
            else:
                Y=empty       
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
            else:
                Y=empty      
            X=glueing(X,Y,0,-1)             
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=equivalence(realisation(M,p))
                else:
                    Y=empty        
                X=glueing(X,Y,i*M,j*M)
        X=reduce(X)
    elif n<buildup:
        empty=np.array([[0],[0],[0],[0]])
        if ber.rvs(p) == 1:
            X=simulation(M,p,n-1,buildup)
        else:
            X=empty
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1,buildup)
            else:
                Y=empty        
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1,buildup)
            else:
                Y=empty     
            X=glueing(X,Y,0,-1)             
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=simulation(M,p,n-1,buildup)
                else:
                    Y=empty       
                X=glueing(X,Y,i,j)
        X=reduce(X)
    else:
        empty=np.array([np.zeros(M**(n-buildup),int),np.zeros(M**(n-buildup),int),np.zeros(M**(n-buildup),int),np.zeros(M**(n-buildup),int)])
        if ber.rvs(p) == 1:
            X=simulation(M,p,n-1,buildup)
        else:
            X=empty
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1,buildup)
            else:
                Y=empty      
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1,buildup)
            else:
                Y=empty     
            X=glueing(X,Y,0,-1) 
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=simulation(M,p,n-1,buildup)
                else:
                    Y=empty      
                X=glueing(X,Y,i*M**(n-buildup),j*M**(n-buildup))         
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



