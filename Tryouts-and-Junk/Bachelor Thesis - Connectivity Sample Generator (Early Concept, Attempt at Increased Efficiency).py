
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
import time

## Sampling

M=2
p=0.78
n=11
N=10

start_time = time.time()
e=steekproef(M,p,n,N)
print("--- %s seconds ---" % (time.time() - start_time))
print(e)


a=(e-0.556)*math.sqrt(N/(e*(1-e)))
d=norm.cdf(a)
print((1-d)*100)

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

def emptyglueing(X1,X2,beginx,beginy):
    if beginy==0:
        X=np.array([X1[0],np.append(X1[1],X2[1]),X2[2],np.append(X1[3],X2[3])])        
    elif beginx==0:
        for i in range (0,len(X2[3])):
            X1[3][i]=0
        X=np.array([np.append(X1[0],X2[0]),X1[1],np.append(X1[2],X2[2]),X1[3]])  
    else:
        for i in range (0,len(X2[3])):
            X1[2][i+beginy]=0
        for i in range (0,len(X2[3])):
            X1[3][i+beginx]=0
        X=X1
    return(X)

def realisation(M,p):   
    X=np.zeros((M,M),int)
    for m in range (0,M):
        for l in range (0,M):
            X[m,l]=ber.rvs(p)
    return(X)

def checker(X):
    path=False
    for i in range (0,3):
        for j in range (0,len(X[0])):
            if X[i][j]!=0:
                for k in range (i+1,4):
                    if X[i][j] in X[k]:
                        path=True
                        break
    return(path)

def simulation(M,p,n):
    if n==2:
        empty=np.array([np.zeros(M,int),np.zeros(M,int),np.zeros(M,int),np.zeros(M,int)])
        if ber.rvs(p) == 1:
            X=equivalence(realisation(M,p))
        else:
            X=empty
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
                X=glueing(X,Y,-1,0) 
            else:
                Y=empty        
                X=emptyglueing(X,Y,-1,0)        
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=equivalence(realisation(M,p))
                X=glueing(X,Y,0,-1) 
            else:
                Y=empty          
                X=emptyglueing(X,Y,0,-1)           
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=equivalence(realisation(M,p))
                    X=glueing(X,Y,i*M,j*M)
                else:
                    Y=empty        
                    X=emptyglueing(X,Y,i*M,j*M)
    else:
        empty=np.array([np.zeros(M**(n-1),int),np.zeros(M**(n-1),int),np.zeros(M**(n-1),int),np.zeros(M**(n-1),int)])
        if ber.rvs(p) == 1:
            X=simulation(M,p,n-1)
            if checker(X) is False:
                X=empty
        else:
            X=empty
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1)
                if checker(Y) is False:
                    Y=empty
                    X=emptyglueing(X,Y,-1,0)
                else:
                    X=glueing(X,Y,-1,0)
            else:
                Y=empty        
                X=emptyglueing(X,Y,-1,0)       
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=simulation(M,p,n-1)
                if checker(Y) is False:
                    Y=empty
                    X=emptyglueing(X,Y,0,-1) 
                else:
                    X=glueing(X,Y,0,-1) 
            else:
                Y=empty
                X=emptyglueing(X,Y,0,-1)                         
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=simulation(M,p,n-1)
                    if checker(Y) is False:
                        Y=empty
                        X=emptyglueing(X,Y,i*M,j*M)
                    else:
                        X=glueing(X,Y,i*M,j*M)    
                else:
                    Y=empty
                    X=emptyglueing(X,Y,i*M,j*M)        
    return(X)

def steekproef(M,p,n,N):
    c=0
    for i in range (0,N):
        if checker(simulation(M,p,n)) is True:
            c+=1
        print(i,c)    
    return(c/N)
        