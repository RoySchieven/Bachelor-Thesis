
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
import time

## Sampling

M=2
p=0.8
n=7
N=1000

bet=np.array([[[0],[0],[0],[0]],[[1],[1],[1],[1]]])
#bet=np.array([[[0,0],[0,0],[0,0],[0,0]],[[1,1],[1,1],[1,1],[1,1]]])
x=np.array([1-p,p])

start_time = time.time()
e,r=drawing(M,p,bet,x,n,N)
print("--- %s seconds ---" % (time.time() - start_time))

print(e,len(r))

## Functions

def drawing(M,p,alphabet,x,n,N):
    count=np.zeros(len(alphabet),int)
    alphanew=np.ndarray.tolist(alphabet)
    for i in range (0,N):
        X=stepdraw(M,p,alphabet,x,n)
        t=np.ndarray.tolist(X)
        if t in alphanew:
            count[alphanew.index(t)]+=1
        else:
            alphanew.append(t)
            count=np.append(count,1)
        if i%100==0:
            print(i)
    return(count,np.array(alphanew))

def stepdraw(M,p,alphabet,x,n):
    if n==1:
        numbers=len(alphabet)  
        X=np.copy(alphabet[np.random.choice(numbers,p=x)])
        for i in range (1,M):
            Y=np.copy(alphabet[np.random.choice(numbers,p=x)])
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            Y=np.copy(alphabet[np.random.choice(numbers,p=x)])
            X=glueing(X,Y,0,-1)
            for i in range (1,M):
                Y=np.copy(alphabet[np.random.choice(numbers,p=x)]) 
                X=glueing(X,Y,i*len(alphabet[1][0]),j*len(alphabet[0][0]))
        X=reduce(X)
    else:
        if ber.rvs(p) == 1:
            X=stepdraw(M,p,alphabet,x,n-1)
        else:
            X=alphabet[0]
        for i in range (1,M):
            if ber.rvs(p) == 1:
                Y=stepdraw(M,p,alphabet,x,n-1)
            else:
                Y=alphabet[0]     
            X=glueing(X,Y,-1,0)            
        for j in range (1,M):
            if ber.rvs(p) == 1:
                Y=stepdraw(M,p,alphabet,x,n-1)
            else:
                Y=alphabet[0]     
            X=glueing(X,Y,0,-1)             
            for i in range (1,M):
                if ber.rvs(p) == 1:
                    Y=stepdraw(M,p,alphabet,x,n-1)
                else:
                    Y=alphabet[0]        
                X=glueing(X,Y,i*len(alphabet[0][1]),j*len(alphabet[0][0]))
        X=reduce(X)
    return(X)

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
        for i in range (0,len(X2[2])):
            X1[2][i+beginy]=X2[2][i]
        for i in range (0,len(X2[3])):
            X1[3][i+beginx]=X2[3][i]
        X=X1
    return(X)

def reduce(X):
    redto1=int(len(X[0])/2)
    redto2=int(len(X[1])/2)
    Y=np.array([np.arange(1,redto1+1),np.arange(redto1+1,redto1+redto2+1),np.arange(redto1+redto2+1,2*redto1+redto2+1),np.arange(2*redto1+redto2+1,2*redto1+2*redto2+1)])
    for i in range (0,4):
        for j in range (0,len(Y[i])):
            flip=False
            options=[X[i][2*j],X[i][2*j+1]]
            options=[-1 if s==0 else s for s in options]
            for k in range (i+1,4):
                for l in range (0,len(Y[k])):
                    if X[k][2*l] in options or X[k][2*l+1] in options:
                        mini=min(Y[i][j],Y[k][l])
                        maxi=max(Y[i][j],Y[k][l])
                        for m in range (0,4):
                           Y[m]=[mini if s==maxi else s for s in Y[m]] 
                        flip=True
            for l in range (j+1,len(Y[k])):
                if X[i,2*l] in options or X[i,2*l+1] in options:
                    mini=min(Y[i][j],Y[i][l])
                    maxi=max(Y[i][j],Y[i][l])
                    for m in range (0,4):
                        Y[m]=[mini if s==maxi else s for s in Y[m]] 
                    flip=True
            if flip==False and Y[i][j]==int(i/2+0.5)*redto1+int(i/2)*redto2+j+1:
                Y[i][j]=0
    return(Y)