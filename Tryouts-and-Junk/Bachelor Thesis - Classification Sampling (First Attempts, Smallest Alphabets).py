
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector
desc = importr("DescTools")

##

import numpy as np
from scipy.stats import bernoulli as ber
import time
import matplotlib.pyplot as plt

##Poging

M=2
p=0.72
n=40
N=10000

start_time = time.time()
e,r,dat,end=tausample(n,p,N,alphabet,0.0001)
print("--- %s seconds ---" % (time.time() - start_time))

print(r)

st=np.arange(end+1)
for i in range (0,14):
    plt.plot(st, dat[i][0:end+1])
plt.show()


## A(2,0)

alphabet=np.array([[[0],[0],[0],[0]],[[1],[1],[0],[0]],[[0],[2],[2],[0]],[[0],[0],[3],[3]],[[1],[0],[0],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[0],[1],[0]],[[0],[2],[0],[2]],[[1],[1],[1],[0]],[[0],[2],[2],[2]],[[1],[0],[1],[1]],[[1],[1],[0],[1]],[[1],[1],[1],[1]]])

classes=np.array([[[1],[2],[3],[4]],[[1],[1],[3],[4]],[[1],[2],[2],[4]],[[1],[2],[3],[3]],[[1],[2],[3],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[2],[1],[4]],[[1],[2],[3],[2]],[[1],[1],[1],[4]],[[1],[2],[2],[2]],[[1],[2],[1],[1]],[[1],[1],[3],[1]],[[1],[1],[1],[1]]])

def tausample(n,p,N,alphabet,alpha):
    conflevel=1
    tau=np.zeros(len(alphabet))
    tau[len(alphabet)-1]=1
    taumin=np.zeros(len(alphabet))
    taumin[0]=1
    steps=np.zeros(n)
    data=np.array([steps,steps,steps,steps,steps,steps,steps,steps,steps,steps,steps,steps,steps,steps])
    for k in range (0,n):
        x=p*tau+(1-p)*taumin
        taunew=np.zeros(len(classes))
        for m in range (0,N):
            dr=draw(alphabet,x)
            taunew[dr]+=1
        print(k,taunew[0],taunew[13])
        tau=confidence(taunew,N,1-alpha)
        #nonnorm=taunew/N
        #tau=nonnorm/np.sum(nonnorm)        
        #print(tau)
        for i in range (0,14):
            data[i,k]=tau[i]
        if tau[0]>1-0.556:
            print('Break at', k)
            conflevel=conflevel*(1-alpha)
            return(tau,conflevel,data,k)
        conflevel=conflevel*(1-alpha)**(len(alphabet)-1)
        print(conflevel)
    return(tau,conflevel,data,n-1)
        
def draw(alphabet,x): 
    numbers=len(alphabet)
    result=np.array([alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)]])   
    X=glueing(result[0],result[1],-1,0)
    X=glueing(X,result[2],0,-1)
    X=glueing(X,result[3],1,1)
    X=reduce(X)
    for i in range (0,len(classes)):
        if np.all(X==classes[i]):
            return(i)
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

def reduce(X):
    Y=np.array([[1],[2],[3],[4]])
    for i in range (0,4):
        for j in range (0,len(X[0])):
            if X[i,j]!=0:
                x=X[i,j]
                for k in range (i+1,4):
                    if x in X[k] and Y[k]!=k+1:
                        if Y[i]!=i+1:
                            Y[Y[i]-1]=Y[k]
                        Y[i]=Y[k]
                        X[k]=[0 if s==x else s for s in X[k]]
                for k in range (i+1,4):        
                    if x in X[k]:
                        Y[k]=Y[i]
                        X[k]=[0 if s==x else s for s in X[k]]
                X[i]=[0 if s==x else s for s in X[i]]
    return(Y)

def confidence(theta,N,alpha):
    lowerbound=np.zeros(len(theta))
    for i in range (0,len(lowerbound)-1):
        z=int(theta[i])
        y=desc.BinomCI(z,int(N),alpha,"l")
        lowerbound[i]=y[1]
    s=sum(lowerbound)
    lowerbound[len(lowerbound)-1]=max(1-s,0)
    return(lowerbound)
    
## Simplest Alphabet

def tausample(n,p,N):
    tau=np.array([0,1])
    alphabet=np.array([0,1])
    taumin=np.array([1,0])
    for k in range (0,n):
        x=p*tau+(1-p)*taumin
        taunew=np.array([0,0])
        for m in range (0,N):
            taunew=taunew+draw(alphabet,x)
        tau=taunew/N
    return(tau)
        
def draw(alphabet,x):        
    result=np.zeros(4,int)        
    for i in range (0,4):
        result[i]=np.random.choice(alphabet,p=x)
    for i in range (0,4):
        if result[i]==1:
            return(np.array([0,1]))
    return(np.array([1,0]))