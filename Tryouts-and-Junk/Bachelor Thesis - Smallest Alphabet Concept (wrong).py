import numpy as np
from scipy.stats import bernoulli as ber
import time

## Sampling
M=2
p=0.75
n=7
N=100000


start_time = time.time()
print(tausample(n,p,N,alphabet))
print("--- %s seconds ---" % (time.time() - start_time))


## Functions

alphabet=np.array([[[0],[0],[0],[0]],[[1],[1],[0],[0]],[[0],[2],[2],[0]],[[0],[0],[3],[3]],[[1],[0],[0],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[1],[1],[0]],[[0],[2],[2],[2]],[[1],[0],[1],[1]],[[1],[1],[0],[1]],[[1],[1],[1],[1]]])

classes=np.array([[[1],[2],[3],[4]],[[1],[1],[3],[4]],[[1],[2],[2],[4]],[[1],[2],[3],[3]],[[1],[2],[3],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[1],[1],[4]],[[1],[2],[2],[2]],[[1],[2],[1],[1]],[[1],[1],[3],[1]],[[1],[1],[1],[1]]])

def tausample (n,p,N,alphabet):
    tau=np.zeros(len(alphabet))
    tau[len(alphabet)-1]=1
    taumin=np.zeros(len(alphabet))
    taumin[0]=1
    for k in range (0,n):
        x=p*tau+(1-p)*taumin
        taunew=np.zeros(len(classes))
        for m in range (0,N):
            taunew[draw(alphabet,x)]+=1
        nonnorm=taunew/N
        tau=nonnorm/np.sum(nonnorm)
        print(tau)
    return(tau)
        
        
def draw (alphabet,x): 
    numbers=np.array(len(alphabet))  
    result=np.array([alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)]])    
    X=glueing(result[0],result[1],-1,0)
    X=glueing(X,result[2],0,-1)
    X=glueing(X,result[3],1,1)
    X=reduce(X)
    for i in range (0,len(classes)):
        if np.all(X==classes[i]):
            return(i)

def equaling(X1,X2,glueside,begin,m):
    if X1[glueside][begin]!=0:
        x=X1[glueside][begin]
        y=X2[glueside-2,0]
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
        X1[3][0]=X2[3]
        X=np.array([np.append(X1[0],X2[0]),X1[1],np.append(X1[2],X2[2]),X1[3]])  
    else:
        X1,X2=equaling(X1,X2,2,beginy,m)
        X1,X2=equaling(X1,X2,3,beginx,m)
        X1[2][beginy]=X2[2]
        X1[3][beginx]=X2[3]
        X=X1
    return(X)

def reduce(X):
    Y=np.array([[1],[2],[3],[4]])
    for i in range (0,4):
        for j in range (0,2):
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


