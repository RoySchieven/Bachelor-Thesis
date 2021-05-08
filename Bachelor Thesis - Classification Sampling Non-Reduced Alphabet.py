"""
This code samples classifications of two-dimensional fractal percolation for a non-reduced alphabet A_(2,k). 
"""

import numpy as np
from scipy.stats import bernoulli as ber
import time
import matplotlib.pyplot as plt
from scipy.stats import chi2

## Sampling
# First run the functions defined below!
# For dropping confidence: enable corresponding lines in 'tausample'

p=0.78
n=100
N=10**5 
alph=10**-4

# Select starting alphabet (only A_(2,k) possible):
bet=np.array([[[0],[0],[0],[0]],[[1],[1],[1],[1]]]) 
#bet=np.array([[[0,0],[0,0],[0,0],[0,0]],[[1,1],[1,1],[1,1],[1,1]]]) 

start_time = time.time()
e,r,dat,end=tausample(n,p,N,bet,alph)
print("--- %s seconds ---" % (time.time() - start_time))

print(p,e[0],r)

st=np.arange(0,end+1)
for i in range (0,2):
    plt.plot(st, dat[i][0:end+1])
plt.show()

## Multiple Samples
# Running the classification sampling several times to check consistency

M=2
p=0.7
n=10
N=10**5
alph=10**-4
trials=5

bet=np.array([[[0,0],[0,0],[0,0],[0,0]],[[1,1],[1,1],[1,1],[1,1]]])
#bet=np.array([[[0],[0],[0],[0]],[[1],[1],[1],[1]]]) 

for i in range(0,trials):
    start_time = time.time()
    e,r,dat,end=tausample(n,p,N,bet,alph)
    print(M,p,e[0],r)
    print("--- %s seconds ---" % (time.time() - start_time))
    for j in range (0,2):
        st=np.arange(0,end+1)
        plt.plot(st, dat[j][0:end+1],label=(i,j))
    print(i)
    
plt.legend()
plt.show()

## Functions

# Function which samples classifications from starting distribution recursively:
def tausample(n,p,N,alphabet,alpha):
    conflevel=1
    tau=np.array([0,1])
    steps=np.zeros(n+1)
    data=np.array([steps,steps])
    for i in range (0,2):
        data[i,0]=tau[i]
    for k in range (0,n):
        x=p*tau
        x[0]+=(1-p)
        taunew,alphabet=draw(alphabet,x,N)
        tau=multiconf(taunew,N,1-alpha)
        #nonnorm=taunew/N
        #tau=nonnorm/np.sum(nonnorm)        
        print(k,tau[0],tau[1])
        print(len(alphabet))
        for i in range (0,2):
            data[i,k+1]=tau[i]
        if tau[0]>1-0.556:
            print('Break after step', k+1)
            conflevel=(1-alpha)**(k+1)
            return(tau,conflevel,data,k+1)
        conflevel=(1-alpha)**(k+1)
        print(conflevel)
    return(tau,conflevel,data,n)

# Draw a sample of classifications from distribution x:     
def draw(alphabet,x,N):
    count=np.zeros(len(alphabet),int)
    alphanew=np.ndarray.tolist(alphabet)
    for i in range (0,N):
        numbers=len(alphabet)  
        result=np.array([alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)],alphabet[np.random.choice(numbers,p=x)]])
        X=glueing(result[0],result[1],-1,0)
        X=glueing(X,result[2],0,-1)
        X=glueing(X,result[3],len(alphabet[0,1]),len(alphabet[0,0]))
        X=reduce(X)
        t=np.ndarray.tolist(X)
        if t in alphanew:
            count[alphanew.index(t)]+=1
        else:
            alphanew.append(t)
            count=np.append(count,1)
    return(count,np.array(alphanew))

# Two functions for glueing letters together into a word:
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

# Reduce a word back into a letter using weak connectivity:
def reduce(X):
    redto=int(len(X[0])/2)
    Y=np.array([np.arange(1,redto+1),np.arange(redto+1,2*redto+1),np.arange(2*redto+1,3*redto+1),np.arange(3*redto+1,4*redto+1)])
    for i in range (0,4):
        for j in range (0,redto):
            flip=False
            options=[X[i,2*j],X[i,2*j+1]]
            options=[-1 if s==0 else s for s in options]
            for k in range (i+1,4):
                for l in range (0,redto):
                    if X[k,2*l] in options or X[k,2*l+1] in options:
                        mini=min(Y[i,j],Y[k,l])
                        maxi=max(Y[i,j],Y[k,l])
                        Y=np.where(Y==maxi, mini, Y)
                        flip=True
            for l in range (j+1,redto):
                if X[i,2*l] in options or X[i,2*l+1] in options:
                    mini=min(Y[i,j],Y[i,l])
                    maxi=max(Y[i,j],Y[i,l])
                    Y=np.where(Y==maxi, mini, Y)    
                    flip=True
            if flip==False and Y[i,j]==i*redto+j+1:
                Y[i,j]=0
    return(Y)

# Determine a lower bound of significance alpha based on the sample theta:
def multiconf(theta,N,alpha):
    lowerbound=np.zeros(len(theta))
    a=chi2.ppf(alpha,len(theta)-1)
    y0=(N-theta[0])/N
    z0=4*theta[0]*y0
    x0=(a+2*theta[0]-np.sqrt(a*(a+z0)))/(2*(N+a))
    lowerbound[0]=x0
    for i in range (2,len(lowerbound)):
        y=(N-theta[i])/N
        z=4*theta[i]*y
        x=(a+2*theta[i]-np.sqrt(a*(a+z)))/(2*(N+a))
        lowerbound[i]=x
    s=sum(lowerbound)
    lowerbound[1]=max(1-s,0)
    return(lowerbound)




