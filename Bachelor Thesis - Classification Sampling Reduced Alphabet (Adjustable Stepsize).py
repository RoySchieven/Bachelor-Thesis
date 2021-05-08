"""
This updated code samples classifications of two-dimensional fractal percolation for a (reduced) alphabet A_(2,k). It is somewhat slower than the code for non-reduced alphabets only, so it is recommended to use this only for reduced alphabets. Adjusting step size is possible.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import bernoulli as ber
import time
from scipy.stats import chi2


## Sampling
# First run the functions defined below!
# For dropping confidence: enable corresponding lines in 'tausample'

M=2
p=0.840
n=1000
N=100000
alph=10**-4
step=5 # Stepsize

# Select starting alphabet (A_(2,k) reduced/non-reduced possible): 
#bet=np.array([[[0,0],[0,0],[0,0,0,0],[0,0,0,0]],[[1,1],[1,1],[1,1,1,1],[1,1,1,1]]])
bet=np.array([[[0],[0],[0,0],[0,0]],[[1],[1],[1,1],[1,1]]])

start_time = time.time()
e,r,dat,end=tausample(n,p,N,bet,alph,step)
print("--- %s seconds ---" % (time.time() - start_time))

print(M,p,e[0],r)

st=np.arange(0,end*step+1,step,int)
const=np.full(end+1,1-0.556)
plt.plot(st, dat[0][0:end+1], label='pmin')
plt.plot(st, dat[1][0:end+1], label='pmax')
plt.plot(st, const, linestyle='--', color='r', label='_nolegend_') 
plt.xticks(st)
plt.xlabel('n')
plt.ylabel('Probability')
plt.legend()
plt.show()

## Functions

# Function which samples classifications from starting distribution recursively:
def tausample(n,p,N,alphabet,alpha,stepsize):
    conflevel=1
    tau=np.array([0,1])
    steps=np.zeros(n+1)
    data=np.array([steps,steps])
    for i in range (0,2):
        data[i,0]=tau[i]
    for k in range (0,n):
        x=p*tau
        x[0]+=(1-p)
        taunew,alphabet=draw2(alphabet,x,N,stepsize)
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
def draw2(alphabet,x,N,stepsize):
    count=np.zeros(len(alphabet),int)
    alphanew=np.ndarray.tolist(alphabet)
    for j in range (0,N):
        X=drawstep(stepsize,x,alphabet)
        t=np.ndarray.tolist(X)
        if t in alphanew:
            count[alphanew.index(t)]+=1
        else:
            alphanew.append(t)
            count=np.append(count,1)     
    return(count,np.array(alphanew))

# Draw a single classification from distribution x after s steps:
def drawstep(s,x,alphabet):
    if s==1:
        numbers=len(alphabet)
        result=np.array([alphabet[np.random.choice(numbers,p=x)],rotate(alphabet[np.random.choice(numbers,p=x)],1),rotate(alphabet[np.random.choice(numbers,p=x)],3),rotate(alphabet[np.random.choice(numbers,p=x)],2)])
        X=glueing(result[0],result[1],-1,0)
        X=glueing(X,result[2],0,-1)
        X=glueing(X,result[3],len(alphabet[0,1]),len(alphabet[0,0]))
        X=reduce2(X)
    else:
        if ber.rvs(p) == 1:
            X=drawstep(s-1,x,alphabet)
        else:
            X=np.copy(alphabet[0])
        if ber.rvs(p) == 1:
            Y=rotate(drawstep(s-1,x,alphabet),1)
        else:
            Y=rotate(alphabet[0],1)     
        X=glueing(X,Y,-1,0)            
        if ber.rvs(p) == 1:
            Y=rotate(drawstep(s-1,x,alphabet),3)
        else:
            Y=rotate(alphabet[0],3)     
        X=glueing(X,Y,0,-1)             
        if ber.rvs(p) == 1:
            Y=rotate(drawstep(s-1,x,alphabet),2)
        else:
            Y=rotate(alphabet[0],2)        
        X=glueing(X,Y,len(alphabet[0,1]),len(alphabet[0,0]))
        X=reduce2(X)
    return(X)

# Two functions for glueing letters together into a word:
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
        X=np.array([np.append(X1[0],X2[0]),X1[1],np.append(X1[2],X2[2]),np.append(X2[3],X1[3][2*len(X2[3]):])])  
    else:
        X1,X2=equaling(X1,X2,2,beginy,m)
        X1,X2=equaling(X1,X2,3,beginx,m)
        X=np.array([X1[0],X1[1],np.append(X1[2][:len(X2[2])],X2[2]),np.append(X1[3][:len(X2[3])],X2[3])],dtype=object)  
    return(X)
 
# Reduce a word back into a letter using weak connectivity:    
def reduce2(X):
    redto=int(len(X[0])/2)
    Y=np.array([list(np.arange(1,redto+1)),list(np.arange(redto+1,2*redto+1)),list(np.arange(2*redto+1,4*redto+1)),list(np.arange(4*redto+1,6*redto+1))])
    for i in range(0,2):
        for j in range (0,redto):
            flip=False
            options=[X[i][2*j],X[i][2*j+1]]
            options=[-1 if s==0 else s for s in options]
            if i==0:
                for l in range (0,redto):
                    if X[1][2*l] in options or X[1][2*l+1] in options:
                        mini=min(Y[0][j],Y[1][l])
                        maxi=max(Y[0][j],Y[1][l])
                        for m in range (0,4):
                            Y[m]=[mini if s==maxi else s for s in Y[m]] 
                            flip=True            
            for k in range (2,4):
                for l in range (0,len(Y[k])):
                    if X[k][l] in options:
                        mini=min(Y[i][j],Y[k][l])
                        maxi=max(Y[i][j],Y[k][l])
                        for m in range (0,4):
                            Y[m]=[mini if s==maxi else s for s in Y[m]] 
                            flip=True  
            for l in range (j+1,redto):
                if X[i,2*l] in options or X[i,2*l+1] in options:
                    mini=min(Y[i][j],Y[i][l])
                    maxi=max(Y[i][j],Y[i][l])
                    for m in range (0,4):
                        Y[m]=[mini if s==maxi else s for s in Y[m]] 
                    flip=True
            if flip==False and Y[i][j]==i*redto+j+1:
                Y[i][j]=0
    for i in range(2,4):
        for j in range (0,2*redto):
            flip=False
            if X[i][j]!=0:
                if i==2:         
                    for l in range (0,len(Y[3])):
                        if X[3][l]==X[2][j]:
                            mini=min(Y[2][j],Y[3][l])
                            maxi=max(Y[2][j],Y[3][l])
                            for m in range (0,4):
                                Y[m]=[mini if s==maxi else s for s in Y[m]] 
                                flip=True  
                for l in range (j+1,2*redto):
                    if X[i][l]==X[i][j]:
                        mini=min(Y[i][j],Y[i][l])
                        maxi=max(Y[i][j],Y[i][l])
                        for m in range (0,4):
                            Y[m]=[mini if s==maxi else s for s in Y[m]] 
                        flip=True
            if flip==False and Y[i][j]==2*redto+2*(i-2)*redto+j+1:
                Y[i][j]=0
    return(Y)

# Rotate letters to orientate them correctly for the process:
def rotate(X,t):
    g=[[],[],[],[]]
    for i in range (0,4):
        g[i]=X[(i-t)%4]
    if t!=3:
        g[1]=g[1][::-1]
        g[3]=g[3][::-1]
    if t!=1:
        g[0]=g[0][::-1]
        g[2]=g[2][::-1]        
    return(np.array(g))

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
        