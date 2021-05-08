count=np.zeros(14)

for i in range(0,14):
    xi=alphabet[i].copy()
    for j in range(0,14):
        xj=alphabet[j].copy()
        for k in range(0,14):
            xk=alphabet[k].copy()
            for m in range(0,14):
                xm=alphabet[m].copy()
                X=glueing(xi,xj,-1,0)
                X=glueing(X,xk,0,-1)
                X=glueing(X,xm,1,1)
                X=reduce(X)
                for l in range (0,14):
                    if np.all(X==classes[l]):
                        count[l]+=1
    print(i)
print(count)

##
import numpy as np

alphabet=np.array([[[0],[0],[0],[0]],[[1],[1],[0],[0]],[[0],[2],[2],[0]],[[0],[0],[3],[3]],[[1],[0],[0],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[0],[1],[0]],[[0],[2],[0],[2]],[[1],[1],[1],[0]],[[0],[2],[2],[2]],[[1],[0],[1],[1]],[[1],[1],[0],[1]],[[1],[1],[1],[1]]])

classes=np.array([[[1],[2],[3],[4]],[[1],[1],[3],[4]],[[1],[2],[2],[4]],[[1],[2],[3],[3]],[[1],[2],[3],[1]],[[1],[1],[3],[3]],[[1],[2],[2],[1]],[[1],[2],[1],[4]],[[1],[2],[3],[2]],[[1],[1],[1],[4]],[[1],[2],[2],[2]],[[1],[2],[1],[1]],[[1],[1],[3],[1]],[[1],[1],[1],[1]]])

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
