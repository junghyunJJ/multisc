#!/usr/bin/python

import sys
from pylmm.input import plink
from pylmm import lmm

import numpy as np

sys.path = ['./src'] + sys.path
from mvLMM import mvLMM


#Y = np.loadtxt("../GAMMA/testData/Y.txt")[:2].T
#K = np.loadtxt("../GAMMA/testData/K.txt")

# geno = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt")
Y = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt")
Y = Y[:,[0,5]]
K = np.loadtxt("/common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.cXX.txt")

X = np.loadtxt("../rs8275764").reshape(-1,1)


# Create the mvLMM objectq
M = mvLMM(Y, K, X0 = X, norm=True)
#M = mvLMM(Y, K, norm=True)
# Perform the opitimization
R = M.getMax()

# The important results are then stored in M.mxCor
# Genetic correlation = M.mxCor[0]
# Environmental correlation = M.mxCor[1]
gcor = M.mxCor[0]
# In [2]: gcor
# Out[2]: -0.26000000000000001

ecor = M.mxCor[1]
# In [3]: ecor
# Out[3]: 0.55000000000000004

Psi,Phi = M.getParameterMatrices(gcor, ecor)
# In [4]: Psi
# Out[4]: 
# array([[ 2.25105485, -0.26768148],
#        [-0.26768148,  0.47087311]])

# In [5]: Phi
# Out[5]: 
# array([[ 0.29493425,  0.24340108],
#        [ 0.24340108,  0.66404027]])

# X = np.loadtxt("../GAMMA/testData/X.txt")
# xx = X[:1].T

# ###########################################################
# ### summary ###############################################
# ###########################################################
# def normPhenos(Y):
#     M = Y.shape[1]
#     for i in range(M): Y[:,i] = (Y[:,i] - Y[:,i].mean())/np.sqrt(Y[:,i].var())
#     return Y

# def cleanPhenos(Y):
#     M = Y.shape[1] 
#     for i in range(M):
#         y = Y[:,i]
#         x = True ^ np.isnan(y)
#         if sum(x) == len(y): continue
#         m = y[x].mean()
#         y[np.isnan(y)] = m
#         Y[:,i] = y
#     return Y

# # K
# # Kva,Kve = linalg.eigh(K)

# # leftTransform
# # self.leftTransform = self.Kve.T


# # np.hstack([M.X0_T, np.dot(M.leftTransform,X)])

from scipy import linalg
Kva,Kve = linalg.eigh(K)
leftTransform = Kve.T

N = K.shape[0]
X0 = np.ones((N,1))
X0_T = np.dot(leftTransform, X0)
X_T = np.dot(leftTransform,X)
new_X = np.hstack([X0_T, X_T])

# Y = cleanPhenos(Y)
# Y = normPhenos(Y)
      
# Ystar = np.dot(leftTransform, Y)


###########################################################
### get beta ##############################################
###########################################################

from scipy import linalg
# Check that they are positive semi-def
Psi_Kva,Psi_Kve = linalg.eigh(Psi)
Phi_Kva,Phi_Kve = linalg.eigh(Phi)

Psi_Kva[Psi_Kva == 0] = 1e-6
Phi_Kva[Phi_Kva == 0] = 1e-6

# Get the D matrix by diagonalizing Psi and Phi 
R = Psi_Kve*np.sqrt(1.0/Psi_Kva) # Now R %*% R.T is = Psi^{-1}
RR = np.dot(np.dot(R.T,Phi),R)
RR_Kva,RR_Kve = linalg.eigh(RR)
D = RR_Kva
rightTransform = np.dot(RR_Kve.T,R.T)


# Get Transformed Y
Yt = np.dot(M.Ystar,rightTransform.T).T.reshape(-1,1)

X = M.X0_T
# else: X = np.hstack([M.X0_T, np.dot(M.leftTransform,X)])


# beta_T,mu,beta_T_stderr,_REML_part = M._getBetaT(X,rightTransform,D,Yt)
Ap = rightTransform
XStar = X_T
P = []
for i in range(M.M): P += (M.Kva + D[i]).tolist()
P = np.array(P)
L = np.kron(Ap,XStar)


A = L.T * 1.0/(P+1.0)
B = np.dot(A,L)
Bi = linalg.inv(B)
beta = np.dot(np.dot(Bi,A),Yt)

_REML_part = np.log(linalg.det(np.dot(L.T,L))) + np.log(linalg.det(B))

mu = np.dot(L,beta)


###########################################################
### get original ##########################################
###########################################################

# Xt = np.kron(Ap,XStar)
# A = Xt.T * 1.0/(P+1.0)
# B = np.dot(A,Xt)
# Bi = linalg.inv(B)
# np.dot(np.dot(Bi,A),Yt)

###########################################################
### get results using new data ############################
###########################################################

Xt = np.kron(Ap,XStar)
new_x = (Xt.T * 1.0/(P+1.0)).T
new_x = np.dot(Xt.T, linalg.inv(np.diag(P))).T

new_y = (Yt.T * 1.0/(P+1.0)).T
new_y = np.dot(Yt.T, linalg.inv(np.diag(P))).T

t1 = linalg.inv(np.dot(new_x.T, new_x))
t2 = np.dot(new_x.T, new_y)

np.dot(t1, t2)