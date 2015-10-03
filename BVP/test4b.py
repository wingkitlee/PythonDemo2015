import numpy as np
import matplotlib.pyplot as plt
from cheb import cheb
from scipy.linalg import eig 

"""
Solve BVP:

	u''+ xu' + (1/gam)*u = B*lam*u for x=[0,5]
	B = F^2 - (gam-1)/gam *x^2
"""

b=5.0
N=20
gam = 5.0/3.0

D, x = cheb(N,0.0,b)
D2 = np.dot(D,D)

A = D2 + np.dot(np.diag(x),D) + np.eye(N+1)/gam
B=np.eye(N+1)-(gam-1.0)/gam*np.diag(x*x)

A=A[1:-1,1:-1]
B=B[1:-1,1:-1]

AA = 

#print "D2=", D2
#print "A =", A
#print "B =", B

lam,V = eig(A,B)
lam = lam.real
ii = np.argsort(-lam)
lam = lam[ii]; V = V[:,ii]

for l in lam:
	print l



