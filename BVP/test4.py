import numpy as np
import matplotlib.pyplot as plt
from cheb import cheb

"""
Solve BVP:

	u''+ xu' + (1/gam)*u = lam*u for x=[0,5]
"""

N=50
gam = 5.0/3.0

D, x = cheb(N,0.0,5.0)
D2 = np.dot(D,D)

A = np.dot(np.diag(x),D)
B = np.eye(N+1)/gam

D2 = D2[1:-1,1:-1]
A=A[1:-1,1:-1]
B=B[1:-1,1:-1]

print "D2=", D2
print "A =", A
print "B =", B

M = D2 - A - B

lam,V = np.linalg.eig(D2)
ii = np.argsort(-lam)
lam = lam[ii]; V = V[:,ii]

print lam



