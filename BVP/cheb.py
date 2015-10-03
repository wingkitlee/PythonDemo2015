import numpy as np

def cheb(N,a=-1.0,b=1.0):
	if N==0:
		D = 0.
		x = 1.0
	else:
		n = np.arange(0,N+1)
		x = np.cos(np.pi*n/float(N)).reshape(N+1,1)
		c = (np.hstack(( [2.], np.ones(N-1), [2.]))*(-1)**n).reshape(N+1,1)
		X = np.tile(x,(1,N+1))
		dX = X - X.T
		D = np.dot(c,1./c.T)/(dX+np.eye(N+1))
		D -= np.diag(np.sum(D.T,axis=0))
		
		x = (b-a)*(x+1.0)*0.5 + a
		D = 2.0/(b-a)*D
	return D, x.reshape(N+1)
