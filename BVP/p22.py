from cheb import *
from numpy import dot,diag,real,argsort,zeros,linspace,polyval,polyfit,where
from scipy.linalg import eig
from scipy.special import airy
import matplotlib.pyplot as plt

plt.figure(figsize=(10,8))
for N in range(12,60,12):
	print N
	D,x = cheb(N); D2 = dot(D,D); D2 = D2[1:N,1:N]
	Lam,V = eig(D2,diag(x[1:N]))
	Lam = real(Lam); ii = where(Lam>0)[0]
	V = real(V[:,ii]); Lam = Lam[ii]
	ii = argsort(Lam); ii=ii[4]; Lam=Lam[ii]
	v = zeros(N+1); v[1:N] = V[:,ii]; v = v/v[N/2]*airy(0.0)[0]
	xx = linspace(-1.0,1.0,200); vv = polyval(polyfit(x,v,N),xx);
	plt.subplot(2,2,N/12); plt.plot(xx,vv)
	plt.title("N = %d   eig = %15.10f"%(N,Lam));

plt.show()
