# solve eigenvalue BVP u_xx = lambda*u, u(-1)=u(1)=0

import numpy as np
import matplotlib.pyplot as plt
from cheb  import cheb

"""
u_xx = lam*u, u(0) = u(1.0) = 0

u = a cos(kx) + b sin(kx), where lam=-k^2
u(0) = a = 0
u(1) = b sin(k) = 0
=> k = n pi
=> lam = -k^2 = -pi^2 * n^2

wavelength = 2pi/k = 2.0/n

"""
A = 0.0
B = 1.0
N = 100
D,x = cheb(N,A,B)
D2 = np.dot(D,D)

# ignore the first and last row/column <=> setting u[0]=u[N]=0
D2 = D2[1:-1,1:-1] 
lam,V = np.linalg.eig(D2)
ii = np.argsort(-lam)  # sort eigenvalues and -vectors
lam = lam[ii]; V = V[:,ii]

fig = plt.figure()
#for j in range(4,30,5):              # plot 6 eigenvectors
for j in range(0,5):
	u = np.hstack(( [0.], V[:,j], [0.] ))
#	fig.add_subplot(7,1,j/5+1)
	fig.add_subplot(7,1,j+1)
	plt.plot(x,u,'.',markersize=5)
	xx = np.linspace(A, B, 200); uu = np.polyval(np.polyfit(x,u,N),xx)
	plt.plot(xx,uu); plt.xlim(-1,1); plt.ylim(-.5,.5); plt.axis('off')
	plt.text(-.4,.4,r'$\lambda_%s$' % str(j+1)+'='+str('%.13f' % (lam[j]/np.pi**2))+r'$\pi^2$') #typo in p15.m
	plt.text( .7,.4,'%.1f' % ( 2*N/(j+1) )+' ppw')

#savefig('p15.png')
plt.show()

