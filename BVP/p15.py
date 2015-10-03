# p15.py - solve eigenvalue BVP u_xx = lambda*u, u(-1)=u(1)=0
#          Python/NumPy translation 12/23/12
from numpy import *
from pylab import *
from cheb  import cheb
N = 36; D,x = cheb(N); D2 = dot(D,D); D2 = D2[1:-1,1:-1] 
lam,V = eig(D2)
ii = argsort(-lam)  # sort eigenvalues and -vectors
lam = lam[ii]; V = V[:,ii]
for j in range(4,30,5):              # plot 6 eigenvectors
	u = hstack(( [0.], V[:,j], [0.] )); subplot(7,1,j/5+1)
	plot(x,u,'.',markersize=5)
	xx = arange(-1,1.01,.01); uu = polyval(polyfit(x,u,N),xx)
	plot(xx,uu); xlim(-1,1); ylim(-.5,.5); axis('off')
	text(-.4,.4,'eig'+str(j+1)+'='+str('%.13f' % (lam[j]*4/pi**2))+'*pi^2/4') #typo in p15.m
	text( .7,.4,str('%.1f' % (4*N/(pi*(j+1))))+' ppw')
savefig('p15.png')
show()

