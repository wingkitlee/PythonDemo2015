import numpy as np
import matplotlib.pyplot as plt
from cheb import cheb

"""
% p33.m - solve linear BVP u_xx = exp(4x), u'(-1)=u(1)=0
N = 16
[D,x] = cheb(N) D2 = D^2
D2(N+1,:) = D(N+1,:)
% Neumann condition at x=-1
D2 = D2(2:N+1,2:N+1)
f = exp(4*x(2:N))
u = D2\ f 0]
u = 0 u]
clf, subplot('position', .1 .4 .8 .5])
plot(x,u,'.','markersize',16)
axis( -1 1 -4 0])
xx = -1:.01:1
uu = polyval(polyfit(x,u,N),xx)
line(xx,uu,'linewidth',.8)
grid on
exact = (exp(4*xx) - 4*exp(-4)*(xx-1) - exp(4))/16
title( 'max err = ' num2str(norm(uu-exact,inf))],'fontsize',12)
"""

N=16
D, x = cheb(N); D2=np.dot(D,D)

# replace the last row (x=-1) of D2 by entries in D
D2[-1,:] = D[-1,:] #Neumann condition at x=-1
D2 = D2[1:,1:] #u(x=1) = 0

f=np.hstack(( np.exp(4.0*x[1:-1]), [0.] )) # RHS of f for the last row is zero
u = np.linalg.solve(D2,f)
u = np.hstack(([0],u))

xx = np.linspace(-1.0,1.0,100)
uu = np.polyval(np.polyfit(x,u,N),xx)

exact = (np.exp(4*xx) - 4*np.exp(-4)*(xx-1) - np.exp(4))/16

plt.figure()
plt.plot(x,u,'ko',markersize=8)
plt.plot(xx,uu,'k-')
plt.grid()
plt.xlim(-1,1)
plt.ylim(-4.0,0.0)
plt.title( 'max err = %e ' % np.linalg.norm(uu-exact,np.inf),fontsize=12)
plt.show()


