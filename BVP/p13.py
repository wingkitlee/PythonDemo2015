# p13.py - solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0
#          Python/NumPy translation 12/23/12
from numpy import *
from pylab import *
from cheb  import cheb
N = 16
D,x = cheb(N)        
D2 = dot(D,D)              
D2 = D2[1:-1,1:-1]                 # boundary conditions
f = exp(4*x[1:-1])           
u = linalg.solve(D2,f)             # Poisson eq. solved here
u = hstack(([0.],u,[0.]))
plot(x,u,'.',markersize=8)
xx = arange(-1,1.01,.01)
uu = polyval(polyfit(x,u,N),xx)    # interpolate grid data
plot(xx,uu)
grid(); xlim(-1,1)
exact = ( exp(4*xx) - sinh(4)*xx - cosh(4) )/16 
title('max err = '+str('%.3e' % norm(uu-exact,inf)),fontsize=12)
show()
