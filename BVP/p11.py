# p11.py - Chebyshev differentation of a smooth function
#          Python/NumPy translation by JR 12/22/12
from numpy import *
from pylab import *
from cheb  import cheb
xx = arange(-1,1.01,.01); uu = exp(xx)*sin(5*xx)
for N in [10, 20]:
	D,x = cheb(N); u = exp(x)*sin(5*x)
	subplot(2,2,2*(N==20)+1)
	plot(x,u,'.',markersize=8); grid()
	plot(xx,uu)
	xlim(-1,1); ylim(-4,2); title('u(x),  N='+str(N))
	error = dot(D,u) - exp(x)*(sin(5*x)+5*cos(5*x))
	subplot(2,2,2*(N==20)+2)
	plot(x,error,'.',markersize=8); grid()
	plot(x,error)
	xlim(-1,1); title("    error in u'(x),  N="+str(N))
show()
