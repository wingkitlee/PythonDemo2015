# Chebyshev differentation of a smooth function

import numpy as np
import matplotlib.pyplot as plt
from cheb  import cheb

a = 0.0
b = 2.0

def f(x):
	return np.cos(2.0*np.pi*x)-1.0

def df(x):
	return -2.0*np.pi*np.sin(2.0*np.pi*x)

xx = np.linspace(a,b,100)
uu = f(xx)

fig = plt.figure()

for N in [10, 20]:
	D,x = cheb(N,a,b)
	u = f(x)
	Du= df(x)
	Du_cheb = np.dot(D,u)
	
	fig.add_subplot(2,2,2*(N==20)+1)
	plt.plot(x,Du,'.',markersize=8); plt.grid()
	plt.plot(xx,df(xx))
	plt.xlim(a,b); plt.ylim(-7,7); plt.title("u'(x),  N="+str(N))
	error = Du_cheb - Du

	fig.add_subplot(2,2,2*(N==20)+2)
	plt.plot(x,error,'o',markersize=8); plt.grid()
	plt.plot(x,error)
	plt.xlim(a,b); plt.title("    error in u'(x),  N="+str(N))

plt.show()
