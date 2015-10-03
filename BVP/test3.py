import numpy as np
from cheb import cheb
import matplotlib.pyplot as plt

def f(x):
	return np.exp(-x)*np.sin(5.0*x)

a = 0.0
b = 5.0
x = np.linspace(a, b, 100)
y = f(x)

N = 100

DD,xx = cheb(N,a,b)
yy = f(xx)

fig = plt.figure()
ax = fig.add_subplot(121)
plt.plot(x,y,'.',markersize=8)
plt.plot(xx,yy)

ax = fig.add_subplot(122)
error = np.dot(DD,yy) - np.exp(-xx)*(-np.sin(5.0*xx)+5.*np.cos(5.*xx))
plt.plot(xx,error,'.',markersize=8)
plt.plot(xx,error)

plt.show()
