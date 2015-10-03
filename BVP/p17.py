# p17.py - Helmholtz eq. u_xx + u_yy + (k^2)u = f
#          on [-1,1]x[-1,1]    (compare p16.py)
#          Python translation 12/25/12
from numpy import *
from scipy.interpolate import interp2d
from pylab import *
from cheb  import cheb
# Set up spectral grid and tensor product Helmholtz operator:
N = 24; D,x = cheb(N); y = x
xx,yy = meshgrid(x[1:-1],y[1:-1])
Nm1 = N-1; xx = xx.T.reshape(Nm1*Nm1); yy = yy.T.reshape(Nm1*Nm1) 
f = exp(-10*((yy-1)**2+(xx-.5)**2));
D2 = dot(D,D); D2 = D2[1:-1,1:-1]; I = eye(Nm1)
k = 9.
L = kron(I,D2) + kron(D2,I) + k**2*eye(Nm1*Nm1)  
# Solve for u, reshape to 2D grid, and plot:
u = linalg.solve(L,f)
uu = zeros((N+1,N+1)); uu[1:-1,1:-1] = u.reshape((N-1,N-1)).T
xx,yy = meshgrid(x,y)
fg = arange(-1.,1.001,.0333); [xxx,yyy] = meshgrid(fg,fg)
f = interp2d(x,y,uu,kind='cubic'); uuu = f(fg,fg) 
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xxx, yyy, uuu, cstride=1,rstride=1,cmap='jet')
ax.azim = -138; ax.elev = 25
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('u')
ax.text(.2,1,.022,'u(0,0) ='+str('%.11f' % uu[N/2,N/2]) )
figure(2); 
contour(xxx,yyy,uuu)
plt.show()


