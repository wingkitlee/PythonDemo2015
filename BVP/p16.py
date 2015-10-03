# p16.py - Poisson eq. on [-1,1]x[-1,1] with u=0 on boundary
#          Python started 12/23-24/12
#          Note: interp2d generates spurious structures near y=-0.75
#          when N is near 24. Not sure why. Have changed N to 32.
from numpy import *
from scipy.interpolate import interp2d
from pylab import *
from cheb  import cheb
from time  import time
# Set up grids and tensor product Laplacian and solve for u:
N = 24 
D,x = cheb(N); y = x
xx,yy = meshgrid(x[1:-1],y[1:-1])
Nm1 = N-1
xx = xx.T.reshape(Nm1*Nm1); yy = yy.T.reshape(Nm1*Nm1) # stretch 2D grids to 1D vectors
                        # Transpose because default storage in numpy arrays is by rows.
f = 10*sin(8*xx*(yy-1))
D2 = dot(D,D); D2 = D2[1:-1,1:-1]; I = eye(Nm1)
L = kron(I,D2) + kron(D2,I)   # Laplacian 
figure(1); spy(L) 
tic=time(); u = linalg.solve(L,f);    # solve problem and watch the clock
print L.shape,'linear solve took',time()-tic,'secs'
# Reshape long 1D results onto 2D grid:
uu = zeros((N+1,N+1)); uu[1:-1,1:-1] = u.reshape((N-1,N-1)).T
xx,yy = meshgrid(x,y)
value = uu[N/4,N/4]
# Interpolate to finer grid and plot:
fg = arange(-1.,1.01,.04); [xxx,yyy] = meshgrid(fg,fg)
f = interp2d(x,y,uu,kind='cubic'); uuu = f(fg,fg) # artifacts generated near y=-0.75 if N=24
#figure(2); imshow( uu,interpolation='nearest')
#figure(3); imshow(uuu,interpolation='nearest')
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(xx,  yy,  uu, cstride=1,rstride=1,cmap='jet')
ax.plot_surface(xxx, yyy, uuu, cstride=1,rstride=1,cmap='jet')
ax.azim = -138; ax.elev = 25
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('u')
ax.text(.4,-.3,-.3,'u(2^-1/2,2^-1/2) ='+str('%.11f' % value) )
plt.show()


