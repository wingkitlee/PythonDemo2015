import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dt = 0.05
t = np.arange(0.0, 20.0, dt)
w = 1.0+0.05*np.sin(t*2.0*np.pi)

x = np.linspace(-1.0,1.0,100)
y = np.linspace(-1.0,1.0,100)
xx, yy = np.meshgrid(x,y)

fig, ax = plt.subplots()

w2 = w[0]**2
z = np.exp(-(xx**2+yy**2)/w2)
im = ax.imshow(z,extent=[-1,1,-1,1])

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)

def init():	
	w2 = w[0]**2
	z = np.exp(-(xx**2+yy**2)/w2)
	im.set_data(z)
	return im

def update(i):

	w2 = w[i]**2
	z = np.exp(-(xx**2+yy**2)/w2)
	im.set_data(z)

	return im

ani = animation.FuncAnimation(fig, update, np.arange(1,len(t)), interval=100,init_func=init,blit=False)
plt.show()
