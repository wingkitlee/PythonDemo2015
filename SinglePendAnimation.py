from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from simplependulum2 import SinglePendulum
#------------------------------------------------------------
# set up initial state and global variables
pendulum = SinglePendulum([5., 0.0], L1=2.0)
dt = 1./30 # 30 fps

#------------------------------------------------------------
# set up figure and animation
#plt.rc('text', usetex=True)
plt.rc('font', size=18)
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-2, 2), ylim=(-3, 1))                

#ax.grid()
ax.axis('off')

line, = ax.plot([], [], 'o-', lw=2)
time_text = ax.text(0.02, 0.90, '', transform=ax.transAxes, size=18)
#energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)

ax.set_xlabel(r'$x [m]$')
ax.set_ylabel(r'$y [m]$')
plt.tight_layout()

def init():
    """initialize animation"""
    line.set_data([], [])
    time_text.set_text('')
#    energy_text.set_text('')

    return line, time_text

def animate(i):
    """perform animation step"""
    global pendulum, dt
    pendulum.step(dt)
    
    theta = pendulum.state[0]*180.0/np.pi
    
    line.set_data(*pendulum.position())
    time_text.set_text('time = %.1f s' % pendulum.time_elapsed)
#    energy_text.set_text('angle (deg) = %.2f' % theta)
    
    return line, time_text

# choose the interval based on dt and the time to animate one step
from time import time
t0 = time()
animate(0)
t1 = time()
interval = 1000 * dt - (t1 - t0)
print "interval = ", interval
#interval = 30

ani = animation.FuncAnimation(fig, animate, frames=300,
                              interval=interval, blit=True, init_func=init)

FFwriter = animation.FFMpegWriter(fps=30)
ani.save('single_pendulum.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('double_pendulum.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#plt.tight_layout()
#plt.show()