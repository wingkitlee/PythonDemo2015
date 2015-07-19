from numpy import *
import pylab
import matplotlib.animation as animation

dt = 0.05
t = arange(0.0, 5.0, dt)

fig = pylab.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
#ax.grid()
ax.set_aspect('equal')

circ = pylab.Circle((0.0,0.0), ls='dashed', radius=1, color='k', fill=False)
ax.add_patch(circ)

cl = linspace(0.1, 1.0, 20)
line1, = ax.plot([], [], 'ro', lw=2)
line2, = ax.plot([], [], 'ro', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def getxy(t):
    omega0 = 2.0*pi # omega = 2pi / T
    rl = linspace(1.0, 2, 20)
    
    # angular rate, omega = omega0/r ; v=r*omega = const.
    xl = array([ r*cos(omega0/r*t) for r in rl])
    yl = array([ r*sin(omega0/r*t) for r in rl])
    
    return xl, yl

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return line1, line2, time_text

def animate(i):
    thisx, thisy = getxy(t[i])

    line1.set_data(thisx, thisy)
    line2.set_data(-thisx, -thisy)
    time_text.set_text(time_template%(i*dt))
    return line1, line2, time_text
    
ani = animation.FuncAnimation(fig, animate, arange(1, len(t)),
    interval=200, blit=True, init_func=init)

FFMpegWriter = animation.writers['ffmpeg']
mywriter = FFMpegWriter()

#mywriter = animation.MencoderWriter()
ani.save('winding.mp4', fps=15, writer=mywriter)
#pylab.show()  