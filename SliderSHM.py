import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

g = 9.8 #m/s^2

plt.rc('font', size=18)
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.2, bottom=0.4)
t = np.arange(0.0, 10.0, 0.001)
x0 = 5.0
v0 = 0.0
l0 = 1.0
omega0 = np.sqrt(g/l0)
P0 = 2.0*np.pi/omega0
s = x0*np.cos(omega0*t) + (v0/omega0)*np.sin(omega0*t)
l, = plt.plot(t,s, lw=2, color='red')
period_text = ax.text(0.02, 0.90, r'$T = %.2f\,{\rm s}$'% P0, transform=ax.transAxes)
plt.axis([0, 10, -10, 10])
plt.title(r'Pendulum (small amplitude)')
plt.xlabel(r'time [s]')
plt.ylabel(r'$\theta$ [deg]')
plt.grid()

axcolor = 'white'
axlength = plt.axes([0.2, 0.2, 0.65, 0.03], axisbg=axcolor)
axx0  = plt.axes([0.2, 0.15, 0.65, 0.03], axisbg=axcolor)
axv0  = plt.axes([0.2, 0.1, 0.65, 0.03], axisbg=axcolor)

slength = Slider(axlength, r'$l$ [m]', 0.1, 10.0, valinit=l0)
sx0 = Slider(axx0, r'$\alpha$ [deg]', 0.1, 10.0, valinit=x0)
sv0 = Slider(axv0, r'$\beta$ [deg/s]', 0.1, 100.0, valinit=v0)

def update(val):
    xx0 = sx0.val
    vv0 = sv0.val
    ll  = slength.val
    omega = np.sqrt(g/ll)
    p = 2.0*np.pi/omega
    s = xx0*np.cos(omega*t) + (vv0/omega)*np.sin(omega*t)
    l.set_ydata(s)
    period_text.set_text(r'$T = %.2f\,{\rm s}$'% p)
    fig.canvas.draw_idle()
    
sx0.on_changed(update)
sv0.on_changed(update)
slength.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    sx0.reset()
    sv0.reset()
    slength.reset()
button.on_clicked(reset)

#rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
#radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
#def colorfunc(label):
#    l.set_color(label)
#    fig.canvas.draw_idle()
#radio.on_clicked(colorfunc)

#plt.tight_layout()
plt.show()