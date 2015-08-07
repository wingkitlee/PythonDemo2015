from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from simplependulum2 import SinglePendulum


if __name__=='__main__':
    #------------------------------------------------------------
    # set up initial state and global variables
    pendulum = SinglePendulum([20., 0.0])
    ymin =-pendulum.L
    ymax =-pendulum.L*cos(pendulum.state[0])
    dt = 1./30 # 30 fps
    
    #------------------------------------------------------------
    # set up figure and animation
    plt.rc('font', size=18)
#    fig = plt.figure(figsize=(6,6))
    fig = plt.figure()
    ax = fig.add_subplot(211, aspect='equal', autoscale_on=False,
                        xlim=(-3, 3), ylim=(-2, 1.5))                
    #ax1= fig.add_subplot(122, xlim=(0.0, 1000*dt))
    ax1 = fig.add_subplot(212, autoscale_on=False,\
                    xlim=(0.0,10.0), ylim=(-0.1, pendulum.energy()+0.5))                
    
    ax.grid()
    #ax1.grid()
    
    line, = ax.plot([], [], 'o-', lw=2)
    time_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
#    energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    lineU, = ax1.plot([], [], 'b-', lw=2)
    lineK, = ax1.plot([], [], 'r-', lw=2)
    lineE, = ax1.plot([], [], 'k--', lw=3)
    ax1.legend([lineU, lineK, lineE], [r'$U$',r'$K$',r'$E$'],loc='upper left',ncol=3,fontsize=14)
        
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel(r'energy [J]')
    
    xlist = []
    ylist = []
    tlist = []
    Ulist = []
    Klist = []
    Elist = []
    
    plt.tight_layout()
    
    def init():
        """initialize animation"""
        line.set_data([], [])
        time_text.set_text('')
#        energy_text.set_text('')
    
        lineU.set_data([],[])
        lineK.set_data([],[])
        lineE.set_data([],[])
        return line, time_text, lineU, lineK, lineE
    
    def animate(i):
        """perform animation step"""
        global pendulum, dt
        pendulum.step(dt)
        
        E = pendulum.energy()
        U = pendulum.potentialenergy()
        K = pendulum.kineticenergy()
        
        line.set_data(*pendulum.position())
        time_text.set_text('time = %.1f s' % pendulum.time_elapsed)
#        energy_text.set_text('energy = %.3f J' % E)
        
        x1 = pendulum.position()[0][1]
        y1 = pendulum.position()[1][1]
        
        xlist.append(x1)
        ylist.append(y1)
        tlist.append(pendulum.time_elapsed)
        Ulist.append(U)
        Klist.append(K)
        E = K+U
        Elist.append(E)
    
        lineU.set_data(tlist, Ulist)
        lineK.set_data(tlist, Klist)
        lineE.set_data(tlist, Elist)
    
        return line, time_text, lineU, lineK, lineE
    
    # choose the interval based on dt and the time to animate one step
    from time import time
    t0 = time()
    animate(0)
    t1 = time()
    interval = 1000 * dt - (t1 - t0)
    
    ani = animation.FuncAnimation(fig, animate, frames=300,
                                interval=interval, blit=True, init_func=init)
    
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #ani.save('double_pendulum.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# Windows version. Need to make sure ffmpeg.exe is in the path.    
    FFwriter = animation.FFMpegWriter(fps=30)
    ani.save('single_pendulum_energy.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)

#    plt.tight_layout()
#    plt.show()