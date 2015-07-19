"""
Solver for Simple Pendulum


----- original text ------
General Numerical Solver for the 1D Time-Dependent Schrodinger's equation.

adapted from code at http://matplotlib.sourceforge.net/examples/animation/double_pendulum_animated.py

Double pendulum formula translated from the C code at
http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
----- end of original text ------
"""

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

class SinglePendulum:
    """Single Pendulum Class

    init_state is [theta1, omega1] in degrees,
    where theta1, omega1 is the angular position and velocity of the first
    pendulum arm.
    
    state[0] = theta
    state[1] = theta-dot
    """
    def __init__(self,
                 init_state = [120, 0],
                 L1=1.0,  # length of pendulum 1 in m
                 M1=1.0,  # mass of pendulum 1 in kg
                 G=9.8,  # acceleration due to gravity, in m/s^2
                 origin=(0, 0)): 
        self.init_state = np.asarray(init_state, dtype='float')
        self.params = (L1, M1, G)
        self.origin = origin
        self.time_elapsed = 0

        self.state = self.init_state * np.pi / 180.
        self.L = L1
    
    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        (L1, M1, G) = self.params

        x = np.cumsum([self.origin[0],
                       L1 * sin(self.state[0])])
        y = np.cumsum([self.origin[1],
                       -L1 * cos(self.state[0])])
        return (x, y)
        
    def energy(self):
        """compute energy of the current state
        
        x = L1*sin(theta)
        y =-L1*cos(theta)
        vx = L1*cos(theta)*theta_dot
        vy = L1*sin(theta)*theta_dot
        
        
        """
        (L1, M1, G) = self.params
        
        theta = self.state[0]
        theta_dot = self.state[1]
        
        x = L1*sin(theta)
        y =-L1*cos(theta)
        vx= L1*cos(theta)*theta_dot
        vy= L1*sin(theta)*theta_dot
        
        U = G * M1 * (y + L1)
        K = 0.5 * M1 * (vx**2 + vy**2)
        
        return U + K
        
    def potentialenergy(self):
        """compute energy of the current state

        x = L1*sin(theta)
        y =-L1*cos(theta)
        """
        (L1, M1, G) = self.params
        
        theta = self.state[0]
        y =-L1*cos(theta)
        
        U = G * M1 * (y+L1)
#        K = 0.5 * M1 * (vx**2 + vy**2)
        
        return U
        
    def kineticenergy(self):
        """compute energy of the current state
        
        vx = L1*cos(theta)*theta_dot
        vy = L1*sin(theta)*theta_dot
        
        vx^2 + vy^2 = L1^2 * theta_dot^2
        """
        (L1, M1, G) = self.params
        
#        theta = self.state[0]
        theta_dot = self.state[1]
        
#        vx = L1*cos(theta)*theta_dot
#        vy = L1*sin(theta)*theta_dot

        v2 = L1**2 * theta_dot**2        
#        U = G * M1 * (y+L1)
#        K = 0.5 * M1 * (vx**2 + vy**2)
        K = 0.5 * M1 * v2
        return K

    def dstate_dt(self, state, t):
        """compute the derivative of the given state
        
        d/dt (theta) = theta_dot = state[1]
        d/dt (theta_dot) = -g/l*sin(theta)
        
        """
        (M1, L1, G) = self.params

        dydx = np.zeros_like(state)
        dydx[0] = state[1]
        dydx[1] = - G/L1*sin(state[0])
        
        return dydx

    def step(self, dt):
        """execute one time step of length dt and update state"""
        self.state = integrate.odeint(self.dstate_dt, self.state, [0, dt])[1]
        self.time_elapsed += dt
        
class SimplePendulum:
    """Simple Pendulum Class

    init_state is [theta1, omega1, theta2, omega2] in degrees,
    where theta1, omega1 is the angular position and velocity of the first
    pendulum arm, and theta2, omega2 is that of the second pendulum arm
    """
    def __init__(self,
                 init_state = [120, 0, -20, 0],
                 L1=1.0,  # length of pendulum 1 in m
                 L2=1.0,  # length of pendulum 2 in m
                 M1=1.0,  # mass of pendulum 1 in kg
                 M2=1.0,  # mass of pendulum 2 in kg
                 G=9.8,  # acceleration due to gravity, in m/s^2
                 origin=(0, 0)): 
        self.init_state = np.asarray(init_state, dtype='float')
        self.params = (L1, L2, M1, M2, G)
        self.origin = origin
        self.time_elapsed = 0

        self.state = self.init_state * np.pi / 180.
    
    def position(self):
        """compute the current x,y positions of the pendulum arms"""
        (L1, L2, M1, M2, G) = self.params

        x = np.cumsum([self.origin[0],
                       L1 * sin(self.state[0]),
                       L2 * sin(self.state[2])])
        y = np.cumsum([self.origin[1],
                       -L1 * cos(self.state[0]),
                       -L2 * cos(self.state[2])])
        return (x, y)

    def energy(self):
        """compute the energy of the current state"""
        (L1, L2, M1, M2, G) = self.params

        x = np.cumsum([L1 * sin(self.state[0]),
                       L2 * sin(self.state[2])])
        y = np.cumsum([-L1 * cos(self.state[0]),
                       -L2 * cos(self.state[2])])
        vx = np.cumsum([L1 * self.state[1] * cos(self.state[0]),
                        L2 * self.state[3] * cos(self.state[2])])
        vy = np.cumsum([L1 * self.state[1] * sin(self.state[0]),
                        L2 * self.state[3] * sin(self.state[2])])

        U = G * (M1 * y[0] + M2 * y[1])
        K = 0.5 * (M1 * np.dot(vx, vx) + M2 * np.dot(vy, vy))

        return U + K

    def dstate_dt(self, state, t):
        """compute the derivative of the given state"""
        (M1, M2, L1, L2, G) = self.params

        dydx = np.zeros_like(state)
        dydx[0] = state[1]
        dydx[2] = state[3]

        cos_delta = cos(state[2] - state[0])
        sin_delta = sin(state[2] - state[0])

        den1 = (M1 + M2) * L1 - M2 * L1 * cos_delta * cos_delta
        dydx[1] = (M2 * L1 * state[1] * state[1] * sin_delta * cos_delta
                   + M2 * G * sin(state[2]) * cos_delta
                   + M2 * L2 * state[3] * state[3] * sin_delta
                   - (M1 + M2) * G * sin(state[0])) / den1

        den2 = (L2 / L1) * den1
        dydx[3] = (-M2 * L2 * state[3] * state[3] * sin_delta * cos_delta
                   + (M1 + M2) * G * sin(state[0]) * cos_delta
                   - (M1 + M2) * L1 * state[1] * state[1] * sin_delta
                   - (M1 + M2) * G * sin(state[2])) / den2
        
        return dydx

    def step(self, dt):
        """execute one time step of length dt and update state"""
        self.state = integrate.odeint(self.dstate_dt, self.state, [0, dt])[1]
        self.time_elapsed += dt
        
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
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(211, aspect='equal', autoscale_on=False,
                        xlim=(-3, 3), ylim=(-2, 1.5))                
    #ax1= fig.add_subplot(122, xlim=(0.0, 1000*dt))
    ax1 = fig.add_subplot(212, autoscale_on=False,\
                    xlim=(0.0,10.0), ylim=(-0.1, pendulum.energy()+0.5))                
    
    ax.grid()
    #ax1.grid()
    
    line, = ax.plot([], [], 'o-', lw=2)
    time_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
#    energy_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    lineU, = ax1.plot([], [], 'b-', lw=2)
    lineK, = ax1.plot([], [], 'r-', lw=2)
    lineE, = ax1.plot([], [], 'k--', lw=3)
    
    ax.set_xlabel(r'$x [m]$')
    ax.set_ylabel(r'$y [m]$')
    ax1.set_xlabel(r'$t [s]$')
    ax1.set_ylabel(r'energy [J]')
    
    xlist = []
    ylist = []
    tlist = []
    Ulist = []
    Klist = []
    Elist = []
    
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
        time_text.set_text('time = %.1f' % pendulum.time_elapsed)
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
    
    plt.tight_layout()
    plt.show()