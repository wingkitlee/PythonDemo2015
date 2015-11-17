import numpy as np
import matplotlib.pyplot as plt

v0 = 200.0 # km/s
r0 = 1.0   # kpc
omegap = 30.0 # km/s/kpc

# circular velocity
r = np.linspace(0.01, 10.0, 200)
vc = v0 * np.sqrt( 1.0 - (r0/r) * np.arctan(r/r0) )

omega = vc / r
kappa2 = v0**2 * ( 2.0*r - r0*np.arctan(r/r0) - r*r0**2/(r0**2 + r**2) )/r**3
kappa = np.sqrt(kappa2)




# plotting
fig = plt.figure()
ax = fig.add_subplot(211)
plt.plot(r, vc, 'k-')
ax = fig.add_subplot(212)
plt.plot(r,omega,'k-',lw=2)
plt.plot(r,omega-0.5*kappa,'k--',lw=2)
plt.plot(r,omega+0.5*kappa,'k--',lw=2)
plt.axhline(omegap, c='b', ls='-')
ax.set_ylim(0.0, 150.0)
plt.show()

