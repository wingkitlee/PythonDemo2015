import numpy as np
import matplotlib.pyplot as plt

v0 = 200.0    # km/s
r0 = 1.0      # kpc
omegap = 30.0 # km/s/kpc
Sigma0 = 80.0   # Msun/pc^2
G = 4.302e-3  # (km/s)^2 /Msun *pc

# effective sound speed
c0 = 10.0     # km/s

# circular velocity
r = np.linspace(0.01, 10.0, 200)
vc = v0 * np.sqrt( 1.0 - (r0/r) * np.arctan(r/r0) )

omega = vc / r
kappa2 = v0**2 * ( 2.0*r - r0*np.arctan(r/r0) - r*r0**2/(r0**2 + r**2) )/r**3
kappa = np.sqrt(kappa2)

s = np.log(5.0)/ 6.0**2
Sigma = Sigma0 * np.exp(-s*r**2)

Q = kappa*c0 / (np.pi * G * Sigma) / 1000.0

# plotting
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(311)
plt.plot(r, vc, 'k-')
plt.plot(r, Sigma, 'b-')
ax = fig.add_subplot(312)
plt.plot(r,omega,'k-',lw=2)
plt.plot(r,omega-0.5*kappa,'k--',lw=2)
plt.plot(r,omega+0.5*kappa,'k--',lw=2)
plt.axhline(omegap, c='b', ls='-')
ax.set_ylim(0.0, 150.0)
ax = fig.add_subplot(313)
plt.plot(r, Q, 'k-', lw=2)
plt.axhline(1.0, c='b', ls='--')
plt.show()

