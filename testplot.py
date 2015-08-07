import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

x = np.linspace(0.0, 10.0)
y = np.sin(x)

p, = ax.plot(x,y)

ax.legend([p],['test'],loc='best')
plt.show()
