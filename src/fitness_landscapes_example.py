## From https://www.kofler.or.at/bioinformatic/python-fitness-surface-quantitativ-trait/
#!/usr/bin/env python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
 
fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, 1, 0.01)
Y = np.arange(0, 1, 0.01)
X, Y = np.meshgrid(X, Y)
 
 
def getfitness(x,y):
    fitnesshash={0:1.0, 1:1.1,2:1.2,3:1.1,4:1.0}
    xgt={0:x**2, 1:2*x*(1-x),2:(1-x)**2}
    ygt={0:y**2, 1:2*y*(1-y),2:(1-y)**2}
    z=0.0
    for k1,fx in xgt.items():
        for k2,fy in ygt.items():
            k=k1+k2
            f=fx*fy
            z+=f*fitnesshash[k]
    return z
 
 
Z = np.array([getfitness(x,y) for (x,y) in zip(X.ravel(), Y.ravel())]).reshape(X.shape)
#Z = np.sqrt(X**2 + Y**2)
 
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.set_zlim(1.0, 1.2)
 
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
 
fig.colorbar(surf, shrink=0.5, aspect=5)
 
plt.show()