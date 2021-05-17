import numpy as np
import matplotlib.pyplot as plt
from Subsystem_design import fuel_required as H2calc
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import meshgrid,cm


#Approximation of the wing podded tank sizes for various volumes and tank radii

def l_func(x,y):
    return (x-(4/3)*np.pi*y**3)/(np.pi*y*y)+2*y

VH2 = H2calc.VH2_calc()/2000
r = np.arange(0.5,1.6,0.01)

X,Y = meshgrid(VH2,r)
print(X, Y)
L = l_func(X,Y)
print(L)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, L, rstride=1, cstride=1,
                      cmap=cm.RdBu,linewidth=0, antialiased=False)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_xlabel('H2 volume per tank [m3]')
ax.set_ylabel('Radius of the tank [m]')
ax.set_zlabel('Length of the tank [m]')
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

plt.plot(r, l_func(max(H2calc.VH2_calc()/2000),r))
plt.show()
