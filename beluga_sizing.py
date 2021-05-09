'''
This piece of code helps to evaluate the feasibility of a beluga configutation
'''

import numpy as np
import matplotlib.pyplot as plt
from fuel_constants import *
import fuel_required as H2calc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from sympy.solvers import solve
from sympy import Symbol

def l_solver(V,r,n=1):
    return (V-(8/3)*np.pi*r**3)/(n*np.pi*r**2)+2*r

VH2 = H2calc.VH2_calc()/1000
r = np.arange(0.56,1.5,0.01)
X,Y = meshgrid(VH2,r)
L = l_solver(X,Y)

print(np.amax(np.array(L)))
print(np.amin(np.array(L)))

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

plt.plot(r,l_solver(max(H2calc.VH2_calc()/1000),r))
plt.title('One tank')
plt.xlabel('Radius [m]')
plt.ylabel('Tank length [m]')
plt.show()

plt.plot(r,l_solver(max(H2calc.VH2_calc()/1000),r,n=2))
plt.title('Two tanks')
plt.xlabel('Radius [m]')
plt.ylabel('Tank length [m]')
plt.show()

plt.plot(r,l_solver(max(H2calc.VH2_calc()/1000),r,n=3))
plt.title('Three tanks')
plt.xlabel('Radius [m]')
plt.ylabel('Tank length [m]')
plt.show()
