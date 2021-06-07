#Victor and valiant
import numpy as np
import matplotlib.pyplot as plt

r = np.arange(1,2.6,0.1)
A = np.pi*r*r
p_cd0 = 0.34
S = 122.6
cd0 = p_cd0*A/S


plt.plot(r,cd0)
plt.show()

V = np.arange(200,1000,1)

rho = 0.4

D = cd0*0.5*rho*V*V

plt.plot(V,D)
plt.show()