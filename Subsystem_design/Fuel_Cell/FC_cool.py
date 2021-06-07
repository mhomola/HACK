import numpy as np
import matplotlib.pyplot as plt

n1 = 0.4
n2 = 0.5
n3 = 0.6
n4 = 0.7
P = 400000

ch2o = 4184
dT = np.arange(10,80,0.01)

m1 = ((1-n1)*P)/(n1*ch2o*dT)
m2 = ((1-n2)*P)/(n2*ch2o*dT)
m3 = ((1-n3)*P)/(n3*ch2o*dT)
m4 = ((1-n4)*P)/(n4*ch2o*dT)

plt.plot(dT, m1, '#8AA1CF', label = 'Efficiency 40%')
plt.plot(dT, m2, '#6f83e3', label = 'Efficiency 50%')
plt.plot(dT, m3, '#3c4a8c', label = 'Efficiency 60%')
plt.plot(dT, m4, '#1d2759', label = 'Efficiency 70%')
plt.plot(60, ((1-n3)*P)/(n3*ch2o*60), 'ro', label = 'Design point')
plt.legend()
plt.title('Mass flow vs delta T')
plt.ylabel("Mass Flow [kg/s]")
plt.xlabel("delta Temperature [K]")
plt.show()

