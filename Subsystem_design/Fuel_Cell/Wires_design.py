import numpy as np
import matplotlib.pyplot as plt

Cu_res20 = 1.68*10**(-8) #Ohm*m
alpha = 0.00404 #1/T
T0 = 20
Tres_AIR_45 = 5.88#36.075 #m*K/W
Tres_XLPE = 3.5 #m*K/W
Tres_tot = Tres_AIR_45+Tres_XLPE
Tamb = 40

def wire_sizing(I,Tamb):
    return 2*np.sqrt((I*I*Cu_res20*(1+alpha*(T-T0))*Tres_tot)/(np.pi*(T-Tamb)))

T = 50
I = np.arange(100,1501)
d1 = 2*np.sqrt((I*I*Cu_res20*(1+alpha*(T-T0))*Tres_tot)/(np.pi*(T-Tamb)))

I01 = 100
T2 = np.arange(45,100)
d01 = 2*np.sqrt((I01*I01*Cu_res20*(1+alpha*(T2-T0))*Tres_tot)/(np.pi*(T2-Tamb)))

I02 = 300
d02 = 2*np.sqrt((I02*I02*Cu_res20*(1+alpha*(T2-T0))*Tres_tot)/(np.pi*(T2-Tamb)))

I1 = 500
d21 = 2*np.sqrt((I1*I1*Cu_res20*(1+alpha*(T2-T0))*Tres_tot)/(np.pi*(T2-Tamb)))

I2 = 700
d22 = 2*np.sqrt((I2*I2*Cu_res20*(1+alpha*(T2-T0))*Tres_tot)/(np.pi*(T2-Tamb)))

I3 = 900
d23 = 2*np.sqrt((I3*I3*Cu_res20*(1+alpha*(T2-T0))*Tres_tot)/(np.pi*(T2-Tamb)))


print('d2')

plt.plot(I,d1,"k")
plt.ylabel("Diameter [m]")
plt.xlabel("Current [A]")
plt.show()

plt.plot(T2, d01, '#BBCFF5', label = '100 A')
plt.plot(T2, d02, '#8AA1CF', label = '300 A')
plt.plot(T2, d21, '#6f83e3', label = '500 A')
plt.plot(T2, d22, '#3c4a8c', label = '700 A')
plt.plot(T2, d23, '#1d2759', label = '900 A')
plt.legend()
plt.title('Diameter vs Temperature; T = 45 deg')
plt.ylabel("Diameter [m]")
plt.xlabel("Temperature [deg]")
plt.show()