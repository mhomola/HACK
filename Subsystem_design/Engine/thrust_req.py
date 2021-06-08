# In constants:
import numpy as np
from Subsystem_design.common_constants import Constants
from math import *
#In constants
'''
self.points = np.array(['Station 1: Idle',
                        'Station 2: Liftoff / Takeoff / start of climb',
                        'Station 3: LG retraction point',
                        'Station 4: End of climb/start of cruise',
                        'Station 5: End of cruise/start of descent',
                        'Station 6: LG extension point',
                        'Station 7: Touchdown',
                        'Station 8: Full Stop'])

self.v0 = np.array([0.015,5,])  # [-] Mach number
self.h = np.array([1, 1, 50, 3000, 11280, 300, 1])  # [m] altitude
self.T0, self.p0, self.rho0, self.a0 = np.array([]), np.array([]), np.array([]), np.array([])
for i in self.h:
    self.ISA_calculator(h_input=i)
    self.T0 = np.append(self.T0, self.T)
    self.p0 = np.append(self.p0, self.p)
    self.rho0 = np.append(self.rho0, self.rho)

'''
class thrust_req(Constants):
    def __init__(self):
        super().__init__()

    def weight(self): #Calculated using
        self.weight_st_1 = self.MTOW_320hack / (self.W3_W2 * self.W2_W1 * self.W1_Wto)
        self.weight_st_2 = self.MTOW_320hack



if __name__ == '__main__':
    t = thrust_req()
    t.weight()

    print('weight at idle', t.weight_st_1)


'''
self.W = np.array([])
nm = 1.852   #nm to km conversion
ft = 0.3048

climb_distance_fraction = 0.037080636
cruise_distance_fraction = 0.932901707
landing_distance_fraction = 0.030017657
harmonic_range = 3200000  #m
cruise_altitude = 11400 #
climb_time = 20 #mins
#climb_angle = np.degrees(np.arcsin(cruise_altitude / (climb_distance_fraction * harmonic_range )))
ROC = cruise_altitude/(climb_time * 60)

T_req_to_climb = ROC * W[i]




E_CO2 = [200,300,300,100]
t=2039
tb=np.arange(2036,2040,1)
ta=np.arange(2035,2039,1)
a1= 0.067
a2=0.1135
a3=0.152
a4=0.097
a5=0.041
tau_2=313.8
tau_3=79.8
tau_4=18.8
tau_5=1.7
for i in range(len(tb)):
    summ = E_CO2[i] * ( a1 + a2 * tau_2 * ( e ** ((tb-t) / tau_2) - e **((ta-t)/tau_2)) +
                 a3 * tau_3 * ( e ** ((tb-t) / tau_3) - e **((ta-t)/tau_3)) +
                 a4 * tau_4 * (e ** ((tb - t) / tau_4) - e ** ((ta - t) / tau_4)) +
                 a5 * tau_5 * (e ** ((tb - t) / tau_5) - e ** ((ta - t) / tau_5)))


print(summ)
'''


