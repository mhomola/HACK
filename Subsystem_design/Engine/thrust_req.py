# In constants:
import numpy as np
from Subsystem_design.common_constants import Constants
from math import *
'''
#In constants:
self.points = np.array(['Station 1: Idle',
                        'Station 12: Start of TO',
                        'Station 2: Liftoff / End of Takeoff / Start of climb',
                        'Station 3: End of climb / start of cruise',
                        'Station 4: End of cruise/start of descent',
                        'Station 5: Touchdown'])


self.v0 = np.array([0.015,5,])  # [-] Mach number
self.h = np.array([1, 1, 5, 11280, 11280, 300, 1])  # [m] altitude
self.T0, self.p0, self.rho0, self.a0 = np.array([]), np.array([]), np.array([]), np.array([])

for m in self.h:
    self.ISA_calculator(h_input=m)
    self.T0 = np.append(self.T0, self.T)
    self.p0 = np.append(self.p0, self.p)
    self.rho0 = np.append(self.rho0, self.rho)

nm = 1.852   #nm to km conversion
ft = 0.3048

climb_distance_fraction = 0.037080636
cruise_distance_fraction = 0.932901707
landing_distance_fraction = 0.030017657
duration_1_to_12 = 7.5 #min                          taxi out phase
duration_12_to_2 = 1 #min                            take off phase
duration_2_to_4_climb_time = 20 #min                 climb phase
duration_4_to_5_cruise_time = 217 #min
duration_5_to_7_descent_time = 18 #min
duration_7_to_8 = 8 #min includes approach + taxi in
'''

harmonic_range = 3200000  #m

#climb_angle = np.degrees(np.arcsin(cruise_altitude / (climb_distance_fraction * harmonic_range )))
#from sara:
mf_fuel = [1.2,5,7,6,3,2]

durations = 60* [7.5,1,20,217,18,8] #taxi out, takeoff, climb, cruise, descent, taxi in
h = [1,1,1,11280,11280,1]

v_tas = []




class thrust_req(Constants):
    def __init__(self):
        super().__init__()


    def weight(self,t):
        w = w0 - (t-t0) * mf_fuel[i]
        return w

    def velocity

    def drag

    def thrust

while t < total time:
    if t <= taxi_out:
        i = 0
        w0 = ini_hack_weight
        t0 = 0

    if t <= take_off_time and t > taxi out:
     i = 1
     w0 = t.weight(taxi_out)
    t0 = taxi_out



    def weight(self):#Calculated using
        self.ROC = self.cruise_altitude / (self.durations[2] * 60)

        fuel_for_phases = sum(mf_fuel * durations)
        self.weight_st_1 = fuel_for_phases + self.OEW_320hack + self.cargo_320hack
        self.weight_st_12 = self.weight_st_1 - mf_fuel[0] * durations[0]
        self.weight_st_2 = self.weight_st_12 - mf_fuel[1] * durations[1]
        self.weight_st_3 = self.weight_st_2 - mf_fuel[2] * durations[2]
        self.weight_st_4 = self.weight_st_3 - mf_fuel[3] * durations[3]
        self.weight_st_5 = self.weight_st_4 - mf_fuel[4] * durations[4]

        self.T_climb = (self.ROC * self.weight_st_2)/84.9
        self.acc =

    def thrust(self):
        self.weight()


        #self.weight_st_4 = self.weight_st_3 - mf_fuel_climb * duration_3_to_4




if __name__ == '__main__':
    t = thrust_req()
    con = Constants()
    t.weight()


    print('weight at idle', t.weight_st_1)



self.W = np.array([])

T_req_to_climb = ROC * W[i]



#Climate effect integration CO2
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



