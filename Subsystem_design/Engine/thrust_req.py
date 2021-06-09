# In constants:
import numpy as np
from Subsystem_design.common_constants import Constants
from math import *
from Subsystem_design.aerodynamic_subsys import AerodynamicCharacteristics
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
#mf_fuel = [1.2,5,7,6,3,2]

#durations = 60* [7.5,1,20,217,18,8] #taxi out, takeoff, climb, cruise, descent, taxi in   16290 s total flight
#h = [0,0,0,10400,11600,0]

#v_tas = []

'''
class thrust_req(Constants):
    def __init__(self):
        super().__init__()
        #self.C_D_0_HACK =
        self.durations = [60* 7.5, 60 * 1, 60 * 20, 60 * 217, 60 * 18, 60 * 8]
        self.mf_fuel = [0.09,1.4,1.4,0.52,0.6,0.05]
        self.Cd0_takeoff = 0.078
        self.AR = 10

    def weight(self):  # Calculated using
        self.ROC = self.cruise_altitude / (self.durations[2] )

        fuel_for_phases=[]
        for i in range(len(self.mf_fuel)):
            fuel_for_phases.append(self.mf_fuel[i] * self.durations[i])
        fuel_weight = sum(fuel_for_phases)
        self.weight_st_1 = fuel_weight + self.OEW_320hack + self.payload_320hack
        self.weight_st_12 = self.weight_st_1 - self.mf_fuel[0] * self.durations[0]
        self.weight_st_2 = self.weight_st_12 - self.mf_fuel[1] * self.durations[1]
        self.weight_st_3 = self.weight_st_2 - self.mf_fuel[2] * self.durations[2]
        self.weight_st_4 = self.weight_st_3 - self.mf_fuel[3] * self.durations[3]
        self.weight_st_5 = self.weight_st_4 - self.mf_fuel[4] * self.durations[4]

        self.T_climb = (self.ROC * self.weight_st_2) / 84.9
        self.acc = 0.12
        self.T_Acc = self.weight_st_2 * self.acc

        Cl_takeoff = (self.weight_st_2 * 2) / (1.225 * 84.9 ** 2 * self.S)
        Cd = self.Cd0_takeoff + Cl_takeoff ** 2 / np.pi * self.AR * self.e
        self.T_drag = Cd * 0.5 * 1.225 * 84.9 ** 2 ** 2 * self.S

        self.total_T_climb = self.T_climb + self.T_Acc + self.T_drag



    def weight(self,t):
        w = w0 - (t-t0) * mf_fuel[i]
        return w
'''
class thrust_req(Constants):
    def __init__(self):
        super().__init__()
        self.durations = 60 * np.array([7.5, 1., 20., 217., 18., 8.])
        self.t_array =np.arange(0, np.sum(self.durations)+30, 30)
        self.h = np.array([0., 0., 0., 10400., 11600, 0.])
        self.taxiout_time = 7.5*60
        self.takeoff_time =8.5*60
        self.cruise_start_time = 28.5 *60
        self.cruise_end_time = 245.5*60
        self.land_time = 263.5*60
        self.stop_time = 271.5*60
        self.liftoffv = 84.9
        self.cruisev= 230
        self.landv = 67.9
        self.takeoff_acc = (self.liftoffv-7) / (self.takeoff_time - self.taxiout_time)
        self.cruise_altitude = 11600
        self.ROC = self.cruise_altitude/ (self.cruise_start_time-self.takeoff_time)
        self.climb_acc = (self.cruisev-self.liftoffv) / (self.cruise_start_time-self.takeoff_time)
        self.descent_acc = (self.landv-self.cruisev)/ (self.land_time-self.cruise_end_time)       #negative
        self.land_decc = (0-self.landv) / (self.stop_time-self.land_time)         #negative

        self.mf_fuel = np.array([0.09,1.4,1.4,0.52,0.6,0.05])
        self.fuel_for_phases = sum(self.mf_fuel * self.durations)

        self.w1 = self.fuel_for_phases + self.OEW_320hack + self.payload_320hack
        self.w2 = self.w1 - self.mf_fuel[0] * self.durations[0]
        self.w3 = self.w2 - self.mf_fuel[1] * self.durations[1]
        self.w4 = self.w3 - self.mf_fuel[2] * self.durations[2]
        self.w5 = self.w4 - self.mf_fuel[3] * self.durations[3]
        self.w6 = self.w5 - self.mf_fuel[4] * self.durations[4]
        self.w7 = self.w6 - self.mf_fuel[5] * self.durations[5]
        #self.warray  = np.array([self.w1,self.w2,self.w3,self.w4,self.w5,self.w6,self.w7])




    def velocity(self):
        v = np.zeros(len(self.t_array))
        v[self.t_array<self.taxiout_time] = 7
        v[(self.t_array>=self.taxiout_time) & (self.t_array<self.takeoff_time)] = np.arange(7,self.liftoffv,
                                                                                        self.takeoff_acc*30)
        v[(self.t_array>=self.takeoff_time) & (self.t_array<self.cruise_start_time)] = np.arange(self.liftoffv,self.cruisev,
                                                                                    self.climb_acc*30)
        v[(self.t_array>=self.cruise_start_time) & (self.t_array<self.cruise_end_time)] = self.cruisev
        v[(self.t_array>=self.cruise_end_time) & (self.t_array<self.land_time)] = np.arange(self.cruisev,self.landv,
                                                                                 self.descent_acc*30)
        v[(self.t_array>=self.land_time) & (self.t_array<self.stop_time)] = np.arange(self.landv,0,self.land_decc*30)


        return (self.t_array, v)

    def mass(self):
        m = np.zeros(len(self.t_array))
        print('hELLO')
        m[self.t_array <self.taxiout_time] = np.linspace(self.w1, self.w2 -self.mf_fuel[0] * 30, len(self.t_array[self.t_array <self.taxiout_time]))
        # print(np.linspace(self.w1, self.w2 -self.mf_fuel[0] * 30, len(self.t_array[self.t_array <self.taxiout_time])))
        # print(np.linspace(self.w2,self.w3-self.mf_fuel[1]*30,len(self.t_array[(self.t_array >= self.taxiout_time) & (self.t_array < self.takeoff_time)])))
        # print(np.linspace(self.w3,self.w4 -self.mf_fuel[2] * 30,len(self.t_array[(self.t_array >= self.takeoff_time) & (self.t_array < self.cruise_start_time)])))
        m[(self.t_array >= self.taxiout_time) & (self.t_array < self.takeoff_time)] = np.linspace(self.w2,
                                                                                                  self.w3-self.mf_fuel[1]*30,
                                                                                                  len(self.t_array[(self.t_array >= self.taxiout_time) & (self.t_array < self.takeoff_time)]))
        m[(self.t_array >= self.takeoff_time) & (self.t_array < self.cruise_start_time)] = np.linspace(self.w3,
                                                                                                       self.w4 -self.mf_fuel[2] * 30,
                                                                                                       len(self.t_array[(self.t_array >= self.takeoff_time) & (self.t_array < self.cruise_start_time)]))
        m[(self.t_array >= self.cruise_start_time) & (self.t_array < self.cruise_end_time)] = np.linspace(self.w4,
                                                                                                        self.w5-self.mf_fuel[3] * 30,
                                                                                                        len(self.t_array[(self.t_array >= self.cruise_start_time) & (self.t_array < self.cruise_end_time)]))
        m[(self.t_array >= self.cruise_end_time) & (self.t_array < self.land_time)] = np.linspace(self.w5,
                                                                                                self.w6-self.mf_fuel[4] * 30,
                                                                                                len(self.t_array[(self.t_array >= self.cruise_end_time) & (self.t_array < self.land_time)]))
        m[(self.t_array >= self.land_time) & (self.t_array <= self.stop_time)] = np.linspace(self.w6,
                                                                                          self.w7-self.mf_fuel[5] * 30,
                                                                                          len(self.t_array[(self.t_array >= self.land_time) & (self.t_array <= self.stop_time)]))
        #print(self.t_array[[(self.t_array >= self.land_time) & (self.t_array <= self.stop_time)]])
        return (m)




        #self.ISA_calculator(h_input=h)


   # def drag(self):
        #Cd0 = 0.025409934074969512
       # CL2 =

   # def thrust



'''
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
        self.acc = 0.12
        self.T_Acc = self.weight_st_2 * self.acc


        Cl_takeoff = (self.weight_st_2 * 2) / (1.225*84.9**2 * hack_wing_area)
        Cd = Cd0_takeoff + Cl_takeoff **2 / np.pi * ae.AR * self.e
        self.T_drag = Cd * 0.5 * rho * 84.9**2 ** 2 * hack_wing_area

        self.total_T_climb = self.T_climb + self.T_Acc + self.T_drag

    def thrust(self):
        self.weight()

'''
if __name__ == '__main__':
    t = thrust_req()
    con = Constants()
    ae = AerodynamicCharacteristics()

    #t.velocity()
    t.mass()





    print('TO T', t.mass(), t.velocity())


'''
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

'''


'''
durations = 60 *np.array([7.5, 1, 20, 217, 18, 8])
print(durations)
t_array = np.arange(0, np.sum(durations), 0.5)
print(t_array)

taxi_out = durations[0]
take_off_time = 8.5
v_array = []
for i in t_array:
    while i <= taxi_out:
        v_array.append(10.3)
# while i>taxi_out and i<=take_off_time:
# velocity_array
'''
