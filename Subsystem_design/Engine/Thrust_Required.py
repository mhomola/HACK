# In constants:
import numpy as np
from Subsystem_design.common_constants import Constants
from math import *
from Subsystem_design.aerodynamic_subsys import cd0clean, wingar
import matplotlib.pyplot as plt


class thrust_req(Constants):
    def __init__(self, cd0clean,wingar):
        super().__init__()
        self.durations = 60 * np.array([7.5, 1., 20., 217., 18., 8.])
        self.t_array =np.arange(0, np.sum(self.durations)+30, 30)
        self.h = np.array([0., 0., 0., 11600., 11600, 0.])
        self.taxiout_time = 7.5*60
        self.takeoff_time =8.5*60
        self.mid_climb_time = 18.5 *60
        self.cruise_start_time = 28.5 *60
        self.cruise_end_time = 245.5*60
        self.mid_descent_time = 254.5 *60
        self.land_time = 263.5*60
        self.stop_time = 271.5*60
        self.taxiout_velocity = 13
        self.liftoffv = 84.9
        self.cruisev= 230
        self.landv = 67.9
        self.takeoff_acc = (self.liftoffv-self.taxiout_velocity) / (self.takeoff_time - self.taxiout_time)
        self.cruise_altitude = 11600
        self.ROC = self.cruise_altitude/ (self.cruise_start_time-self.takeoff_time)
        self.ROD = (0-self.cruise_altitude)/(self.land_time-self.cruise_end_time)
        self.climb_acc = (self.cruisev-self.liftoffv) / (self.cruise_start_time-self.takeoff_time)
        self.descent_acc = (self.landv-self.cruisev)/ (self.land_time-self.cruise_end_time)       #negative
        self.land_decc = (0-self.landv) / (self.stop_time-self.land_time)         #negative

        self.mf_fuel = np.array([0.3,1.4,1.4,0.52,0.6,0.3])
        self.fuel_for_phases = sum(self.mf_fuel * self.durations)

        self.w1 = self.fuel_for_phases + self.OEW_320hack + self.payload_320hack
        self.w2 = self.w1 - self.mf_fuel[0] * self.durations[0]
        self.w3 = self.w2 - self.mf_fuel[1] * self.durations[1]
        self.w4 = self.w3 - self.mf_fuel[2] * self.durations[2]
        self.w5 = self.w4 - self.mf_fuel[3] * self.durations[3]
        self.w6 = self.w5 - self.mf_fuel[4] * self.durations[4]
        self.w7 = self.w6 - self.mf_fuel[5] * self.durations[5]
        self.warray = np.array([self.w1,self.w2,self.w3,self.w4,self.w5,self.w6,self.w7])

        self.cd0array = np.concatenate([(np.repeat(cd0clean[0],15)),(np.repeat(cd0clean[1],2)), (np.repeat(cd0clean[2],474)),
                                                                     (np.repeat(cd0clean[4], 36)),(np.repeat(cd0clean[0],17))])  #decided that landing is halfway through descent
        self.acceleration = np.concatenate([(np.repeat(0,15)),(np.repeat(self.takeoff_acc,2)), (np.repeat(self.climb_acc,40)),
                                                                     (np.repeat(0,434)),(np.repeat(self.descent_acc,36)), np.repeat(self.land_decc,17)])
        self.rates_of_climb_and_descent = np.concatenate([(np.repeat(0,15)), (np.repeat(self.ROC,42)),
                                                                     (np.repeat(0,434)),(np.repeat(self.ROD,36)), np.repeat(0,17)])

        self.wingar = wingar
        self.rolling_fric_coefficient = 0.05

    def velocity(self):
        self.v = np.zeros(len(self.t_array))
        self.v[self.t_array<=self.taxiout_time] = self.taxiout_velocity
        self.v[(self.t_array>self.taxiout_time) & (self.t_array<=self.takeoff_time)] = np.arange(self.taxiout_velocity,self.liftoffv,
                                                                                        self.takeoff_acc*30)
        self.v[(self.t_array>self.takeoff_time) & (self.t_array<=self.cruise_start_time)] = np.arange(self.liftoffv,self.cruisev,
                                                                                    self.climb_acc*30)
        self.v[(self.t_array>self.cruise_start_time) & (self.t_array<=self.cruise_end_time)] = self.cruisev
        self.v[(self.t_array>self.cruise_end_time) & (self.t_array<=self.land_time)] = np.arange(self.cruisev,self.landv,
                                                                                 self.descent_acc*30)
        self.v[(self.t_array>self.land_time) & (self.t_array<=self.stop_time)] = np.arange(self.landv,0,self.land_decc*30)


        #return self.v[(self.t_array>=self.taxiout_time) & (self.t_array<self.taxiout_time+30)]

    def mass(self):
        self.m = np.zeros(len(self.t_array))
        print('hELLO')
        self.m[self.t_array <=self.taxiout_time] = np.linspace(self.w1, self.w2 -self.mf_fuel[0] * 30, len(self.t_array[self.t_array <=self.taxiout_time]))
        # print(np.linspace(self.w1, self.w2 -self.mf_fuel[0] * 30, len(self.t_array[self.t_array <self.taxiout_time])))
        # print(np.linspace(self.w2,self.w3-self.mf_fuel[1]*30,len(self.t_array[(self.t_array >= self.taxiout_time) & (self.t_array < self.takeoff_time)])))
        # print(np.linspace(self.w3,self.w4 -self.mf_fuel[2] * 30,len(self.t_array[(self.t_array >= self.takeoff_time) & (self.t_array < self.cruise_start_time)])))
        self.m[(self.t_array > self.taxiout_time) & (self.t_array <= self.takeoff_time)] = np.linspace(self.w2,
                                                                                                  self.w3-self.mf_fuel[1]*30,
                                                                                                  len(self.t_array[(self.t_array > self.taxiout_time) & (self.t_array <= self.takeoff_time)]))
        self.m[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)] = np.linspace(self.w3,
                                                                                                       self.w4 -self.mf_fuel[2] * 30,
                                                                                                       len(self.t_array[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)]))
        self.m[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)] = np.linspace(self.w4,
                                                                                                        self.w5-self.mf_fuel[3] * 30,
                                                                                                        len(self.t_array[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)]))
        self.m[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)] = np.linspace(self.w5,
                                                                                                self.w6-self.mf_fuel[4] * 30,
                                                                                                len(self.t_array[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)]))
        self.m[(self.t_array > self.land_time) & (self.t_array <= self.stop_time)] = np.linspace(self.w6,
                                                                                          self.w7-self.mf_fuel[5] * 30,
                                                                                          len(self.t_array[(self.t_array > self.land_time) & (self.t_array <= self.stop_time)]))
        #print(self.t_array[[(self.t_array >= self.land_time) & (self.t_array <= self.stop_time)]])
        self.testt = self.m[(self.t_array >= self.takeoff_time) & (self.t_array<self.takeoff_time+30)]

    def dens(self):
        self.alt = np.zeros(len(self.t_array))
        self.alt[np.where(self.t_array <= self.takeoff_time)] = 0
        self.alt[np.where((self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time))] = np.arange(0,self.cruise_altitude, self.ROC*30)
        self.alt[np.where((self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time))] = self.cruise_altitude
        self.alt[np.where((self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time))] = np.arange(self.cruise_altitude, 0, self.ROD*30)
        self.alt[np.where((self.t_array > self.land_time) & (self.t_array <= self.stop_time))] = 0

        self.density_array = np.ones(len(self.alt))
        self.temp_array = np.ones(len(self.alt))

        for alt in range(np.shape(self.alt)[0]):
            self.ISA_calculator(h_input=self.alt[alt])
            self.density_array[alt] = self.rho
            self.temp_array[alt] = self.T
        #return self.density_array, len(self.alt[np.where((self.t_array >= self.land_time) & (self.t_array <= self.stop_time))])
        #return self.alt[np.where((self.t_array >= 263.5*60))], self.temp_array[np.where((self.t_array >= 263.5*60))]
        plt.plot(self.t_array, self.alt)
        plt.show()

    def drag(self):
        self.dens()
        self.mass()
        self.velocity()
        self.CL = np.ones(len(self.t_array))

        self.CL[self.t_array < self.taxiout_time] = 0

        self.CL[(self.t_array >= self.taxiout_time) & (self.t_array <= self.takeoff_time)] = 1.6
            #(self.m[(self.t_array > self.taxiout_time) & (self.t_array <= self.takeoff_time)]* self.g_0 *2) /(self.density_array[(self.t_array > self.taxiout_time) & (self.t_array <= self.takeoff_time)] *(84.9)**2 * self.S)

        self.CL[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)] =(self.m[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)]* self.g_0 *2) /(self.density_array[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)] *(self.v[(self.t_array > self.takeoff_time) & (self.t_array <= self.cruise_start_time)])**2 * self.S)

        self.CL[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)] = (self.m[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)]
                          * self.g_0 *2) /(self.density_array[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)]*
                         (self.v[(self.t_array > self.cruise_start_time) & (self.t_array <= self.cruise_end_time)])**2 * self.S)


        self.CL[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)] = (self.m[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)]
                          * self.g_0 *2) /(self.density_array[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)]*
                         (self.v[(self.t_array > self.cruise_end_time) & (self.t_array <= self.land_time)])**2 * self.S)

        self.CL[(self.t_array > self.land_time) & (self.t_array <= self.stop_time)] = 0

        self.Lift = np.zeros(len(self.t_array))
        self.Lift[(self.t_array >= self.taxiout_time) & (self.t_array <= self.takeoff_time)] = \
            self.CL[(self.t_array >= self.taxiout_time) & (self.t_array <= self.takeoff_time)] * 0.5 * self.density_array[(self.t_array >= self.taxiout_time) & (self.t_array <= self.takeoff_time)] \
            * self.v[(self.t_array >= self.taxiout_time) & (self.t_array <= self.takeoff_time)] **2 * self.S


        self.cdarray = self.cd0array + np.square(self.CL) / (np.pi * self.wingar * self.e)

        self.T_to_overcome_drag = np.array(self.cdarray*0.5*self.density_array*self.v**2 * self.S)
        self.T_to_Acc = self.m * self.acceleration
        self.T_to_climb_and_descend = (self.rates_of_climb_and_descent * self.m * self.g_0)/self.v
        self.T_to_overcome_friction = self.rolling_fric_coefficient * (self.m * self.g_0 - self.Lift)

        self.thrust_required = self.T_to_Acc + self.T_to_overcome_drag + self.T_to_climb_and_descend + self.T_to_overcome_friction


        plt.plot(self.t_array, self.T_to_overcome_drag, label='A320-HACK')
        #ax1.plot(C_L_range, C_D_range_neo, label='A320neo', linestyle='--')

        plt.show()
        # > self.takeoff_time) & (self.t_array <= self.takeoff_time+30))]

        self.ThrustReq_TO = self.thrust_required[np.where((self.t_array == self.takeoff_time))]
        self.ThrustReq_TaxiOut = self.thrust_required[np.where((self.t_array >= self.taxiout_time/2) & (self.t_array < (self.taxiout_time/2)+30))]
        self.ThrustReq_Climb =  self.thrust_required[np.where((self.t_array >= self.mid_climb_time) & (self.t_array < self.mid_climb_time+30))]
        self.ThrustReq_Cruise = self.thrust_required[np.where((self.t_array >= self.cruise_start_time) & (self.t_array < self.cruise_start_time+30))]
        self.ThrustReq_Descent = self.thrust_required[np.where((self.t_array >= self.mid_descent_time) & (self.t_array < self.mid_descent_time+30))]
        self.ThrustReq_TaxiIn = self.thrust_required[np.where((self.t_array >= self.land_time + (self.stop_time-self.land_time)/4) & (self.t_array < self.land_time + (self.stop_time-self.land_time)/4 + 30))]


        #self.testing = self.CL[(self.t_array > self.taxiout_time) & (self.t_array <= self.takeoff_time)]
        self.masstest = self.T_to_climb_and_descend
        #return self.testing

        #return self.thrust_required
        #return self.T_to_overcome_drag
        #
        #return self.Cl_Array
        #return self.Cl_Array





        #self.CD = self.CD0 + (self.CL) **2 / np.pi * self.wingar * self.e

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
    t = thrust_req(cd0clean,wingar)
    con = Constants()

    t.drag()
    t.dens()
    #print(t.drag())

    print('\n Thrust Req TaxiOut = ', t.ThrustReq_TaxiOut/1000 , '[kN]'
          '\n Thrust Req TO = ', t.ThrustReq_TO/1000 , '[kN]'
          '\n Thrust Req Climb = ', t.ThrustReq_Climb/1000,'[kN]'
          '\n Thrust Req Cruise =', t.ThrustReq_Cruise/1000, '[kN]'
          '\n Thrust Req Descend =', t.ThrustReq_Descent/1000, '[kN]'
          '\n Thrust Req Taxi In =', t.ThrustReq_TaxiIn/1000, '[kN]')

    #t.mass()

    #print(t.drag())





    #print('TO T', t.mass(), t.velocity(),len(t.cd0array), t.drag() )
