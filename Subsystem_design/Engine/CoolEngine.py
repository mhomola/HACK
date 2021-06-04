# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:53:39 2021

@author: sarar
"""

"""
    SOME NOTES:

    mr = mass ratio
    mf = mass flow
    _cc = during the combustion, so only the reactants
    _mix = after the cooling air is added

"""
# TODO - adapt this code to classes and to be able to do the iteration

from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
import numpy as np
import matlab.engine
eng = matlab.engine.start_matlab()

test = eng.reactor1('kerosene')

class Engine_Cool(Constants):
    def __init__(self):
        super().__init__()

    ''' DEFINITION OF THE FUNCTION '''
    def cp_regression(self, data, T):

        self.T_before = data[ np.where(data[:,0] < T)[0][-1] ][0]
        self.T_after = data[ np.where(data[:,0] > T)[0][0] ][0]
        self.cp_before = data[ np.where(data[:,0] < T)[0][-1] ][1]
        self.cp_after = data[ np.where(data[:,0] > T)[0][0] ][1]

        self.slope = (self.cp_after - self.cp_before) / (self.T_after - self.T_before)
        self.cp = self.slope * T + (self.cp_after - self.slope * self.T_after)
        self.index_after = np.where(data[:,0] > T)[0][0]


    def cp_temperature(self, data, T):
        self.cp = data[ np.where( data[:,0] == T ) ][0][1]
        self.index_after = np.where( data[:,0] == T )[0][0]

    def cp_between(self, data, i0, T, T_max):
        # input T is only to initialise loop
        for i in range(i0, len(data)):
            if T < T_max:
                self.cp_array = np.append(self.cp_array, data[i][1])
                self.T_array = np.append(self.T_array, data[i][0])
                T = data[i][0]
            else:
                break

    def cp_first_last(self, data, T):

        if not(T in data[:][0]):
            self.cp_regression(data, T)
        else:
            self.cp_temperature(data, T)

        self.cp_array = np.append(self.cp_array, self.cp)
        self.T_array = np.append(self.T_array, T)


    def integral(self, data, T0, Tmax):
        self.cp_array, self.T_array = np.array([]), np.array([])
        # Find initial data to retrieve
        if data[0][0] < 300:
            T0 = T0
        elif data[0][0] == 300:
            T0 = 300 # [K]

        self.cp_first_last(data, T0)
        # Find all other cp's and temperatures
        self.cp_between(data, self.index_after, T0, Tmax)
        # Find final data
        if self.T_array[-1] != Tmax and data[-1][1] >= Tmax:
            self.cp_first_last(data, Tmax)

        self.cp_integral = np.array([])
        for i in range(len(self.cp_array)-1):
            self.deltaT = self.T_array[i+1] - self.T_array[i]
            self.cp_avg = (self.cp_array[i+1] - self.cp_array[i] ) / 2
            self.cp_integral = np.sum( np.append(self.cp_integral, self.deltaT * self.cp_avg) )


    def enthalpy(self, h0):
        self.h = h0 + self.cp_integral

    def SZ_air(self, Tpz, mf_hot, mf_fuel, T03, T04):

        self.mr_SZair = ( (mf_hot + mf_fuel) * self.cp_gas * (T04 - Tpz)) / ( mf_hot * ( self.cp_gas*(T04-Tpz) + self.cp_air*(T03-T04) ) )




aircraft = input('Would you like to analyse the engine of A320neo or A320-HACK? Answer neo or hack: ')
cycle, const, cool = Engine_Cycle(), Constants(), Engine_Cool()

class Emissions:
    def __init__(self):
        self.Tpz = 2000 # [K]  ---> Take from Ivan's code. This is just so that this code can run for now


emiss = Emissions()


save_list = list()
for i in range(len(const.phases)):
    cycle.cycle_analysis(aircraft, i)
    if aircraft == 'neo':
        const.engine_data_neo()
    elif aircraft == 'hack':
        const.engine_data_hack()

    # INITIALIZE WHILE LOOP
    cool.SZ_air(emiss.Tpz, cycle.mf_hot, cycle.mf_fuel, cycle.T03, const.T04)
    e = 0.1
    delta_mr_air = 2*e

    while abs(delta_mr_air) > e * cool.mr_SZair:
        cool.SZ_air(emiss.Tpz, cycle.mf_hot, cycle.mf_fuel, cycle.T03, const.T04)
        delta_mr_air = ( (1 - const.ratio_air_cc[i]) - cool.mr_SZair ) / cool.mr_SZair
        # print('WHILE LOOP\nMass ratio of air needed to be injected on secondary zone:', cool.mr_SZair)
        # print('Mass ratio of air actually injected on secondary zone:', 1 - const.ratio_air_cc[i])
        const.ratio_air_cc[i] = 1-cool.mr_SZair
        np.savetxt("mr_cc_" + aircraft + ".dat", const.ratio_air_cc)

    if __name__ == "__main__":
        print('\nPHASE:', const.phases[i])
        print('\nFINAL\nMass ratio of air needed to be injected on secondary zone:', cool.mr_SZair)
        print('Mass ratio of air actually injected on secondary zone:', 1 - const.ratio_air_cc[i],'\n')
        print('Mass ratio of air needed to be injected to the cc:', 1 - cool.mr_SZair)

    save_list.append(cool.mr_SZair)

np.savetxt("mr_cc_"+aircraft+".dat", save_list)

# if cool.mr_SZair < const.ratio_air_cc:
#     const.ratio_air_cc -= 0.1
# elif cool.mr_SZair > const.ratio_air_cc:
#     const.ratio_air_cc += 0.1
# Check if mr_air_SZ == 1 - const.ratio_air_cc, if not iterate