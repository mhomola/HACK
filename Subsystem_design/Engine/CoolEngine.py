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

    def SZ_air(self, Tpz, mf_hot, mf_h2, mf_ker, T03, T04):

        # Contribution of air that enters the Primary Zone
        self.integral(self.N2_cp_data, T04, Tpz)
        # self.A = - self.cp_integral * mf_hot
        self.A = self.cp_gas * mf_hot * (T04 - Tpz)

        # Contribution of hydrogen
        self.integral(self.h2_cp_data, T04, Tpz)
        # self.B = - self.cp_integral * mf_h2
        self.B = self.cp_gas * mf_h2 * (T04 - Tpz)

        # Contribution of kerosene
        self.integral(self.C12H26_cp_data, T04, Tpz)
        # self.C = - self.cp_integral * mf_ker
        self.C = self.cp_gas * mf_ker * (T04 - Tpz)

        # (Partial) contribution of air that enters the Secondary Zone
        self.integral(self.N2_cp_data, T03, T04)
        # self.D = self.cp_integral * mf_hot
        self.D = self.cp_air * mf_hot * (T04 - T03)

        self.mr_SZair = (self.A + self.B + self.C) / (self.A + self.D)



if __name__ == "__main__":
    cycle = Engine_Cycle()
    const = Constants()
    cool = Engine_Cool()

    class Emissions:
        def __init__(self):
            self.Tpz = 2000 # [K]  ---> Take from Ivan's code. This is just so that this code can run for now


    emiss = Emissions()
    cycle.cycle_analysis('neo', -3)
    const.engine_data_neo()

    cool.SZ_air(emiss.Tpz, cycle.mf_hot, cycle.mf_h2, cycle.mf_ker, cycle.T03, const.T04)

    print('Mass ratio of air needed to be injected on secondary zone:', cool.mr_SZair)
    print('Mass ratio of air actually injected on secondary zone:', const.ratio_air_cc)

    # if cool.mr_SZair < const.ratio_air_cc:
    #     const.ratio_air_cc -= 0.1
    # elif cool.mr_SZair > const.ratio_air_cc:
    #     const.ratio_air_cc += 0.1
    # Check if mr_air_SZ == 1 - const.ratio_air_cc, if not iterate