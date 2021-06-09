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

from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EngineCycle import Engine_Cycle
import numpy as np

import matlab.engine
eng = matlab.engine.start_matlab()

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

    # Get TPZ from Ivan's code. Inputs are ( gas, P, T, phi )
    Tpz = eng.reactor1('kerosene', float(cycle.p03), float(cycle.T03), 1) # [K]
    cycle.cycle_analysis('neo', -3)
    const.engine_data_neo()

    cool.SZ_air(Tpz, cycle.mf_hot, cycle.mf_h2, cycle.mf_ker, cycle.T03, const.T04)

    print('Mass ratio of air needed to be injected on secondary zone:', cool.mr_SZair)
    print('Mass ratio of air actually injected on secondary zone:', const.ratio_air_cc)

    # if cool.mr_SZair < const.ratio_air_cc:
    #     const.ratio_air_cc -= 0.1
    # elif cool.mr_SZair > const.ratio_air_cc:
    #     const.ratio_air_cc += 0.1
    # Check if mr_air_SZ == 1 - const.ratio_air_cc, if not iterate



# ========= PARTS OF THE OLD VERSION ============ #
#             break
#
#     return cp_array, T_array
#
# def cp_first_last(data, T, cp_array, T_array):
#
#     if not(T in data[:][0]):
#         cp, index_after = cp_regression(data, T)
#     else:
#         cp, index_after = cp_temperature(data, T)
#
#     cp_array = np.append(cp_array, cp)
#     T_array = np.append(T_array, T)
#
#     return cp_array, T_array, index_after
#
#
# def h(h0, data, T0, Tmax):
#     cp_array, T_array = np.array([]), np.array([])
#     # Find initial data to retrieve
#     cp_array, T_array, index_after = cp_first_last(data, T0, cp_array, T_array)
#     # Find all other cp's and temperatures
#     cp_array, T_array = cp_between(data, index_after, cp_array, T_array, T0, Tmax)
#     # Find final data
#     if T_array[-1] != Tmax and data[-1][1] >= Tmax:
#         cp_array, T_array, index_after = cp_first_last(data, Tmax, cp_array, T_array)
#
#     cp_integral = np.array([])
#     for i in range(len(cp_array)-1):
#         deltaT = T_array[i+1] - T_array[i]
#         cp_avg = (cp_array[i+1] - cp_array[i] ) / 2
#         cp_integral = np.append(cp_integral, deltaT*cp_avg)
#
#     h = h0 + np.sum(cp_integral)
#
#     return h
#
#
# ''' BEGINNING OF CODE '''
# ''' Data '''
# N2_cp_data = np.array(np.genfromtxt('N2_cp.dat')) # T[K]; cp[kJ/(kg*K)]
# molarmass_N2 = 28.01340 # [g/mol]
#
# h2_cp_data = np.array(np.genfromtxt('h2_cp.dat')) # T[K]; cp[kJ/(kg*K)]
# molarmass_h2 = 2.01588 # [g/mol]
#
# # https://webbook.nist.gov/cgi/cbook.cgi?ID=C124185&Units=SI&Mask=7
# # C10H22_cp_data = np.array(np.genfromtxt('C10H22_cp.dat'))  # T[K]; cp[J/(mol*K)]
# # h0_C10H22 = -249.7 # [kJ/mol]
# # https://webbook.nist.gov/cgi/cbook.cgi?ID=C95636&Units=SI&Mask=7
# # C9H12_cp_data = np.array(np.genfromtxt('C9H12_cp.dat'))  # T[K]; cp[J/(mol*K)]
# # h0_C9H12 = -13.9 # [kJ/mol]
#
# #
# C12H26_cp_data = np.array(np.genfromtxt('C12H26_cp.dat'))  # T[K]; cp[J/(mol*K)]
# # https://www.chemeo.com/cid/34-125-5/n-Dodecane
# h0_C12H26 = -290.90 # [kJ/mol]
# molarmass_C12H26 = 170.3348 # [g/mol]
#
#
# Tair = 800 # [K] # Temperature of air when injected (both in PZ and SZ)
# Tcc = 2000 # [K] # Temperature of fuel after the combustion, i.e. maximum temperature
# T0 = 300 # [K] 15 deg C = 298 K ~ 300 K # To be used in the formulae
#
# # Mass flows
# # Before combustion
# mf_h2_cc = 1
# mf_ker_cc = 1
# mf_air_cc = 1
# mf_cc = mf_h2_cc + mf_ker_cc + mf_air_cc
# # After cooling air is injected
# mf_cool = 1
# mf_end = mf_cc + mf_cool
#
# # Mass ratios
# # Before combustion
# mr_h2_cc = mf_h2_cc / mf_cc
# mr_ker_cc = mf_ker_cc / mf_cc
# mr_air_cc = mf_air_cc / mf_cc
# # After cooling air is injected
# mr_h2_mix = mf_h2_cc / mf_end
# mr_ker_mix = mf_ker_cc / mf_end
# mr_air_mix = (mf_air_cc + mf_cool) / mf_end
# # When we mix reactants to air to cool engine
# mr_cc = mf_cc / mf_end
# mr_cool = mf_cool/ mf_end
#
#
# '''Find h_cc: find h_ker, h_h2, h_air [J/kg]'''
# # h_decane = h(h0_C10H22*1000, C10H22_cp_data, T0, Tcc) # [J/mol] C10H22
# # h_benzene = h(h0_C9H12*1000, C9H12_cp_data, T0, Tcc) # [J/mol] C9H12
# # h_hexane: couldn't find values
# # h_ker_cc = 0.74/(0.74+0.15)*h_decane + 0.15/(0.74+0.15)*h_benzene # [J/mol] --> Cp not found for C8H19
#                                         # can I use this formula? Or maybe not dividing by the sum
#
# h_ker_cc = h(h0_C12H26*1000, C12H26_cp_data, T0, Tcc) # [J/mol]
#
# # Find h_h2_cc [J/kg]
# h_h2_cc = h(0, h2_cp_data, T0, Tcc) * molarmass_h2 # [J/mol]
#
# # Find h_air_cc [J/kg]
# h_air_cc = h(0, h2_cp_data, T0, Tcc) * molarmass_N2 # [J/mol]
#
# # Total enthalpy in combustion
# h_cc = h_ker_cc*mr_ker_cc + h_h2_cc*mr_h2_cc + h_air_cc*mr_air_cc
#
# # Find h_cool [kJ/kg/K] * [g/mol]
# h_cool = h(0, h2_cp_data, T0, Tair) * molarmass_h2 # [J/mol]
#
# # Find h_mix
# h_mix = mr_cc*h_cc + mr_cool*h_cool # [J/mol]
#
# # Find delta_h_mix (only contribution comes from kerosene, but now mr_ker changes because of the cooling air added)
# # r_h0_C9H12 = 1
# # r_C10H22 = 1
# # delta_h_mix = mr_ker_mix*r_h0_C9H12*h0_C10H22 + mr_ker_mix*r_C10H22*h0_C9H12
# delta_h_mix = mr_ker_mix*h0_C12H26*1000 # [J/mol]
#
# # Find cp_mix
# # Assume a temperature at the end of the cc: T_assumed
# # Find cp_mix for that temperature: C_p_mix(T_assumed) = SUM( c_p_k(T_assumed)*mr_k )
#
# T = Tcc
# e = 50 # [K]
#
# while T >= T0:
#     cp_air_mix = cp_first_last(N2_cp_data, T, np.array([]), np.array([]))[0][0]
#     cp_h2_mix = cp_first_last(h2_cp_data, T, np.array([]), np.array([]))[0][0]
#     cp_ker_mix = cp_first_last(C12H26_cp_data, T, np.array([]), np.array([]))[0][0]
#     cp_mix = mr_air_mix*cp_air_mix + mr_h2_mix*cp_h2_mix + mr_ker_mix*cp_ker_mix
#
#     # Find the temperature at the end of the combustion chamber
#     T_end = T0 + (h_mix - delta_h_mix) / cp_mix
#     print('T_end = ',T_end,' for starting T = ',T)
#     # Check if T_end ~ T_assumed, if not iterate T_assumed
#
#     if T_end < T + e and T_end > T  - e:
#         break
#
#     else:
#         T -= 25 # [K] temperature steps
#
# print('Temperature at the end of the combustion chamber [K]: ', T_end)
#
# # Find the temperature in case delta Cp is sufficiently small
# T_linear = Tcc*mr_cc + Tair*mr_cool
# >>>>>>>>> Temporary merge branch 2
