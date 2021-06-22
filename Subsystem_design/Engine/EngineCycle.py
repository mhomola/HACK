# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:39:04 2021

@author: sarar
"""

# TODO
# Find estimates for all constants --> Aero Engine Technology course

# DONE
# Find altitude and speed of each phase
# Determine LHV_fuel according to mass ratio of fuels
# Determine atmospheric conditions
# Determine stoichiometric ratio --> find equivalence ratio

import numpy as np
from DataEngine import DataFrame
from Subsystem_design.common_constants import Constants
import matlab.engine
import math as m
import matplotlib.pyplot as plt

class Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()

    def data(self, aircraft, phase):
        data, i = self.get_dataframe(aircraft, phase)

        self.M0 = float(data[0])
        self.h = float(data[1])
        self.A_eff = float(data[3]) * np.pi * (float(data[28]) * 0.0254) ** 2 / 4

        self.ISA_calculator(h_input=self.h) # gives self.T0, self.p0, self.rho0, self.a0
        self.T0, self.p0, self.rho0, self.a0 = self.T, self.p, self.rho, self.a
        self.v0 = self.M0 * np.sqrt(self.cp_air * self.R * self.T0)

        self.Thrust = float(data[2])
        self.eta_inlet = float(data[4])
        self.PR_fan = float(data[5])
        self.eta_fan = float(data[6])
        self.BPR = float(data[7])
        self.eta_LPC = float(data[8])
        self.eta_HPC = float(data[9])
        self.PR_LPC = float(data[10])
        self.PR_HPC = float(data[11])
        self.eta_mech = float(data[12])
        self.eta_cc = float(data[13])
        self.PR_cc = float(data[14])
        self.T04 = float(data[15])
        self.eta_LPT = float(data[16])
        self.eta_HPT = float(data[17])
        self.PR_LPT = float(data[18])
        self.PR_HPT = float(data[19])
        self.eta_nozzle = float(data[20])
        self.PR_noz_core = float(data[21])
        self.PR_noz_fan = float(data[22])
        self.mr_h2 = float(data[23])
        self.mr_ker = float(data[24])
        self.ER_h2 = float(data[25])
        self.ER_ker = float(data[26])
        self.LHV_f = float(data[27])

        if aircraft == 'neo':
            self.mr_cc = 0.37
        else:
            if self.mr_h2 == 0:         # When we run HACK only with kerosene - use valves
                self.mr_cc = 0.37
            else:
                self.mr_cc = 0.73

    def get_dataframe(self, aircraft, phs):

        if aircraft == 'neo':
            d = DataFrame().neo

        elif aircraft == 'hack':
            d = DataFrame().hack

        ###################################
        if phs == 'taxi_out':
            data, i = d.taxi_out, 0

        elif phs == 'take_off':
            data, i = d.take_off, 1

        elif phs == 'climb':
            data, i = d.climb, 2

        elif phs == 'cruise':
            data, i = d.cruise, 3

        elif phs == 'approach':
            data, i = d.approach, 4

        elif phs == 'taxi_in':
            data, i = d.taxi_in, 5

        elif phs == 'idle':
            data, i = d.idle, 6

        return data, i

    def cycle_analysis(self, aircraft, phase, flag, alph):          # i = phase
        self.data(aircraft, phase)

        self.v0 = self.M0 * np.sqrt(self.k_air * self.R * self.T0)
        self.rho0 = self.p0 / (self.R*self.T0)

        self.mf_air_init = self.rho0 * self.A_eff * self.v0

        # Total temperature and pressure at inlet
        self.T00 = self.T0 * (1 + (self.k_air-1)/2 * self.M0**2)
        self.p00 = self.p0 * (1 + (self.k_air-1)/2 * self.M0**2) ** (self.k_air / (self.k_air-1))

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0 * (1 + self.eta_inlet * (self.k_air-1)/2 * self.M0**2) ** (self.k_air / (self.k_air-1))

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + (self.T02/self.eta_fan) * (self.PR_fan ** ((self.k_air-1)/self.k_air) - 1)
        self.p021 = self.p02 * self.PR_fan

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init / (self.BPR+1)
        self.mf_cold = self.mf_air_init * self.BPR / (self.BPR+1)

        # Further on the bypass duct
        self.T016 = self.T021
        self.p016 = self.p021 * self.PR_noz_fan

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + (self.T021/self.eta_LPC) * (self.PR_LPC ** ((self.k_air-1)/self.k_air) - 1)
        self.p025 = self.p021 * self.PR_LPC

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + (self.T025/self.eta_HPC) * (self.PR_HPC ** ((self.k_air-1)/self.k_air) - 1)
        self.p03 = self.p025 * self.PR_HPC
        self.OPR = self.p03 / self.p02         # Overall Pressure Ratio

        ''' FIND NECESSARY FUEL FLOW '''
        # START WHILE LOOP
        if self.mr_h2 == 0:                                                 # only kerosene
                self.mf_fuel = (1/14.79) * 0.6 * self.mf_hot * self.mr_cc
        elif self.mr_ker == 0:                                              # only H2
            self.mf_fuel = (1 / 33.99) * 0.25 * self.mf_hot * self.mr_cc
        else:                                                               # H2 + ker
            self.mf_fuel = 0.048086 * 0.35 * self.mf_hot * self.mr_cc

        self.T04 = (self.mf_hot * self.cp_gas * self.T03 + self.mf_fuel * self.LHV_f * 10 ** 6 * self.eta_cc) / (
                    self.cp_gas * (self.mf_fuel + self.mf_hot))
        self.cycle_after_cc()

        if phase in ['taxi_out', 'taxi_in', 'idle', 'approach', 'take_off']:
            # when T_core is a real number and larger than 0, the loop stops
            while m.isnan(self.T_core) or self.T_core < 0:
                self.T04 += 0.5       # if we do not meet thrust condition, increase T04 by 5K, meaning we put more fuel
                self.mf_fuel = (self.mf_hot * self.cp_gas * (self.T04 - self.T03)) / (self.LHV_f * 10 ** 6 *
                                                                        self.eta_cc - self.cp_gas * self.T04)
                self.cycle_after_cc()    # when loop ends, we have our final T04 and mf_fuel
                print('T_total', self.T_total / 1000, 'T_core', self.T_core / 1000)
                print('T04', self.T04, 'fuel flow', self.mf_fuel)

        else:
            while m.isnan(self.T_total) or self.T_total < self.Thrust or self.T_total > self.Thrust + 1000 or self.T_core < 0 or m.isnan(self.T_core):
                # if we do not meet thrust condition, increase/decrease T04 by 5K, meaning we put more/less fuel
                if self.T_total > self.Thrust + 1000 and not self.T_core < 0:
                    self.T04 -= 0.5
                else:
                    self.T04 += 0.5

                self.mf_fuel = (self.mf_hot * self.cp_gas * (self.T04 - self.T03)) / (self.LHV_f * 10 ** 6 *
                                                                        self.eta_cc - self.cp_gas * self.T04)
                self.cycle_after_cc()    # when loop ends, we have our final T04 and mf_fuel

                # if phase == 'climb':
                #     print('T_total', self.T_total/1000, 'T_core', self.T_core/1000)
                #     print('T04', self.T04, 'fuel flow', self.mf_fuel)

        if self.flag_noz == 'yes':
            print('The core nozzle is chocked')
        else:
            print('The core nozzle is NOT chocked')
        if self.flag_fan == 'yes':
            print('The duct nozzle is chocked')
        else:
            print('The duct nozzle is NOT chocked')

        print('[kN] Required Thrust:', round(self.Thrust / 1000, 3), 'Actual Thrust:', round(self.T_total / 1000, 3))
        # print('[kN] Core Thrust:', round(self.T_core / 1000, 3), 'Fan Thrust:', round(self.T_fan / 1000, 3))
        # print('T04 [K]', self.T04, 'Fuel flow [kg/s]:', round(self.mf_fuel, 3))

        # self.T_total = T_total
        self.TSFC_m = (self.mf_fuel*10**3) / (self.T_total*10**(-3))        # [g/kN/s]
        self.TSFC_e = self.mf_fuel * self.LHV_f / (self.T_total * 10 ** (-3))    # [MJ/kN/s]

        ''' USE EQR FROM CoolEngine.py ON THE FIRST ITERATION OF IVAN'S CODE '''
        self.mole_rate()            # gives stoichiometric ratio
        self.air_cool(aircraft)     # gives eqr_PZ

        if flag == True:
            self.plot_TS(aircraft, phase, alph)

    def cycle_after_cc(self):
        self.mf_airfuel = self.mf_hot + self.mf_fuel  # at the end of the cc
        self.mf_h2 = self.mf_fuel * self.ER_h2
        self.mf_ker = self.mf_fuel * self.ER_ker

        # T04 is given
        self.p04 = self.p03 * self.PR_cc

        # Power to drive fan, LPC, HPC, HPT, LPT [W]
        self.W_fan = self.mf_air_init * self.cp_air * (self.T021-self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025-self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03-self.T025)
        self.W_HPT = self.W_HPC / self.eta_mech
        self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech

        # Exit of HPT - Entrance of LPT
        self.T045 = self.T04 - self.W_HPT / (self.mf_airfuel * self.cp_gas)
        self.p045 = self.p04 * (1 - (1 - self.T045/self.T04) / self.eta_HPT) ** (self.k_gas / (self.k_gas-1))

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * (1 - (1 - self.T05/self.T045) / self.eta_LPT) ** (self.k_gas / (self.k_gas-1))

        # Nozzle
        self.T07 = self.T05
        self.p07 = self.p05 * self.PR_noz_core

        # Is the nozzle chocked?
        self.PR_cr_noz_core = 1 / ((1 - (self.k_gas-1)/(self.k_gas+1)/self.eta_nozzle) ** (self.k_gas /
                                                                                           (self.k_gas - 1)))

        # stations = [0, 2, 2.1, 2.5, 3, 4, 4.5, 5, 7]
        # pressures = np.array([self.p00, self.p02, self.p021, self.p025, self.p03, self.p04, self.p045, self.p05, self.p07])/self.p0
        # plt.plot(stations, pressures, label=phase)
        # plt.scatter(1.6, self.p016/self.p00, label=phase)
        # plt.hlines(self.p0/self.p0, 0, 8)
        # plt.legend()
        # plt.show()

        # Exit of the nozzle
        if self.p07/self.p0 > self.PR_cr_noz_core:
            self.flag_noz = 'yes'
            # print('The nozzle is chocked')
            self.TR_cr_noz_core = (self.k_gas + 1) / 2
            self.T8 = self.T07 / self.TR_cr_noz_core
            self.p8 = self.p07 / self.PR_cr_noz_core
            self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
            self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
            self.T_core = self.mf_airfuel * (self.v8 - self.v0) + self.A8 * (self.p8 - self.p0)  # [N]

        elif self.p07/self.p0 <= self.PR_cr_noz_core:
            self.flag_noz = 'no'
            # print('The nozzle is NOT chocked')
            self.p8 = self.p0
            self.T8 = self.T07 * (1 - self.eta_nozzle * (1 - (self.p8/self.p07) ** ((self.k_gas-1)/self.k_gas)))
            self.v8 = np.sqrt(2 * self.cp_gas * (self.T07 - self.T8))
            self.T_core = self.mf_airfuel * (self.v8 - self.v0)  # [N]

        # Is the fan chocked?
        self.PR_cr_fan = 1 / ((1 - (self.k_air-1)/(self.eta_nozzle*(self.k_air+1))) ** (self.k_air / (self.k_air-1)))

        # Exit of bypassed air
        if self.p016/self.p0 > self.PR_cr_fan:
            self.flag_fan = 'yes'
            # print('The fan is chocked')
            self.TR_cr_bypassed = (self.k_air + 1) / 2
            self.T18 = self.T016 / self.TR_cr_bypassed
            self.p18 = self.p016 / self.PR_cr_fan
            self.v18 = np.sqrt(self.k_air * self.R * self.T18)
            self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
            self.T_fan = self.mf_cold * (self.v18 - self.v0) + self.A18 * (self.p18 - self.p0)          # [N]

        elif self.p016/self.p0 <= self.PR_cr_fan:
            self.flag_fan = 'no'
            # print('The fan is NOT chocked')
            self.p18 = self.p0
            self.T18 = self.T016 - self.T016 * self.eta_fan * (1 - (self.p18/self.p016) ** ((self.k_air-1)/self.k_air))
            self.v18 = np.sqrt(2 * self.cp_air * (self.T016 - self.T18))
            self.T_fan = self.mf_cold * (self.v18 - self.v0)        # [N]

        self.T_total = self.T_core + self.T_fan             # [N]


    def air_cool(self, aircraft):
        mr_cool = 1-self.mr_cc
        self.eqr_PZ = (self.mf_fuel / (self.mr_cc * self.mf_hot)) / self.stoichiometric_ratio
        self.eqr_overall = (self.mf_fuel / self.mf_hot) / self.stoichiometric_ratio

        TPZ1 = self.T04 - (mr_cool * self.mf_hot * self.cp_air * (self.T03 - self.T04)) / ((self.mr_cc * self.mf_hot +
                                                                                           self.mf_fuel) * self.cp_gas)
        TPZ2 = self.T03 + (self.mf_fuel * self.eta_cc * self.LHV_f * 10**6) / (self.cp_gas * (self.mr_cc * self.mf_hot +
                                                                                              self.mf_fuel))

        self.TPZ = (TPZ1 + TPZ2) / 2            # Total temperature at the end of the PZ [K]


    def mole_rate(self):
        #Gives me mole rate
        self.n_h2 = self.mf_h2 / (self.molarmass_h2 * 10**-3)
        self.n_ker = self.mf_ker / (self.molar_mass_kerosene * 10**-3)
        self.n_O2 = self.n_h2 * 0.5 + self.n_ker * 14.76
        self.n_N2 = self.n_h2 * 1.88 + self.n_ker * 55.45
        self.m_O2 = self.n_O2 * 32 * 10**-3
        self.m_N2 = self.n_N2 * self.molarmass_N2 * 10**-3
        self.stoichiometric_ratio = (self.mf_h2 + self.mf_ker) / (self.m_O2 + self.m_N2)


    def plot_TS(self, aircraft, phase, alp):

        n = 10              # number of intermediate points between stations 3 and 4

        TR_arr1 = np.array([self.T00/self.T0, self.T02/self.T00, self.T021/self.T02, self.T025/self.T021,
                            self.T03/self.T025])
        TR_arr3 = np.array([self.T045/self.T04, self.T05/self.T045, self.T07/self.T05, self.T07/self.T8])

        PR_arr1 = np.array([self.p00/self.p0, self.p02/self.p00, self.p021/self.p02, self.p025/self.p021,
                            self.p03/self.p025])
        PR_arr3 = np.array([self.p045/self.p04, self.p05/self.p045, self.p07/self.p05, self.p07/self.p8])

        TR_arr4 = np.array([self.T02/self.T00, self.T021/self.T02, self.T016/self.T021, self.T016/self.T18])
        PR_arr4 = np.array([self.p02/self.p00, self.p021/self.p02, self.p016/self.p021, self.p016/self.p18])

        ds1 = self.cp_gas * np.log(TR_arr1) - self.R * np.log(PR_arr1)
        ds3 = self.cp_gas * np.log(TR_arr3) - self.R * np.log(PR_arr3)
        ds2 = np.array([(self.cp_gas * np.log(self.T04/self.T03) - self.R * np.log(self.p04/self.p03))/(n+1)] * (n+1))
        ds4 = self.cp_gas * np.log(TR_arr4) - self.R * np.log(PR_arr4)

        ds = np.append(ds1, ds2)
        ds = np.append(ds, ds3)
        s = np.cumsum(ds)
        s4 = s[0] + np.cumsum(ds4)

        T1 = [self.T00, self.T02, self.T021, self.T025, self.T03]
        T3 = [self.T045, self.T05, self.T07, self.T8]
        T4 = [self.T02, self.T021, self.T016, self.T18]
        T2, PL_arr = list(), list()
        PL4 = 1 - self.PR_cc

        for i in range(1, n+2):
            PR = 1 - (PL4 * i / (n+1))
            T_cc = self.T03 * np.exp((i * ds2[0] + self.R * np.log(PR)) / self.cp_gas)
            # T_cc = self.T03 * 2 ** (i/(n+1))
            T2.append(T_cc)

        print('T in cc', T2)
        print('ds2', ds2)

        T = np.array(T1+T2+T3)
        plt.gcf().canvas.set_window_title(phase)
        plt.plot(s, T, 'g', label=aircraft, alpha=alp)
        plt.plot(s[0:5], T[0:5], 'go', alpha=alp)
        plt.plot(s[-5:], T[-5:], 'go', alpha=alp)
        plt.plot(s4, T4, 'r-o', alpha=alp)
        plt.xlabel('Entropy [J/kg/K]', fontsize=15)
        plt.ylabel('Temperature [K]', fontsize=15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=15)

''' FORMULAE

Stations:
    0 - ambient
    13 - bypassed air
    2 - entrance of the fan
    21 - entrance of LPC
    25 - entrance of HPC
    3 - entrance of cc
    4 - exit of cc
    45 - exit of HPT
    5 - exit of LPT
    18 - exit of bypassed air
    8 - exit of nozzle from core
'''


if __name__ == '__main__':
    ec = Engine_Cycle()
    aircraft = ['neo', 'hack']
    phases = ['taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in', 'idle']
    #phases = ['taxi_out', 'take_off1', 'take_off2', 'climb1', 'climb2', 'cruise1', 'cruise2', 'approach1', 'approach2', 'taxi_in', 'idle']
    # aircraft = ['neo']
    # phases = ['cruise']

    TO_arr, cruise_arr = list(), list()

    for a in aircraft:
        print("\n= = = = Analysis for A320", a, "= = = =")
        for p in phases:
            print("\n", p)
            ec.cycle_analysis(a, p, flag=False, alph=0)         # since flag = False, alph doesn't matter


            # print('\nInlet: T0 = ', round(ec.T0, 3), '[K]; p0 = ', round(ec.p0, 3), '[Pa]; v0 = ', round(ec.v0, 3),
            #       '[m/s]')
            # print('T00 = ', round(ec.T00, 3), '[K]; p00 = ', round(ec.p00, 3), '[Pa]')
            # print('Entrance of fan: T02 = ', round(ec.T02, 3), '[K]; p02 = ', round(ec.p02,3), '[Pa]')
            # print('Entrance of LPC: T021 = ', round(ec.T021, 3), '[K]; p021 = ', round(ec.p021,3), '[Pa]')
            #print('Mass flow of air: Total = ', round(ec.mf_air_init, 3), '[kg/s]; Core = ', round(ec.mf_hot, 3),
                  #'[kg/s]; Bypassed = ', round(ec.mf_cold, 3), '[kg/s]')
            print('Entrance of HPC: T025 = ', round(ec.T025, 3), '[K]; p025 = ', round(ec.p025, 3), '[Pa]')
            print('Entrance of CC: T03 = ', round(ec.T03, 3), '[K]; p03 = ', round(ec.p03, 3), '[Pa]; OPR = ',
                   round(ec.OPR, 3))
            print('Mass flow CC: Fuel = ', round(ec.mf_fuel, 3), '[kg/s]; mf_h2 = ', round(ec.mf_h2,3),
                  '[kg/s]; mf_ker = ', round(ec.mf_ker, 3), '[kg/s]')
            print('Moles reacting per second of kerosene:', round(ec.n_ker, 3), 'of H2:', round(ec.n_h2, 3), 'of O2:', round(ec.n_O2, 3), 'of N2:', round(ec.n_N2, 3))
            # print('LHV fuel = ', round(ec.LHV_f, 3))
            print('TPZ = ', round(ec.TPZ, 3))
            # print('Power: Fan = ', round(ec.W_fan, 3), '[W]; LPC = ', round(ec.W_LPC, 3), '[W]; HPC = ',
            #       round(ec.W_HPC, 3), '[W]')
            # print('LPT = ', round(ec.W_LPT, 3), '[W]; HPT = ', round(ec.W_HPT, 3), '[W]')
            print('Entrance of HPT: T04 = ', round(ec.T04, 3), '[K]; p04 = ', round(ec.p04, 3), '[Pa]')
            # print('Entrance of LPT: T045 = ', round(ec.T045, 3), '[K]; p045 = ', round(ec.p045, 3), '[Pa]')
            # print('Entrance of nozzle: T05 = ', round(ec.T05, 3), '[K]; p05 = ', round(ec.p05, 3), '[Pa]')
            # print('Exit of nozzle: T07 = ', round(ec.T07, 3), '[K]; p07 = ', round(ec.p07, 3), '[Pa]; PR_cr_noz = ',
            #       round(ec.PR_cr_noz_core, 3))
            # print('Exit of nozzle: T8 = ', round(ec.T8, 3), '[K]; p8 = ', round(ec.p8, 3), '[Pa]; v8 = ',
            #       round(ec.v8, 3), '[m/s]')
            # print('Exit of fan: T016 = ', round(ec.T016, 3), '[K]; p016 = ', round(ec.p016, 3), '[Pa]; PR_cr_fan = ',
            #       round(ec.PR_cr_fan, 3))
            # print('Exit of fan: T18 = ', round(ec.T18, 3), '[K]; p18 = ', round(ec.p18, 3), '[Pa]; v18 = ',
            #       round(ec.v18,3), '[m/s]')
            print('Provided Thrust: Fan = ', round(ec.T_fan, 3), '[N]; Core = ', round(ec.T_core, 3), '[N]; Total = ',
                  round(ec.T_total, 3), '[N]')
            # print('Thrust SFC = ', round(ec.TSFC_m, 5), '[g/kN/s];', round(ec.TSFC_e, 5), '[MJ/kN/s]')
            print('Eqr PZ:', round(ec.eqr_PZ, 3), 'Eqr Overall:', round(ec.eqr_overall, 3))

            amb = [['T0', round(ec.T0, 3), 'K'], ['p0', round(ec.p0, 3), 'Pa'], ['v0', round(ec.v0, 3), 'm/s']]
            air = [['m_intake', round(ec.mf_air_init, 3), 'kg/s'], ['m_hot', round(ec.mf_hot, 3), 'kg/s'],
                   ['m_cold', round(ec.mf_cold, 3), 'kg/s']]
            st0 = [['T00', round(ec.T00, 3), 'K'], ['p00', round(ec.p00, 3), 'Pa']]
            st2 = [['T02', round(ec.T02, 3), 'K'], ['p02', round(ec.p02, 3), 'Pa']]
            BPR = ['BPR', round(ec.BPR, 3), '-']
            st21 = [['T021', round(ec.T021, 3), 'K'], ['p021', round(ec.p021, 3), 'Pa']]
            st25 = [['T025', round(ec.T025, 3), 'K'], ['p025', round(ec.p025, 3), 'Pa']]
            st3 = [['T03', round(ec.T03, 3), 'K'], ['p03', round(ec.p03, 3), 'Pa']]
            st4 = [['T04', round(ec.T04, 3), 'K'], ['p04', round(ec.p04, 3), 'Pa']]
            fuel = [['m_fuel', round(ec.mf_fuel, 3), 'kg/s'], ['m_h2', round(ec.mf_h2, 3), 'kg/s'],
                    ['m_ker', round(ec.mf_ker, 3), 'kg/s']]
            st45 = [['T045', round(ec.T045, 3), 'K'], ['p045', round(ec.p045, 3), 'Pa']]
            st5 = [['T05', round(ec.T05, 3), 'K'], ['p05', round(ec.p05, 3), 'Pa']]
            st7 = [['T07', round(ec.T07, 3), 'K'], ['p07', round(ec.p07, 3), 'Pa']]
            st8 = [['T8', round(ec.T8, 3), 'K'], ['p8', round(ec.p8, 3), 'Pa'], ['v8', round(ec.v8, 3), 'm/s']]
            st16 = [['T016', round(ec.T016, 3), 'K'], ['p016', round(ec.p016, 3), 'Pa']]
            st18 = [['T18', round(ec.T18, 3), 'K'], ['p18', round(ec.p18, 3), 'Pa'], ['v18', round(ec.v18, 3), 'm/s']]
            Thr = [['T_fan', round(ec.T_fan, 3), 'N'], ['T_core', round(ec.T_core, 3), 'N'],
                   ['T_tot', round(ec.T_total, 3), 'N'], ['TSCF', round(ec.TSFC_m, 5), 'g/kN/s'],
                   ['TSCF', round(ec.TSFC_e, 5), 'MJ/kN/s']]
            OPR = ['OPR', round(ec.OPR, 3), '-']
            eqr = [['PZ eqr', round(ec.eqr_PZ, 3), '-'], ['Overall eqr', round(ec.eqr_overall, 3), '-']]

            save_txt = amb + air + st0 + st2 + [BPR] + st21 + st25 + st3 + st4 + fuel + st45 + st5 + st7 + st8 + \
                       st16 + st18 + Thr + [OPR] + eqr
            name = a+'_'+p+'.txt'

            F = open(name, 'w')
            for i in range(len(save_txt)):
                for j in range(0, 3):
                    F.write(str(save_txt[i][j]) + '\t')
                F.write('\n')

            F.close()

    def plot_m(mark, l, T, SFC):
        plt.plot(T, SFC, mark, label=l)
        plt.xlabel('Net Thrust [kN]', fontsize=15)
        plt.ylabel('TSFC [kg/kN/s]', fontsize=15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=15)

    def plot_e(mark, l, T, SFC):
        plt.plot(T, SFC, mark, label=l)
        plt.xlabel('Net Thrust [kN]', fontsize=15)
        plt.ylabel('TSFC [MJ/kN/s]', fontsize=15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=15)

    ec.cycle_analysis('neo', 'take_off', flag=False, alph=1)
    Th_nto, SFC_m_nto, SFC_e_nto = ec.T_total/1000, ec.TSFC_m, ec.TSFC_e
    ec.cycle_analysis('hack', 'take_off', flag=False, alph=1)
    Th_hto, SFC_m_hto, SFC_e_hto = ec.T_total / 1000, ec.TSFC_m, ec.TSFC_e
    ec.cycle_analysis('neo', 'cruise', flag=False, alph=1)
    Th_nc, SFC_m_nc, SFC_e_nc = ec.T_total / 1000, ec.TSFC_m, ec.TSFC_e
    ec.cycle_analysis('hack', 'cruise', flag=False, alph=1)
    Th_hc, SFC_m_hc, SFC_e_hc = ec.T_total / 1000, ec.TSFC_m, ec.TSFC_e

    plt.figure()
    plot_e('bx', 'A320neo - Cruise', Th_nc, SFC_e_nc)
    plot_e('gx', 'A320-HACK - Cruise', Th_hc, SFC_e_hc)
    plot_e('bo', 'A320neo - Take-off', Th_nto, SFC_e_nto)
    plot_e('go', 'A320-HACK - Take-off', Th_hto, SFC_e_hto)
    plt.show()

    plt.figure()
    plot_m('bx', 'A320neo - Cruise', Th_nc, SFC_m_nc)
    plot_m('gx', 'A320-HACK - Cruise', Th_hc, SFC_m_hc)
    plot_m('bo', 'A320neo - Take-off', Th_nto, SFC_m_nto)
    plot_m('go', 'A320-HACK - Take-off', Th_hto, SFC_m_hto)
    plt.show()

    plt.figure()
    print('\nNEO CRUISE - TS')
    ec.cycle_analysis('neo', 'cruise', flag=True, alph=0.2)
    print('\nHACK CRUISE - TS')
    ec.cycle_analysis('hack', 'cruise', flag=True, alph=1)
    plt.show()

    plt.figure()
    ec.cycle_analysis('neo', 'take_off', flag=True, alph=0.2)
    ec.cycle_analysis('hack', 'take_off', flag=True, alph=1)
    plt.show()


