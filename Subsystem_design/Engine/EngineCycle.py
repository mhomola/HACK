# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:39:04 2021

@author: sarar
"""

# TODO
# Find estimates for all constants
# Find altitude and speed of each phase
# Determine stoichiometric ratio --> find equivalence ratio
# Fuel flow of each fuel

# DONE
# Determine LHV_fuel according to mass ratio of fuels
# Determine atmospheric conditions


from Subsystem_design.common_constants import Constants

import math as m
import numpy as np
import pandas as pd


class Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()

    def data(self, aircraft):
        if aircraft == 'neo':
            self.engine_data_neo()
        elif aircraft == 'hack':
            self.engine_data_hack()

    def cycle_analysis(self, aircraft, i): # i = phase
        self.data(aircraft)

        self.mf_air_init = self.rho0 * self.A_fan * self.v0
        # Total temperature and pressure at inlet
        self.T00 = self.T0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 )
        self.p00 = self.p0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )
        # print('\nInlet:\nT0 = ', T0[i], '[K]; p0 = ', p0[i], ',[Pa]; T00 = ', T00, '[K]; p00 = ', p00, '[Pa]')

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0[i] * (1 + self.eta_inlet * (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )
        # print('\nEntrance of fan:\nT02 = ', self.T02, '[K]; p02 = ', self.p02, '[Pa]')

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + ( self.T02/self.eta_fan ) * ( self.PR_fan ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p021 = self.p02 * self.PR_fan
        # print('\nEntrance of LPC:\nT021 = ', self.T021, '[K]; p021 = ', self.p021, '[Pa]')

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init[i] / (self.BR+1)
        self.mf_cold = self.mf_air_init[i] * self.BR / (self.BR+1)
        # print('\nMass flow of air:\nTotal = ', self.mf_air_init[i], '[kg/s]; Core = ', self.mf_hot, '[kg/s]; Bypassed = ', self.mf_cold,'[kg/s]')

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + ( self.T021/self.eta_LPC ) * ( self.PR_LPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p025 = self.p021 * self.PR_LPC
        # print('\nEntrance of HPC:\nT025 = ', self.T025, '[K]; p025 = ', self.p025, '[Pa]')

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + ( self.T025/self.eta_HPC ) * ( self.PR_HPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p03 = self.p025 * self.PR_HPC
        # print('\nEntrance of CC:\nT03 = ', self.T03, '[K]; p03 = ', self.p03, '[Pa]')

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.ratio_air_cc * self.cp_gas * (self.T04-self.T03)) / (self.LHV_f[i]*10**6 * self.eta_cc)
        self.mf_airfuel = self.mf_hot + self.mf_fuel # at the end of the cc
        # print('\nMass flow CC:\nFuel = ', mf_fuel, '[kg/s]; air CC = ', mf_hot*ratio_air_cc, '[kg/s]; Total end of CC = ', mf_airfuel,'[kg/s]')

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * self.PR_cc
        # print('\nEntrance of HPT:\nT04 = ', self.T04, '[K]; p04 = ', self.p04, '[Pa]')

        # Power to drive fan, LPC, HPC, HPT, LPT [W]
        self.W_fan = self.mf_air_init[i] * self.cp_air * (self.T021-self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025-self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03-self.T025)
        self.W_HPT = self.W_HPC / self.eta_mech
        self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech

        # Exit of HPT - Entrance of LPT
        self.T045 = self.T04 - self.W_HPT / ( self.mf_airfuel * self.cp_gas )
        self.p045 = self.p04 * ( 1 - ( 1 - self.T045/self.T04 ) / self.eta_HPT ) ** ( self.k_gas / (self.k_gas-1) )
        # print('\nEntrance of LPT:\nT045 = ', T045, '[K]; p045 = ', p045, '[Pa]')

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * ( 1 - ( 1 - self.T05/self.T045 ) / self.eta_LPT ) ** ( self.k_gas / (self.k_gas-1) )
        # print('\nEntrance of nozzle:\nT05 = ', T05, '[K]; p05 = ', p05, '[Pa]')

        # Is the nozzle chocked?
        self.PR_cr_nozzle = 1 / ( ( 1 - (self.k_gas-1)/(self.eta_nozzle*(self.k_gas+1)) ) ** (self.k_gas / (self.k_gas-1)) )
        if self.p05/self.p0[i] > self.PR_cr_nozzle:
            print('The nozzle is chocked')
        else:
            print('The nozzle is NOT chocked')

        # Is the fan chocked?
        self.PR_cr_fan = 1 / ( ( 1 - (self.k_air-1)/(self.eta_nozzle*(self.k_air+1)) ) ** (self.k_air / (self.k_air-1)) )
        if self.p021/self.p0[i] > self.PR_cr_fan:
            print('The fan is chocked')
        else:
            print('The fan is NOT chocked')

        # Exit of the nozzle
        self.TR_cr_nozzle = ( self.k_gas + 1 ) / 2
        self.T8 = self.T05 / self.TR_cr_nozzle
        self.p8 = self.p05 / self.PR_cr_nozzle
        self.v8 = np.sqrt( self.k_gas * self.R * self.T8 )

        # Exit of bypassed air
        self.TR_cr_bypassed = ( self.k_air + 1 ) / 2
        self.T18 = self.T021 / self.TR_cr_bypassed
        self.p18 = self.p021 / self.PR_fan
        self.v18 = np.sqrt( self.k_air * self.R * self.T18 )

        # Find thrust of fan, core, and total
        self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
        self.T_fan = self.mf_cold * (self.v18 - self.v0[i]) + self.A18 * (self.p18 - self.p0[i]) # [N]
        self.A8 = (self.mf_airfuel * R * T8) / (p8 * v8)
        self.T_core = self.mf_airfuel * (self.v8 - self.v0[i]) + self.A8 * (self.p8 - self.p0[i]) # [N]
        self.T_total = self.T_fan + self.T_core # [N]
        # print('\nProvided Trhust:\nFan = ', self.T_fan, '[N]; Core = ', self.T_core, '[N]; Total = ', self.T_total, '[N]')


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
    #def cycle_analysis(self):



if __name__ == '__main__':
    aircraft = input('Would you like to do the analysis for the engine of A320neo or A320-HACK? Answer neo or HACK: ')

    ec = Engine_Cycle()
    c = Constants()

    for i in range(len(c.phases)):
        print('** Analysis for ', c.phases[i], ' **')
        ec.cycle_analysis(aircraft, i)

        print('\nInlet:\nT0 = ', c.T0[i], '[K]; p0 = ', c.p0[i], ',[Pa]; T00 = ', ec.T00, '[K]; p00 = ', ec.p00, '[Pa]')
        print('\nEntrance of fan:\nT02 = ', ec.T02, '[K]; p02 = ', ec.p02, '[Pa]')
        print('\nEntrance of LPC:\nT021 = ', ec.T021, '[K]; p021 = ', ec.p021, '[Pa]')
        print('\nEntrance of LPC:\nT021 = ', ec.T021, '[K]; p021 = ', ec.p021, '[Pa]')
        print('\nEntrance of LPC:\nT021 = ', ec.T021, '[K]; p021 = ', ec.p021, '[Pa]')
        print('\nMass flow of air:\nTotal = ', ec.mf_air_init[i], '[kg/s]; Core = ', ec.mf_hot, '[kg/s]; Bypassed = ', ec.mf_cold,'[kg/s]')
        print('\nEntrance of HPC:\nT025 = ', ec.T025, '[K]; p025 = ', ec.p025, '[Pa]')
        print('\nEntrance of CC:\nT03 = ', ec.T03, '[K]; p03 = ', ec.p03, '[Pa]')
        print('\nMass flow CC:\nFuel = ', ec.mf_fuel, '[kg/s]; air CC = ', ec.mf_hot*c.ratio_air_cc, '[kg/s]; Total end of CC = ', ec.mf_airfuel,'[kg/s]')
        print('\nEntrance of HPT:\nT04 = ', ec.T04, '[K]; p04 = ', ec.p04, '[Pa]')
        print('\nEntrance of LPT:\nT045 = ', ec.T045, '[K]; p045 = ', ec.p045, '[Pa]')
        print('\nEntrance of nozzle:\nT05 = ', ec.T05, '[K]; p05 = ', ec.p05, '[Pa]')
        print('\nProvided Trhust:\nFan = ', ec.T_fan, '[N]; Core = ', ec.T_core, '[N]; Total = ', ec.T_total, '[N]')


