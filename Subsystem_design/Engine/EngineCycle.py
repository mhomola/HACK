# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:39:04 2021

@author: sarar
"""

# TODO
# Find estimates for all constants --> Aero Engine Technology course
# Find altitude and speed of each phase --> Elena
# Check why when the nozzle is not chocked T8 > T05 (hypothesis: error with p05)

# DONE
# Determine LHV_fuel according to mass ratio of fuels
# Determine atmospheric conditions
# Fuel flow of each fuel
# Determine stoichiometric ratio --> find equivalence ratio --> Elena (check what is used rn for stoichiometric ratio)
# Verification: ran code with values from P&P past paper and results matched


from Subsystem_design.common_constants import Constants
import numpy as np
import pandas as pd

class Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()

    def data(self, aircraft):
        if aircraft == 'neo':
            self.engine_data_neo()
            print("Efficiency cc: ",self.eta_cc)
        elif aircraft == 'hack':
            self.engine_data_hack()

    def cycle_analysis(self, aircraft, i): # i = phase
        self.data(aircraft)
        print( 'Bleed: ', self.mf_bleed )

        self.mf_air_init = self.rho0 * self.A_fan * self.v0
        self.mf_air_init[0] = 298 # [kg/s]
        # Total temperature and pressure at inlet
        self.T00 = self.T0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 )
        self.p00 = self.p0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0[i] * (1 + self.eta_inlet * (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + ( self.T02/self.eta_fan ) * ( self.PR_fan ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p021 = self.p02 * self.PR_fan

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init[i] / (self.BR+1)
        self.mf_cold = self.mf_air_init[i] * self.BR / (self.BR+1)

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + ( self.T021/self.eta_LPC ) * ( self.PR_LPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p025 = self.p021 * self.PR_LPC

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + ( self.T025/self.eta_HPC ) * ( self.PR_HPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p03 = self.p025 * self.PR_HPC

        # Power to drive fan, LPC, HPC[W]
        self.W_fan = self.mf_air_init[i] * self.cp_air * (self.T021-self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025-self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03-self.T025)

        # Bleed air
        self.mf_hot = self.mf_hot - self.mf_bleed

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.ratio_air_cc[i] * self.cp_gas * (self.T04-self.T03)) / (self.LHV_f[i]*10**6 * self.eta_cc)
        self.mf_airfuel = self.mf_hot + self.mf_fuel # at the end of the cc
        self.mf_h2 = self.mf_fuel * self.ER_h2[i]   #energy ratio H2
        self.mf_ker = self.mf_fuel * self.ER_ker[i]

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * self.PR_cc

        # Power to drive HPT, LPT [W]
        self.W_HPT = self.W_HPC / self.eta_mech_H
        self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech_L

        # Exit of HPT - Entrance of LPT
        self.T045 = self.T04 - self.W_HPT / ( self.mf_airfuel * self.cp_gas )
        self.p045 = self.p04 * ( 1 - ( 1 - self.T045/self.T04 ) / self.eta_HPT ) ** ( self.k_gas / (self.k_gas-1) )


        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * ( 1 - ( 1 - self.T05/self.T045 ) / self.eta_LPT ) ** ( self.k_gas / (self.k_gas-1) )
        self.OPR = self.p03 / self.p02

        # Is the nozzle chocked?
        self.PR_cr_nozzle = 1 / ( ( 1 - (self.k_gas-1)/(self.k_gas+1)/self.eta_nozzle) ** (self.k_gas / (self.k_gas-1)) )

        # Exit of the nozzle
        if self.p05/self.p0[i] > self.PR_cr_nozzle:
            print('The nozzle is chocked')
            self.TR_cr_nozzle = (self.k_gas + 1) / 2
            self.T8 = self.T05 / self.TR_cr_nozzle
            self.p8 = self.p05 / self.PR_cr_nozzle
            self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
            self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
            self.T_core = self.mf_airfuel * (self.v8 - self.v0[i]) + self.A8 * (self.p8 - self.p0[i])  # [N]

        elif self.p05/self.p0[i] < self.PR_cr_nozzle:
            print('The nozzle is NOT chocked')
            self.p8 = self.p0[i]
            self.T8 = self.T05 * ( 1- ( self.eta_nozzle * ( 1- ( (self.p8/self.p05)**((self.k_gas-1)/self.k_gas) ) ) ) )
            self.v8 = np.sqrt( 2 * self.cp_gas * (self.T05 - self.T8) )
            self.T_core = self.mf_airfuel * ( self.v8 - self.v0[i] )

        # Is the fan chocked?
        self.PR_cr_fan = 1 / ( ( 1 - (self.k_air-1)/(self.eta_nozzle*(self.k_air+1)) ) ** (self.k_air / (self.k_air-1)) )

        # Exit of bypassed air
        if self.p021/self.p0[i] > self.PR_cr_fan:
            print('The fan is chocked')
            self.TR_cr_bypassed = (self.k_air + 1) / 2
            self.T18 = self.T021 / self.TR_cr_bypassed
            self.p18 = self.p021 / self.PR_cr_fan
            self.v18 = np.sqrt(self.k_air * self.R * self.T18)
            self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
            self.T_fan = self.mf_cold * (self.v18 - self.v0[i]) + self.A18 * (self.p18 - self.p0[i])  # [N]

        elif self.p021/self.p0[i] < self.PR_cr_fan:
            print('The fan is NOT chocked')
            self.p18 = self.p0[i]
            self.T18 = self.T021 - self.T021 * self.eta_fan * ( 1 - (self.p18/self.p021) ** ( (self.k_air-1)/self.k_air ) )
            self.v18 = np.sqrt( 2 * self.cp_air * (self.T021 - self.T18) )
            self.T_fan = self.mf_cold * ( self.v18 - self.v0[i] )


        self.T_total = self.T_fan + self.T_core # [N]

        # self.stoichiometric_ratio = self.mr_h2[i] * self.stoich_ratio_h2 + self.mr_ker[i] * self.stoich_ratio_ker
        # self.equivalence_ratio = (self.mf_fuel/(self.mf_hot * self.ratio_air_cc)) / self.stoichiometric_ratio #TBD what mf_air to use


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
    aircraft = input('Would you like to analyse the engine of A320neo or A320-HACK? Answer neo or HACK: ')

    ec = Engine_Cycle()
    c = Constants()
    if aircraft == 'neo':
        c.engine_data_neo()
    elif aircraft == 'hack':
        c.engine_data_hack()


    for i in range(len(c.phases)):
        print('\n** Analysis for', c.phases[i], ' **')
        ec.cycle_analysis(aircraft=aircraft, i=i)

        print('\nInlet: v0 = ', ec.v0[i], '[m/s]; T0 = ', c.T0[i], '[K]; p0 = ', c.p0[i], ',[Pa]; T00 = ', ec.T00, '[K]; p00 = ', ec.p00, '[Pa]')
        print('Entrance of fan: T02 = ', ec.T02, '[K]; p02 = ', ec.p02, '[Pa]')
        print('Entrance of LPC: T021 = ', ec.T021, '[K]; p021 = ', ec.p021, '[Pa]')
        print('Mass flow of air: Total = ', ec.mf_air_init[i], '[kg/s]; Core = ', ec.mf_hot, '[kg/s]; Bypassed = ', ec.mf_cold,'[kg/s]')
        print('Entrance of HPC: T025 = ', ec.T025, '[K]; p025 = ', ec.p025, '[Pa]')
        print('Entrance of CC: T03 = ', ec.T03, '[K]; p03 = ', ec.p03, '[Pa]')
        print('\nMass flow CC: Fuel = ', ec.mf_fuel, '[kg/s]; air CC = ', ec.mf_hot*c.ratio_air_cc[i], '[kg/s]; Total end of CC = ', ec.mf_airfuel,'[kg/s]')

        print('Hydrogen = ', ec.mf_h2, '[kg/s]; Kerosene = ', ec.mf_ker,'[kg/s]')

        print('\nEntrance of HPT: T04 = ', ec.T04, '[K]; p04 = ', ec.p04, '[Pa]')
        print('Entrance of LPT: T045 = ', ec.T045, '[K]; p045 = ', ec.p045, '[Pa]')
        print('Entrance of nozzle: T05 = ', ec.T05, '[K]; p05 = ', ec.p05, '[Pa]')
        print('Overall pressure ratio, OPR = ', ec.OPR)
        print('Exit of nozzle: T8 = ', ec.T8, '[K]; p8 = ', ec.p8, '[Pa]; v8 = ', ec.v8, '[m/s]')
        print('Exit of fan: T18 = ', ec.T18, '[K]; p18 = ', ec.p18, '[Pa]; v18 = ', ec.v18, '[m/s]')

        print('\nW_fan = ', ec.W_fan, 'W_LPC = ', ec.W_LPC, 'W_HPC = ', ec.W_HPC, 'W_LPT = ', ec.W_LPT, 'W_HPT = ', ec.W_HPT)
        print('Provided Thrust: Fan = ', ec.T_fan, '[N]; Core = ', ec.T_core, '[N]; Total = ', ec.T_total, '[N]')

        #print('Equivalence Ratio = ', ec.equivalence_ratio)

