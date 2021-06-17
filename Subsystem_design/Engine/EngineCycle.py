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

class Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()

    def data(self, aircraft, phase):
        data, i = self.get_dataframe(aircraft, phase)
        self.M0 = float(data[0])
        self.h = float(data[1])
        self.ISA_calculator(h_input=self.h) # gives self.T0, self.p0, self.rho0, self.a0
        self.T0, self.p0, self.rho0, self.a0 = self.T, self.p, self.rho, self.a
        self.v0 = self.M0 * np.sqrt(self.cp_air * self.R * self.T0)
        self.Thrust = float(data[2])
        self.A_eff = float(data[3]) * np.pi * (float(data[28]) * 0.0254)**2 / 4
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

        # percentage of core air that is used in combustion
        # if aircraft == 'neo':
        #     self.mr_air_cc = np.array(np.genfromtxt('mr_cc_neo.dat'))[i]
        # elif aircraft == 'hack':
        #     self.mr_air_cc = np.array(np.genfromtxt('mr_cc_hack.dat'))[i]

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


    def cycle_analysis(self, aircraft, phase): # i = phase
        self.data(aircraft, phase)

        self.v0 = self.M0 * np.sqrt(self.k_air * self.R * self.T0)
        self.rho0 = self.p0 / (self.R*self.T0)

        self.mf_air_init = self.rho0 * self.A_eff * self.v0

        # Total temperature and pressure at inlet
        self.T00 = self.T0 * ( 1 + (self.k_air-1)/2 * self.M0**2 )
        self.p00 = self.p0 * ( 1 + (self.k_air-1)/2 * self.M0**2 ) ** (self.k_air / (self.k_air-1) )

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0 * (1 + self.eta_inlet * (self.k_air-1)/2 * self.M0**2 ) ** (self.k_air / (self.k_air-1) )

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + ( self.T02/self.eta_fan ) * ( self.PR_fan ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p021 = self.p02 * self.PR_fan

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init / (self.BPR+1)
        self.mf_cold = self.mf_air_init * self.BPR / (self.BPR+1)

        # Further on the bypass duct
        self.T016 = self.T021
        self.p016 = self.p021 * self.PR_noz_fan

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + ( self.T021/self.eta_LPC ) * ( self.PR_LPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p025 = self.p021 * self.PR_LPC

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + ( self.T025/self.eta_HPC ) * ( self.PR_HPC ** ( (self.k_air-1)/self.k_air ) - 1 )
        self.p03 = self.p025 * self.PR_HPC
        self.OPR = self.p03 / self.p02         # Overall Pressure Ratio

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.cp_gas * (self.T04-self.T03)) / (self.LHV_f*10**6 * self.eta_cc)
        self.mf_airfuel = self.mf_hot + self.mf_fuel # at the end of the cc
        self.mf_h2 = self.mf_fuel * self.ER_h2
        self.mf_ker = self.mf_fuel * self.ER_ker

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * self.PR_cc

        # Power to drive fan, LPC, HPC, HPT, LPT [W]
        self.W_fan = self.mf_air_init * self.cp_air * (self.T021-self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025-self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03-self.T025)
        self.W_HPT = self.W_HPC / self.eta_mech
        self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech

        # Exit of HPT - Entrance of LPT
        self.T045 = self.T04 - self.W_HPT / ( self.mf_airfuel * self.cp_gas )
        self.p045 = self.p04 * ( 1 - ( 1 - self.T045/self.T04 ) / self.eta_HPT ) ** ( self.k_gas / (self.k_gas-1) )

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * ( 1 - ( 1 - self.T05/self.T045 ) / self.eta_LPT ) ** ( self.k_gas / (self.k_gas-1) )

        # Nozzle
        self.T07 = self.T05
        self.p07 = self.p05 * self.PR_noz_core

        # Is the nozzle chocked?
        self.PR_cr_noz_core = 1 / ( ( 1 - (self.k_gas-1)/(self.k_gas+1)/self.eta_nozzle) ** (self.k_gas / (self.k_gas-1)) )

        # Exit of the nozzle
        if self.p07/self.p0 > self.PR_cr_noz_core:
            print('The nozzle is chocked')
            self.TR_cr_noz_core = (self.k_gas + 1) / 2
            self.T8 = self.T07 / self.TR_cr_noz_core
            self.p8 = self.p07 / self.PR_cr_noz_core
            self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
            self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
            self.T_core = self.mf_airfuel * (self.v8 - self.v0) + self.A8 * (self.p8 - self.p0)  # [N]

        elif self.p07/self.p0 <= self.PR_cr_noz_core:
            print('The nozzle is NOT chocked')
            self.p8 = self.p0
            self.T8 = self.T07 * ( 1 - self.eta_nozzle * ( 1 - (self.p8/self.p07) ** ( (self.k_gas-1)/self.k_gas ) ) )
            self.v8 = np.sqrt( 2 * self.cp_gas * (self.T07 - self.T8) )
            self.T_core = self.mf_airfuel * ( self.v8 - self.v0 )  # [N]

        # Is the fan chocked?
        self.PR_cr_fan = 1 / ( ( 1 - (self.k_air-1)/(self.eta_nozzle*(self.k_air+1)) ) ** (self.k_air / (self.k_air-1)) )

        # Exit of bypassed air
        if self.p016/self.p0 > self.PR_cr_fan:
            print('The fan is chocked')
            self.TR_cr_bypassed = (self.k_air + 1) / 2
            self.T18 = self.T016 / self.TR_cr_bypassed
            self.p18 = self.p016 / self.PR_cr_fan
            self.v18 = np.sqrt(self.k_air * self.R * self.T18)
            self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
            self.T_fan = self.mf_cold * (self.v18 - self.v0) + self.A18 * (self.p18 - self.p0)  # [N]

        elif self.p016/self.p0 <= self.PR_cr_fan:
            print('The fan is NOT chocked')
            self.p18 = self.p0
            self.T18 = self.T016 - self.T016 * self.eta_fan * ( 1 - (self.p18/self.p016) ** ( (self.k_air-1)/self.k_air ) )
            self.v18 = np.sqrt( 2 * self.cp_air * (self.T016 - self.T18) )
            self.T_fan = self.mf_cold * ( self.v18 - self.v0 ) # [N]


        self.T_total = self.T_fan + self.T_core # [N]
        self.TSFC = self.mf_fuel / (self.T_total*10**(-3)) # [g/kN/s]

        ''' USE EQR FROM CoolEngine.py ON THE FIRST ITERATION OF IVAN'S CODE '''
        self.stoichiometric_ratio = self.stoich_ratio_ker_h2
        #self.stoichiometric_ratio = self.mr_h2 * self.stoich_ratio_h2 + self.mr_ker * self.stoich_ratio_ker # UPDATE THIS, SOFIA
        #self.mf_air_combustion = self.mf_hot * self.mr_air_cc
        #self.equivalence_ratio = (self.mf_fuel / (self.mf_air_combustion)) / \
                                  #self.stoichiometric_ratio

        self.air_cool()
        self.mole_rate()



    def air_cool(self):

        self.mr_SZair_simpl1 = (self.mf_airfuel*self.cp_gas*(self.T04-self.T03) - self.mf_fuel*self.eta_cc*self.LHV_f*10**6) /\
                               (self.mf_hot * (self.T04-self.T03)*(self.cp_gas-self.cp_air))
        self.TPZ = self.T03 + ( self.mf_fuel*self.eta_cc*self.LHV_f*10**6 ) / (self.cp_gas*( (1-self.mr_SZair_simpl1)*self.mf_hot+self.mf_fuel ))

        self.mr_SZair_simpl = (self.mf_airfuel * self.cp_gas * (self.TPZ - self.T04)) / (
                    self.mf_hot * (self.cp_air * (self.T04 - self.T03) + self.cp_gas * (self.TPZ - self.T04))) # just to check but should be the same as simpl1

    def mole_rate(self):
        #Gives me mole rate
        self.n_h2 = self.mf_h2/(self.molarmass_h2 * 10**-3)
        self.n_ker = self.mf_ker/(self.molar_mass_kerosene*10**-3)
        self.n_O2 = self.n_h2 * 0.5 + self.n_ker * 14.76
        self.n_N2 = self.n_h2 * 1.88 + self.n_ker * 55.45
        self.m_O2 = self.n_O2* 32 *10**-3
        self.m_N2 = self.n_N2 * self.molarmass_N2 *10**-3


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

    phases = ['taxi_out', 'take_off', 'climb', 'cruise', 'approach', 'taxi_in']

    for a in aircraft:
        print("\n= = = = Analysis for A320",a,"= = = =")
        for p in phases:
            print("\n",p)
            ec.cycle_analysis(a, p)

            print('\nInlet: T0 = ', round(ec.T0,3), '[K]; p0 = ', round(ec.p0,3), '[Pa]; v0 = ', round(ec.v0,3), '[m/s]')
            print('T00 = ', round(ec.T00,3), '[K]; p00 = ', round(ec.p00,3), '[Pa]')
            print('Entrance of fan: T02 = ', round(ec.T02,3), '[K]; p02 = ', round(ec.p02,3), '[Pa]')
            print('Entrance of LPC: T021 = ', round(ec.T021,3), '[K]; p021 = ', round(ec.p021,3), '[Pa]')
            print('Mass flow of air: Total = ', round(ec.mf_air_init,3), '[kg/s]; Core = ', round(ec.mf_hot,3), '[kg/s]; Bypassed = ', round(ec.mf_cold,3),'[kg/s]')
            print('Entrance of HPC: T025 = ', round(ec.T025,3), '[K]; p025 = ', round(ec.p025,3), '[Pa]')
            print('Entrance of CC: T03 = ', round(ec.T03,3), '[K]; p03 = ', round(ec.p03,3), '[Pa]; OPR = ', round(ec.OPR,3))
            print('Mass flow CC: Fuel = ', round(ec.mf_fuel,3), '[kg/s]; air CC = ', round(ec.mf_hot,3), '[kg/s]; Total end of CC = ', round(ec.mf_airfuel,3),'[kg/s]')
            print('LHV fuel = ',round(ec.LHV_f,3),'m air to cool / m air core', round(ec.mr_SZair_simpl,4), round(ec.mr_SZair_simpl1, 4) )
            print('Initial estimate TPZ = ', round(ec.TPZ,3))
            print('Power: Fan = ', round(ec.W_fan,3), '[W]; LPC = ', round(ec.W_LPC,3), '[W]; HPC = ', round(ec.W_HPC,3), '[W]')
            print('LPT = ', round(ec.W_LPT,3), '[W]; HPT = ', round(ec.W_HPT,3), '[W]')
            print('Entrance of HPT: T04 = ', round(ec.T04,3), '[K]; p04 = ', round(ec.p04,3), '[Pa]')
            print('Entrance of LPT: T045 = ', round(ec.T045,3), '[K]; p045 = ', round(ec.p045,3), '[Pa]')
            print('Entrance of nozzle: T05 = ', round(ec.T05,3), '[K]; p05 = ', round(ec.p05,3), '[Pa]')
            print('Exit of nozzle: T07 = ', round(ec.T07,3), '[K]; p07 = ', round(ec.p07,3), '[Pa]; PR_cr_noz = ', ec.PR_cr_noz_core)
            print('Exit of nozzle: T8 = ', round(ec.T8,3), '[K]; p8 = ', round(ec.p8,3), '[Pa]; v8 = ', round(ec.v8,3), '[m/s]')
            print('Exit of fan: T016 = ', round(ec.T016,3), '[K]; p016 = ', round(ec.p016,3), '[Pa]; PR_cr_fan = ', ec.PR_cr_fan)
            print('Exit of fan: T18 = ', round(ec.T18,3), '[K]; p18 = ', round(ec.p18,3), '[Pa]; v18 = ', round(ec.v18,3), '[m/s]')
            print('Provided Thrust: Fan = ', round(ec.T_fan,3), '[N]; Core = ', round(ec.T_core,3), '[N]; Total = ', round(ec.T_total,3), '[N]')
            print('Thrust SFC = ', round(ec.TSFC,5), '[g/kN/s]; Equivalence ratio = ', round(ec.equivalence_ratio,4))

            amb = [['T0', round(ec.T0,3), 'K'], ['p0', round(ec.p0,3), 'Pa'], ['v0', round(ec.v0,3), 'm/s']]
            air = [['m_intake', round(ec.mf_air_init,3), 'kg/s'], ['m_hot', round(ec.mf_hot,3), 'kg/s'], ['m_cold', round(ec.mf_cold,3), 'kg/s']]
            st0 = [['T00', round(ec.T00,3), 'K'], ['p00', round(ec.p00,3), 'Pa']]
            st2 = [['T02', round(ec.T02,3), 'K'], ['p02', round(ec.p02,3), 'Pa']]
            BPR = ['BPR', round(ec.BPR, 3), '-']
            st21 = [['T021', round(ec.T021,3), 'K'], ['p021', round(ec.p021,3), 'Pa']]
            st25 = [['T025', round(ec.T025,3), 'K'], ['p025', round(ec.p025,3), 'Pa']]
            st3 = [['T03', round(ec.T03,3), 'K'], ['p03', round(ec.p03,3), 'Pa']]
            st4 = [['T04', round(ec.T04,3), 'K'], ['p04', round(ec.p04,3), 'Pa']]
            fuel = [['m_fuel', round(ec.mf_fuel,3), 'kg/s'], ['m_h2', round(ec.mf_h2,3), 'kg/s'], ['m_ker', round(ec.mf_ker,3), 'kg/s']]
            st45 = [['T045', round(ec.T045,3), 'K'], ['p045', round(ec.p045,3), 'Pa']]
            st5 = [['T05', round(ec.T05,3), 'K'], ['p05', round(ec.p05,3), 'Pa']]
            st7 = [['T07', round(ec.T07,3), 'K'], ['p07', round(ec.p07,3), 'Pa']]
            st8 = [['T8', round(ec.T8,3), 'K'], ['p8', round(ec.p8,3), 'Pa'], ['v8', round(ec.v8,3), 'm/s']]
            st16 = [['T016', round(ec.T016,3), 'K'], ['p016', round(ec.p016,3), 'Pa']]
            st18 = [['T18', round(ec.T18,3), 'K'], ['p18', round(ec.p18,3), 'Pa'], ['v18', round(ec.v18,3), 'm/s']]
            Thr = [['T_fan', round(ec.T_fan,3), 'N'], ['T_core', round(ec.T_core,3), 'N'], ['T_tot', round(ec.T_total,3), 'N'], ['TSCF', round(ec.TSFC,5), 'g/kN/s']]
            OPR = ['OPR',round(ec.OPR,3), '-']

            save_txt = amb + air + st0 + st2 + [BPR] + st21 + st25 + st3 + st4 + fuel + st45 + st5 + st7 + st8 + st16 + st18 + Thr + [OPR]
            name = a+'_'+p+'.txt'

            F = open(name,'w')
            for i in range(len(save_txt)):
                for j in range(0,3):
                    F.write(str(save_txt[i][j]) + '\t')
                F.write('\n')

            F.close()
