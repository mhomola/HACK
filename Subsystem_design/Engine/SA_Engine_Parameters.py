# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 09:39:04 2021

@author: sarar
"""

# TODO
# Find estimates for all constants --> Aero Engine Technology course
# Find altitude and speed of each phase
# Check why when the nozzle is not chocked T8 > T05 (hypothesis: error with p05)
# Fuel flow of each fuel

# DONE
# Determine LHV_fuel according to mass ratio of fuels
# Determine atmospheric conditions
# Determine stoichiometric ratio --> find equivalence ratio


from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EnergySplit import LHV_hack, ER_h2, ER_ker
import numpy as np


class SA_Engine_Cycle(Constants):
    def __init__(self, i):
        super().__init__()

        self.data(self.phases[i])
        """ I want to perform a sensitivity analysis on the following variables """
        var_orig = [self.eta_inlet, self.PR_fan, self.eta_fan, self.BPR, self.eta_LPC, self.eta_HPC, self.PR_LPC, self.PR_HPC,
                             self.PR_cc, self.eta_LPT, self.eta_HPT, self.PR_LPT, self.PR_HPT, self.eta_nozzle,
                             self.PR_noz_core]  #, self.D_fan_eff ])

        # PR_fan, BPR, PR_LPC, PR_HPC

        var_names = [ 'eta_inlet','PR_fan','eta_fan','BPR','eta_LPC','eta_HPC','PR_LPC','PR_HPC','PR_cc','eta_LPT', 'eta_HPT','PR_LPT','PR_HPT','eta_nozzle','PR_noz_core' ]

        if self.phases[i] == 'takeoff':
            var_orig[0] = self.PR_inlet
            var_names[0] = 'PR_inlet'

        # var = np.array(var_orig)
        for k in range(len(var_orig)):
            var = np.array(var_orig)
        ## ---------- GET ARRAY OF DELTAS ------------ ##
            if "eta" in var_names[k]:
                vals = np.arange(var_orig[k], 1, 0.005)
                deltas = vals - var_orig[k]
            else:
                deltas = np.arange(0.0, 0.31, 0.01)
        ## ---------------------------------------------- ##

            store_var, store_mf, store_TSFC = np.array([]), np.array([]), np.array([])
            # var = np.array(var_orig)

            for j in deltas:
                var[k] = var_orig[k] * (1 + j)
                store_var = np.append(store_var, var[k])

                # Total temperature and pressure at inlet
                self.T00 = self.T0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 )
                self.p00 = self.p0[i] * ( 1 + (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )

                # Entrance of the fan
                self.T02 = self.T00

                if self.M0[i] == 0.:
                    self.p02 = self.p0[i] * self.PR_inlet
                else:
                    self.p02 = self.p0[i] * (1 + var[0] * (self.k_air-1)/2 * self.M0[i]**2 ) ** (self.k_air / (self.k_air-1) )

                # Exit of the fan - Entrance of LPC
                self.T021 = self.T02 + ( self.T02/var[2] ) * ( var[1] ** ( (self.k_air-1)/self.k_air ) - 1 )
                self.p021 = self.p02 * var[1]

                # Hot and cold mass flow of air
                self.mf_hot = self.mf_air_init / (var[3]+1)
                self.mf_cold = self.mf_air_init * var[3] / (var[3]+1)

                # Exit of LPC - Entrance of HPC
                self.T025 = self.T021 + ( self.T021/var[4] ) * ( var[6] ** ( (self.k_air-1)/self.k_air ) - 1 )
                self.p025 = self.p021 * var[6]

                # Exit of HPC - Entrance of cc
                self.T03 = self.T025 + ( self.T025/var[5] ) * ( var[7] ** ( (self.k_air-1)/self.k_air ) - 1 )
                self.p03 = self.p025 * var[7]
                self.OPR = self.p03 / self.p02         # Overall Pressure Ratio

                # Exit of cc - Entrance of HPT
                self.mf_fuel = (self.mf_hot * self.cp_gas * (self.T04-self.T03)) / (self.LHV_f[i]*10**6 * self.eta_cc)
                self.mf_airfuel = self.mf_hot + self.mf_fuel # at the end of the cc
                store_mf = np.append(store_mf, self.mf_fuel)

                self.mf_h2 = self.mf_fuel * ER_h2[i]   #energy ratio H2
                self.mf_ker = self.mf_fuel * ER_ker[i]

                # T04 = 1500 [K], is given
                self.p04 = self.p03 * var[8]

                # Power to drive fan, LPC, HPC, HPT, LPT [W]
                # self.W_fan = self.mf_air_init[i] * self.cp_air * (self.T021-self.T00)
                # self.W_LPC = self.mf_hot * self.cp_air * (self.T025-self.T021)
                # self.W_HPC = self.mf_hot * self.cp_air * (self.T03-self.T025)
                # self.W_HPT = self.W_HPC / self.eta_mech
                # self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech

                # Exit of HPT - Entrance of LPT
                self.T045 = var[8] + (var[8]/var[10]) * ( (1/var[12])**((self.k_gas-1)/self.k_gas) - 1 )
                self.p045 = self.p04 / var[12]
                # self.T045 = self.T04 - self.W_HPT / ( self.mf_airfuel * self.cp_gas )
                # self.p045 = self.p04 * ( 1 - ( 1 - self.T045/self.T04 ) / self.eta_HPT ) ** ( self.k_gas / (self.k_gas-1) )


                # Exit of LPT - Entrance of nozzle
                self.T05 = self.T045 + (self.T045/var[9]) * ( (1/var[11])**((self.k_gas-1)/self.k_gas) - 1 )
                self.p05 = self.p045 / var[11]
                # self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
                # self.p05 = self.p045 * ( 1 - ( 1 - self.T05/self.T045 ) / self.eta_LPT ) ** ( self.k_gas / (self.k_gas-1) )

                # Nozzle
                self.T07 = self.T05
                self.p07 = self.p05 * var[14]

                # Is the nozzle chocked?
                self.PR_cr_nozzle =  ( 1 - (self.k_gas-1)/((self.k_gas+1)*self.eta_nozzle)) ** (-self.k_gas / (self.k_gas-1))

                # Exit of the nozzle
                if self.p07/self.p0[i] > self.PR_cr_noz_core:
                    self.TR_cr_noz_core = (self.k_gas + 1) / 2
                    self.T8 = self.T07 / self.TR_cr_noz_core
                    self.p8 = self.p07 / self.PR_cr_noz_core
                    self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
                    self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
                    self.T_core = self.mf_airfuel * (self.v8 - self.v0[i]) + self.A8 * (self.p8 - self.p0[i])  # [N]

                elif self.p07/self.p0[i] <= self.PR_cr_noz_core:
                    self.p8 = self.p0[i]
                    self.T8 = self.T07 * ( 1 - var[13] * ( 1 - (self.p8/self.p05) ** ( (self.k_gas-1)/self.k_gas ) ) )
                    self.v8 = np.sqrt( 2 * self.cp_gas * (self.T07 - self.T8) )
                    self.T_core = self.mf_airfuel * ( self.v8 - self.v0[i] )  # [N]

                # Is the fan chocked?
                self.T016 = self.T021
                self.p016 = self.p021 * self. PR_noz_fan
                self.PR_cr_fan = 1 / ( ( 1 - (self.k_air-1)/(var[13]*(self.k_air+1)) ) ** (self.k_air / (self.k_air-1)) )

                # Exit of bypassed air
                if self.p016/self.p0[i] > self.PR_cr_fan:
                    self.TR_cr_bypassed = (self.k_air + 1) / 2
                    self.T18 = self.T016 / self.TR_cr_bypassed
                    self.p18 = self.p016 / self.PR_cr_fan
                    self.v18 = np.sqrt(self.k_air * self.R * self.T18)
                    self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
                    self.T_fan = self.mf_cold * (self.v18 - self.v0[i]) + self.A18 * (self.p18 - self.p0[i])  # [N]

                elif self.p016/self.p0[i] <= self.PR_cr_fan:
                    self.p18 = self.p0[i]
                    self.T18 = self.T016 - self.T016 * var[2] * ( 1 - (self.p18/self.p021) ** ( (self.k_air-1)/self.k_air ) )
                    self.v18 = np.sqrt( 2 * self.cp_air * (self.T016 - self.T18) )
                    self.T_fan = self.mf_cold * ( self.v18 - self.v0[i] ) # [N]


                self.T_total = self.T_fan + self.T_core # [N]
                self.TSFC = self.mf_fuel/(self.T_total * 10**(-3))
                store_TSFC = np.append(store_TSFC, self.TSFC ) # [kg/s/kN]

            a = np.where( store_TSFC == np.min(store_TSFC) )[0][0]
            # print('\n',var_names, '\n', store_var,'\n',store_mf)
            incr = round( (store_var[a] - var_orig[k]) / var_orig[k] * 100, 1)
            print(var_names[k], "=", store_var[a], " (", incr ,"% higher) leads to the min TSFC, equal to ", store_TSFC[a], '; Thrust = ', self.T_total )

    def data(self, phase):
        self.LHV_f = LHV_hack

        if phase == 'cruise':
            self.eta_inlet = 0.9208608681597723  # calculated
            self.PR_fan = 1.4206
            self.eta_fan = 0.90445
            self.BPR = 11.24426
            self.eta_LPC = 0.90019
            self.eta_HPC = 0.95469  # calculated # 0.91449 given
            self.PR_LPC = 2.69419
            self.PR_HPC = 9.73784
            self.eta_mech = 0.9  # not used
            self.eta_cc = 0.995  # that of Leap-1B, assumed
            self.PR_cc = 0.9395309126896629
            self.T04 = 1459.30433 + 20  # [K], increases 10K per year: https://www.researchgate.net/publication/221919089_Future_Aero_Engine_Designs_An_Evolving_Vision
            self.eta_LPT = 0.9405
            self.eta_HPT = 0.9328  # computed # 0.91898 #(given)
            self.PR_LPT = 7.9204
            self.PR_HPT = 3.81933
            self.eta_nozzle = 0.981797  # 1.0737340755627587 (computed) # previous assumption: 0.98
            self.PR_noz_core = 0.985443  # Between stations 5 and 7
            self.PR_cr_noz_core = 1.653828  # computed
            self.PR_noz_fan = 0.987444  # Between stations 21 and 16

            self.D_fan_eff = 1.787  # Effective fan diameter for air ingestion      [m]
            self.A_fan_eff = np.pi * self.D_fan_eff ** 2 / 4  # Area of the fan                               [m2]
            self.mf_air_init = self.rho0[4] * self.A_fan_eff * self.v0[4]

        elif phase == 'takeoff':
            self.eta_inlet = 0.92  # assumed
            self.PR_inlet = 0.9699975  # calculated
            self.PR_fan = 1.4
            self.eta_fan = 0.93
            self.BPR = 11.1
            self.eta_LPC = 0.92
            self.eta_HPC = 0.97735  # calculated # 0.92 given
            self.PR_LPC = 1.99241
            self.PR_HPC = 11.92998521
            self.eta_mech = 0.9  # not used
            self.eta_cc = 0.995  # that of Leap-1B, assumed
            self.PR_cc = 0.94
            self.T04 = 1630.39416  # [K]
            self.eta_LPT = 0.94
            self.eta_HPT = 0.94  # computed # 0.91898 #(given)
            self.PR_LPT = 6.36217496  # 4.835766
            self.PR_HPT = 3.82808
            self.eta_nozzle = 0.98  # assumed
            self.PR_noz_core = 0.99  # Between stations 5 and 7
            self.PR_cr_noz_core = 1.2380755  # computed
            self.PR_noz_fan = 0.99  # Between stations 21 and 16
            self.mf_air_init = 510  # [kg/s]


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
    7 - a bit more further in the nozzle
    18 - exit of bypassed air
    8 - exit of nozzle from core
'''


if __name__ == '__main__':
    c = Constants()
    # index = [2,4] # takeoff and cruise
    # for i in index:
    #     ec = SA_Engine_Cycle(i=i)
    ec = SA_Engine_Cycle(i=4)



        # print('\nInlet: T0 = ', c.T0[i], '[K]; p0 = ', c.p0[i], '[Pa]; v0 = ', c.v0[i], '[m/s]')
        # print('T00 = ', ec.T00, '[K]; p00 = ', ec.p00, '[Pa]')
        # print('Entrance of fan: T02 = ', ec.T02, '[K]; p02 = ', ec.p02, '[Pa]')
        # print('Entrance of LPC: T021 = ', ec.T021, '[K]; p021 = ', ec.p021, '[Pa]')
        # print('Mass flow of air: Total = ', ec.mf_air_init, '[kg/s]; Core = ', ec.mf_hot, '[kg/s]; Bypassed = ', ec.mf_cold,'[kg/s]')
        # print('Entrance of HPC: T025 = ', ec.T025, '[K]; p025 = ', ec.p025, '[Pa]')
        # print('Entrance of CC: T03 = ', ec.T03, '[K]; p03 = ', ec.p03, '[Pa]')
        # print('Mass flow CC: Fuel = ', ec.mf_fuel, '[kg/s]; air CC = ', ec.mf_hot, '[kg/s]; Total end of CC = ', ec.mf_airfuel,'[kg/s]')
        #
        #
        # print('Hydrogen = ', ec.mf_h2, '[kg/s]; Kerosene = ', ec.mf_ker,'[kg/s]')
        #
        # print('Entrance of HPT: T04 = ', ec.T04, '[K]; p04 = ', ec.p04, '[Pa]')
        # print('Entrance of LPT: T045 = ', ec.T045, '[K]; p045 = ', ec.p045, '[Pa]')
        # print('Entrance of nozzle: T05 = ', ec.T05, '[K]; p05 = ', ec.p05, '[Pa]')
        # print('Exit of nozzle: T07 = ', ec.T07, '[K]; p07 = ', ec.p07, '[Pa]')
        # print('Exit of nozzle: T8 = ', ec.T8, '[K]; p8 = ', ec.p8, '[Pa]; v8 = ', ec.v8, '[m/s]')
        # print('Exit of fan: T016 = ', ec.T016, '[K]; p016 = ', ec.p016, '[Pa]')
        # print('Exit of fan: T18 = ', ec.T18, '[K]; p18 = ', ec.p18, '[Pa]; v18 = ', ec.v18, '[m/s]')
        # print('Provided Trhust: Fan = ', ec.T_fan, '[N]; Core = ', ec.T_core, '[N]; Total = ', ec.T_total, '[N]')
