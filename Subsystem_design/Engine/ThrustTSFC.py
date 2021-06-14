from Subsystem_design.common_constants import Constants
from DataEngine import DataFrame
from Subsystem_design.Engine.EnergySplit import LHV_hack, ER_h2, ER_ker
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import ndimage


class SA_Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()


        PR_HPC_OG, T04_OG = self.get_data()
        print(PR_HPC_OG, T04_OG)
        plt.figure()
        self.cycle(PR_HPC_OG, T04_OG)
        plt.plot( self.T_core/1000, self.TSFC, 'go')

        # PLOT PR LINES
        self.new_PR = np.arange(7, 17.1, 1)
        # self.deltas_T = np.arange(-0.2, 1.05, 0.05)
        self.new_T = np.arange(1100,1651,50)
        for PR_new in self.new_PR: # increase one variable at a time, the remaining stay the original value
            PR_HPC = PR_new
            save_TSFC, save_var, save_T4, save_thrust = np.array([]), np.array([]), np.array([]), np.array([])

            for T_new in self.new_T: # increase each time by 1%, 2%, etc
                T04 = T_new

                self.cycle(PR_HPC, T04)
                if not m.isnan(self.TSFC):
                    save_TSFC = np.append(save_TSFC, self.TSFC)
                    save_thrust = np.append(save_thrust, self.T_core/1000) # [kN]
                    save_T4 = np.append(save_T4, T04)

            # plt.plot(save_thrust, save_TSFC, color='black', label='HPC PR ='+str(round(PR_HPC,2)))
            plt.plot(save_thrust, save_TSFC, color='black')
            plt.xlabel('Net thrust from core [kN]', fontsize=16)
            plt.ylabel('Thrust Specific Fuel Consumptino [g/kN/s]', fontsize=16)
            # plt.legend()

        # PLOT T04 LINES
        # self.deltas_PR = np.arange(-0.2, 1.05, 0.05)
        # self.deltas_T = np.arange(-0.2, 1.05, 0.2)
        self.new_PR = np.arange(7, 17.1, 1)
        self.new_T = np.arange(1100,1651,50)
        for T_new in self.new_T:  # increase one variable at a time, the remaining stay the original value
            T04 = T_new
            save_TSFC, save_var, save_PR, save_thrust = np.array([]), np.array([]), np.array([]), np.array([])

            for PR_new in self.new_PR:  # increase each time by 1%, 2%, etc
                PR_HPC = PR_new

                self.cycle(PR_HPC, T04)
                if not m.isnan(self.TSFC):
                    save_TSFC = np.append(save_TSFC, self.TSFC)
                    save_thrust = np.append(save_thrust, self.T_core / 1000)  # [kN]
                    save_PR = np.append(save_PR, PR_HPC)


            # plt.plot(save_thrust, save_TSFC, color='red', label='T04 =' + str(round(T04, 2)))
            plt.plot(save_thrust, save_TSFC, color='red')
            # plt.legend()

        plt.show()



            # # print(save_TSFC)
            # v_opt = save_var[ np.where( save_TSFC == np.min( save_TSFC ) )[0][0] ]
            # incr = self.deltas[ np.where( save_TSFC == np.min( save_TSFC ) )[0][0] ] * 100 # [%]
            # delta_TSFC = ( np.min( save_TSFC ) - save_TSFC[0] ) / save_TSFC[0] * 100
            # self.opt_vars.append( [ self.names[k+4], round(v_opt,3), round(incr, 3), 'Delta TSFC:', round(delta_TSFC,3) ] )
            # print(self.opt_vars[k])



    def get_data(self):
        data = np.array(DataFrame().neo.cruise)
        # self.names = np.array(DataFrame().neo.parameter)
        # self.names = np.delete(self.names, np.where(self.names == 'PR_LPT')[0][0])
        # self.names = np.delete(self.names, np.where(self.names == 'PR_HPT')[0][0])

        self.M0 = float(data[0])
        self.h = float(data[1])
        self.ISA_calculator(h_input=self.h) # gives self.T0, self.p0, self.rho0, self.a0
        self.T0, self.p0, self.rho0, self.a0 = self.T, self.p, self.rho, self.a
        self.v0 = self.M0 * np.sqrt(self.cp_air * self.R * self.T0)
        self.Thrust = float(data[2])
        self.A_eff = float(data[3])*self.A_fan

        self.eta_inlet = float(data[4])
        self.PR_fan = float(data[5])
        self.eta_fan = float(data[6])
        self.BPR = float(data[7])
        self.eta_LPC = float(data[8])
        self.eta_HPC = float(data[9])
        self.PR_LPC = float(data[10])
        PR_HPC = float(data[11])
        self.eta_mech = float(data[12])
        self.eta_cc = float(data[13])
        self.PR_cc = float(data[14])
        T04 = float(data[15])
        self.eta_LPT = float(data[16])
        self.eta_HPT = float(data[17])
        self.eta_nozzle = float(data[20])
        self.PR_noz_core = float(data[21])
        self.PR_noz_fan = float(data[22])
        self.mr_h2 = float(data[23])
        self.m_ker = float(data[24])
        self.ER_h2 = float(data[25])
        self.ER_ker = float(data[26])
        self.LHV_f = float(data[27])

        self.var_OG = np.array([PR_HPC, T04])
        self.data = data
        return PR_HPC, T04


    # def get_vals(self, low, up):
    #     self.vals = np.linspace(low, up, 100)


    def cycle(self, PR_HPC, T04):

        '''
        INDEX
        0: eta_inlet
        1: PR_fan
        2: eta_fan
        3: BPR
        4: eta_LPC
        5: eta_HPC
        6: PR_LPC
        7: PR_HPC
        8: eta_mech
        9: eta_cc
        10: PR_cc
        11: T04
        12: eta_LPT
        13: eta_HPT
        16-->14: eta_nozzle
        17-->15: PR_noz_core
        18-->16: PR_noz_fan
        '''

        self.v0 = self.M0 * np.sqrt(self.k_air * self.R * self.T0)
        self.rho0 = self.p0 / (self.R * self.T0)

        self.mf_air_init = self.rho0 * self.A_eff * self.v0

        # Total temperature and pressure at inlet
        self.T00 = self.T0 * (1 + (self.k_air - 1) / 2 * self.M0 ** 2)
        self.p00 = self.p0 * (1 + (self.k_air - 1) / 2 * self.M0 ** 2) ** (self.k_air / (self.k_air - 1))

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0 * (1 + self.eta_inlet * (self.k_air - 1) / 2 * self.M0 ** 2) ** (
                    self.k_air / (self.k_air - 1))

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + (self.T02 / self.eta_fan) * (self.PR_fan ** ((self.k_air - 1) / self.k_air) - 1)
        self.p021 = self.p02 * self.PR_fan

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init / (self.BPR + 1)
        self.mf_cold = self.mf_air_init * self.BPR / (self.BPR + 1)

        # Further on the bypass duct
        self.T016 = self.T021
        self.p016 = self.p021  # * self.PR_noz_fan

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + (self.T021 / self.eta_LPC) * (self.PR_LPC ** ((self.k_air - 1) / self.k_air) - 1)
        self.p025 = self.p021 * self.PR_LPC

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + (self.T02 / self.eta_HPC) * (PR_HPC ** ((self.k_air - 1) / self.k_air) - 1)
        self.p03 = self.p025 * PR_HPC
        self.OPR = self.p03 / self.p02  # Overall Pressure Ratio

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.cp_gas * (T04 - self.T03)) / (self.LHV_f * 10 ** 6 * self.eta_cc)
        self.mf_airfuel = self.mf_hot + self.mf_fuel  # at the end of the cc
        self.mf_h2 = self.mf_fuel * self.ER_h2
        self.mf_ker = self.mf_fuel * self.ER_ker

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * self.PR_cc

        # Power to drive fan, LPC, HPC, HPT, LPT [W]
        self.W_fan = self.mf_air_init * self.cp_air * (self.T021 - self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025 - self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03 - self.T025)
        self.W_HPT = self.W_HPC / self.eta_mech
        self.W_LPT = (self.W_fan + self.W_LPC) / self.eta_mech

        # Exit of HPT - Entrance of LPT
        self.T045 = T04 - self.W_HPT / (self.mf_airfuel * self.cp_gas)
        self.p045 = self.p04 * (1 - (1 - self.T045 / T04) / self.eta_HPT) ** (self.k_gas / (self.k_gas - 1))

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * (1 - (1 - self.T05 / self.T045) / self.eta_LPT) ** (self.k_gas / (self.k_gas - 1))

        # Nozzle
        self.T07 = self.T05
        self.p07 = self.p05  # * self.PR_noz_core

        # Is the nozzle chocked?
        self.PR_cr_noz_core = 1 / (
                    (1 - (self.k_gas - 1) / (self.k_gas + 1) / self.eta_nozzle) ** (self.k_gas / (self.k_gas - 1)))

        # Exit of the nozzle
        if self.p07 / self.p0 > self.PR_cr_noz_core:
            self.TR_cr_noz_core = (self.k_gas + 1) / 2
            self.T8 = self.T07 / self.TR_cr_noz_core
            self.p8 = self.p07 / self.PR_cr_noz_core
            self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
            self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
            self.T_core = self.mf_airfuel * (self.v8 - self.v0) + self.A8 * (self.p8 - self.p0)  # [N]

        elif self.p07 / self.p0 <= self.PR_cr_noz_core:
            self.p8 = self.p0
            self.T8 = self.T07 * (1 - self.eta_nozzle * (1 - (self.p8 / self.p07) ** ((self.k_gas - 1) / self.k_gas)))
            self.v8 = np.sqrt(2 * self.cp_gas * (self.T07 - self.T8))
            self.T_core = self.mf_airfuel * (self.v8 - self.v0)  # [N]

        # Is the fan chocked?
        self.PR_cr_fan = 1 / (
                    (1 - (self.k_air - 1) / (self.eta_nozzle * (self.k_air + 1))) ** (self.k_air / (self.k_air - 1)))

        # Exit of bypassed air
        if self.p016 / self.p0 > self.PR_cr_fan:
            self.TR_cr_bypassed = (self.k_air + 1) / 2
            self.T18 = self.T016 / self.TR_cr_bypassed
            self.p18 = self.p016 / self.PR_cr_fan
            self.v18 = np.sqrt(self.k_air * self.R * self.T18)
            self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
            self.T_fan = self.mf_cold * (self.v18 - self.v0) + self.A18 * (self.p18 - self.p0)  # [N]

        elif self.p016 / self.p0 <= self.PR_cr_fan:
            self.p18 = self.p0
            self.T18 = self.T016 - self.T016 * self.eta_fan * (
                        1 - (self.p18 / self.p016) ** ((self.k_air - 1) / self.k_air))
            self.v18 = np.sqrt(2 * self.cp_air * (self.T016 - self.T18))
            self.T_fan = self.mf_cold * (self.v18 - self.v0)  # [N]

        self.T_total = self.T_fan + self.T_core  # [N]
        self.TSFC = self.mf_fuel / (self.T_total) * 10 ** 3  # [g/kN/s] (kg*10**3)/s / (N*10**-3)


if __name__=='__main__':
    c = Constants()
    sa = SA_Engine_Cycle()
