from Subsystem_design.common_constants import Constants
from DataEngine import DataFrame
from Subsystem_design.Engine.EnergySplit import LHV_hack, ER_h2, ER_ker
import numpy as np
import math as m
import matplotlib.pyplot as plt


class SA_Engine_Cycle(Constants):
    def __init__(self):
        super().__init__()
        self.opt_vars = list()
        # var_OG = np.array([self.PR_fan, self.BPR, self.PR_LPC, self.PR_HPC, self.PR_LPT, self.PR_HPT])
        # print(var_OG)
        # self.var_name = np.array(['PR_fan', 'BPR', 'PR_LPC', 'PR_HPC', 'PR_LPT', 'PR_HPT'])
        # self.up_bound = np.array([ 2, 15.808, 7, 15, 0.3, 0.7 ])
        self.data()

        for k in range(len(self.var_OG)): # increase one variable at a time, the remaining stay the original value
            var_loop = self.var_OG.copy()
            # print('\nvar:', self.var_name[k])
            # self.get_vals(self.var_OG[k], self.up_bound[k])
            save_TSFC, save_var = np.array([]), np.array([])
            save_T45, save_T5 = np.array([]), np.array([])
            save_thrust = np.array([])


            if 'eta' in self.names[k+4] or 'PR_noz' in self.names[k+4] or 'PR_cc' in self.names[k+4]:
                self.deltas = (np.linspace(0.9*self.var_OG[k], 1, 1000) - self.var_OG[k]) / self.var_OG[k]
            else:
                self.deltas = np.arange(-0.2, 1.05, 0.05)

            for d in self.deltas: # increase each time by 1%, 2%, etc
                var_loop[k] = self.var_OG[k].copy() * (d + 1)
                save_var = np.append(save_var, var_loop[k])

                self.cycle(var_loop)
                if not m.isnan(self.TSFC):
                    save_TSFC = np.append(save_TSFC, self.TSFC)
                else:
                    save_TSFC = np.append(save_TSFC, 1000) # add a dummy number that will never be the minimum of the TSFC array

                if self.names[k+4] == 'eta_HPC' or  self.names[k+4] == 'eta_LPC':
                    save_T45 = np.append(save_T45, self.p045)
                    save_T5 = np.append(save_T5, self.p05)
                    save_thrust = np.append(save_thrust, self.T_core)

            # print(save_TSFC)
            v_opt = save_var[ np.where( save_TSFC == np.min( save_TSFC ) )[0][0] ]
            incr = self.deltas[ np.where( save_TSFC == np.min( save_TSFC ) )[0][0] ] * 100 # [%]
            delta_TSFC = ( np.min( save_TSFC ) - save_TSFC[0] ) / save_TSFC[0] * 100
            self.opt_vars.append( [ self.names[k+4], round(v_opt,3), round(incr, 3), 'Delta TSFC:', round(delta_TSFC,3) ] )
            print(self.opt_vars[k])

            self.plot(save_var, save_TSFC, k, flag=False)
            if self.names[k + 4] == 'eta_HPC' or self.names[k + 4] == 'eta_LPC':
                self.plot(save_var, save_T45, k, flag='p045')
                self.plot(save_var, save_thrust, k, flag='Thrust core')



    def data(self):
        data = np.array(DataFrame().neo.cruise)
        self.names = np.array(DataFrame().neo.parameter)
        self.names = np.delete(self.names, np.where(self.names == 'PR_LPT')[0][0])
        self.names = np.delete(self.names, np.where(self.names == 'PR_HPT')[0][0])

        self.M0 = float(data[0])
        self.h = float(data[1])
        self.ISA_calculator(h_input=self.h) # gives self.T0, self.p0, self.rho0, self.a0
        self.T0, self.p0, self.rho0, self.a0 = self.T, self.p, self.rho, self.a
        self.v0 = self.M0 * np.sqrt(self.cp_air * self.R * self.T0)
        self.Thrust = float(data[2])
        self.A_eff = float(data[3])*self.A_fan

        eta_inlet_OG = float(data[4])
        PR_fan_OG = float(data[5])
        eta_fan_OG = float(data[6])
        BPR_OG = float(data[7])
        eta_LPC_OG = float(data[8])
        eta_HPC_OG = float(data[9])
        PR_LPC_OG = float(data[10])
        PR_HPC_OG = float(data[11])
        eta_mech_OG = float(data[12])
        eta_cc_OG = float(data[13])
        PR_cc_OG = float(data[14])
        T04_OG = float(data[15])
        eta_LPT_OG = float(data[16])
        eta_HPT_OG = float(data[17])
        eta_nozzle_OG = float(data[20])
        PR_noz_core_OG = float(data[21])
        PR_noz_fan_OG = float(data[22])
        self.mr_h2 = float(data[23])
        self.m_ker = float(data[24])
        self.ER_h2 = float(data[25])
        self.ER_ker = float(data[26])
        self.LHV_f = float(data[27])

        self.var_OG = np.array([eta_inlet_OG, PR_fan_OG, eta_fan_OG, BPR_OG,
         eta_LPC_OG, eta_HPC_OG, PR_LPC_OG, PR_HPC_OG, eta_mech_OG, eta_cc_OG,
         PR_cc_OG, T04_OG, eta_LPT_OG, eta_HPT_OG,
         eta_nozzle_OG, PR_noz_core_OG, PR_noz_fan_OG])
        self.data = data


    # def get_vals(self, low, up):
    #     self.vals = np.linspace(low, up, 100)


    def cycle(self, var_loop):

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

        # Total temperature and pressure at inlet
        self.v0 = self.M0 * np.sqrt(self.k_air * self.R * self.T0)
        self.rho0 = self.p0 / (self.R * self.T0)

        self.mf_air_init = self.rho0 * self.A_eff * self.v0

        # Total temperature and pressure at inlet
        self.T00 = self.T0 * (1 + (self.k_air - 1) / 2 * self.M0 ** 2)
        self.p00 = self.p0 * (1 + (self.k_air - 1) / 2 * self.M0 ** 2) ** (self.k_air / (self.k_air - 1))

        # Entrance of the fan
        self.T02 = self.T00
        self.p02 = self.p0 * (1 + var_loop[0] * (self.k_air - 1) / 2 * self.M0 ** 2) ** (
                    self.k_air / (self.k_air - 1))

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + (self.T02 / var_loop[2]) * (var_loop[1] ** ((self.k_air - 1) / self.k_air) - 1)
        self.p021 = self.p02 * var_loop[1]

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init / (var_loop[3] + 1)
        self.mf_cold = self.mf_air_init * var_loop[3] / (var_loop[3] + 1)

        # Further on the bypass duct
        self.T016 = self.T021
        self.p016 = self.p021  * var_loop[16]

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + (self.T021 / var_loop[4]) * (var_loop[6] ** ((self.k_air - 1) / self.k_air) - 1)
        self.p025 = self.p021 * var_loop[6]

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + (self.T025 / var_loop[5]) * (var_loop[7] ** ((self.k_air - 1) / self.k_air) - 1)
        self.p03 = self.p025 * var_loop[7]
        self.OPR = self.p03 / self.p02  # Overall Pressure Ratio

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.cp_gas * (var_loop[11] - self.T03)) / (self.LHV_f * 10 ** 6 * var_loop[9])
        self.mf_airfuel = self.mf_hot + self.mf_fuel  # at the end of the cc
        self.mf_h2 = self.mf_fuel * self.ER_h2
        self.mf_ker = self.mf_fuel * self.ER_ker

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * var_loop[10]

        # Power to drive fan, LPC, HPC, HPT, LPT [W]
        self.W_fan = self.mf_air_init * self.cp_air * (self.T021 - self.T00)
        self.W_LPC = self.mf_hot * self.cp_air * (self.T025 - self.T021)
        self.W_HPC = self.mf_hot * self.cp_air * (self.T03 - self.T025)
        self.W_HPT = self.W_HPC / var_loop[8]
        self.W_LPT = (self.W_fan + self.W_LPC) / var_loop[8]

        # Exit of HPT - Entrance of LPT
        self.T045 = var_loop[11] - self.W_HPT / (self.mf_airfuel * self.cp_gas)
        self.p045 = self.p04 * (1 - (1 - self.T045 / var_loop[11]) / var_loop[13]) ** (self.k_gas / (self.k_gas - 1))

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T045 - self.W_LPT / (self.mf_airfuel * self.cp_gas)
        self.p05 = self.p045 * (1 - (1 - self.T05 / self.T045) / var_loop[12]) ** (self.k_gas / (self.k_gas - 1))

        # Nozzle
        self.T07 = self.T05
        self.p07 = self.p05  * var_loop[15]

        # Is the nozzle chocked?
        self.PR_cr_noz_core = 1 / (
                    (1 - (self.k_gas - 1) / (self.k_gas + 1) / var_loop[14]) ** (self.k_gas / (self.k_gas - 1)))

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
            self.T8 = self.T07 * (1 - var_loop[14] * (1 - (self.p8 / self.p07) ** ((self.k_gas - 1) / self.k_gas)))
            self.v8 = np.sqrt(2 * self.cp_gas * (self.T07 - self.T8))
            self.T_core = self.mf_airfuel * (self.v8 - self.v0)  # [N]

        # Is the fan chocked?
        self.PR_cr_fan = 1 / (
                    (1 - (self.k_air - 1) / (var_loop[14] * (self.k_air + 1))) ** (self.k_air / (self.k_air - 1)))

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
            self.T18 = self.T016 - self.T016 * var_loop[14] * (
                        1 - (self.p18 / self.p016) ** ((self.k_air - 1) / self.k_air))
            self.v18 = np.sqrt(2 * self.cp_air * (self.T016 - self.T18))
            self.T_fan = self.mf_cold * (self.v18 - self.v0)  # [N]

        self.T_total = self.T_fan + self.T_core  # [N]
        self.TSFC = self.mf_fuel / (self.T_total * 10 ** (-3))  # [g/kN/s]


    def plot(self, save_var, save_TSFC, k, flag):
        i = np.where(save_TSFC == 1000)
        if len(i[0]) != 0:
            save_var = np.delete(save_var, i)
            save_TSFC = np.delete(save_TSFC, i)

        plt.figure()
        plt.gcf().canvas.set_window_title(self.names[k+4])
        plt.plot(save_var, save_TSFC)
        plt.xlabel(self.names[k + 4], fontsize=20)
        if flag == False:
            plt.ylabel('TSFC [g/kN/s]', fontsize=20)
        else:
            plt.ylabel(flag+'[K]', fontsize=20)

        plt.savefig('C:\\Users\\sarar\\OneDrive\\Ambiente de Trabalho\\Folders\\DELFT\\3rd year\\DSE\\HACK\\Subsystem_design\\Engine\\Parameters\\'+self.names[k + 4])
        plt.show()

if __name__=='__main__':
    c = Constants()
    sa = SA_Engine_Cycle()

