from Subsystem_design.common_constants import Constants
from Subsystem_design.Engine.EnergySplit import LHV_hack, ER_h2, ER_ker
import numpy as np


class SA_Engine_Cycle(Constants):
    def __init__(self, i):
        super().__init__()
        self.opt_vars = list()
        self.data_OG(self.phases[i])
        var_OG = np.array([self.PR_fan, self.BPR, self.PR_LPC, self.PR_HPC, self.PR_LPT, self.PR_HPT])
        print(var_OG)
        self.var_name = np.array(['PR_fan', 'BPR', 'PR_LPC', 'PR_HPC', 'PR_LPT', 'PR_HPT'])
        self.up_bound = np.array([ 2, 15.808, 7, 15, 0.3, 0.7 ])

        for k in range(len(var_OG)):
            print('\nvar:', self.var_name[k])
            var_loop = var_OG.copy()
            self.get_vals(var_OG[k], self.up_bound[k])
            save_TSFC = list()
            print(var_OG)

            for v in self.vals:
                var_loop[k] = v
                print(var_OG)

                PR_fan, BPR, PR_LPC, PR_HPC, PR_LPT, PR_HPT = \
                    var_loop[0], var_loop[1], var_loop[2], var_loop[3], var_loop[4], var_loop[5]

                self.cycle(i, PR_fan, BPR, PR_LPC, PR_HPC, PR_LPT, PR_HPT)
                save_TSFC.append(self.TSFC)

            TSFC_array = np.array(save_TSFC)
            # print(save_TSFC)
            v_opt = self.vals[ np.where( TSFC_array == np.min( TSFC_array ) )[0][0] ]
            incr = (v_opt - var_OG[k]) / var_OG[k] * 100
            self.opt_vars.append( [ self.var_name[k], round(v_opt,3), round(incr, 3) ] )









    def data_OG(self, phase):
        if phase == 'cruise':
            self.eta_inlet = 0.94  # 0.92086 to 0.94
            self.PR_fan = 1.4206
            self.eta_fan = 0.915 # 0.904 to 0.915
            self.BPR = 11.24426
            self.eta_LPC = 0.92 # 0.90 to 0.92
            self.eta_HPC = 0.975  # 0.955 to 0.975
            self.PR_LPC = 2.69419
            self.PR_HPC = 9.73784
            self.eta_cc = 0.9995  # 0.995 to 0.9995
            self.PR_cc = 0.95 # 0.9395 to 0.95
            self.T04 = 1459.30433 + 20  # [K], increases 10K per year: https://www.researchgate.net/publication/221919089_Future_Aero_Engine_Designs_An_Evolving_Vision
            self.eta_LPT = 0.96 # 0.9405 to 0.96
            self.eta_HPT = 0.95 # 0.9328 to 0.945
            self.PR_LPT = 1/7.9204
            self.PR_HPT = 1/3.81933
            self.eta_nozzle = 0.981797  # 0.981 to 0.985
            self.PR_noz_core = 0.985443  # 0.985 to 0.99
            self.PR_noz_fan = 0.987444  # 0.987 to 0.99
            self.mf_air_init = self.rho0 * self.A_fan_eff * self.v0

        elif phase == 'takeoff':
            self.eta_inlet = 0.92
            self.PR_inlet = 0.9699975
            self.PR_fan = 1.4
            self.eta_fan = 0.93
            self.BPR = 11.1
            self.eta_LPC = 0.92
            self.eta_HPC = 0.97735
            self.PR_LPC = 1.99241
            self.PR_HPC = 11.92998521
            self.eta_cc = 0.9995  # from 0.995 to 0.9995
            self.PR_cc = 0.96 # from 0.94 to 0.96
            self.T04 = 1630.39416 + 20  # [K]
            self.eta_LPT = 0.94
            self.eta_HPT = 0.94  # computed # 0.91898 #(given)
            self.PR_LPT = 1/6.36217496  # 4.835766
            self.PR_HPT = 1/3.82808
            self.eta_nozzle = 0.98  # assumed
            self.PR_noz_core = 0.99  # Between stations 5 and 7
            self.PR_noz_fan = 0.99  # Between stations 21 and 16
            self.mf_air_init = 510  # [kg/s] from Arvind's document


        self.LHV_f = np.array([self.LHV_ker] * len(self.phases))  # [MJ/kg]

    def get_vals(self, low, up):
        self.vals = np.linspace(low, up, 100)


    def cycle(self, i, PR_fan, BPR, PR_LPC, PR_HPC, PR_LPT, PR_HPT):

        # Total temperature and pressure at inlet
        self.T00 = self.T0[i] * (1 + (self.k_air - 1) / 2 * self.M0[i] ** 2)
        self.p00 = self.p0[i] * (1 + (self.k_air - 1) / 2 * self.M0[i] ** 2) ** (self.k_air / (self.k_air - 1))

        # Entrance of the fan
        self.T02 = self.T00

        if self.M0[i] == 0.:
            self.p02 = self.p0[i] * self.PR_inlet
        else:
            self.p02 = self.p0[i] * (1 + self.eta_inlet * (self.k_air - 1) / 2 * self.M0[i] ** 2) ** (
                        self.k_air / (self.k_air - 1))

        # Exit of the fan - Entrance of LPC
        self.T021 = self.T02 + (self.T02 / self.eta_fan) * (PR_fan ** ((self.k_air - 1) / self.k_air) - 1)
        self.p021 = self.p02 * PR_fan

        # Hot and cold mass flow of air
        self.mf_hot = self.mf_air_init / (BPR + 1)
        self.mf_cold = self.mf_air_init * BPR / (BPR + 1)

        # Further on the bypass duct
        self.T016 = self.T021
        self.p016 = self.p021 * self.PR_noz_fan

        """ USE THIS WHEN WE HAVE THE REQUIRED THRUST SETTINGS """
        # if self.v0[i] == 0.:
        #     self.T_fan_tot = (104340.3935 / 12401.48002) * T_tot # T_tot from Elena's program
        #     self.p18 = self.p0[i]
        #     self.T18 = self.T016 - self.T016 * self.eta_fan * ( 1 - (self.p18/self.p016) ** ( (self.k_air-1)/self.k_air ) )
        #     self.v18 = np.sqrt( 2 * self.cp_air * (self.T016 - self.T18) )
        #     self.mf_cold = self.T_fan / self.v18
        #     self.mf_air_init = (self.BPR+1)/self.BPR * self.mf_cold
        #     self.mf_hot = self.mf_air_init / (self.BPR + 1)

        # Exit of LPC - Entrance of HPC
        self.T025 = self.T021 + (self.T021 / self.eta_LPC) * (PR_LPC ** ((self.k_air - 1) / self.k_air) - 1)
        self.p025 = self.p021 * PR_LPC

        # Exit of HPC - Entrance of cc
        self.T03 = self.T025 + (self.T025 / self.eta_HPC) * (PR_HPC ** ((self.k_air - 1) / self.k_air) - 1)
        self.p03 = self.p025 * PR_HPC
        self.OPR = self.p03 / self.p02  # Overall Pressure Ratio

        # Exit of cc - Entrance of HPT
        self.mf_fuel = (self.mf_hot * self.cp_gas * (self.T04 - self.T03)) / (self.LHV_f[i] * 10 ** 6 * self.eta_cc)
        self.mf_airfuel = self.mf_hot + self.mf_fuel  # at the end of the cc

        # T04 = 1500 [K], is given
        self.p04 = self.p03 * self.PR_cc

        # Exit of HPT - Entrance of LPT
        self.T045 = self.T04 + (self.T04 / self.eta_HPT) * ((PR_HPT) ** ((self.k_gas - 1) / self.k_gas) - 1)
        self.p045 = self.p04 * PR_HPT

        # Exit of LPT - Entrance of nozzle
        self.T05 = self.T04 + (self.T04 / self.eta_LPT) * ((PR_LPT) ** ((self.k_gas - 1) / self.k_gas) - 1)
        self.p05 = self.p045 * PR_LPT

        # Nozzle
        self.T07 = self.T05
        self.p07 = self.p05 * self.PR_noz_core

        # Is the nozzle chocked?
        self.PR_cr_noz_core = 1 / (
                        (1 - (self.k_gas - 1) / (self.k_gas + 1) / self.eta_nozzle) ** (self.k_gas / (self.k_gas - 1)))

        # Exit of the nozzle
        if self.p07 / self.p0[i] > self.PR_cr_noz_core:
            # print('The nozzle is chocked')
            self.TR_cr_noz_core = (self.k_gas + 1) / 2
            self.T8 = self.T07 / self.TR_cr_noz_core
            self.p8 = self.p07 / self.PR_cr_noz_core
            self.v8 = np.sqrt(self.k_gas * self.R * self.T8)
            self.A8 = (self.mf_airfuel * self.R * self.T8) / (self.p8 * self.v8)
            self.T_core = self.mf_airfuel * (self.v8 - self.v0[i]) + self.A8 * (self.p8 - self.p0[i])  # [N]

        elif self.p07 / self.p0[i] <= self.PR_cr_noz_core:
            # print('The nozzle is NOT chocked')
            self.p8 = self.p0[i]
            self.T8 = self.T07 * (1 - self.eta_nozzle * (1 - (self.p8 / self.p07) ** ((self.k_gas - 1) / self.k_gas)))
            self.v8 = np.sqrt(2 * self.cp_gas * (self.T07 - self.T8))
            self.T_core = self.mf_airfuel * (self.v8 - self.v0[i])  # [N]

        # Is the fan chocked?
        self.PR_cr_fan = 1 / (
                    (1 - (self.k_air - 1) / (self.eta_nozzle * (self.k_air + 1))) ** (self.k_air / (self.k_air - 1)))

        # Exit of bypassed air
        if self.p016 / self.p0[i] > self.PR_cr_fan:
            # print('The fan is chocked')
            self.TR_cr_bypassed = (self.k_air + 1) / 2
            self.T18 = self.T016 / self.TR_cr_bypassed
            self.p18 = self.p016 / self.PR_cr_fan
            self.v18 = np.sqrt(self.k_air * self.R * self.T18)
            self.A18 = (self.mf_cold * self.R * self.T18) / (self.p18 * self.v18)
            self.T_fan = self.mf_cold * (self.v18 - self.v0[i]) + self.A18 * (self.p18 - self.p0[i])  # [N]

        elif self.p016 / self.p0[i] <= self.PR_cr_fan:
            # print('The fan is NOT chocked')
            self.p18 = self.p0[i]
            self.T18 = self.T016 - self.T016 * self.eta_fan * (
                        1 - (self.p18 / self.p016) ** ((self.k_air - 1) / self.k_air))
            self.v18 = np.sqrt(2 * self.cp_air * (self.T016 - self.T18))
            self.T_fan = self.mf_cold * (self.v18 - self.v0[i])  # [N]

        self.T_total = self.T_fan + self.T_core  # [N]
        self.TSFC = self.mf_fuel / self.T_total
        # print("TSFC:", self.TSFC)


if __name__=='__main__':
    c = Constants()
    for i in [2]:
        sa = SA_Engine_Cycle(i)
        print('Phase:', c.phases[i])
        print(sa.opt_vars)
