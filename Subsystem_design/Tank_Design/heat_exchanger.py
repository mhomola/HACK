from Subsystem_design.common_constants import Constants
from Fuel_Masses_Estimation import kff, h2ff
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint
import pandas as pd

class HeatExchanger(Constants):

    def __init__(self):
        super().__init__()
        self.k_therm = 50  # Thermal conductivity of Nickel Chrome Alloys [W/mK]
        self.T_exhaust = 1060 + 273.15  # Minimum exhaust temperature of the engine [K]
        self.visc_air_eng = 4.25 * 10**(-5)  # https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
        self.thermal_conduct_air = 0.07  # https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
        self.cp_air = 1150  # https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
        self.T_H2_0 = 30  # Initial temperature of hydrogen [k]
        self.T_H2_final = 700  # Final temperature of hydrogen


    def heat_flux(self, r_i, m_flow_core, D_o, D_T0, D_T1):

        # Fuel flow of hydrogen into the engine:

        F_flow_H2 = [h2ff[0], h2ff[2], h2ff[4], h2ff[6], h2ff[8], h2ff[10]]
        flight_phases = ['Idle', 'Taxi out', 'Climb', 'Cruise', 'Descent', 'Taxi in']

        plt.bar(flight_phases, F_flow_H2)
        plt.show()

        max_ff = max(F_flow_H2)
        print(max_ff)

        # Specific heat of hydrogen

        T_H2 = np.array([20, 30, 40, 45, 50, 60, 70, 80, 100, 120, 140, 200, 240, 264])
        cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4, 3.45]) * 4184  # J/(kg*K) https://www.bnl.gov/magnets/Staff/Gupta/cryogenic-data-handbook/Section3.pdf
        visc_H2 = np.array([2500, 1000, 600, 500, 400, 375, 380, 400, 450, 500, 550, 700, 800, 820]) * 10**(-5) * 0.001  # Pa*s https://pubs-acs-org.tudelft.idm.oclc.org/doi/pdf/10.1021/i160007a014
        k_H2 = np.array([0.119, 0.11, 0.11, 0.1, 0.1, 0.09, 0.09, 0.09, 0.08, 0.09, 0.1, 0.14, 0.165, 0.170])  # https://www.engineeringtoolbox.com/hydrogen-H2-thermal-conductivity-temperature-pressure-d_2106.html

        plt.plot(T_H2, cp_H2)
        plt.ylabel(r'$C_p$ [ J / $kg \cdot K$ ]')
        plt.xlabel(r'T [$K$]')
        plt.show()

        q = max_ff * spint.simps(y=cp_H2, x=T_H2)
        LMTD = (D_T0 - D_T1) / np.log(D_T0 / D_T1)

        Re_D = m_flow_core / (np.pi * r_i**2 * self.visc_air_eng)
        Pr = self.visc_air_eng * self.cp_air / self.thermal_conduct_air
        h_i = self.thermal_conduct_air * 0.023 * Re_D**0.8 * Pr**0.3 / (2 * r_i)

        print('q required = ', q, 'W')
        print('\n Re_D = ', Re_D)
        print('\n Pr = ', Pr)
        print('\n NuD = ', 0.023 * Re_D**0.8 * Pr**0.3)
        print('\n h_i = ', h_i)

        def h_o_calc(D_o):
            Dh = D_o - 2 * r_i
            Re_Dh = max_ff / (np.pi * (Dh/4)**2 * visc_H2)
            print('\n Re_Dh = ', Re_Dh)
            Pr_h = visc_H2 * cp_H2 / k_H2
            h_out = k_H2 * 0.023 * Re_Dh**0.8 * Pr_h**0.4 / (Dh)
            # h_out = k_H2 * 4.86 / Dh
            h_o_avg = np.average(h_out)
            plt.plot(T_H2, h_out)
            plt.show()

            return h_o_avg

        h_o = h_o_calc(D_o)

        Ui = (1/h_i + r_i/self.k_therm * np.log(D_o/2 / r_i)
              + r_i / (h_o * D_o/2))**(-1)

        L = q / (Ui * 2 * np.pi * r_i * LMTD)

        print('Length of the heat exchanger = ', L, ' m')



if __name__ == "__main__":

    print('\n \n ========= Heat Exchanger ==================================== \n')

    dat = pd.read_csv('../Engine/hack_take_off.txt', sep='\t', header=None, names=['variable', 'value', 'unit', 'nan'])
    print(dat)
    T_comp = dat.loc[dat['variable'] == 'T03']['value'].values[0]
    P_comp = dat.loc[dat['variable'] == 'p03']['value'].values[0]
    m_flow_h2 = dat.loc[dat['variable'] == 'm_h2']['value'].values[0]
    m_flow_hot = dat.loc[dat['variable'] == 'm_hot']['value'].values[0]

    T_H2 = np.array([20, 30, 40, 45, 50, 60, 70, 80, 100, 120, 140, 200, 240, 264])
    cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4, 3.45]) * 4184  # J/(kg*K) https://www.bnl.gov/magnets/Staff/Gupta/cryogenic-data-handbook/Section3.pdf
    visc_H2 = np.array([2500, 1000, 600, 500, 400, 375, 380, 400, 450, 500, 550, 700, 800, 820]) * 10**(-5) * 0.001  # Pa*s https://pubs-acs-org.tudelft.idm.oclc.org/doi/pdf/10.1021/i160007a014
    k_H2 = np.array([0.119, 0.11, 0.11, 0.1, 0.1, 0.09, 0.09, 0.09, 0.08, 0.09, 0.1, 0.14, 0.165, 0.170])  # https://www.engineeringtoolbox.com/hydrogen-H2-thermal-conductivity-temperature-pressure-d_2106.html

    # plt.plot(T_H2, cp_H2)
    # plt.ylabel(r'$C_p$ [ J / $kg \cdot K$ ]')
    # plt.xlabel(r'T [$K$]')
    # plt.show()

    he = HeatExchanger()

    # Turbine Heat Exchanger
    T_turbine_cooling = 300
    T_H2_out = 150
    print('\n The pressure of the air coming out of the compressor is = ', P_comp*10**(-5), ' bar')

    Ch = m_flow_hot * he.cp_air
    Cc = m_flow_h2 * np.average(cp_H2[np.where(T_H2 < T_H2_out)[0]])

    print('C_c = ', Cc, ' and C_h = ', Ch)


    # Oil Heat Exchanger
    # D_

    # Nozzle Heat exchanger
    # D_T0 = he.T_exhaust - he.T_LH2in
    # D_T1 = he.T_exhaust - he.T_H2out
    # he.heat_flux(r_i=0.5, D_o=2*0.5+0.005, m_flow_core=43, D_T0=D_T0, D_T1=D_T1)