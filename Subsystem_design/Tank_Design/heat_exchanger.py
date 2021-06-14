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
        self.visc_air_eng = 4.25 * 10 ** (
            -5)  # https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
        self.thermal_conduct_air = 0.07  # https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
        self.cp_air = 1150  # https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
        self.T_H2_0 = 30  # Initial temperature of hydrogen [k]
        self.T_H2_final = 700  # Final temperature of hydrogen

    def heat_flux(self, r_i, m_flow_core, D_o, D_T0, D_T1, q):
        # Fuel flow of hydrogen into the engine:

        F_flow_H2 = [h2ff[0], h2ff[2], h2ff[4], h2ff[6], h2ff[8], h2ff[10]]
        flight_phases = ['Idle', 'Taxi out', 'Climb', 'Cruise', 'Descent', 'Taxi in']

        plt.bar(flight_phases, F_flow_H2)
        plt.show()

        max_ff = max(F_flow_H2)
        print(max_ff)

        # Specific heat of hydrogen

        T_H2 = np.array([20, 30, 40, 45, 50, 60, 70, 80, 100, 120, 140, 200, 240, 264])
        cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4,
                          3.45]) * 4184  # J/(kg*K) https://www.bnl.gov/magnets/Staff/Gupta/cryogenic-data-handbook/Section3.pdf
        visc_H2 = np.array([2500, 1000, 600, 500, 400, 375, 380, 400, 450, 500, 550, 700, 800, 820]) * 10 ** (
            -5) * 0.001  # Pa*s https://pubs-acs-org.tudelft.idm.oclc.org/doi/pdf/10.1021/i160007a014
        k_H2 = np.array([0.119, 0.11, 0.11, 0.1, 0.1, 0.09, 0.09, 0.09, 0.08, 0.09, 0.1, 0.14, 0.165,
                         0.170])  # https://www.engineeringtoolbox.com/hydrogen-H2-thermal-conductivity-temperature-pressure-d_2106.html

        plt.plot(T_H2, cp_H2)
        plt.ylabel(r'$C_p$ [ J / $kg \cdot K$ ]')
        plt.xlabel(r'T [$K$]')
        plt.show()

        LMTD = (D_T0 - D_T1) / np.log(D_T0 / D_T1)

        Re_D = m_flow_core / (np.pi * r_i ** 2 * self.visc_air_eng)
        Pr = self.visc_air_eng * self.cp_air / self.thermal_conduct_air
        h_i = self.thermal_conduct_air * 0.023 * Re_D ** 0.8 * Pr ** 0.3 / (2 * r_i)

        print('q required = ', q, 'W')
        print('\n Re_D = ', Re_D)
        print('\n Pr = ', Pr)
        print('\n NuD = ', 0.023 * Re_D ** 0.8 * Pr ** 0.3)
        print('\n h_i = ', h_i)

        def h_o_calc(D_o):
            Dh = D_o - 2 * r_i
            Re_Dh = max_ff / (np.pi * (Dh / 4) ** 2 * visc_H2)
            # print('\n Re_Dh = ', Re_Dh)
            Pr_h = visc_H2 * cp_H2 / k_H2
            h_out = k_H2 * 0.023 * Re_Dh ** 0.8 * Pr_h ** 0.4 / (Dh)
            h_o_avg = np.average(h_out)
            # plt.plot(T_H2, h_out)
            # plt.show()

            return h_o_avg

        h_o = h_o_calc(D_o)

        Ui = (1 / h_i + r_i / self.k_therm * np.log(D_o / 2 / r_i)
              + r_i / (h_o * D_o / 2)) ** (-1)

        L = q / (Ui * 2 * np.pi * r_i * LMTD)

        print('Length of the heat exchanger = ', L, ' m')


if __name__ == "__main__":
    print('\n \n ========= Heat Exchangers ==================================== \n')

    F_flow_H2 = [h2ff[0], h2ff[2], h2ff[4], h2ff[6], h2ff[8], h2ff[10]]
    flight_phases = ['Idle', 'Taxi out', 'Climb', 'Cruise', 'Descent', 'Taxi in']

    # plt.bar(flight_phases, F_flow_H2)
    # plt.show()
    #
    # max_ff = max(F_flow_H2)

    dat = pd.read_csv('../Engine/hack_take_off.txt', sep='\t', header=None, names=['variable', 'value', 'unit', 'nan'])
    print(dat)
    T_comp = dat.loc[dat['variable'] == 'T03']['value'].values[0]
    T_nozz = dat.loc[dat['variable'] == 'T07']['value'].values[0]
    P_comp = dat.loc[dat['variable'] == 'p03']['value'].values[0]
    m_flow_h2 = dat.loc[dat['variable'] == 'm_h2']['value'].values[0]
    m_flow_hot = dat.loc[dat['variable'] == 'm_hot']['value'].values[0]

    T_H2_0 = 30  # Initial temperature of hydrogen [k]
    T_H2_final = 600
    T_H2 = np.array([20, 30, 40, 45, 50, 60, 70, 80, 100, 120, 140, 200, 240, 264])
    cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4,
                      3.45]) * 4184  # J/(kg*K) https://www.bnl.gov/magnets/Staff/Gupta/cryogenic-data-handbook/Section3.pdf
    visc_H2 = np.array([2500, 1000, 600, 500, 400, 375, 380, 400, 450, 500, 550, 700, 800, 820]) * 10 ** (
        -5) * 0.001  # Pa*s https://pubs-acs-org.tudelft.idm.oclc.org/doi/pdf/10.1021/i160007a014
    k_H2 = np.array([0.119, 0.11, 0.11, 0.1, 0.1, 0.09, 0.09, 0.09, 0.08, 0.09, 0.1, 0.14, 0.165,
                     0.170])  # https://www.engineeringtoolbox.com/hydrogen-H2-thermal-conductivity-temperature-pressure-d_2106.html

    # plt.plot(T_H2, cp_H2)
    # plt.ylabel(r'$C_p$ [ J / $kg \cdot K$ ]')
    # plt.xlabel(r'T [$K$]')
    # plt.show()

    # ============================================
    # ========== Turbine Heat Exchanger ==========
    # ============================================

    print('\n \n ========== Turbine Heat Exchanger ==========')

    U_turb = 180
    T_turbine_cooling = 300
    cp_air = 1150
    cp_H2_turb = np.average(cp_H2[np.where(T_H2 < 150)[0]])
    m_flow_air_turb = 49.9 / 60
    print('\n The pressure of the air coming out of the compressor is = ', P_comp * 10 ** (-5), ' bar')

    q_air_turb = m_flow_air_turb * cp_air * (T_comp - T_turbine_cooling)
    T_H2_out_turb = T_H2_0 + q_air_turb / (m_flow_h2 * cp_H2_turb)

    print('\n T_H2_out = ', T_H2_out_turb)

    Ch = m_flow_hot * cp_air
    Cc = m_flow_h2 * np.average(cp_H2[np.where(T_H2 < T_H2_out_turb)[0]])

    print('\n C_c = ', Cc, ' and C_h = ', Ch)
    C_min = min(Cc, Ch)
    C_max = max(Cc, Ch)

    print(' C_min/C_max = ', C_min / C_max)
    print(' eps_c = ', (T_H2_out_turb - T_H2_0) / (T_comp - T_H2_0))

    NTU = 0.1
    A_turb_heat = NTU * C_min / U_turb

    m2_m3 = 700  # meters squared of surface to heat per meter cubed

    print('\n The required area for the cross flow heat exchanger is = ', A_turb_heat, 'm^2')
    print(' The required volume for the cross flow heat exchanger is = ', A_turb_heat / 700, 'm^3')

    # ========================================
    # ========== Oil Heat Exchanger ==========
    # ========================================

    print('\n \n ========== Oil Heat Exchanger ==========')

    T_oil_in, T_oil_out = 378, 343
    cp_H2_oil = 3.3 * 4184
    m_flow_oil = 57.6 / 60
    cp_oil = 2000

    q_oil = m_flow_oil * cp_oil * (T_oil_in - T_oil_out)
    T_H2_out_oil = T_H2_out_turb + q_oil / (cp_H2_oil * m_flow_h2)

    Ch = m_flow_oil * cp_oil
    Cc = m_flow_h2 * cp_H2_oil

    print('\n C_c = ', Cc, ' and C_h = ', Ch)
    C_min = min(Cc, Ch)
    C_max = max(Cc, Ch)

    print(' C_min/C_max = ', C_min / C_max)
    print(' eps_h = ', (T_oil_in - T_oil_out) / (T_oil_in - T_H2_0))

    NTU = 0.1
    A_oil_heat = NTU * C_min / U_turb

    print('\n The required area for the cross flow heat exchanger is = ', A_oil_heat, 'm^2')
    print(' The required volume for the cross flow heat exchanger is = ', A_oil_heat / 700, 'm^3')

    # ===========================================
    # ========== Nozzle Heat Exchanger ==========
    # ===========================================

    print('\n \n ========== Nozzle Heat Exchanger ========== \n')

    visc_air = 4.25 * 10 ** (-5)
    cp_H2_nozz = 3.5 * 4184
    k_air = 0.075
    r_i = 0.1  # np.linspace(0.15, 0.6, 100)
    dr = 0.005

    r_o = r_i + dr
    D_o = 2 * r_o
    T_nozz_arr = np.linspace(T_H2_out_oil, T_H2_final, 100)
    H2_visc = 208 * 10 ** (-5) * (T_nozz_arr / 33.3) ** 0.65 * 0.001
    H2_k = (0.45 - 0.3) / (630 - 270) * (T_nozz_arr - 270) + 0.3

    q_nozz = m_flow_h2 * cp_H2_nozz * (T_H2_final - T_H2_out_oil)
    T_air_nozz_out = T_nozz - q_nozz / (cp_air * m_flow_hot)

    print('Heat transfer q = ', q_nozz)
    print('T_air_nozz_out = ', T_air_nozz_out)

    D_T0 = T_nozz - T_H2_out_oil
    D_T1 = T_air_nozz_out - T_H2_final

    LMTD = (D_T0 - D_T1) / np.log(D_T0 / D_T1)

    Re_D = m_flow_hot / (np.pi * r_i ** 2 * visc_air)
    Pr = visc_air * cp_air / k_air
    h_i = k_air * 0.023 * Re_D ** 0.8 * Pr ** 0.3 / (2 * r_i)

    print('\n Re_D = ', Re_D)
    print(' Pr = ', Pr)
    print(' NuD = ', 0.023 * Re_D ** 0.8 * Pr ** 0.3)
    print(' h_i = ', h_i)

    # plt.plot(T_nozz_arr, H2_visc)
    # plt.show()
    # plt.plot(T_nozz_arr, H2_k)
    # plt.show()
    Dh = D_o - 2 * r_i
    Re_Dh = m_flow_h2 / ((Dh / 2) ** 2 * np.average(H2_visc))
    Pr_h = np.average(H2_visc) * cp_H2_nozz / np.average(H2_k)
    h_out = np.average(H2_k) * 0.023 * Re_Dh ** 0.8 * Pr_h ** 0.4 / Dh

    print('\n Re_Dh = ', Re_Dh)
    print(' Pr_h = ', Pr_h)
    print('\n h_o = ', h_out)

    Ui = (1 / h_i + r_i / (h_out * D_o / 2)) ** (-1)

    L = q_nozz / (Ui * 2 * np.pi * r_i * LMTD)

    print('\n Length of the heat exchanger = ', L, ' m')

