from Subsystem_design.common_constants import Constants
from Fuel_Masses_Estimation import kff, h2ff
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

class HeatExchanger(Constants):

    def __init__(self):
        super().__init__()
        self.A_i = 2 * np.pi * self.r_i_nozzle * 0.5
        self.k_therm = 12  # Thermal conductivity of Nickel Chrome Alloys [W/mK]
        self.T_exhaust = 750 + 273.15  # Minimum exhaust temperature of the engine [K]
        self.visc_air_eng = 4.25 * 10**(-5)  # https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm
        self.thermal_conduct_air = 0.07  # https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
        self.cp_air = 1150  # https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html


    def heat_flux(self):

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
        D_T0 = self.T_exhaust - self.T_LH2in
        D_T1 = self.T_exhaust - self.T_H2out
        LMTD = (D_T0 - D_T1) / np.log(D_T0 / D_T1)

        Ui = q / (self.A_i * LMTD)
        print('\n Ui = ', Ui)

        m_flow_core = 43

        Re_D = 4 * m_flow_core / (np.pi * 2 * self.r_i_nozzle * self.visc_air_eng)
        Pr = self.visc_air_eng * self.cp_air / self.thermal_conduct_air
        h_i = self.thermal_conduct_air * 0.0023 * Re_D**0.8 * Pr**0.3 / (2 * self.r_i_nozzle)

        # print(h_i)

        def h_o_calc(D_o):
            Dh = D_o - 2 * self.r_i_nozzle
            Re_Dh = 4 * max_ff / (np.pi * (D_o + 2 * self.r_i_nozzle) * visc_H2)
            Pr_h = visc_H2 * cp_H2 / k_H2
            h_out = k_H2 * 0.023 * Re_Dh**0.8 * Pr_h**0.4 / (D_o)
            h_o_avg = np.average(h_out)

            return h_o_avg

        D_o_arr = np.arange(2*self.r_i_nozzle, 2, 0.001)
        h_o_arr = np.zeros(len(D_o_arr))
        for i, d in enumerate(D_o_arr):
            h_o_arr[i] = h_o_calc(d)

        print(h_o_arr)

        print(Ui)
        solver = (1/h_i + self.r_i_nozzle/self.k_therm * np.log(D_o_arr/2 / self.r_i_nozzle)
                  + self.r_i_nozzle / (h_o_arr * D_o_arr/2))**(-1)

        plt.plot(D_o_arr, solver)
        plt.show()



if __name__ == "__main__":
    he = HeatExchanger()
    he.heat_flux()