from Subsystem_design.common_constants import Constants
from Fuel_Masses_Estimation import kff, h2ff
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

class HeatExchanger(Constants):

    def __init__(self):
        super().__init__()
        self.A_i = 2 * np.pi * self.r_i_nozzle * 0.352
        self.k = 12  # Thermal conductivity of Nickel Chrome Alloys [W/mK]
        self.T_exhaust = 750 - 273.15  # Minimum exhaust temperature of the engine [K]

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
        cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4, 3.45]) * 4184  # J/(kg*K)

        plt.plot(T_H2, cp_H2)
        plt.ylabel(r'$C_p$ [ J / $kg \cdot K$ ]')
        plt.xlabel(r'T [$K$]')
        plt.show()

        q = max_ff * spint.simps(y=cp_H2, x=T_H2)
        D_T0 = self.T_exhaust - self.T_LH2in
        D_T1 = self.T_exhaust - self.T_H2out
        LMTD = (D_T0 - D_T1) / np.ln(D_T0 / D_T1)

        Ui = q / (self.A_i * LMTD)




if __name__ == "__main__":
    he = HeatExchanger()
    he.heat_flux()