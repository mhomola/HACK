from Subsystem_design.common_constants import Constants
from Fuel_Masses_Estimation import kff, h2ff
import numpy as np
import matplotlib.pyplot as plt

class HeatExchanger(Constants):

    def __init__(self):
        super().__init__()
        self.A_i = 2 * np.pi * self.r_i_nozzle * 0.352

    def heat_flux(self):

        # Fuel flow of hydrogen into the engine:

        F_flow_H2 = [h2ff[0], h2ff[2], h2ff[4], h2ff[6], h2ff[8], h2ff[10]]
        flight_phases = ['Idle', 'Taxi out', 'Climb', 'Cruise', 'Descent', 'Taxi in']

        plt.bar(flight_phases, F_flow_H2)
        plt.show()

        max_ff = max(F_flow_H2)

        # Specific heat of hydrogen

        T_H2 = np.array([20, 30, 40, 45, 50, 60, 70, 80, 100, 120, 140, 200, 240, 264])
        cp_H2 = np.array([2, 3.2, 4.8, 5.4, 4.8, 3.9, 3.4, 3.2, 3.05, 3.1, 3.2, 3.35, 3.4, 3.45]) * 4184  # J/(g*K)

        plt.plot(T_H2, cp_H2)
        plt.ylabel(r'$C_p$ [$J/(g \cdot K)$]')

        Delta_T = T_H2[1:] - T_H2[:-1]
        cp = (cp_H2[1:] + cp_H2[:-1]) / 2