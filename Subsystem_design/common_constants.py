import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

"""
This file contains one class only, which is meant to contain the variables which are common to the entire subsystem 
design. It may also contain some simple functions to compute constants derived from other constants (e.g. ISA).
"""

class Constants():
    def __init__(self):
        self.rho_0= 1.225  # Sea level density [kg/m^3]
        self.p_0 = 101325  # Sea level pressure [Pa]
        self.T_0 = 288.15  # Sea level temperature [K]
        self.g_0 = 9.80665  # Gravity at sea level [m/s^2]
        self.R = 287.0  # Specific gas constant of air [J/(kg*K)]
        self.gamma = 1.4  # Heat capacity ratio of air [-]
        self.a_0 = 340.294  # Sea level speed of sound [m/s]

        self.S = 122.6  #  Wing surface area [m^2]
        self.l_f = 37.57  # Fuselage length in [m] :todo: Change this value to account for H2
        self.d_f = 4.14  # Fuselage maximum diameter in [m]
        self.l_cockpit = 5.04  # Length of the cockpit [m]
        self.l_cabin = 29.53 - self.l_cockpit  # Length of the cabin [m] :todo: Change this value to account for H2
        self.l_tail = self.l_f - 29.53  # Length of the tail [m] :todo: Change this value to account for H2
        self.S_b_fus = np.pi * 0.3/2 * 0.45/2  # Base surface area [m^2]
        self.sweep = 27  # Wing sweep [deg]

    def speed_of_sound(self, T):
        """
        Compute the speed of sound in air for any temperature
        :param T: Temperature [K]
        :return: the speed of sound in air [m/s]
        """
        a = np.sqrt(self.gamma * self.R * T)
        return a

    def ISA_temp(self, gr, T_init, h_init, h):
        """
        International Standard Atmosphere equation for temperature
        :param gr: Temperature gradient [K/m]
        :param T_init: Temperature at the beginning of the atmospheric layer [K]
        :param h_init: Altitude at the beginning of the atmospheric layer [m]
        :param h: Altitude of interest [m]
        :return: Temperature at altitude h according to the ISA [K]
        """
        ISA_T = T_init + gr * (h - h_init)
        return ISA_T

    def ISA_grad_pressure(self, p1_grad, T1_grad, ISA_T, gr):
        """
        International Standard Atmosphere equation for pressure in a layer with changing temperature
        :param p1_grad: Pressure at the beginning of the atmospheric layer [Pa]
        :param T1_grad: Temperature at the beginning of the atmospheric layer [K]
        :param ISA_T: Temperature at the desired altitude [K]
        :param gr: Temperature gradient [K/m]
        :return: Pressure at altitude h according to the ISA [Pa]
        """
        ISA_p = p1_grad * (ISA_T / T1_grad)**(-self.g_0 / (gr * self.R))
        return ISA_p

    def ISA_iso_pressure(self, p1_iso, T1_iso, h1_iso, h):
        """
        International Standard Atmosphere equation for pressure in a layer with constant temperature
        :param p1_iso: Pressure at the beginning of the atmospheric layer [Pa]
        :param T1_iso: Temperature at the beginning of the atmospheric layer [K]
        :param h1_iso: Altitude at the beginning of the atmospheric layer [m]
        :param h: Altitude of interest [m]
        :return: Pressure at altitude h according to the ISA [Pa]
        """
        ISA_p = p1_iso * np.exp(-self.g_0 / (self.R * T1_iso) * (h - h1_iso))
        return ISA_p

    def ISA_density(self, ISA_p, ISA_T):
        """
        Ideal gas law for air
        :param ISA_p: Pressure at altitude h according to the ISA [Pa]
        :param ISA_T: Temperature at altitude h according to the ISA [K]
        :return: Air density at altitude h according to the ISA [kg/m^3]
        """
        rho = ISA_p / (ISA_T * self.R)
        return rho

    def ISA_calculator(self, h_input):
        """
        Compute the atmospheric conditions at altitude h_input
        :param h_input: Altitude of interest [m]
        :return: Nothing, it only defines attributes of the class
        """

        h_atmos = np.array([0., 11000., 20000., 32000., 47000., 51000., 71000., 86000])
        p_atmos = np.array([101325., 22625.79, 5471.9, 867.25, 110.766, 66.848, 3.94898, 0.301618])
        T_atmos = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 184.65])
        grad_atmos = np.array([-0.0065, 0., 0.001, 0.0028, 0., -0.0028, -0.002, 0.])

        loc_atmos = np.where(h_atmos <= h_input)[0][-1]
        h1 = h_atmos[loc_atmos]
        p1 = p_atmos[loc_atmos]
        T1 = T_atmos[loc_atmos]
        grad = grad_atmos[loc_atmos]

        self.T = self.ISA_temp(grad, T1, h1, h_input)
        self.a = self.speed_of_sound(self.T)

        if grad != 0.:
            self.p = self.ISA_grad_pressure(p1, T1, self.T, grad)

        else:
            self.p = self.ISA_iso_pressure(p1, T1, h1, h_input)

        self.rho = self.ISA_density(self.p, self.T)



# Try out the class

if __name__ == '__main__':
    c = Constants()
    c.ISA_calculator(h_input=53000)
    print('\n T = ', c.T, ' K',
          '\n P = ', c.p, ' Pa',
          '\n rho = ', c.rho, ' kg/m^3')



