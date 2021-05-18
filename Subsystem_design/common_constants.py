import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from math import pi
from Subsystem_design.fuel_required import V_H2, V_k
"""
This file contains one class only, which is meant to contain the variables which are common to the entire subsystem 
design. It may also contain some simple functions to compute constants derived from other constants (e.g. ISA).
"""

class Constants():
    def __init__(self):
        self.rho_0 = 1.225                                          # Sea level density                         [kg/m^3]
        self.p_0 = 101325                                           # Sea level pressure                           [Pa]
        self.T_0 = 288.15                                           # Sea level temperature                        [K]
        self.g_0 = 9.80665                                          # Gravity at sea level                       [m/s^2]
        self.R = 287.0                                              # Specific gas constant of air            [J/(kg*K)]
        self.gamma = 1.4                                            # Heat capacity ratio of air                   [-]
        self.a_0 = 340.294                                          # Sea level speed of sound                     [m/s]

        '''Performance'''
        self.M = 0.78
        self.V_cruise = 230
        self.cruise_altitude = 11280

        ''' Dimensions of A320HACK'''
        self.S = 122.6                                              # Wing surface area                            [m^2]
        self.l_f = 37.57 + 2.8                                      # Fuselage length                               [m]
        self.height_f = 4.14                                        # Fuselage height                               [m]
        self.width_f = 3.95                                         # Fuselage width                                [m]
        self.l_cockpit = 5.04                                       # Length of the cockpit                         [m]
        self.l_cabin = 29.53 + 3.7 - self.l_cockpit                 # Length of the cabin + H2 tank                 [m]
        self.l_tail = self.l_f - self.l_cabin - self.l_cockpit      # Length of the tail                            [m]
        self.S_b_fus = np.pi * 0.3/2 * 0.45/2                       # Base surface area                            [m^2]
        self.sweep_LE = 27                                          # Wing sweep                                   [deg]

        """Fuel constant A320HACK"""
        # self.V_H2 = 37.893
        self.V_H2 = V_H2                                            # Volume required of hydrogen                  [m^3]
        self.V_H2 = V_k
        # self.V_k = 14.316                                           # Volume required of kerosene                  [m^3]

        """Weights of A320neo"""
        self.MTOW_320neo = 73500                                    # Maximum Take-Off weight of A320neo            [kg]
        self.MLW_320neo = 66300                                     # Maximum Landing weight of A320neo             [kg]
        self.MRW_320neo = 73900                                     # Maximum Ramp weight of A320neo                [kg]
        self.MZFW_320neo = 62800                                    # Maximu Zero fuel weight of A320neo            [kg]
        self.MPLW_320neo = 18240                                    # Maximum Payload weight of A320neo             [kg]
        self.OEW_320neo = 44300                                     # Operational Empty weight of A320neo           [kg]

        """Weights of A321neo"""
        self.MTOW_321neo = 89000                                    # Maximum Take-Off weight of A321neo            [kg]
        self.MLW_321neo = 77300                                     # Maximum Landing weight of A321neo             [kg]
        self.MRW_321neo = 89400                                     # Maximum Ramp weight of A321neo                [kg]
        self.MZFW_321neo = 73300                                    # Maximu Zero fuel weight of A321neo            [kg]
        self.MPLW_321neo = 22910                                    # Maximum Payload weight of A321neo             [kg]
        self.OEW_321neo = self.MZFW_321neo - self.MPLW_321neo       # Operational Empty weight of A320neo           [kg]

        """Dimensions of A320neo and A321neo"""
        self.l_f_321neo= 44.51                                      # Fuselage length of A321neo                    [m]
        self.l_f_320neo = 37.57                                     # Fuselage length of A320neo                    [m]
        self.l_cockpit_320neo = 5.04                                # Length of the cockpit of A320neo              [m]
        self.l_cabin_320neo = 29.53 - self.l_cockpit_320neo         # Length of the cabin of A320neo                [m]
        self.l_tail_320neo = self.l_f_320neo - 29.53                # Length of the tail of A320neo                 [m]

        self.c_j_kerosene = 16.68 * 10**(-6)                        # Specific cruise fuel consumption of neo   [kg/N*s]

    def fuselage_length(self,vol_eff, vol_fus):
        """

        :param vol_eff: Volumetric efficiency of integral tanks (Ratio of usable tank volume-to-volume occupied
                        in the fuselage                                         [-]
        :param vol_fus: Volume available for tanks in the center wingbox        [m^3]
        :return:
        """
        self.V_H2_center_w = vol_eff * vol_fus                           # Volume of hydrogen stored on center wing box [m^3]
        self.V_H2_ext_fus = self.V_H2 - self.V_H2_center_w               # Volume of hydrogen stored on extended fuselage[m^3]
        self.l_f = self.V_H2_ext_fus/(pi * self.width_f/2 * self.height_f/2) #

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



