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
        self.visc = 1.458 * 10**(-5)                                # Air viscosity                            [N*s/m^2]


        """Properties of H2"""

        self.H2_ed = 33.5                                   # Energy density of hydrogen                         [kWh/kg]

        """Properties of kerosene"""
        self.k_ed = 12.0                                    # Energy density of kerosene                         [kWh/kg]

        '''Performance'''
        self.cruise_altitude = 11280
        self.ISA_calculator(h_input=self.cruise_altitude)

        self.rho_c = self.rho
        self.a_c = self.a
        self.T_c = self.T
        self.p_c = self.p

        self.M = 0.78
        self.V_cruise = self.M * self.a_c

        '''Aerodynamics'''
        self.e = 0.992                                              # Oswald efficiency factor
        self.C_D_0_TO_neo = 0.078                                   # Zero-lift drag coefficient of A320neo - TO
        self.C_D_0_clean_neo = 0.023                                # Zero-lift drag coefficient of A320neo - cruise
        self.C_L_max_clean = 1.8                                    # Max lift coefficient clean
        self.C_L_max_TO = 2.2                                       # Max lift coefficient during TO
        self.C_L_max_land = 3.0                                     # Max lift coefficient during landing

        self.b_in = 2 * 6.4                                         # Span of inner wing trapezoid                  [m]
        self.b_out = 2 * 10.616                                     # Span of outer wing trapezoid                  [m]
        self.b = self.b_in + self.b_out                             # Span of the entire wing
        self.c_root = 7.0465                                        # Root chord including the kink                 [m]
        self.c_kink_out = 3.72                                      # Chord at the point where the kink ends        [m]
        self.c_tip = 1.488                                          # Chord at the tip of the wing excl sharklet    [m]
        self.taper_in = self.c_kink_out / self.c_root               # Taper ratio of the inner trapezoid
        self.taper_out = self.c_tip / self.c_kink_out               # Taper ratio of the outer trapezoid
        self.h_sharklet = 2.43                                      # Height of the wing sharklets                  [m]
        self.sweep_025 = 25                                         # Sweep at the quarter chord                   [deg]
        self.sweep_05 = 22.4                                        # Sweep at the half chord                      [deg]

        self.b_h = 2 * 6.12                                         # Span of the horizontal tail                   [m]
        self.c_r_h = 3.814                                          # Root chord of the horizontal tail             [m]
        self.c_t_h = 1.186                                          # Tip chord of the horizontal tail              [m]
        self.taper_h = self.c_t_h / self.c_r_h                      # Taper ratio of the horizontal tail
        self.sweep_LE_h = 33                                        # Sweep of the LE of the horizontal tail       [deg]

        ''' Dimensions of A320-HACK'''
        self.S = 122.6                                              # Wing surface area                            [m^2]
        self.l_f = 37.57 + 3.15                                     # Fuselage length                               [m]
        self.height_f = 4.14                                        # Fuselage height                               [m]
        self.width_f = 3.95                                         # Fuselage width                                [m]
        self.l_cockpit = 5.04                                       # Length of the cockpit                         [m]
        self.l_cabin = 29.53 + 3.7 - self.l_cockpit                 # Length of the cabin + H2 tank                 [m]
        self.l_tail = self.l_f - self.l_cabin - self.l_cockpit      # Length of the tail                            [m]
        self.S_b_fus = np.pi * 0.3/2 * 0.45/2                       # Base surface area                            [m^2]
        self.sweep_LE = 27                                          # Wing sweep                                   [deg]

        """Fuel constant A320-HACK"""
        # self.V_H2 = 37.893
        self.V_H2 = V_H2                                            # Volume required of hydrogen                  [m^3]
        self.V_H2 = V_k
        # self.V_k = 14.316                                           # Volume required of kerosene                  [m^3]
        self.W1_Wto = 0.990
        self.W2_W1 = 0.990
        self.W3_W2 = 0.995
        self.W4_W3 = 0.980
        self.W6_W5 = 0.990
        self.W7_W6 = 0.992
        self.c_j_kerosene = 16.68 * 10 ** (-6)                     # Specific cruise fuel consumption of neo    [kg/N*s]
        self.c_j_k_H2_cruise = 11.83 * 10 ** (-6)                  # Specific cruise fuel consumption of HACK   [kg/N*s]

        """Fuel Contants A320neo"""
        self.fuel_capacity_a320neo_FUTURE = 23859*0.9
        self.fuel_capacity_320neo = 23.859                          # Maximum Fuel capacity of A320neo              [m^3]
        self.k_d = 810.0                                            # Mass density of kerosene                           [kg/m^3]

        """Tank design constants""" #Plsss automate these, for design changes
        self.center_tank_mass = 390.6                               # Mass of center tanks in total (2 tanks)       [kg]
        self.fuselage_tank_mass = 286.6                             # Mass of aft tank (1 tank)                     [kg]

        """Weights of HACK"""
        self.Fuel_idel_taxi_take_off_HACK = 262.88                # Fuel for before take -off                     [kg]

        """Weights of A320neo"""
        self.MTOW_320neo = 73500                                    # Maximum Take-Off weight of A320neo            [kg]
        self.MLW_320neo = 66300                                     # Maximum Landing weight of A320neo             [kg]
        self.MRW_320neo = 73900                                     # Maximum Ramp weight of A320neo                [kg]
        self.MZFW_320neo = 62800                                    # Maximu Zero fuel weight of A320neo            [kg]
        self.MPLW_320neo = 18240                                    # Maximum Payload weight of A320neo             [kg]
        self.OEW_320neo = 44560                                     # Operational Empty weight of A320neo           [kg]
        self.Fuel_idel_taxi_take_off_320neo = 400                   # Fuel for before take -off                     [kg]
        self.Max_fuel_mass_capacity_320neo = self.fuel_capacity_320neo * self.k_d   #Maximum kerosene mass of A320neo [kg]

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

    def fuselage_length(self,vol_eff, vol_fus):
        """

        :param vol_eff: Volumetric efficiency of integral tanks (Ratio of usable tank volume-to-volume occupied
                        in the fuselage                                         [-]
        :param vol_fus: Volume available for tanks in the center wingbox        [m^3]
        :return:
        """
        self.V_H2_center_w = vol_eff * vol_fus  # Volume of hydrogen stored on center wing box [m^3]
        self.V_H2_ext_fus = self.V_H2 - self.V_H2_center_w  # Volume of hydrogen stored on extended fuselage[m^3]
        self.l_f = self.V_H2_ext_fus/(pi * self.width_f/2 * self.height_f/2)

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

#Weigth_Centre_Tanks

if __name__ == '__main__':
    c = Constants()
    c.ISA_calculator(h_input=53000)
    print('\n T = ', c.T, ' K',
          '\n P = ', c.p, ' Pa',
          '\n rho = ', c.rho, ' kg/m^3')



