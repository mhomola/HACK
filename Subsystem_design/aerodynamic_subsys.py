from common_constants import Constants
import numpy as np


class AerodynamicCharacteristics(Constants):
    def __init__(self):
        super().__init__()
        self.M = 0.78  # Cruising Mach number
        self.cruise_alt = 11280  # Initial cruising altitude [m]
        self.visc = 1.458 * 10**(-5)  # Air viscosity [N*s/m^2]

        self.b_in = 2 * 6.4  # Span of inner wing trapezoid [m]
        self.b_out = 2 * 10.616  # Span of outer wing trapezoid [m]
        self.b = self.b_in + self.b_out  # Span of the entire wing
        self.c_root = 7.0465  # Root chord including the kink [m]
        self.c_kink_out = 3.72  # Chord at the point where the kink ends [m]
        self.c_tip = 1.488  # Chord at the tip of the wing excluding sharklet [m]
        self.taper_in = self.c_kink_out / self.c_root  # Taper ratio of the inner trapezoid
        self.taper_out = self.c_tip / self.c_kink_out  # Taper ratio of the outer trapezoid
        self.h_sharklet = 2.43  # Height of the wing sharklets [m]

        self.b_h = 2 * 6.12  # Span of the horizontal tail [m]
        self.c_r_h = 3.814  # Root chord of the horizontal tail [m]
        self.c_t_h = 1.186  # Tip chord of the horizontal tail [m]
        self.taper_h = self.c_t_h / self.c_r_h  # Taper ratio of the horizontal tail
        self.sweep_LE_h = 33  # Sweep of the LE of the horizontal tail [deg]

        self.MTOW_neo = 79000  # Maximum take off weight of A320neo [kg]
        self.RoC_neo = 12.7  # Max rate of climb of the A320neo [m/s]

        self.C_D_0_TO_neo = 0.078  # Zero-lift drag coefficient of A320 during take-off
        # self.
        
    def wing_MAC(self):
        """
        The wing MAC is computed using the ADSEE-II slides. The wing surface is composed out of two trapezoids.
        First the mac of these two trapezoids is computed separately and after the global MAC is computed.
        :return: The length and position of the wing's MAC
        """
        # Inner trapezoid
        mac_in = (2/3) * self.c_root * (1 + self.taper_in + self.taper_in**2) / (1 + self.taper_in)  # MAC [m]
        y_mac_in = (self.b_in / 6) * (1 + 2 * self.taper_in) / (1 + self.taper_in)  # y position of MAC [m]
        x_mac_in = y_mac_in * np.tan(self.sweep_LE * np.pi / 180)  # x location of LEMAC from start of root chord [m]
        S_in = self.b_in * (self.c_root + self.c_kink_out) / 2  # Surface area [m^2]

        # Outer trapezoid
        mac_out = (2/3) * self.c_kink_out * (1 + self.taper_out + self.taper_out**2) / (1 + self.taper_out)
        y_mac_out = (self.b_out / 6) * (1 + 2 * self.taper_out) / (1 + self.taper_out)
        x_mac_out = y_mac_out * np.tan(self.sweep_LE * np.pi / 180)
        S_out = self.b_out * (self.c_kink_out + self.c_tip) / 2

        # Complete wing
        self.mac = (S_out * mac_out + S_in * mac_in) / (S_out + S_in)  # MAC of the complete wing [m]
        self.y_mac = (y_mac_in * S_in + (0.5 * self.b_in + y_mac_out) * S_out) / (S_out + S_in)  # y position of MAC [m]
        self.x_mac = self.y_mac * np.tan(self.sweep_LE * np.pi / 180)  # x location of LEMAC from root chord [m]

        # V&V - Unit tests
        print('\n MAC of the inner trapezoid = ', mac_in, ' m',
              '\n with y position of the MAC = ', y_mac_in, ' m',
              '\n and with surface area = ', S_in, ' m^2')

        print('\n MAC of the outer trapezoid = ', mac_out, ' m',
              '\n with y position of the MAC = ', 0.5 * self.b_in + y_mac_out, ' m',
              '\n and with surface area = ', S_out, ' m^2')

        print('\n MAC of the wing = ', self.mac, ' m',
              '\n y position of the MAC = ', self.y_mac, ' m',
              '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac, ' m')

    def wing_AR(self):
        """
        Compute the aspect ratio of the wing including sharklets with ADSEE-II folrmula
        :return: The AR of the main wing
        """
        AR_no_sharklets =  self.b**2 / self.S
        self.AR = AR_no_sharklets + 1.9 * (self.h_sharklet / self.b) * AR_no_sharklets

    def h_tail_MAC(self):
        """
        The MAC is computed using the ADSEE-II slides.
        :return: The length and position of the horizontal wing's MAC
        """
        self.mac_h = (2/3) * self.c_r_h * (1 + self.taper_h + self.taper_h**2) / (1 + self.taper_h)
        self.y_mac_h = (self.b_h / 6) * (1 + 2 * self.taper_h) / (1 + self.taper_h)
        self.x_mac_h = self.y_mac_h * np.tan(self.sweep_LE_h * np.pi / 180)

    def Roskam_drag_prediction_cruise(self):
        self.ISA_calculator(h_input=self.cruise_alt)
        self.R_n_fus = self.rho * self.M * self.a * self.l_f / self.visc
        print('The Fuselage Reynolds Number R_f_fus is: ', self.R_n_fus, ' [-]')
        print('The Mach number M is: ', self.M)
        self.R_wf = 1.015
        self.C_f_fus = 0.002
        ld = self.l_f / self.d_f
        self.C_D_o_fus_exc_base = self.R_wf * self.C_f_fus * (1 + 60 / ld**3 + 0.0025 * ld) * self.S_wet_fus / self.S
        bf = np.sqrt(4 / np.pi * self.S_b_fus) / self.d_f
        self.C_D_b_fus = (0.029 * bf**3 / ((self.C_D_o_fus_exc_base * self.S / self.S_fus)**0.5)) * self.S_fus / self.S
        self.C_D_o_fus = self.C_D_o_fus_exc_base + self.C_D_b_fus




# Try out the class

if __name__ == '__main__':
    Ae = AerodynamicCharacteristics()
    Ae.wing_MAC()
