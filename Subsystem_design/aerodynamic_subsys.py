from Subsystem_design.common_constants import Constants
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
        self.sweep_025 = 25  # Sweep at the quarter chord [deg]
        self.sweep_05 = 22.4  # Sweep at the half chord [deg]

        self.b_h = 2 * 6.12  # Span of the horizontal tail [m]
        self.c_r_h = 3.814  # Root chord of the horizontal tail [m]
        self.c_t_h = 1.186  # Tip chord of the horizontal tail [m]
        self.taper_h = self.c_t_h / self.c_r_h  # Taper ratio of the horizontal tail
        self.sweep_LE_h = 33  # Sweep of the LE of the horizontal tail [deg]

        self.MTOW_neo = 79000  # Maximum take off weight of A320neo [kg]
        self.RoC_neo = 12.7  # Max rate of climb of the A320neo [m/s]

        self.C_D_0_TO_neo = 0.078  # Zero-lift drag coefficient of A320 during take-off
        self.C_D_0_clean_neo = 0.023  # Zero-lift drag coefficient of A320 during cruise

        self.R_neo = 4800  # Harmonic range of A320neo [km]
        self.e = 0.992  # Oswald efficiency factor
        
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
        AR_no_sharklets = self.b**2 / self.S
        self.AR = AR_no_sharklets + 1.9 * (self.h_sharklet / self.b) * AR_no_sharklets

    def h_tail_MAC(self):
        """
        The MAC is computed using the ADSEE-II slides.
        :return: The length and position of the horizontal wing's MAC
        """
        self.mac_h = (2/3) * self.c_r_h * (1 + self.taper_h + self.taper_h**2) / (1 + self.taper_h)
        self.y_mac_h = (self.b_h / 6) * (1 + 2 * self.taper_h) / (1 + self.taper_h)
        self.x_mac_h = self.y_mac_h * np.tan(self.sweep_LE_h * np.pi / 180)

        print('\n MAC of the horizontal tail = ', self.mac_h, ' m',
              '\n y position of the MAC = ', self.y_mac_h, ' m',
              '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac_h, ' m')

    def Roskam_drag_prediction_cruise(self, rho, u1, l_f, l_cockpit, l_cabin, l_tail, AoA):
        """
        This function computes the drag from the fuselage according to the procedure given by Roskam at transonic
        conditions.
        :param rho: Air density [kg/m^3]
        :param u1: Airspeed [m/s]
        :param l_f: Length of the entire fuselage [m]
        :param l_cockpit: Length of the cockpit [m]
        :param l_cabin: Length of the cabin [m]
        :param l_tail: Length of the tail [m]
        :param AoA: Angle of attack of the aircraft [deg]
        :return: The zero lift drag coefficient of the fuselage and the drag coefficient of the fusselage due to lift
        """
        # Wetted area as computed in ADSEE-II
        S_fus = np.pi * 0.25 * self.height_f * self.width_f
        self.d_f = np.sqrt(4 * S_fus / np.pi)
        S_wet_fus = np.pi * self.d_f / 4 * \
                    (1 / (3 * l_cockpit**2) * ((4 * l_cockpit**2 + self.d_f**2 / 4) - self.d_f**3 / 8)
                     - self.d_f + 4 * l_cabin + 2 * np.sqrt(l_tail**2 + self.d_f**2 / 4))


        # Compute zero lift drag for M = 0.6 for fuselage exclusive of base
        R_n_fus = rho * u1 * l_f / self.visc
        # print('The Fuselage Reynolds Number R_f_fus is: ', R_n_fus, ' [-]')
        # print('The Mach number M is: ', self.M)
        R_wf = 1.015  # The wing/fuselage iterference factor from Figure 4.1 in Roskam-VI
        C_f_fus = 0.0016  # The turbulent flat plate skin friction coefficient from Figure 4.3 in Roskam-VI

        ld = self.l_f / self.d_f
        C_D_o_fus_exc_base = R_wf * C_f_fus * (1 + 60 / ld**3 + 0.0025 * ld) * S_wet_fus / self.S

        # Compute the fuselage base drag coefficient
        bf = np.sqrt(4 / np.pi * self.S_b_fus) / self.d_f
        C_D_b_fus = 0.9 * bf**2

        # The zero lift drag coefficient of the fuselage becomes:
        C_D_0_fus = C_D_o_fus_exc_base + C_D_b_fus

        # The lift induced drag
        C_D_L_fus = (AoA * np.pi / 180)**2 * self.S_b_fus / self.S

        return C_D_0_fus, C_D_L_fus

    def drag_increase_cruise(self, AoA_cruise):
        """
        Compute the drag increase due to the new fuselage length
        :param AoA_cruise: Angle of attack during cruise for the lift induced drag [deg]
        :return: The most important parameter which is computed is C_D_0_HACK
        """
        # Compute density at cruise altitude
        self.ISA_calculator(h_input=self.cruise_alt)

        # A320neo
        self.C_D_0_fus_neo, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a, l_f=self.l_f_320neo,
                                                                   l_cockpit=self.l_cockpit_320neo,
                                                                   l_cabin=self.l_cabin_320neo,
                                                                   l_tail=self.l_tail_320neo, AoA=AoA_cruise)

        # A320-HACK
        self.C_D_0_fus_HACK, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a, l_f=self.l_f,
                                                                    l_cockpit=self.l_cockpit, l_cabin=self.l_cabin,
                                                                    l_tail=self.l_tail, AoA=AoA_cruise)

        self.C_D_0_HACK = self.C_D_0_clean_neo - self.C_D_0_fus_neo + self.C_D_0_fus_HACK




    def L_over_D_cruise(self):
        self.ISA_calculator(h_input=self.cruise_alt)
        self.wing_AR()
        self.drag_increase_cruise(2)

        W_start_cruise = self.MTOW_320neo * (0.995 * 0.98)
        V = self.M * self.a

        self.C_L_start_cruise = W_start_cruise * self.g_0 / (0.5 * self.rho * V**2 * self.S)

        self.C_D_start_cruise_neo = self.C_D_0_clean_neo + self.C_L_start_cruise**2 / (np.pi * self.AR * self.e)
        self.C_D_start_cruise_HACK = self.C_D_0_HACK + self.C_L_start_cruise**2 / (np.pi * self.AR * self.e)

        self.D_start_cruise_HACK = self.C_D_start_cruise_HACK * 0.5 * self.rho * V**2 * self.S
        print(self.C_D_start_cruise_HACK, '*0.5*', self.rho, '*', V, '**2 * ', self.S)

        self.L_D_ratio_neo = self.C_L_start_cruise / self.C_D_start_cruise_neo
        self.L_D_ratio_HACK = self.C_L_start_cruise / self.C_D_start_cruise_HACK

    def lift_gradient(self, M):
        self.wing_AR()
        beta = np.sqrt(1 - M**2)
        self.CL_alpha = 2 * np.pi * self.AR / (2 + np.sqrt(4 + (self.AR * beta / 0.95)**2 *
                                                      (1 + (np.tan(self.sweep_05) / beta)**2)))

    def plot_lift_drag_characteristics(self):
        C_L_range = np.linspace(-5, 15, 500)
        C_D_range = self.C_D_0_HACK + C_L_range**2 / (np.pi * self.AR * self.e)
        CL_CD = C_L_range / C_D_range



# Try out the class

if __name__ == '__main__':
    ae = AerodynamicCharacteristics()

    ae.wing_MAC()
    ae.h_tail_MAC()
    ae.wing_AR()
    print('\n Wing AR = ', ae.AR)

    ae.drag_increase_cruise(AoA_cruise=2)

    print('\n The zero-lift drag coefficient of the fuselage of the A320neo = ', ae.C_D_0_fus_neo,
          '\n For the A32-HACK it is = ', ae.C_D_0_fus_HACK)

    print('\n Assuming the C_D_0 of the A320neo during cruise is = ', ae.C_D_0_clean_neo,
          '\n The C_D_0 of the A320-HACK now becomes = ', ae.C_D_0_HACK,
          '\n That is a ', (ae.C_D_0_HACK / ae.C_D_0_clean_neo - 1) * 100, '% increase')

    ae.L_over_D_cruise()
    print('\n The lift coefficient is = ', ae.C_L_start_cruise,
          '\n The drag coefficient of the A320-HACK during cruise is  = ', ae.C_D_start_cruise_HACK,
          '\n The drag then is = ', ae.D_start_cruise_HACK,
          '\n The L/D ratio is = ', ae.L_D_ratio_HACK)
