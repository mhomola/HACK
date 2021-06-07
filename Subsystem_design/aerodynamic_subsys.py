from Subsystem_design.common_constants import Constants
import numpy as np
import matplotlib.pyplot as plt
from Subsystem_design.Tank_Design.Main_PreliminaryTank import d_wing_pod, l_wing_pod


class AerodynamicCharacteristics(Constants):
    def __init__(self):
        super().__init__()

    def aero_functions(self, AoA_cruise):
        self.wing_MAC()
        self.wing_AR()
        self.h_tail_MAC()
        self.v_tail_MAC()
        self.drag_increase_cruise(AoA_cruise=AoA_cruise)
        self.lift_gradient_res()
        self.L_over_D_cruise()

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
        # print('\n MAC of the inner trapezoid = ', mac_in, ' m',
        #       '\n with y position of the MAC = ', y_mac_in, ' m',
        #       '\n and with surface area = ', S_in, ' m^2')
        #
        # print('\n MAC of the outer trapezoid = ', mac_out, ' m',
        #       '\n with y position of the MAC = ', 0.5 * self.b_in + y_mac_out, ' m',
        #       '\n and with surface area = ', S_out, ' m^2')
        #
        # print('\n MAC of the wing = ', self.mac, ' m',
        #       '\n y position of the MAC = ', self.y_mac, ' m',
        #       '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac, ' m')

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

        # print('\n MAC of the horizontal tail = ', self.mac_h, ' m',
        #       '\n y position of the MAC = ', self.y_mac_h, ' m',
        #       '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac_h, ' m')

    def v_tail_MAC(self):
        """
        The MAC is computed using the ADSEE-II slides.
        :return: The length and position of the horizontal wing's MAC
        """
        self.mac_v = (2/3) * self.c_r_v * (1 + self.taper_v + self.taper_v**2) / (1 + self.taper_v)
        self.y_mac_v = (self.b_v / 6) * (1 + 2 * self.taper_v) / (1 + self.taper_v)
        self.x_mac_v = self.y_mac_v * np.tan(self.sweep_LE_v * np.pi / 180)

        # print('\n MAC of the vertical tail = ', self.mac_v, ' m',
        #       '\n y position of the MAC = ', self.y_mac_v, ' m',
        #       '\n x position of the LEMAC measured from the start of the root chord = ', self.x_mac_v, ' m')

    def Roskam_drag_prediction_cruise(self, rho, u1, l_cockpit, l_cabin, l_tail, AoA, S_b_fus, height_f, width_f):
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
        S_fus = np.pi * 0.25 * height_f * width_f
        self.d_f = np.sqrt(4 * S_fus / np.pi)
        S_wet_fus = np.pi * self.d_f / 4 * \
                    (1 / (3 * l_cockpit**2) * ((4 * l_cockpit**2 + self.d_f**2 / 4) - self.d_f**3 / 8)
                     - self.d_f + 4 * l_cabin + 2 * np.sqrt(l_tail**2 + self.d_f**2 / 4))

        # Compute zero lift drag for M = 0.6 for fuselage exclusive of base
        l_f = l_cockpit + l_cabin + l_tail
        R_n_fus = rho * u1 * l_f / self.visc
        print('The Fuselage Reynolds Number R_f_fus is: ', R_n_fus, ' [-]')
        # print('The Mach number M is: ', self.M)
        if l_f > 20:
            R_wf = 1.015  # The wing/fuselage iterference factor from Figure 4.1 in Roskam-VI
            C_f_fus = 0.0016  # The turbulent flat plate skin friction coefficient from Figure 4.3 in Roskam-VI
        else:
            R_wf = 1
            C_f_fus = 0.0018


        ld = l_f / self.d_f
        C_D_o_fus_exc_base = R_wf * C_f_fus * (1 + 60 / ld**3 + 0.0025 * ld) * S_wet_fus / self.S
        # Compute the fuselage base drag coefficient
        bf = np.sqrt(4 / np.pi * S_b_fus) / self.d_f
        C_D_b_fus = 0.09 * bf**2

        # The zero lift drag coefficient of the fuselage becomes:
        C_D_0_fus = C_D_o_fus_exc_base + C_D_b_fus

        # The lift induced drag
        C_D_L_fus = (AoA * np.pi / 180)**2 * S_b_fus / self.S

        return C_D_0_fus, C_D_L_fus

    def DATCOM_podded_tanks(self, c_loc, x_fwd, d_tank, alpha):
        """
        With this function you can compute the aerodynamic effects of wing podded tanks on the aircraft.
        First it computes the lift increment (Delta CL), and then it computes the added drag.
        :return:
        """

        # Lift increment
        S_w_ab = (c_loc - x_fwd) * d_tank * (39.3701**2 / 1000)  # Pylon vertically projected area [1000 in.^2]
        L_R = -30 / (16.2 - 4.1) * S_w_ab + 5  # Incremental lift effect due to pylon [ft^2]
        K_sweep = 1  # Wing sweep factor
        K_H = 0  # Pylon height factor todo: Check if this value is correct with the plots from DATCOM
        K_X = 0  # Store-placement factor todo: Check if this value is correct with the plots from DATCOM
        K_W = 0  # Store-installation-width factor todo: Check if this value is correct with the plots from DATCOM
        K_P = 1 + K_H * K_W * K_X
        L_alpha_ws = 0.4  # Incremental wing-stores lift effect due to angle of attack [ft^2/deg]

        self.D_CL_tank = 1 / (self.S * 3.281**2) * (L_R * (K_sweep * K_P) + L_alpha_ws * (alpha - 4))

    def drag_engines(self):
        """
        In this function you you can compute the drag from the engines. This is used to know what the loads are
        which act on the fuselage.
        :return:
        """


    def drag_increase_cruise(self, AoA_cruise):
        """
        Compute the drag increase due to the new fuselage length
        :param AoA_cruise: Angle of attack during cruise for the lift induced drag [deg]
        :return: The most important parameter which is computed is C_D_0_HACK
        """
        # Compute density at cruise altitude
        self.ISA_calculator(h_input=self.cruise_altitude)

        # A320neo
        self.C_D_0_fus_neo, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a,
                                                                   l_cockpit=self.l_cockpit_320neo,
                                                                   l_cabin=self.l_cabin_320neo,
                                                                   l_tail=self.l_tail_320neo, AoA=AoA_cruise,
                                                                   S_b_fus=self.S_b_fus, height_f=self.height_f,
                                                                   width_f=self.width_f)

        # A320-HACK
        self.l_tank_nose, self.l_tank_body, self.l_tank_tail = d_wing_pod / 2, l_wing_pod - d_wing_pod, 2*d_wing_pod
        self.d_tank = d_wing_pod
        self.C_D_0_tank_HACK, _ = self.Roskam_drag_prediction_cruise(rho=self.rho, u1=self.M*self.a,
                                                                     l_cockpit=self.l_tank_nose,
                                                                     l_cabin=self.l_tank_body,
                                                                     l_tail=self.l_tank_tail, AoA=AoA_cruise,
                                                                     S_b_fus=0, height_f=self.d_tank,
                                                                     width_f=self.d_tank)

        self.C_D_0_tank_sys_HACK = 2 * 1.3 * self.C_D_0_tank_HACK
        self.C_D_0_HACK = self.C_D_0_clean_neo + self.C_D_0_tank_sys_HACK

    def L_over_D_cruise(self):
        self.ISA_calculator(h_input=self.cruise_altitude)
        self.wing_AR()
        self.drag_increase_cruise(2)
        self.DATCOM_podded_tanks(c_loc=4, x_fwd=1, d_tank=self.d_tank, alpha=2)

        W_start_cruise = self.MTOW_320neo * (0.995 * 0.98)
        V = self.M * self.a

        self.C_L_start_cruise = W_start_cruise * self.g_0 / (0.5 * self.rho * V**2 * self.S)

        self.C_D_start_cruise_neo = self.C_D_0_clean_neo + self.C_L_start_cruise**2 / (np.pi * self.AR * self.e)
        self.C_D_start_cruise_HACK = self.C_D_0_HACK + self.C_L_start_cruise**2 / (np.pi * self.AR * self.e)

        self.D_start_cruise_HACK = self.C_D_start_cruise_HACK * 0.5 * self.rho * V**2 * self.S

        self.L_D_ratio_neo = self.C_L_start_cruise / self.C_D_start_cruise_neo
        self.L_D_ratio_HACK = self.C_L_start_cruise / self.C_D_start_cruise_HACK

    def lift_gradient_calc(self, M, AR, sweep_05):
        beta = np.sqrt(1 - M**2)
        return 2 * np.pi * AR / (2 + np.sqrt(4 + (AR * beta / 0.95)**2 *
                                            (1 + (np.tan(sweep_05) / beta)**2)))

    def lift_gradient_res(self):
        self.CL_alpha_w = self.lift_gradient_calc(M=self.M, AR=self.AR, sweep_05=self.sweep_05*np.pi / 180)
        sweep_05_h_rad = np.arctan(np.tan(self.sweep_LE_h * np.pi / 180) +
                                   self.c_r_h / 2 / (self.b_h / 2) * (self.taper_h - 1))
        self.CL_alpha_h = self.lift_gradient_calc(M=self.M*self.Vh_V, AR=self.AR_h, sweep_05=sweep_05_h_rad)

    def plot_lift_drag_characteristics(self):
        C_L_range = np.linspace(-0.3, 1.5, 500)
        C_D_range_neo = self.C_D_0_clean_neo + C_L_range**2 / (np.pi * self.AR * self.e)
        C_D_range_HACK = self.C_D_0_HACK + C_L_range**2 / (np.pi * self.AR * self.e)
        CL_CD_neo = C_L_range / C_D_range_neo
        CL_CD_HACK = C_L_range / C_D_range_HACK

        fig1, ax1 = plt.subplots(1,1)
        ax1.plot(C_L_range, C_D_range_HACK, label='A320-HACK')
        ax1.plot(C_L_range, C_D_range_neo, label='A320neo', linestyle='--')
        ax1.legend(loc='best')
        ax1.set_xlabel(r'$C_L$', size=15)
        ax1.set_ylabel(r'$C_D$',size=15)
        # plt.show()

        fig2, ax2 = plt.subplots(1,1)
        ax2.plot(C_L_range, CL_CD_HACK, label='A320-HACK')
        ax2.plot(C_L_range, CL_CD_neo, label='A320neo', linestyle='--')
        ax2.legend(loc='best')
        ax2.set_xlabel(r'$C_L$', size=15)
        ax2.set_ylabel(r'$\frac{C_L}{C_D}$', size=15)
        plt.show()



# Try out the class

if __name__ == '__main__':
    ae = AerodynamicCharacteristics()

    ae.aero_functions(AoA_cruise=2)

    print('\n Wing AR = ', ae.AR)

    print('\n The zero-lift drag coefficient of the fuselage of the A320neo = ', ae.C_D_0_fus_neo,
          '\n For the A320-HACK the tanks produce = ', ae.C_D_0_tank_sys_HACK)

    print('\n Assuming the C_D_0 of the A320neo during cruise is = ', ae.C_D_0_clean_neo,
          '\n The C_D_0 of the A320-HACK now becomes = ', ae.C_D_0_HACK,
          '\n That is a ', (ae.C_D_0_HACK / ae.C_D_0_clean_neo - 1) * 100, '% increase')

    print('\n The lift coefficient is = ', ae.C_L_start_cruise,
          '\n The drag coefficient of the A320-HACK during cruise is  = ', ae.C_D_start_cruise_HACK,
          '\n The drag then is = ', ae.D_start_cruise_HACK,
          '\n The L/D ratio is = ', ae.L_D_ratio_HACK,
          '\n The C_L_alpha of the wing is = ', ae.CL_alpha_w,
          '\n The C_L_alpha of the horizontal tail is = ', ae.CL_alpha_h)


    print('\n rho = ', ae.rho,
          '\n V = ', ae.M * ae.a)

    print(ae.C_L_start_cruise)

    ae.plot_lift_drag_characteristics()

